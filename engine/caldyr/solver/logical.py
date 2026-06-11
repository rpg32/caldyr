"""Flowsheet-level logical operations: **Set** and **Adjust** (the HYSYS
SET/ADJUST equivalents). These are not unit ops — they live on the flowsheet
(``Flowsheet.logical``, persisted in `.flow` under the ``"logical"`` key) and
are interpreted by the solve entry point.

* **Set** — ``{"type": "set", "target": [unit, param], "source": [unit2,
  param2], "multiplier": k, "offset": b}``: before *every* solve, the target
  parameter is overwritten with ``k * source + b`` (k defaults to 1, b to 0).

* **Adjust** — ``{"type": "adjust", "vary": [unit, param], "bounds": [lo, hi],
  "spec": <metric>, "value": v, "tolerance": tol}``: vary the parameter inside
  its bounds until the metric hits ``value`` within ``tolerance``. ``<metric>``
  is a declarative spec (the same shape the optimizer/API use):

    {"type": "duty",           "stream": energy_stream_id}
    {"type": "T"|"P"|"molar_flow", "stream": stream_id}
    {"type": "component_rate", "stream": stream_id, "component": component_id}

  One Adjust is solved with a scalar Brent root find around the full flowsheet
  solve; several simultaneous Adjusts become a vector root find (scipy hybr).
  Iteration counts and outcomes are appended to ``SolveReport.messages``.

**Order of application** (documented contract): Sets are applied immediately
before *each* inner flowsheet solve — so an Adjust that varies the *source*
parameter of a Set drags the Set's target along on every iteration. Adjusts
wrap the whole solve in the outer root find. Because Adjust only wraps
``solve()``, it is backend-agnostic (sequential or equation-oriented).
"""
from __future__ import annotations

from typing import Any, Callable

import numpy as np
from scipy.optimize import brentq, root

from .sequential import SolveReport


class LogicalOpError(ValueError):
    """Malformed or unresolvable Set/Adjust specification."""


class AdjustError(LogicalOpError):
    """An Adjust could not be solved (no bracket, no convergence, ...)."""


# -- declarative metrics ------------------------------------------------------
def eval_metric(fs, report, spec: dict) -> float:
    """Evaluate a declarative metric spec against a solved flowsheet (see the
    module docstring for the supported shapes)."""
    mtype = spec.get("type")
    sid = spec.get("stream")
    if mtype == "duty":
        if sid not in report.duties:
            raise LogicalOpError(
                f"metric: no duty for energy stream {sid!r}; known duties: "
                f"{sorted(report.duties)}"
            )
        return float(report.duties[sid])
    stream = fs.streams.get(sid)
    if stream is None:
        raise LogicalOpError(
            f"metric: stream {sid!r} is not defined (known: {sorted(fs.streams)})"
        )
    if mtype == "component_rate":
        comp = spec.get("component")
        if comp not in stream.components:
            raise LogicalOpError(
                f"metric: component {comp!r} not in stream {sid!r} "
                f"(components: {stream.components})"
            )
        return float((stream.molar_flow or 0.0) * stream.normalized_z().get(comp, 0.0))
    if mtype in ("T", "P", "molar_flow"):
        value = getattr(stream, mtype)
        if value is None:
            raise LogicalOpError(f"metric: stream {sid!r} has no resolved {mtype}")
        return float(value)
    raise LogicalOpError(
        f"unknown metric type {mtype!r}; expected one of "
        f"'duty', 'T', 'P', 'molar_flow', 'component_rate'"
    )


# -- Set ----------------------------------------------------------------------
def _unit_param(fs, ref, role: str) -> tuple[str, str]:
    try:
        uid, param = ref
    except (TypeError, ValueError):
        raise LogicalOpError(
            f"logical op {role} must be a [unit_id, param] pair; got {ref!r}"
        ) from None
    if uid not in fs.units:
        raise LogicalOpError(
            f"logical op {role}: unknown unit {uid!r} (known: {sorted(fs.units)})"
        )
    return uid, param


def apply_sets(fs) -> list[str]:
    """Apply every Set in ``fs.logical``; returns one message per Set."""
    messages: list[str] = []
    for op in fs.logical:
        if op.get("type") != "set":
            continue
        tgt_u, tgt_p = _unit_param(fs, op.get("target"), "'set' target")
        src_u, src_p = _unit_param(fs, op.get("source"), "'set' source")
        if src_p not in fs.units[src_u].params:
            raise LogicalOpError(
                f"'set' source param {src_u}.{src_p} is not defined "
                f"(params: {sorted(fs.units[src_u].params)})"
            )
        k = float(op.get("multiplier", 1.0))
        b = float(op.get("offset", 0.0))
        value = k * float(fs.units[src_u].params[src_p]) + b
        fs.units[tgt_u].params[tgt_p] = value
        messages.append(f"set: {tgt_u}.{tgt_p} = {k:g}*{src_u}.{src_p} + {b:g} = {value:g}")
    return messages


# -- Adjust -------------------------------------------------------------------
def _read_adjust(fs, op: dict) -> tuple[str, str, float, float, float, float, dict]:
    uid, param = _unit_param(fs, op.get("vary"), "'adjust' vary")
    bounds = op.get("bounds")
    if not bounds or len(bounds) != 2 or not bounds[0] < bounds[1]:
        raise LogicalOpError(
            f"'adjust' on {uid}.{param}: bounds must be [lo, hi] with lo < hi; "
            f"got {bounds!r}"
        )
    spec = op.get("spec")
    if not isinstance(spec, dict):
        raise LogicalOpError(f"'adjust' on {uid}.{param}: missing metric 'spec'")
    if "value" not in op:
        raise LogicalOpError(f"'adjust' on {uid}.{param}: missing target 'value'")
    return (uid, param, float(bounds[0]), float(bounds[1]), float(op["value"]),
            float(op.get("tolerance", 1e-6)), spec)


def _solve_adjusts(fs, adjusts: list[dict], inner: Callable[[], SolveReport]) -> SolveReport:
    """Outer root find over the Adjust variables, wrapping the inner solve."""
    parsed = [_read_adjust(fs, op) for op in adjusts]
    counter = {"n": 0}
    last: dict[str, Any] = {}

    def solve_at(values: list[float]) -> SolveReport:
        for (uid, param, *_), v in zip(parsed, values):
            fs.units[uid].params[param] = float(v)
        report = inner()
        counter["n"] += 1
        last["report"] = report
        return report

    def residuals(report: SolveReport) -> list[float]:
        return [eval_metric(fs, report, spec) - target
                for (_, _, _, _, target, _, spec) in parsed]

    if len(parsed) == 1:
        uid, param, lo, hi, target, tol, spec = parsed[0]

        def f(x: float) -> float:
            return residuals(solve_at([x]))[0]

        f_lo, f_hi = f(lo), f(hi)
        if f_lo == 0.0:
            x_star = lo
        elif f_hi == 0.0:
            x_star = hi
        elif f_lo * f_hi > 0.0:
            raise AdjustError(
                f"adjust on {uid}.{param}: spec residual does not change sign over "
                f"bounds [{lo:g}, {hi:g}] (residual {f_lo:.4g} at lo, {f_hi:.4g} at "
                f"hi); the target value {target:g} is not reachable in these bounds"
            )
        else:
            x_star = float(brentq(f, lo, hi, xtol=1e-10 * max(abs(lo), abs(hi), 1.0)))
        resid = f(x_star)                      # final solve leaves fs at the root
        values = [x_star]
        resids, tols = [resid], [tol]
    else:
        los = np.array([p[2] for p in parsed])
        his = np.array([p[3] for p in parsed])
        x0 = np.array([_start(fs, p) for p in parsed])

        def vec(x: np.ndarray) -> np.ndarray:
            clipped = np.clip(x, los, his)
            rep = solve_at(list(clipped))
            scale = np.maximum(np.abs([p[4] for p in parsed]), 1.0)
            return np.asarray(residuals(rep)) / scale

        sol = root(vec, x0, method="hybr")
        values = [float(v) for v in np.clip(sol.x, los, his)]
        resids = residuals(solve_at(values))   # final solve at the clipped root
        tols = [p[5] for p in parsed]

    report: SolveReport = last["report"]
    ok = True
    for (uid, param, *_), v, r, tol in zip(parsed, values, resids, tols):
        met = abs(r) <= tol
        ok = ok and met
        report.messages.append(
            f"adjust: {uid}.{param} = {v:.8g}; |spec residual| = {abs(r):.3g} "
            f"({'within' if met else 'EXCEEDS'} tolerance {tol:g})"
        )
    report.messages.append(f"adjust: {counter['n']} flowsheet solves")
    if not ok:
        report.converged = False
    return report


def _start(fs, parsed_op) -> float:
    """Initial guess for an Adjust variable: the current param value if it lies
    in bounds, else the midpoint."""
    uid, param, lo, hi, *_ = parsed_op
    current = fs.units[uid].params.get(param)
    try:
        x = float(current)
    except (TypeError, ValueError):
        return 0.5 * (lo + hi)
    return x if lo <= x <= hi else 0.5 * (lo + hi)


# -- solve entry point ----------------------------------------------------------
def solve_flowsheet(fs, backend: str = "sequential", *, tol: float = 1e-6,
                    max_iter: int = 200, method: str = "wegstein",
                    on_iteration: Callable[[int, float], None] | None = None) -> SolveReport:
    """Solve a flowsheet, honoring its logical ops: apply Sets before each
    inner solve, then resolve Adjusts by an outer root find around the solve.
    The final report (with Set/Adjust messages merged in) is also stashed on
    ``fs.last_report`` for :func:`caldyr.solver.balance_report`."""

    def inner() -> SolveReport:
        set_msgs = apply_sets(fs)
        report = _backend_solve(fs, backend, tol=tol, max_iter=max_iter, method=method,
                                on_iteration=on_iteration)
        report.messages.extend(set_msgs)
        fs.last_report = report
        return report

    adjusts = [op for op in fs.logical if op.get("type") == "adjust"]
    unknown = [op for op in fs.logical if op.get("type") not in ("set", "adjust")]
    if unknown:
        raise LogicalOpError(
            f"unknown logical op type(s) {[op.get('type') for op in unknown]}; "
            f"expected 'set' or 'adjust'"
        )
    if not adjusts:
        return inner()
    return _solve_adjusts(fs, adjusts, inner)


def _backend_solve(fs, backend: str, *, tol: float, max_iter: int, method: str,
                   on_iteration: Callable[[int, float], None] | None = None) -> SolveReport:
    from ..thermo import make_package

    pp = make_package(fs.property_package, fs.component_ids)
    if backend == "sequential":
        from .sequential import SequentialModularSolver
        return SequentialModularSolver(tol=tol, max_iter=max_iter, method=method,
                                       on_iteration=on_iteration).solve(fs, pp)
    if backend == "equation_oriented":
        from .equation_oriented import EquationOrientedSolver
        return EquationOrientedSolver(tol=tol).solve(fs, pp)
    if backend == "pyomo":
        from .pyomo_backend import PyomoEOSolver
        return PyomoEOSolver(tol=tol).solve(fs, pp)
    raise ValueError(f"unknown solve backend {backend!r}")
