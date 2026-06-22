"""Multi-stream (LNG / plate-fin) heat exchanger.

A single core that exchanges heat among an arbitrary number of passes — the
brazed-aluminium plate-fin exchanger of LNG liquefaction trains, cryogenic
air-separation cold boxes and gas-plant turbo-expander recovery (Hameed 2025,
*Chemical Process Simulations using Aspen HYSYS*, sec. 9.5.2, the "LNG"
operation). It is the N-stream generalisation of
:class:`caldyr.unitops.heat_exchanger.HeatExchanger`.

**Passes, not "hot"/"cold".** The unit takes ``passes`` — a list of dicts, one
per stream passing through the core — and exposes a pair of ports
``pass{i}_in`` / ``pass{i}_out`` for each. Whether a pass is *hot* (releases
heat) or *cold* (absorbs heat) is decided by the **sign of its duty in the
solution**, not by a label (HYSYS does the same: "if the designated hot pass is
actually cold, the operation will still proceed properly"). The engine sign
convention applies: a pass's duty ``Q = n·(H_out - H_in)`` is positive when the
pass is heated (cold pass) and negative when cooled (hot pass).

**Specs and degrees of freedom.** With ``N`` passes there are ``N`` unknown
outlet enthalpies. The overall energy balance ``Σ Q = -heat_loss`` is always one
equation, so ``N-1`` further specifications are needed. Each pass may carry one
of:

* ``T_out`` — a target outlet temperature, or
* ``duty`` — a signed enthalpy change (W; negative = cooled),

and optionally a ``dP`` pressure drop. Exactly **one** pass is left free and
closed by the energy balance (the "single unknown → direct from an energy
balance" case). A single **global** spec may additionally be given —

* ``min_approach`` — the minimum internal temperature approach (MITA, K), or
* ``ua`` — the conductance ``Σ ΔQ/ΔT_lm`` (W/K) —

in which case **two** passes are left free and the system is solved iteratively
for the energy balance *and* the constraint (the "multiple unknowns → iterative"
case). This mirrors the HYSYS LNG solver.

**Weighted (zone) method.** Once every pass duty is known, the hot and cold
**composite curves** are built by sampling each pass's temperature against its
cumulative duty with rigorous PH flashes (so phase changes — the kink where a
stream starts to condense or boil — are captured), then summing the passes over
shared temperature intervals (Kemp, *Pinch Analysis and Process Integration*,
2e, ch. 2; GPSA *Engineering Data Book* zone analysis). From the aligned
counter-current composites the unit reports:

* ``min_approach`` — the minimum internal temperature difference (MITA); a
  non-positive value is a second-law-infeasible temperature **cross** and raises
  :class:`MultiStreamExchangerError` rather than returning a wrong answer;
* ``UA`` — the required conductance, ``Σ ΔQ / LMTD`` over the composite
  intervals (the multi-stream generalisation of ``Q/LMTD``).

The full composite curves and per-pass duties/outlet temperatures live on
``unit.design`` so the cooling curve can be plotted.

The two-pass form reproduces :class:`HeatExchanger` exactly: same duty, and a
required ``UA`` equal to ``Q/LMTD``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register

#: Enthalpy samples per pass when building the final composite curve (captures
#: the nonlinear T(H) through a phase change).
_N_SEG = 30

#: Coarser sampling used while a global spec is being bracketed/root-found — the
#: accepted solution's curves are rebuilt once at full ``_N_SEG`` resolution.
_N_SEG_SCAN = 10

#: Duty (W) below which a pass is treated as carrying no heat (numerical floor).
_MIN_DUTY = 1e-6


class MultiStreamExchangerError(ValueError):
    """A multi-stream exchanger could not be configured or solved (bad spec
    count, an infeasible temperature cross, or a global constraint that cannot
    be met). Raised with diagnostics — never a silent wrong answer."""


@dataclass
class _Pass:
    """One resolved pass: inlet state plus its solved duty/outlet."""
    idx: int
    n: float                      # molar flow, mol/s
    z: dict[str, float]
    components: list[str]
    P_out: float
    H_in: float                   # J/mol
    Q: float = 0.0                # W, signed (n·(H_out - H_in)); + = heated
    # memoised (key, curve) so fixed passes are not re-flashed during a scan
    _curve_cache: tuple | None = None

    @property
    def H_out(self) -> float:
        return self.H_in + self.Q / self.n


@register("MultiStreamExchanger")
class MultiStreamExchanger(UnitOp):
    """N-pass LNG / plate-fin heat exchanger. See the module docstring."""

    def _pass_specs(self) -> list[dict]:
        passes = self.params.get("passes")
        if not isinstance(passes, list) or len(passes) < 2:
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: 'passes' must be a list of at "
                f"least two pass-spec dicts (got {passes!r})"
            )
        out = []
        for i, spec in enumerate(passes):
            if not isinstance(spec, dict):
                raise MultiStreamExchangerError(
                    f"MultiStreamExchanger {self.id!r}: pass {i + 1} spec must be a "
                    f"dict (got {spec!r})"
                )
            out.append(spec)
        return out

    def define_ports(self) -> list[Port]:
        # A genuinely absent 'passes' shows a 2-pass skeleton so the static
        # palette can render the unit; an explicitly-given but invalid 'passes'
        # is a real user error and is raised here as well as in solve().
        if self.params.get("passes") is None:
            n_passes = 2
        else:
            n_passes = len(self._pass_specs())
        ports: list[Port] = []
        for i in range(n_passes):
            ports.append(Port(f"pass{i + 1}_in", "inlet"))
            ports.append(Port(f"pass{i + 1}_out", "outlet"))
        return ports

    # -- main solve ------------------------------------------------------------
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        specs = self._pass_specs()
        heat_loss = float(self.params.get("heat_loss", 0.0))

        # Global constraint (consumes one extra degree of freedom).
        global_spec = None
        if self.params.get("min_approach") is not None:
            global_spec = ("min_approach", float(self.params["min_approach"]))
        if self.params.get("ua") is not None:
            if global_spec is not None:
                raise MultiStreamExchangerError(
                    f"MultiStreamExchanger {self.id!r}: give at most one global spec "
                    f"('min_approach' or 'ua'), not both"
                )
            global_spec = ("ua", float(self.params["ua"]))

        # Resolve inlets and per-pass fixed duties.
        passes: list[_Pass] = []
        fixed: dict[int, float] = {}     # pass index -> signed duty (W)
        free: list[int] = []
        for i, spec in enumerate(specs):
            s = inlets.get(f"pass{i + 1}_in")
            if s is None or not s.molar_flow:
                raise MultiStreamExchangerError(
                    f"MultiStreamExchanger {self.id!r}: missing/empty inlet "
                    f"'pass{i + 1}_in'"
                )
            T, P, n = s.require_state()
            z = s.normalized_z()
            P_out = P - float(spec.get("dP", 0.0))
            H_in = s.H if s.H is not None else pp.enthalpy(T, P, z)
            p = _Pass(idx=i, n=n, z=z, components=list(s.components),
                      P_out=P_out, H_in=H_in)
            passes.append(p)

            has_t = spec.get("T_out") is not None
            has_q = spec.get("duty") is not None
            if has_t and has_q:
                raise MultiStreamExchangerError(
                    f"MultiStreamExchanger {self.id!r}: pass {i + 1} has both "
                    f"'T_out' and 'duty' — give at most one"
                )
            if has_t:
                fixed[i] = n * (pp.enthalpy(float(spec["T_out"]), P_out, z) - H_in)
            elif has_q:
                fixed[i] = float(spec["duty"])
            else:
                free.append(i)

        n_free_needed = 2 if global_spec is not None else 1
        if len(free) != n_free_needed:
            kind = ("one global spec is set, so exactly two passes"
                    if global_spec is not None
                    else "no global spec is set, so exactly one pass")
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: {kind} must be left free "
                f"(without 'T_out'/'duty'); got {len(free)} free of "
                f"{len(specs)} passes. Each pass needs a spec except the free "
                f"one(s); add a 'min_approach' or 'ua' global spec to free a second."
            )

        for i, q in fixed.items():
            passes[i].Q = q
        sum_fixed = sum(fixed.values())

        if global_spec is None:
            # Single unknown: close the energy balance directly.
            passes[free[0]].Q = -heat_loss - sum_fixed
        else:
            self._solve_global(pp, passes, free, sum_fixed, heat_loss, global_spec)

        # Build composites and the performance numbers (also validates MITA > 0).
        perf = self._performance(pp, passes)

        # Outlet streams from rigorous PH flashes.
        outs: dict[str, PortStream] = {}
        pass_T_out: list[float] = []
        for p in passes:
            res = self._safe_flash(pp, p, p.H_out)
            pass_T_out.append(res.T)
            outs[f"pass{p.idx + 1}_out"] = Stream(
                id=f"{self.id}.pass{p.idx + 1}_out", components=p.components,
                T=res.T, P=res.P, molar_flow=p.n, z=p.z,
                H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
            )

        self.design = {
            "n_passes": len(passes),
            "pass_duties": [p.Q for p in passes],
            "pass_T_out": pass_T_out,
            "hot_duty": perf["hot_duty"],
            "cold_duty": perf["cold_duty"],
            "min_approach": perf["min_approach"],
            "UA": perf["UA"],
            "LMTD": perf["LMTD"],
            "hot_composite": perf["hot_composite"],
            "cold_composite": perf["cold_composite"],
            "heat_loss": heat_loss,
        }
        return outs

    # -- iterative global-spec close ------------------------------------------
    def _solve_global(self, pp, passes, free, sum_fixed, heat_loss, global_spec):
        """Two passes free, closed by the energy balance + a global constraint.

        The energy balance ties the two free duties: ``Q_b = -heat_loss -
        sum_fixed - Q_a``. So there is one free scalar ``Q_a``; vary it to meet
        the global constraint by a 1-D bracketed root-find (the MITA decreases
        and the UA increases monotonically as more heat is exchanged)."""
        a, b = free
        name, target = global_spec
        S = -heat_loss - sum_fixed     # Q_a + Q_b must equal S

        def metric(qa: float) -> float:
            passes[a].Q = qa
            passes[b].Q = S - qa
            perf = self._performance(pp, passes, allow_cross=True,
                                     n_seg=_N_SEG_SCAN)
            return perf["min_approach"] if name == "min_approach" else perf["UA"]

        def residual(qa: float) -> float:
            return metric(qa) - target

        # Bracket Q_a: scan a generous duty span (the largest load a single pass
        # could move between the coldest and hottest inlet in the bundle) for a
        # sign change in the residual.
        span = self._duty_span(pp, passes)
        lo, hi = -span, span
        n_scan = 31
        qs = [lo + (hi - lo) * k / (n_scan - 1) for k in range(n_scan)]
        prev_q = qs[0]
        try:
            prev_r = residual(prev_q)
        except Exception:                              # noqa: BLE001 - infeasible probe
            prev_r = None
        root = None
        for qa in qs[1:]:
            try:
                r = residual(qa)
            except Exception:                          # noqa: BLE001 - infeasible probe
                r = None
            if prev_r is not None and r is not None and prev_r == 0.0:
                root = prev_q
                break
            if prev_r is not None and r is not None and prev_r * r <= 0.0:
                root = self._brent(residual, prev_q, qa, prev_r, r)
                break
            prev_q, prev_r = qa, r

        if root is None:
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: could not satisfy the global "
                f"'{name}' = {target:g} spec within the feasible duty range. The "
                f"target may be unreachable (e.g. a minimum approach tighter than "
                f"the streams allow, or a UA larger than a pinch permits)."
            )
        # Lock in the solution and require feasibility (MITA > 0).
        passes[a].Q = root
        passes[b].Q = S - root

    @staticmethod
    def _brent(f, x0, x1, f0, f1, tol=1e-6, maxit=80):
        """Plain bisection/secant hybrid on a sign-changing bracket."""
        a, b, fa, fb = x0, x1, f0, f1
        for _ in range(maxit):
            # secant step, fall back to bisection if it leaves the bracket
            if fb != fa:
                m = b - fb * (b - a) / (fb - fa)
            else:
                m = 0.5 * (a + b)
            if not (min(a, b) < m < max(a, b)):
                m = 0.5 * (a + b)
            fm = f(m)
            if abs(fm) <= tol or abs(b - a) <= tol * max(1.0, abs(b)):
                return m
            if fa * fm <= 0.0:
                b, fb = m, fm
            else:
                a, fa = m, fm
        return 0.5 * (a + b)

    def _duty_span(self, pp, passes) -> float:
        """A duty magnitude that comfortably brackets any single pass's possible
        heat load: each pass taken between the coldest and hottest inlet in the
        bundle."""
        # Sample each pass inlet temperature.
        T_ins = []
        for p in passes:
            res = self._safe_flash(pp, p, p.H_in)
            T_ins.append(res.T)
        T_lo, T_hi = min(T_ins), max(T_ins)
        span = 0.0
        for p, T0 in zip(passes, T_ins):
            h_lo = pp.enthalpy(T_lo, p.P_out, p.z)
            h_hi = pp.enthalpy(T_hi, p.P_out, p.z)
            span = max(span, abs(p.n * (h_hi - h_lo)))
        return max(span, 1.0)

    # -- weighted (zone) composite-curve analysis -----------------------------
    def _performance(self, pp, passes, allow_cross: bool = False,
                     n_seg: int = _N_SEG) -> dict:
        """Build hot/cold composite curves by zone analysis and return the
        minimum approach, required UA and LMTD. Raises on a temperature cross
        unless ``allow_cross`` (used inside the global-spec root scan).

        ``n_seg`` is the per-pass enthalpy sampling; the global-spec scan passes
        a coarse value and the accepted solution is rebuilt at full resolution."""
        hot_samples, cold_samples = [], []
        hot_duty = cold_duty = 0.0
        for p in passes:
            if abs(p.Q) <= _MIN_DUTY:
                continue
            curve = self._pass_curve(pp, p, n_seg)   # [(T, cumulative |duty|)]
            if p.Q < 0.0:                        # cooled -> hot stream
                hot_samples.append(curve)
                hot_duty += -p.Q
            else:                                # heated -> cold stream
                cold_samples.append(curve)
                cold_duty += p.Q

        if not hot_samples or not cold_samples:
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: the exchanger needs at least "
                f"one heated and one cooled pass (got hot_duty={hot_duty:.4g} W, "
                f"cold_duty={cold_duty:.4g} W). Check the specs/energy balance."
            )

        hot_comp = _composite(hot_samples)       # [(H_cum from cold end, T)]
        cold_comp = _composite(cold_samples)
        approach, ua, lmtd = _zone_analysis(hot_comp, cold_comp)

        if approach <= 0.0 and not allow_cross:
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: infeasible — the hot and cold "
                f"composite curves cross (minimum internal approach "
                f"{approach:.2f} K <= 0). The duty cannot be transferred without a "
                f"temperature cross; relax an outlet spec or add a min_approach."
            )
        return {
            "hot_duty": hot_duty, "cold_duty": cold_duty,
            "min_approach": approach, "UA": ua, "LMTD": lmtd,
            "hot_composite": hot_comp, "cold_composite": cold_comp,
        }

    def _pass_curve(self, pp, p: _Pass, n_seg: int = _N_SEG) -> list[tuple[float, float]]:
        """Sample one pass: (T, cumulative |duty| released/absorbed from its
        inlet), monotonic in duty. PH flashes capture phase change. The result is
        memoised on the pass keyed by (duty, n_seg) so the fixed passes are not
        re-flashed at every step of the global-spec scan."""
        cached = p._curve_cache
        if cached is not None and cached[0] == (p.Q, n_seg):
            return cached[1]
        pts = []
        for k in range(n_seg + 1):
            f = k / n_seg
            H = p.H_in + f * (p.Q / p.n)
            res = self._safe_flash(pp, p, H)
            q_cum = abs(p.n * (H - p.H_in))
            pts.append((res.T, q_cum))
        p._curve_cache = ((p.Q, n_seg), pts)
        return pts

    def _safe_flash(self, pp, p: _Pass, H: float):
        """PH flash that turns a property-package convergence failure (thermo can
        fall over near a phase boundary or far off its grid for an extreme duty)
        into a typed :class:`MultiStreamExchangerError` — never a raw traceback."""
        try:
            return pp.flash_ph(p.P_out, H, p.z)
        except Exception as exc:                       # noqa: BLE001 - re-typed
            raise MultiStreamExchangerError(
                f"MultiStreamExchanger {self.id!r}: PH flash failed for pass "
                f"{p.idx + 1} at H={H:.1f} J/mol, P={p.P_out:.0f} Pa — the duty "
                f"likely drives the pass to an unphysical state, or the property "
                f"package cannot flash near its phase boundary ({type(exc).__name__})"
            ) from exc


# -- module-level composite / zone helpers (pure, unit-testable) -------------
def _interp_q_at_T(curve: list[tuple[float, float]], T: float) -> float:
    """Cumulative duty (from the pass inlet) at temperature ``T`` along a pass
    curve ``[(T, q)]``, clamped to the pass's temperature span."""
    Ts = [t for t, _ in curve]
    qs = [q for _, q in curve]
    lo, hi = min(Ts), max(Ts)
    if T <= lo:
        # q at the cold end of this pass
        return qs[Ts.index(lo)]
    if T >= hi:
        return qs[Ts.index(hi)]
    # piecewise-linear in T; curve T is monotonic but may ascend or descend.
    order = sorted(range(len(curve)), key=lambda i: Ts[i])
    Ts_s = [Ts[i] for i in order]
    qs_s = [qs[i] for i in order]
    for a, b in zip(range(len(Ts_s) - 1), range(1, len(Ts_s))):
        if Ts_s[a] <= T <= Ts_s[b]:
            if Ts_s[b] == Ts_s[a]:
                return qs_s[b]
            w = (T - Ts_s[a]) / (Ts_s[b] - Ts_s[a])
            return qs_s[a] + w * (qs_s[b] - qs_s[a])
    return qs_s[-1]


def _composite(samples: list[list[tuple[float, float]]]) -> list[tuple[float, float]]:
    """Composite curve for one side: cumulative duty (from the cold end) vs T.

    ``samples`` is a list of per-pass curves ``[(T, q_from_inlet)]``. Over each
    temperature interval, every pass overlapping it contributes the duty it
    exchanges within the interval; the cumulative sum from the coldest
    temperature upward is the composite curve, returned as ``[(H_cum, T)]``."""
    temps = sorted({round(t, 9) for curve in samples for t, _ in curve})
    points = [(0.0, temps[0])]
    h = 0.0
    for lo, hi in zip(temps[:-1], temps[1:]):
        dq = 0.0
        for curve in samples:
            Ts = [t for t, _ in curve]
            if min(Ts) >= hi or max(Ts) <= lo:
                continue
            dq += abs(_interp_q_at_T(curve, hi) - _interp_q_at_T(curve, lo))
        h += dq
        points.append((h, hi))
    return points


def _T_at_H(comp: list[tuple[float, float]], H: float) -> float:
    """Temperature at cumulative duty ``H`` along a composite ``[(H, T)]``."""
    Hs = [p[0] for p in comp]
    Ts = [p[1] for p in comp]
    if H <= Hs[0]:
        return Ts[0]
    if H >= Hs[-1]:
        return Ts[-1]
    for a, b in zip(range(len(Hs) - 1), range(1, len(Hs))):
        if Hs[a] <= H <= Hs[b]:
            if Hs[b] == Hs[a]:
                return Ts[b]
            w = (H - Hs[a]) / (Hs[b] - Hs[a])
            return Ts[a] + w * (Ts[b] - Ts[a])
    return Ts[-1]


def _zone_analysis(hot_comp, cold_comp, n_grid: int = 200):
    """Minimum approach (MITA), required UA and overall LMTD from the aligned
    counter-current composite curves.

    Both composites carry the same total duty (the exchanger is balanced). They
    are aligned from the cold end (H = 0) and compared on a shared duty grid;
    ``dT(H) = T_hot(H) - T_cold(H)`` is the local approach. ``UA = Σ ΔQ/ΔT_lm``
    over the grid intervals, the multi-stream generalisation of ``Q/LMTD``."""
    Q_total = min(hot_comp[-1][0], cold_comp[-1][0])
    grid = [Q_total * k / n_grid for k in range(n_grid + 1)]
    dT = [_T_at_H(hot_comp, q) - _T_at_H(cold_comp, q) for q in grid]
    approach = min(dT)

    ua = 0.0
    for k in range(n_grid):
        dq = grid[k + 1] - grid[k]
        d1, d2 = dT[k], dT[k + 1]
        if d1 <= 0.0 or d2 <= 0.0:
            # crossed interval: report the conductance as effectively unbounded
            # (the caller raises on approach <= 0 unless scanning).
            ua = math.inf
            continue
        if abs(d1 - d2) < 1e-9:
            lm = d1
        else:
            lm = (d1 - d2) / math.log(d1 / d2)
        ua += dq / lm

    lmtd = Q_total / ua if ua not in (0.0, math.inf) else float("nan")
    return approach, ua, lmtd
