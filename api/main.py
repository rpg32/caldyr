"""Caldyr engine over HTTP. The GUI is a thin client of these endpoints.

    GET  /health
    GET  /unit-types           - registry + port definitions (the palette)
    GET  /property-packages    - selectable property packages
    POST /solve                - solve a .flow document, return resolved streams
    POST /cost                 - full techno-economic analysis of a flowsheet
    POST /optimize             - minimize/maximize a metric over unit parameters
    POST /flow/roundtrip       - parse + re-serialize a .flow document (validate)
    POST /property-table       - stream properties over a (T, P) grid
    POST /relief               - pressure-relief valve sizing (API 520/526)
    POST /pinch                - heat-integration pinch targets + composite curves
"""
from __future__ import annotations

import api  # noqa: F401  - path bootstrap so `caldyr` imports

import asyncio
import re

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware

from dataclasses import asdict, fields as dc_fields, replace as dc_replace

from caldyr.core import Flowsheet
from caldyr.economics import TEAConfig, analyze, monte_carlo, tornado
from caldyr.economics import data as econ_data
from caldyr.economics.data import CostFactors
from caldyr.economics.sizing import SizingOptions
from caldyr.io import from_dict, to_dict
from caldyr.solver import DesignVar, optimize
from caldyr.unitops import REGISTRY

from .param_schemas import PARAM_SCHEMAS
from .models import (
    CostRequest,
    EnvelopeRequest,
    MetricSpec,
    ObjectiveSpec,
    OptimizeRequest,
    PinchRequest,
    PropertyTableRequest,
    ReliefRequest,
    SolveRequest,
    SolveResponse,
)

app = FastAPI(title="Caldyr API", version="0.10.0",
              description="Open, scriptable process simulation + techno-economics.")

# The web dev server (Vite) runs on a different origin.
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"http://(localhost|127\.0\.0\.1):\d+",
    allow_methods=["*"],
    allow_headers=["*"],
)


# -- serialization ----------------------------------------------------------
def _stream_dict(s) -> dict:
    return {
        "id": s.id, "T": s.T, "P": s.P, "molar_flow": s.molar_flow,
        "z": s.z, "H": s.H, "phase": s.phase, "vapor_fraction": s.vapor_fraction,
    }


def _molar_mass_map(fs) -> dict:
    """Per-component molar mass (kg/mol) for the flowsheet's components, so the
    web client can render mass flow / mass fractions from molar data. Covers
    every species in the flowsheet (not just the curated catalog); skips any id
    the chemicals db can't resolve rather than failing the solve."""
    from caldyr.core.components_db import UnknownComponentError, molar_mass

    out: dict = {}
    for cid in fs.component_ids:
        try:
            out[cid] = molar_mass(cid)
        except UnknownComponentError:
            continue
    return out


def _designs_dict(fs) -> dict:
    """JSON-safe per-unit design results (column profiles, FUG numbers, ...)."""
    import json as _json

    out: dict = {}
    for uid, unit in fs.units.items():
        design = getattr(unit, "design", None)
        if not isinstance(design, dict) or not design:
            continue
        try:  # round-trip through json to coerce numpy scalars/arrays safely
            out[uid] = _json.loads(_json.dumps(design, default=float))
        except (TypeError, ValueError):
            continue
    return out


def _report_dict(r) -> dict:
    return {
        "converged": r.converged, "iterations": r.iterations, "residual": r.residual,
        "tol": r.tol, "method": r.method, "order": r.order,
        "tear_streams": r.tear_streams, "duties": r.duties, "messages": r.messages,
        "history": list(getattr(r, "history", []) or []),
    }


def _validated(flow: dict):
    """Structural gate (pydantic) + semantic load (engine io)."""
    from .models import FlowDocModel

    try:
        FlowDocModel.model_validate(flow)
    except Exception as exc:
        raise HTTPException(status_code=422, detail=f"malformed .flow: {exc}") from exc
    try:
        return from_dict(flow)
    except Exception as exc:
        raise HTTPException(status_code=422, detail=f"invalid .flow: {exc}") from exc


def _solve(flow: dict, backend: str, **kw):
    fs = _validated(flow)
    try:
        report = fs.solve(backend=backend, **kw)
    except Exception as exc:  # solver / physics failure
        raise HTTPException(status_code=400, detail=f"solve failed: {exc}") from exc
    return fs, report


# -- metadata ---------------------------------------------------------------
@app.get("/health")
def health() -> dict:
    return {"status": "ok", "service": app.title, "version": app.version}


def _first_paragraph(doc: str | None) -> str:
    """First blank-line-delimited paragraph of a docstring, whitespace-collapsed —
    a complete summary for the palette tooltip (not a mid-sentence first line)."""
    text = (doc or "").strip()
    if not text:
        return ""
    para = " ".join(text.split("\n\n")[0].split())
    return para[:600]


@app.get("/unit-types")
def unit_types() -> list[dict]:
    out = []
    for name, cls in sorted(REGISTRY.items()):
        try:
            ports = [{"name": p.name, "direction": p.direction, "kind": p.kind}
                     for p in cls(f"_{name}").ports]
        except Exception:
            ports = []
        out.append({"type": name, "doc": (cls.__doc__ or "").strip().split("\n")[0],
                    "description": _first_paragraph(cls.__doc__),
                    "ports": ports, "params": PARAM_SCHEMAS.get(name, [])})
    return out


@app.get("/components")
def components() -> list[dict]:
    """Curated catalog of common components (for the UI's autocomplete)."""
    from caldyr.core.components_db import COMMON_COMPONENTS

    return list(COMMON_COMPONENTS)


@app.get("/prices")
def prices() -> dict:
    """Default raw-material/product ($/kg) and utility ($/GJ) prices so the UI can
    pre-fill the costing assumptions (all overridable per flowsheet)."""
    return {
        "prices_per_kg": dict(econ_data.PRICES_PER_KG),
        "utility_prices": {n: u.price_per_GJ for n, u in econ_data.UTILITIES.items()},
        "prices_source": econ_data.PRICES_SOURCE,
    }


@app.post("/ports")
def unit_ports(body: dict) -> list[dict]:
    """Per-instance port list for a unit type with given params.

    Needed because some units derive ports from params (e.g. multi-feed
    columns); the static /unit-types palette shows only the defaults.
    """
    utype = body.get("type")
    if not isinstance(utype, str) or utype not in REGISTRY:
        raise HTTPException(422, f"unknown unit type {utype!r}")
    try:
        unit = REGISTRY[utype]("_ports", body.get("params") or {})
        return [{"name": p.name, "direction": p.direction, "kind": p.kind}
                for p in unit.ports]
    except Exception as exc:
        raise HTTPException(400, f"could not derive ports: {exc}") from exc


@app.get("/property-packages")
def property_packages() -> list[dict]:
    return [
        {"id": "thermo:PR", "name": "Peng-Robinson (cubic EOS)", "use": "non-polar systems"},
        {"id": "thermo:SRK", "name": "Soave-Redlich-Kwong (cubic EOS)", "use": "non-polar"},
        {"id": "thermo:NRTL", "name": "NRTL gamma-phi", "use": "polar mixtures / azeotropes"},
        {"id": "thermo:UNIFAC", "name": "Modified UNIFAC (Dortmund) gamma-phi",
         "use": "predictive VLE+LLE / heteroazeotropes (3-phase distillation)"},
        {"id": "thermo:UNIFAC-LLE", "name": "UNIFAC-LLE gamma-phi",
         "use": "liquid-liquid (Magnussen LLE table)"},
        {"id": "coolprop:Water", "name": "Steam tables (IAPWS-95 via CoolProp)",
         "use": "pure-water / steam systems"},
    ]


# -- core actions -----------------------------------------------------------
@app.post("/solve", response_model=SolveResponse)
def solve(req: SolveRequest) -> SolveResponse:
    kw = ({"tol": req.tol, "max_iter": req.max_iter, "method": req.method}
          if req.backend == "sequential" else {"tol": req.tol})
    fs, report = _solve(req.flow, req.backend, **kw)
    return SolveResponse.model_validate({
        "report": _report_dict(report),
        "streams": {sid: _stream_dict(s) for sid, s in fs.streams.items()},
        "designs": _designs_dict(fs),
        "molar_mass": _molar_mass_map(fs),
    })


def _apply_overrides(defaults, overrides: dict | None):
    """Return a copy of a dataclass instance with only its known fields overridden."""
    if not overrides:
        return defaults
    valid = {k: v for k, v in overrides.items()
             if k in {f.name for f in dc_fields(defaults)}}
    return dc_replace(defaults, **valid) if valid else defaults


_COST_CITATIONS = [
    {"topic": "Plant cost index (CEPCI escalation)", "source": econ_data.CEPCI_SOURCE},
    {"topic": "Raw-material / product prices", "source": econ_data.PRICES_SOURCE},
    {"topic": "Manufacturing cost COM_d + capital roll-up factors",
     "source": "Turton et al., Analysis Synthesis & Design of Chemical Processes 4e, Ch. 7-8"},
    {"topic": "Equipment purchased-cost correlations + bare-module factors",
     "source": "Turton et al. 4e, Appendix A / Ch. 7"},
    {"topic": "Sizing heuristics (U-values, residence, flooding, capacity)",
     "source": "Turton 4e; Seader, Separation Process Principles 3e; GPSA Engineering Data Book 13e"},
]


@app.get("/cost-defaults")
def cost_defaults() -> dict:
    """Default techno-economic assumptions (financial, sizing, factors) so the UI
    can offer a Settings editor before any cost has run. Prices live in /prices."""
    cfg = TEAConfig()
    return {
        "config": {
            "year": cfg.year, "operating_hours": cfg.operating_hours,
            "discount_rate": cfg.discount_rate, "project_years": cfg.project_years,
            "product_min_fraction": cfg.product_min_fraction,
        },
        "sizing": asdict(SizingOptions()),
        "factors": asdict(CostFactors()),
        "citations": _COST_CITATIONS,
    }


def _assumptions(cfg: TEAConfig, fs, res) -> dict:
    """The numbers + correlations that drove this cost result, with citations."""
    eff_prices = {**econ_data.PRICES_PER_KG, **(cfg.prices_per_kg or {})}
    eff_util = {n: u.price_per_GJ for n, u in econ_data.UTILITIES.items()}
    eff_util.update(cfg.utility_prices or {})
    comps = list(fs.component_ids)
    return {
        "config": {
            "year": cfg.year, "operating_hours": cfg.operating_hours,
            "discount_rate": cfg.discount_rate, "project_years": cfg.project_years,
            "product_component": cfg.product_component,
            "product_min_fraction": cfg.product_min_fraction,
        },
        # prices that can apply to this flowsheet's components (raw materials/product)
        "prices_per_kg": {c: eff_prices[c] for c in comps if c in eff_prices},
        # utility prices for the utilities actually selected during sizing
        "utility_prices": {s.utility: eff_util.get(s.utility)
                           for s in res.sizes if s.utility},
        "sizing": asdict(cfg.sizing or SizingOptions()),
        "factors": asdict(cfg.factors or CostFactors()),
        "citations": _COST_CITATIONS,
    }


@app.post("/cost")
def cost(req: CostRequest) -> dict:
    fs, report = _solve(req.flow, req.backend,
                        **({"tol": 1e-7} if req.backend == "sequential" else {}))
    cfg = TEAConfig(
        year=req.config.year, operating_hours=req.config.operating_hours,
        discount_rate=req.config.discount_rate, project_years=req.config.project_years,
        product_component=req.config.product_component,
        product_min_fraction=req.config.product_min_fraction,
        prices_per_kg=req.config.prices_per_kg,
        utility_prices=req.config.utility_prices,
        sizing=_apply_overrides(SizingOptions(), req.config.sizing),
        factors=_apply_overrides(CostFactors(), req.config.factors),
    )
    try:
        res = analyze(fs, report, cfg)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"costing failed: {exc}") from exc

    out: dict = {
        "report": _report_dict(report),
        "annual_production_kg": res.annual_production_kg,
        "annual_revenue": res.annual_revenue,
        "capital": {
            "isbl": res.capital.isbl, "osbl": res.capital.osbl,
            "grassroots": res.capital.grassroots,
            "working_capital": res.capital.working_capital, "tci": res.capital.tci,
        },
        "opex": {
            "raw_materials": res.opex.raw_materials, "utilities": res.opex.utilities,
            "fixed": res.opex.fixed, "total": res.opex.total,
        },
        "profitability": {
            "lcop": res.profitability.lcop, "npv": res.profitability.npv,
            "irr": res.profitability.irr, "payback_years": res.profitability.payback_years,
        },
        "equipment": [
            {"unit_id": c.unit_id, "type": c.equipment_type, "attribute": s.attribute,
             "attribute_name": s.attribute_name, "bare_module": c.bare_module,
             "utility": s.utility}
            for s, c in zip(res.sizes, res.costs)
        ],
        "assumptions": _assumptions(cfg, fs, res),
    }
    if req.tornado:
        out["tornado"] = [
            {"variable": b.variable, "low_lcop": b.low_lcop,
             "high_lcop": b.high_lcop, "swing": b.swing}
            for b in tornado(fs, res.sizes, res.config)
        ]
    if req.monte_carlo > 0:
        try:
            mc = monte_carlo(fs, res.sizes, res.config, n=req.monte_carlo)
        except Exception as exc:
            raise HTTPException(400, f"monte carlo failed: {exc}") from exc
        out["monte_carlo"] = {
            "n": req.monte_carlo,
            "lcop": mc.lcop,
            "npv": mc.npv,
            "lcop_samples": [round(float(x), 6) for x in
                             (mc.lcop_samples if mc.lcop_samples is not None else [])],
        }
    return out


@app.post("/envelope")
def envelope(req: EnvelopeRequest) -> dict:
    """Bubble/dew temperatures vs pressure for one stream's composition.

    Thin orchestration only: loops the engine's ``bubble_dew`` over a pressure
    grid; points where the flash fails are skipped.
    """
    try:
        fs = from_dict(req.flow)
    except Exception as exc:
        raise HTTPException(422, f"invalid .flow: {exc}") from exc
    s = fs.streams.get(req.stream)
    if s is None or not s.z or s.P is None:
        raise HTTPException(422,
            f"stream {req.stream!r} not found or has no composition/pressure "
            "(solve the flowsheet first)")
    from caldyr.thermo import make_package
    pp = make_package(fs.property_package, fs.component_ids)
    z = s.normalized_z()
    p_lo = req.p_min if req.p_min is not None else 0.2 * s.P
    p_hi = req.p_max if req.p_max is not None else 5.0 * s.P
    if not (0 < p_lo < p_hi):
        raise HTTPException(422, "need 0 < p_min < p_max")

    import math
    points = []
    for i in range(req.n):
        p = math.exp(math.log(p_lo) + (math.log(p_hi) - math.log(p_lo)) * i / (req.n - 1))
        try:
            t_bub, t_dew = pp.bubble_dew(p, z)
            points.append({"P": p, "T_bubble": t_bub, "T_dew": t_dew})
        except Exception:
            continue  # outside the envelope's feasible region
    if not points:
        raise HTTPException(400, "no envelope points could be computed for this stream")
    return {"stream": req.stream, "z": z, "points": points}


def _metric(fs, report, spec: MetricSpec) -> float:
    if spec.type == "duty":
        return report.duties[spec.stream]
    s = fs.streams[spec.stream]
    if spec.type == "flow":
        return s.molar_flow
    if spec.type == "component_rate":
        if not spec.component:
            raise HTTPException(422, "component_rate metric needs a component")
        return s.molar_flow * s.normalized_z().get(spec.component, 0.0)
    raise HTTPException(422, f"unknown metric type {spec.type!r}")


def _objective(spec: ObjectiveSpec):
    sign = 1.0 if spec.sense == "min" else -1.0
    return lambda fs, rep: sign * _metric(fs, rep, spec.metric)


@app.post("/optimize")
def optimize_endpoint(req: OptimizeRequest) -> dict:
    try:
        fs = from_dict(req.flow)
    except Exception as exc:
        raise HTTPException(422, f"invalid .flow: {exc}") from exc

    design = [DesignVar(d.unit_id, d.param, d.lower, d.upper, d.initial)
              for d in req.design_vars]
    cons = []
    for c in req.constraints:
        if c.op == ">=":
            cons.append(lambda fs, rep, c=c: _metric(fs, rep, c.metric) - c.value)
        else:
            cons.append(lambda fs, rep, c=c: c.value - _metric(fs, rep, c.metric))

    try:
        res = optimize(fs, _objective(req.objective), design, constraints=cons,
                       backend=req.backend, solve_kwargs={"tol": req.tol})
        report = fs.solve(backend=req.backend, tol=req.tol)
    except Exception as exc:
        raise HTTPException(400, f"optimization failed: {exc}") from exc

    return {
        "success": res.success, "objective": res.objective, "design": res.design,
        "n_solves": res.n_solves, "message": res.message,
        "report": _report_dict(report),
        "streams": {sid: _stream_dict(s) for sid, s in fs.streams.items()},
    }


@app.post("/flow/roundtrip")
def flow_roundtrip(flow: dict) -> dict:
    try:
        fs: Flowsheet = from_dict(flow)
        return to_dict(fs)
    except Exception as exc:
        raise HTTPException(422, f"invalid .flow: {exc}") from exc


# -- balance diagnostics ------------------------------------------------------
@app.post("/balance")
def balance(req: SolveRequest) -> dict:
    """Solve, then report per-unit and overall mass/energy closure."""
    from caldyr.solver import balance_report

    fs, report = _solve(req.flow, req.backend,
                        **({"tol": req.tol} if req.backend != "sequential"
                           else {"tol": req.tol, "max_iter": req.max_iter,
                                 "method": req.method}))
    try:
        bal = balance_report(fs)
    except Exception as exc:
        raise HTTPException(400, f"balance report failed: {exc}") from exc
    return {"report": _report_dict(report), "balance": bal}


# -- analysis tools: property table / relief / pinch -------------------------
@app.post("/property-table")
def property_table_endpoint(req: PropertyTableRequest) -> dict:
    """Stream properties over a (T, P) grid (HYSYS Property Table; §2.1.4).

    Composition from explicit ``z`` or a named ``stream`` in the .flow doc.
    Returns plot-ready 2-D arrays (NaN where a flash failed) + the failure list.
    """
    from caldyr.analysis.property_table import DEFAULT_PROPS, property_table
    from caldyr.thermo import make_package

    fs = _validated(req.flow)
    if req.z:
        z = req.z
    elif req.stream:
        s = fs.streams.get(req.stream)
        if s is None or not s.z:
            raise HTTPException(422,
                f"stream {req.stream!r} not found or has no composition "
                "(provide z, or solve the flowsheet first)")
        z = s.normalized_z()
    else:
        raise HTTPException(422, "provide either a stream name or an explicit z")

    pp = make_package(fs.property_package, fs.component_ids)
    props = req.props or list(DEFAULT_PROPS)
    try:
        grid = property_table(pp, z, T=req.T, P=req.P, props=props)
    except Exception as exc:
        raise HTTPException(400, f"property table failed: {exc}") from exc

    # numpy -> JSON-safe (None for NaN so the client gets nulls, not the string)
    def _clean(a) -> list:
        return [None if v != v else float(v) for v in a.tolist()]

    out: dict = {"T": grid["T"].tolist(), "P": grid["P"].tolist(),
                 "props": props, "z": z,
                 "failures": [{"T": t, "P": p, "error": m}
                              for (t, p, m) in grid["failures"]]}
    out["values"] = {name: [_clean(row) for row in grid[name]] for name in props}
    return out


@app.post("/relief")
def relief_endpoint(req: ReliefRequest) -> dict:
    """Pressure-relief valve sizing (API 520 area + API 526 orifice select)."""
    from caldyr.analysis.relief import relief_liquid, relief_vapor

    try:
        if req.phase == "vapor":
            M = req.M
            if M is None and req.flow and req.stream:
                fs = _validated(req.flow)
                s = fs.streams.get(req.stream)
                if s is not None and s.z:
                    from caldyr.core.components_db import molar_mass
                    M = sum(f * molar_mass(c) for c, f in s.normalized_z().items())
            if req.T is None or M is None or req.k is None:
                raise HTTPException(422,
                    "vapor relief needs T, k, and M (or flow+stream to derive M)")
            res = relief_vapor(
                W=req.W, T=req.T, M=M, Z=req.Z, k=req.k, P1=req.P1,
                backpressure=req.backpressure,
                **({"Kd": req.Kd} if req.Kd is not None else {}),
                Kb=req.Kb, Kc=req.Kc)
        else:
            if req.rho is None or req.P2 is None:
                raise HTTPException(422, "liquid relief needs rho and P2")
            res = relief_liquid(
                W=req.W, rho=req.rho, P1=req.P1, P2=req.P2,
                **({"Kd": req.Kd} if req.Kd is not None else {}),
                Kw=req.Kw, Kc=req.Kc, Kv=req.Kv)
    except HTTPException:
        raise
    except Exception as exc:  # ReliefSizingError + anything else
        raise HTTPException(400, f"relief sizing failed: {exc}") from exc

    return {
        "area_m2": res.area_m2, "area_cm2": res.area_m2 * 1e4,
        "orifice": res.orifice, "orifice_area_m2": res.orifice_area_m2,
        "capacity_used": res.capacity_used, "phase": res.phase,
        "critical": res.critical, "details": res.details, "notes": res.notes,
    }


@app.post("/pinch")
def pinch_endpoint(req: PinchRequest) -> dict:
    """Heat-integration pinch targeting on the solved flowsheet (§heat-integ)."""
    from caldyr.analysis.pinch import pinch_analysis

    fs, report = _solve(req.flow, req.backend, tol=req.tol)
    try:
        r = pinch_analysis(fs, report, dt_min=req.dt_min)
    except Exception as exc:
        raise HTTPException(400, f"pinch analysis failed: {exc}") from exc
    return {
        "report": _report_dict(report),
        "dt_min": r.dt_min,
        "qh_min": r.qh_min, "qc_min": r.qc_min,
        "pinch_T_hot": r.pinch_T_hot, "pinch_T_cold": r.pinch_T_cold,
        "pinch_T_shifted": r.pinch_T_shifted,
        "current_hot_utility": r.current_hot_utility,
        "current_cold_utility": r.current_cold_utility,
        "heat_recovery_potential": r.heat_recovery_potential,
        "hot_composite": [list(p) for p in r.hot_composite],
        "cold_composite": [list(p) for p in r.cold_composite],
        "streams": [{"id": s.unit_id, "T_in": s.T_in, "T_out": s.T_out,
                     "Q": s.Q, "kind": s.kind} for s in r.streams],
    }


# -- AI: direct tool dispatch (no LLM — fallback + one-click buttons) ---------
@app.post("/ai/tool")
def ai_tool(body: dict) -> dict:
    """Run one caldyr.ai tool statelessly against a .flow document.

    body: {"name": str, "args": {}, "flow": {...}, "solve": bool}
    If solve=true (default), the flowsheet is solved first so report-dependent
    tools (explain_convergence) have something to look at.
    """
    from caldyr.ai.session import AgentSession
    from caldyr.ai.tools import dispatch

    name = body.get("name")
    if not isinstance(name, str):
        raise HTTPException(422, "need a tool 'name'")
    session = AgentSession()
    flow = body.get("flow")
    if flow is not None:
        try:
            session.flowsheet = from_dict(flow)
        except Exception as exc:
            raise HTTPException(422, f"invalid .flow: {exc}") from exc
        if body.get("solve", True):
            solve_out = dispatch(session, "solve", {})
            if not solve_out.get("ok"):
                # explain_convergence on a failed solve is still meaningful;
                # other tools get the error surfaced.
                if name != "explain_convergence":
                    raise HTTPException(400, f"pre-solve failed: {solve_out.get('error')}")
    out = dispatch(session, name, body.get("args") or {})
    if not out.get("ok"):
        raise HTTPException(400, str(out.get("error")))
    return out


# Browsers do not apply CORS to WebSocket handshakes, so a malicious page could
# open ws://localhost:<port>/... while the API runs (cross-site WS hijacking).
# Browser connections must come from a local origin; non-browser clients send
# no Origin header and are allowed.
_LOCAL_ORIGIN = re.compile(r"^https?://(localhost|127\.0\.0\.1)(:\d+)?$")


async def _accept_local(ws: WebSocket) -> bool:
    origin = ws.headers.get("origin")
    if origin is not None and not _LOCAL_ORIGIN.match(origin):
        await ws.close(code=1008)  # policy violation
        return False
    await ws.accept()
    return True


# -- WebSocket: solve with live iteration progress ----------------------------
@app.websocket("/ws/solve")
async def ws_solve(ws: WebSocket) -> None:
    """One solve per message: client sends {flow, backend?, tol?}; server
    streams {type:"iteration"|"result"|"error", ...} and awaits the next job."""
    if not await _accept_local(ws):
        return
    loop = asyncio.get_running_loop()
    try:
        while True:
            req = await ws.receive_json()
            queue: asyncio.Queue = asyncio.Queue()

            def on_iter(i: int, residual: float) -> None:
                loop.call_soon_threadsafe(
                    queue.put_nowait, {"type": "iteration", "iteration": i,
                                       "residual": residual})

            def job() -> dict:
                fs = from_dict(req["flow"])
                report = fs.solve(backend=req.get("backend", "sequential"),
                                  tol=float(req.get("tol", 1e-6)),
                                  on_iteration=on_iter)
                return {"type": "result",
                        "report": _report_dict(report),
                        "streams": {sid: _stream_dict(s)
                                    for sid, s in fs.streams.items()},
                        "designs": _designs_dict(fs),
                        "molar_mass": _molar_mass_map(fs)}

            task = loop.run_in_executor(None, job)
            while True:
                done = task.done()
                while not queue.empty():
                    await ws.send_json(queue.get_nowait())
                if done:
                    break
                await asyncio.sleep(0.05)
            try:
                await ws.send_json(task.result())
            except Exception as exc:
                await ws.send_json({"type": "error", "detail": f"solve failed: {exc}"})
    except WebSocketDisconnect:
        return


# -- WebSocket: AI chat over the canvas flowsheet -----------------------------
@app.websocket("/ws/chat")
async def ws_chat(ws: WebSocket) -> None:
    """Persistent chat. Client messages: {text, flow?, provider?, model?}.
    Server streams ChatAgent events plus {"type":"done", text, flow} /
    {"type":"error", detail}. The agent (and its conversation) lives as long as
    the socket."""
    from caldyr.ai.chat import ChatAgent

    if not await _accept_local(ws):
        return
    loop = asyncio.get_running_loop()
    agent: ChatAgent | None = None
    try:
        while True:
            msg = await ws.receive_json()
            if agent is None:
                try:
                    agent = ChatAgent(provider=msg.get("provider"),
                                      model=msg.get("model"))
                except Exception as exc:
                    await ws.send_json({"type": "error",
                                        "detail": f"no LLM backend: {exc}"})
                    continue
            if msg.get("flow") is not None:
                try:
                    agent.load_flow(msg["flow"])
                except Exception as exc:
                    await ws.send_json({"type": "error",
                                        "detail": f"invalid .flow: {exc}"})
                    continue

            queue: asyncio.Queue = asyncio.Queue()
            emit = lambda e: loop.call_soon_threadsafe(queue.put_nowait, e)  # noqa: E731
            task = loop.run_in_executor(
                None, lambda: agent.send(str(msg.get("text", "")), on_event=emit))
            while True:
                done = task.done()
                while not queue.empty():
                    await ws.send_json(queue.get_nowait())
                if done:
                    break
                await asyncio.sleep(0.05)
            try:
                out = task.result()
                await ws.send_json({"type": "done", "text": out["text"],
                                    "flow": out["flow"],
                                    "tool_calls": out["tool_calls"]})
            except Exception as exc:
                await ws.send_json({"type": "error", "detail": f"chat failed: {exc}"})
    except WebSocketDisconnect:
        return
