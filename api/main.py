"""Caldyr engine over HTTP. The GUI is a thin client of these endpoints.

    GET  /health
    GET  /unit-types           - registry + port definitions (the palette)
    GET  /property-packages    - selectable property packages
    POST /solve                - solve a .flow document, return resolved streams
    POST /cost                 - full techno-economic analysis of a flowsheet
    POST /optimize             - minimize/maximize a metric over unit parameters
    POST /flow/roundtrip       - parse + re-serialize a .flow document (validate)
"""
from __future__ import annotations

import api  # noqa: F401  - path bootstrap so `caldyr` imports

import asyncio

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware

from caldyr.core import Flowsheet
from caldyr.economics import TEAConfig, analyze, monte_carlo, tornado
from caldyr.io import from_dict, to_dict
from caldyr.solver import DesignVar, optimize
from caldyr.unitops import REGISTRY

from .models import (
    CostRequest,
    EnvelopeRequest,
    MetricSpec,
    ObjectiveSpec,
    OptimizeRequest,
    SolveRequest,
    SolveResponse,
)

app = FastAPI(title="Caldyr API", version="0.9.0",
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
                    "ports": ports})
    return out


@app.get("/components")
def components() -> list[dict]:
    """Curated catalog of common components (for the UI's autocomplete)."""
    from caldyr.core.components_db import COMMON_COMPONENTS

    return list(COMMON_COMPONENTS)


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
    })


@app.post("/cost")
def cost(req: CostRequest) -> dict:
    fs, report = _solve(req.flow, req.backend,
                        **({"tol": 1e-7} if req.backend == "sequential" else {}))
    cfg = TEAConfig(
        year=req.config.year, operating_hours=req.config.operating_hours,
        discount_rate=req.config.discount_rate, project_years=req.config.project_years,
        product_component=req.config.product_component,
        prices_per_kg=req.config.prices_per_kg,
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


# -- WebSocket: solve with live iteration progress ----------------------------
@app.websocket("/ws/solve")
async def ws_solve(ws: WebSocket) -> None:
    """One solve per message: client sends {flow, backend?, tol?}; server
    streams {type:"iteration"|"result"|"error", ...} and awaits the next job."""
    await ws.accept()
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
                        "designs": _designs_dict(fs)}

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

    await ws.accept()
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
