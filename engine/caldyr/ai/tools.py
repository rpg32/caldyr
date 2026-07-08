"""Typed tools over the engine — the AI-native interface.

Each tool has a JSON-schema input (so it drops straight into Anthropic tool-use
or MCP) and a handler that mutates / reads an :class:`AgentSession`. Handlers
return a dict with a human-readable ``summary`` plus structured data, so the
model gets both a sentence to reason over and machine-usable fields.

`dispatch(session, name, args)` runs one tool; `anthropic_tools()` exports the
schemas in the shape the Anthropic SDK expects.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any, Callable

from ..core import Component
from ..core import Flowsheet as _Flowsheet
from ..economics import TEAConfig, analyze
from ..io import from_dict, to_dict
from ..solver import DesignVar, optimize
from ..thermo import AVAILABLE_PACKAGES
from ..unitops import REGISTRY
from .session import AgentSession

Handler = Callable[[AgentSession, dict], dict]


@dataclass(frozen=True)
class Tool:
    name: str
    description: str
    input_schema: dict
    handler: Handler


# -- small helpers ----------------------------------------------------------
_OBJ = {"type": "object"}


def _check_port(fs, endpoint) -> None:
    """Validate a 'UNIT:port' endpoint against the flowsheet (None = boundary)."""
    if endpoint is None:
        return
    if ":" not in endpoint:
        raise ValueError(f"endpoint must be 'UNIT:port' or null, got {endpoint!r}")
    uid, port = endpoint.split(":", 1)
    unit = fs.units.get(uid)
    if unit is None:
        raise ValueError(f"no unit {uid!r}; units so far: {list(fs.units)}")
    names = [p.name for p in unit.ports]
    if port not in names:
        raise ValueError(f"{type(unit).__name__} {uid!r} has no port {port!r}; "
                         f"its ports are {names}")


def _stream_rows(fs) -> list[dict]:
    rows = []
    for sid, s in fs.streams.items():
        if s.molar_flow is None:
            continue
        rows.append({"id": sid, "T": s.T, "P": s.P, "molar_flow": s.molar_flow,
                     "phase": s.phase, "z": s.z})
    return rows


def _stream_table_text(fs) -> str:
    lines = [f"{'stream':>10} {'T/K':>8} {'n/(mol/s)':>10} {'phase':>6}"]
    for r in _stream_rows(fs):
        lines.append(f"{r['id']:>10} {r['T']:>8.2f} {r['molar_flow']:>10.3f} {r['phase'] or '':>6}")
    return "\n".join(lines)


# -- handlers ---------------------------------------------------------------
def _ports_of(cls, name: str) -> list[dict]:
    return [{"name": p.name, "direction": p.direction, "kind": p.kind}
            for p in cls(f"_{name}").ports]


def _list_unit_types(_s: AgentSession, a: dict) -> dict:
    utype = a.get("type")
    if utype:                              # one type, full parameter documentation
        cls = REGISTRY.get(utype)
        if cls is None:
            raise ValueError(f"unknown unit type {utype!r}; available: {sorted(REGISTRY)}")
        return {"summary": f"Full documentation for {utype}.", "type": utype,
                "doc": (cls.__doc__ or "").strip(), "ports": _ports_of(cls, utype)}
    items = [{"type": name, "doc": (cls.__doc__ or "").strip().split("\n")[0],
              "ports": _ports_of(cls, name)}
             for name, cls in sorted(REGISTRY.items())]
    return {"summary": f"{len(items)} unit types available.", "unit_types": items}


def _list_property_packages(_s: AgentSession, _a: dict) -> dict:
    pkgs = [dict(p) for p in AVAILABLE_PACKAGES]
    return {"summary": f"{len(pkgs)} property packages.", "property_packages": pkgs}


def _invalidate(s: AgentSession) -> None:
    """A structural edit makes any prior solve/cost stale."""
    s.report = None
    s.tea = None


def _new_flowsheet(s: AgentSession, a: dict) -> dict:
    comps = [Component(c) for c in a["components"]]
    s.flowsheet = _Flowsheet(components=comps,
                             property_package=a.get("property_package", "thermo:PR"))
    s.report = None
    s.tea = None
    return {"summary": f"New flowsheet with components {a['components']} "
                       f"and package {s.flowsheet.property_package}."}


def _add_unit(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    utype = a["type"]
    if utype not in REGISTRY:
        raise ValueError(f"unknown unit type {utype!r}; available: {sorted(REGISTRY)}")
    # own the params: the caller's dict must not alias the unit's live state
    unit = fs.add(REGISTRY[utype](a["id"], copy.deepcopy(a.get("params") or {})))
    _invalidate(s)
    ports = [{"name": p.name, "direction": p.direction, "kind": p.kind} for p in unit.ports]
    return {"summary": f"Added {utype} '{a['id']}' with ports "
                       f"{[p['name'] for p in ports]}.", "ports": ports}


def _add_feed(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    _check_port(fs, a["to"])
    z = a["z"]
    total = sum(z.values())
    if total <= 0:
        raise ValueError("feed composition sums to <= 0")
    z = {k: v / total for k, v in z.items()}          # tolerate un-normalized input
    fs.feed(a["id"], a["to"], T=a["T"], P=a["P"], molar_flow=a["molar_flow"], z=z)
    _invalidate(s)
    return {"summary": f"Feed '{a['id']}' -> {a['to']} "
                       f"({a['molar_flow']} mol/s at {a['T']} K, {a['P']} Pa). "
                       f"(add_feed already wires the feed; do not also connect it.)"}


def _connect(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    _check_port(fs, a.get("from"))
    _check_port(fs, a.get("to"))
    fs.connect(a["id"], a.get("from"), a.get("to"))
    _invalidate(s)
    return {"summary": f"Stream '{a['id']}': {a.get('from')} -> {a.get('to')}."}


def _solve(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    backend = a.get("backend", "sequential")
    kw = {"tol": a.get("tol", 1e-6)}
    if backend == "sequential":
        kw["max_iter"] = a.get("max_iter", 200)
    report = fs.solve(backend=backend, **kw)
    s.report = report
    summary = (f"{'Converged' if report.converged else 'DID NOT converge'} "
               f"via {report.method} in {report.iterations} iterations"
               f"{f' (tear: {report.tear_streams})' if report.tear_streams else ''}.")
    return {"summary": summary, "converged": report.converged,
            "iterations": report.iterations, "residual": report.residual,
            "duties": report.duties, "streams": _stream_rows(fs),
            "stream_table": _stream_table_text(fs)}


def _cost(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    if s.report is None:
        _solve(s, {"tol": 1e-7})
    cfg = TEAConfig(product_component=a.get("product_component", "ammonia"),
                    operating_hours=a.get("operating_hours", 8000.0),
                    discount_rate=a.get("discount_rate", 0.10),
                    prices_per_kg=a.get("prices_per_kg"))
    res = analyze(fs, s.report, cfg)
    s.tea = res
    p = res.profitability
    return {
        "summary": (f"LCOP ${p.lcop:.3f}/kg, TCI ${res.capital.tci:,.0f}, "
                    f"OPEX ${res.opex.total:,.0f}/yr, NPV ${p.npv:,.0f}."),
        "lcop": p.lcop, "npv": p.npv, "irr": p.irr, "payback_years": p.payback_years,
        "capital": {"isbl": res.capital.isbl, "tci": res.capital.tci},
        "opex": {"total": res.opex.total, "raw_materials": res.opex.raw_materials,
                 "utilities": res.opex.utilities, "fixed": res.opex.fixed},
        "annual_production_kg": res.annual_production_kg,
        "equipment": [{"unit_id": c.unit_id, "type": c.equipment_type,
                       "bare_module": c.bare_module} for c in res.costs],
    }


def _metric(fs, rep, spec: dict) -> float:
    """Evaluate one observable: {type:'duty'|'flow'|'component_rate', stream, component?}."""
    if spec["type"] == "duty":
        return rep.duties[spec["stream"]]
    st = fs.streams[spec["stream"]]
    if spec["type"] == "flow":
        return st.molar_flow
    return st.molar_flow * st.normalized_z().get(spec["component"], 0.0)


def _optimize(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    design = [DesignVar(d["unit_id"], d["param"], d["lower"], d["upper"], d.get("initial"))
              for d in a["design_vars"]]

    obj_spec = a["objective"]
    sign = 1.0 if obj_spec.get("sense", "min") == "min" else -1.0
    cons = []
    for c in a.get("constraints", []):
        if c["op"] == ">=":
            cons.append(lambda fs_, rep, c=c: _metric(fs_, rep, c["metric"]) - c["value"])
        else:
            cons.append(lambda fs_, rep, c=c: c["value"] - _metric(fs_, rep, c["metric"]))

    res = optimize(fs, lambda fs_, rep: sign * _metric(fs_, rep, obj_spec["metric"]),
                   design, constraints=cons, backend=a.get("backend", "sequential"))
    s.report = fs.solve(backend=a.get("backend", "sequential"))
    return {"summary": f"Optimization {'succeeded' if res.success else 'failed'}; "
                       f"design {res.design}.", "success": res.success,
            "objective": res.objective, "design": res.design, "n_solves": res.n_solves}


def _stream_table(s: AgentSession, _a: dict) -> dict:
    fs = s.require_flowsheet()
    s.require_solved()
    return {"summary": "Current stream table.", "stream_table": _stream_table_text(fs),
            "streams": _stream_rows(fs)}


def _export_flow(s: AgentSession, _a: dict) -> dict:
    fs = s.require_flowsheet()
    return {"summary": "Flowsheet as a .flow document.", "flow": to_dict(fs)}


def _remove_unit(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    removed = fs.remove_unit(a["id"])
    _invalidate(s)
    tail = f" and its streams {removed}" if removed else ""
    return {"summary": f"Removed unit '{a['id']}'{tail}.", "removed_streams": removed}


def _remove_stream(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    fs.remove_stream(a["id"])
    _invalidate(s)
    return {"summary": f"Removed stream '{a['id']}'."}


def _set_param(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    unit = fs.units.get(a["unit_id"])
    if unit is None:
        raise ValueError(f"no unit {a['unit_id']!r}; units so far: {list(fs.units)}")
    key, value = a["param"], a.get("value")
    old = unit.params.get(key)
    if value is None:
        unit.params.pop(key, None)
        summary = f"Cleared {a['unit_id']}.{key} (was {old!r}; back to its default)."
    else:
        unit.params[key] = value
        summary = f"Set {a['unit_id']}.{key} = {value!r} (was {old!r})."
    _invalidate(s)
    return {"summary": summary, "params": unit.params}


def _pinch_analysis(s: AgentSession, a: dict) -> dict:
    from ..analysis.pinch import pinch_analysis
    fs = s.require_flowsheet()
    if s.report is None:
        _solve(s, {"tol": 1e-7})
    r = pinch_analysis(fs, s.report, dt_min=a.get("dt_min", 10.0))
    pinch_txt = (f"pinch at {r.pinch_T_hot:.1f} K hot / {r.pinch_T_cold:.1f} K cold"
                 if r.pinch_T_hot is not None else "no pinch (threshold problem)")
    return {
        "summary": (f"dt_min {r.dt_min:g} K: minimum utilities {r.qh_min:,.0f} W hot / "
                    f"{r.qc_min:,.0f} W cold; {pinch_txt}. Current use "
                    f"{r.current_hot_utility:,.0f} W hot / {r.current_cold_utility:,.0f} W "
                    f"cold -> heat-recovery potential {r.heat_recovery_potential:,.0f} W."),
        "dt_min": r.dt_min, "qh_min": r.qh_min, "qc_min": r.qc_min,
        "pinch_T_hot": r.pinch_T_hot, "pinch_T_cold": r.pinch_T_cold,
        "pinch_T_shifted": r.pinch_T_shifted,
        "current_hot_utility": r.current_hot_utility,
        "current_cold_utility": r.current_cold_utility,
        "heat_recovery_potential": r.heat_recovery_potential,
        "streams": [{"id": t.unit_id, "T_in": t.T_in, "T_out": t.T_out,
                     "Q": t.Q, "kind": t.kind} for t in r.streams],
    }


def _sweep_parameter(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    uid, param, spec = a["unit_id"], a["param"], a["metric"]
    if uid not in fs.units:
        raise ValueError(f"no unit {uid!r}; units so far: {list(fs.units)}")
    steps = int(a.get("steps", 8))
    if not 2 <= steps <= 50:
        raise ValueError(f"steps must be between 2 and 50, got {steps}")
    lo, hi = float(a["start"]), float(a["stop"])
    if lo == hi:
        raise ValueError("start and stop must differ")
    if spec.get("type") != "duty":       # catch a typo'd stream up front, not as
        known = {c.stream_id for c in fs.connections}     # an all-points-failed sweep
        if spec.get("stream") not in known:
            raise ValueError(f"metric stream {spec.get('stream')!r} not found; "
                             f"streams: {sorted(known)}")

    doc = to_dict(fs)                # each point solves a fresh deep copy: the
    points: list[dict] = []          # working flowsheet is never touched
    failures = 0                     # (to_dict/from_dict share nested dicts)
    for i in range(steps):
        x = lo + (hi - lo) * i / (steps - 1)
        trial = from_dict(copy.deepcopy(doc))
        trial.units[uid].params[param] = x
        y = None
        converged = False
        try:
            rep = trial.solve(backend=a.get("backend", "sequential"))
            converged = bool(rep.converged)
            if converged:
                y = _metric(trial, rep, spec)
        except Exception:
            pass
        if y is None:
            failures += 1
        points.append({"value": x, "metric": y, "converged": converged})

    ys = [p["metric"] for p in points if p["metric"] is not None]
    span = (f"metric ranged {min(ys):.6g} .. {max(ys):.6g}" if ys
            else "no point solved")
    tail = f" ({failures} of {steps} points failed)" if failures else ""
    return {"summary": f"Swept {uid}.{param} from {lo:g} to {hi:g} "
                       f"in {steps} steps: {span}{tail}.",
            "unit_id": uid, "param": param, "metric_spec": spec, "points": points}


def _property_table(s: AgentSession, a: dict) -> dict:
    from ..analysis.property_table import DEFAULT_PROPS, PROPERTIES, property_table
    from ..thermo import make_package

    fs = s.require_flowsheet()
    if a.get("z"):
        z = a["z"]
    elif a.get("stream"):
        st = fs.streams.get(a["stream"])
        if st is None or not st.z:
            raise ValueError(f"stream {a['stream']!r} not found or has no "
                             "composition (provide z, or solve first)")
        z = st.normalized_z()
    else:
        raise ValueError("provide either a stream id or an explicit composition z")
    unknown = sorted(set(z) - set(fs.component_ids))
    if unknown:
        raise ValueError(f"components {unknown} are not in the flowsheet "
                         f"(components: {fs.component_ids})")

    T, P = a["T"], a["P"]
    n_pts = (len(T) if isinstance(T, list) else 1) * (len(P) if isinstance(P, list) else 1)
    if n_pts > 200:
        raise ValueError(f"grid too large ({n_pts} points); keep len(T) x len(P) <= 200")
    props = a.get("props") or list(DEFAULT_PROPS)
    bad = [p for p in props if p not in PROPERTIES]
    if bad:
        raise ValueError(f"unknown props {bad}; available: {sorted(PROPERTIES)}")

    pp = make_package(fs.property_package, fs.component_ids)
    grid = property_table(pp, z, T=T, P=P, props=props)

    def _clean(row) -> list:         # NaN -> None so the model sees nulls
        return [None if v != v else float(v) for v in row.tolist()]

    fail_txt = f" {len(grid['failures'])} points failed." if grid["failures"] else ""
    return {"summary": f"{', '.join(props)} over {len(grid['T'])} T x {len(grid['P'])} P "
                       f"points with {fs.property_package}; values[prop][i_T][i_P]."
                       f"{fail_txt}",
            "T": grid["T"].tolist(), "P": grid["P"].tolist(), "z": z, "props": props,
            "values": {name: [_clean(row) for row in grid[name]] for name in props},
            "failures": [{"T": t, "P": p, "error": m} for t, p, m in grid["failures"]]}


def _describe_flowsheet(s: AgentSession, _a: dict) -> dict:
    from .diagnostics import describe_flowsheet
    return describe_flowsheet(s.require_flowsheet(), s.report)


def _explain_convergence(s: AgentSession, _a: dict) -> dict:
    from .diagnostics import explain_convergence
    return explain_convergence(s.require_flowsheet(), s.report)


# -- registry ---------------------------------------------------------------
TOOLS: list[Tool] = [
    Tool("list_unit_types",
         "List the available unit-operation types and their ports. Pass {type} "
         "to get one type's full parameter documentation instead.",
         {"type": "object", "properties": {
             "type": {"type": "string",
                      "description": "a unit type to document in full (params, ports)"}}},
         _list_unit_types),
    Tool("list_property_packages", "List the available thermodynamic property packages.",
         _OBJ, _list_property_packages),
    Tool("new_flowsheet",
         "Start a fresh flowsheet with a component list and a property package "
         "(e.g. 'thermo:PR' for non-polar, 'thermo:NRTL' for polar/azeotropic).",
         {"type": "object", "required": ["components"], "properties": {
             "components": {"type": "array", "items": {"type": "string"},
                            "description": "component ids, e.g. ['nitrogen','hydrogen','ammonia']"},
             "property_package": {"type": "string", "default": "thermo:PR"}}},
         _new_flowsheet),
    Tool("add_unit",
         "Add a unit operation (re-using an existing id replaces that unit). "
         "params depend on the type: Heater {T_out,dP}; Mixer {dP}; Splitter "
         "{split}; Flash {T,P}; EquilibriumReactor {reaction:{stoich,key},T}; "
         "ConversionReactor {reaction,conversion,T_out}; Pump/Compressor "
         "{P_out,eta}; Valve {P_out}. For any other type, call "
         "list_unit_types {type} for its full parameter documentation.",
         {"type": "object", "required": ["id", "type"], "properties": {
             "id": {"type": "string"}, "type": {"type": "string"},
             "params": {"type": "object"}}},
         _add_unit),
    Tool("add_feed",
         "Add a boundary feed stream into a unit port (UNIT:port) with full spec.",
         {"type": "object", "required": ["id", "to", "T", "P", "molar_flow", "z"],
          "properties": {
              "id": {"type": "string"}, "to": {"type": "string"},
              "T": {"type": "number", "description": "K"},
              "P": {"type": "number", "description": "Pa"},
              "molar_flow": {"type": "number", "description": "mol/s"},
              "z": {"type": "object", "description": "{component: mole fraction}"}}},
         _add_feed),
    Tool("connect",
         "Connect a stream from one unit port to another. Use null for 'to' to "
         "make a product/boundary outlet. Ports are 'UNIT:port' strings.",
         {"type": "object", "required": ["id"], "properties": {
             "id": {"type": "string"},
             "from": {"type": ["string", "null"]}, "to": {"type": ["string", "null"]}}},
         _connect),
    Tool("solve",
         "Solve the flowsheet. backend 'sequential' (default) or 'equation_oriented'.",
         {"type": "object", "properties": {
             "backend": {"type": "string", "enum": ["sequential", "equation_oriented"]},
             "tol": {"type": "number"}, "max_iter": {"type": "integer"}}},
         _solve),
    Tool("cost",
         "Run a techno-economic analysis (capital, opex, LCOP, NPV) of the solved "
         "flowsheet. product_component is the species sold (for LCOP).",
         {"type": "object", "properties": {
             "product_component": {"type": "string"},
             "operating_hours": {"type": "number"}, "discount_rate": {"type": "number"},
             "prices_per_kg": {"type": "object"}}},
         _cost),
    Tool("optimize",
         "Minimize/maximize a metric over unit parameters subject to constraints. "
         "metric: {type:'duty'|'flow'|'component_rate', stream, component?}.",
         {"type": "object", "required": ["objective", "design_vars"], "properties": {
             "objective": {"type": "object"}, "design_vars": {"type": "array"},
             "constraints": {"type": "array"},
             "backend": {"type": "string"}}},
         _optimize),
    Tool("stream_table", "Return the current solved stream table.", _OBJ, _stream_table),
    Tool("export_flow", "Export the current flowsheet as a .flow JSON document.",
         _OBJ, _export_flow),
    Tool("describe_flowsheet",
         "Narrate the flowsheet: units, connectivity, feeds/products, logical ops, "
         "and (if solved) the key product flows. Use before explaining a flowsheet.",
         _OBJ, _describe_flowsheet),
    Tool("explain_convergence",
         "Diagnose the last solve: convergence status, residual trend over the "
         "tear iterations, and concrete advice (tear guesses, backend choice, "
         "spec problems). Requires a prior solve.",
         _OBJ, _explain_convergence),
    Tool("remove_unit",
         "Remove a unit and every stream attached to it.",
         {"type": "object", "required": ["id"], "properties": {
             "id": {"type": "string"}}},
         _remove_unit),
    Tool("remove_stream",
         "Remove one stream (a feed, an internal connection, or a product).",
         {"type": "object", "required": ["id"], "properties": {
             "id": {"type": "string"}}},
         _remove_stream),
    Tool("set_param",
         "Set one parameter on an existing unit — the way to edit a single "
         "setting without re-adding the unit. value null clears the param back "
         "to its default.",
         {"type": "object", "required": ["unit_id", "param", "value"], "properties": {
             "unit_id": {"type": "string"}, "param": {"type": "string"},
             "value": {"description": "new value (number/string/bool/object); "
                                      "null clears the param"}}},
         _set_param),
    Tool("pinch_analysis",
         "Heat-integration pinch targeting on the solved flowsheet: minimum "
         "hot/cold utility, pinch temperature, and the heat-recovery potential "
         "vs the current utility use. Solves first if needed.",
         {"type": "object", "properties": {
             "dt_min": {"type": "number",
                        "description": "minimum approach dT in K (default 10)"}}},
         _pinch_analysis),
    Tool("sweep_parameter",
         "Sensitivity study: solve at `steps` evenly spaced values of one unit "
         "parameter and report a metric at each point. Runs on a copy — the "
         "working flowsheet is untouched. metric: {type:'duty'|'flow'|"
         "'component_rate', stream, component?}.",
         {"type": "object", "required": ["unit_id", "param", "start", "stop", "metric"],
          "properties": {
              "unit_id": {"type": "string"}, "param": {"type": "string"},
              "start": {"type": "number"}, "stop": {"type": "number"},
              "steps": {"type": "integer", "minimum": 2, "maximum": 50, "default": 8},
              "metric": {"type": "object"},
              "backend": {"type": "string", "enum": ["sequential", "equation_oriented"]}}},
         _sweep_parameter),
    Tool("property_table",
         "Thermophysical properties of a composition over a T x P grid (K, Pa) "
         "using the flowsheet's property package. Give a stream id (solved, or a "
         "feed) OR an explicit z. props: mass_density, molar_volume, enthalpy, "
         "entropy, vapor_fraction.",
         {"type": "object", "required": ["T", "P"], "properties": {
             "stream": {"type": "string"}, "z": {"type": "object"},
             "T": {"type": ["number", "array"], "items": {"type": "number"},
                   "description": "K; scalar or list"},
             "P": {"type": ["number", "array"], "items": {"type": "number"},
                   "description": "Pa; scalar or list"},
             "props": {"type": "array", "items": {"type": "string"}}}},
         _property_table),
]

_BY_NAME = {t.name: t for t in TOOLS}


def dispatch(session: AgentSession, name: str, args: dict[str, Any]) -> dict:
    """Run one tool. Errors are returned as ``{ok: False, error: ...}`` so an
    agent can read and recover from them rather than crashing the loop."""
    tool = _BY_NAME.get(name)
    if tool is None:
        return {"ok": False, "error": f"unknown tool {name!r}"}
    try:
        out = tool.handler(session, args or {})
        out.setdefault("ok", True)
        return out
    except Exception as exc:  # surface engine/validation errors to the model
        return {"ok": False, "error": f"{type(exc).__name__}: {exc}"}


def anthropic_tools() -> list[dict]:
    """Tool schemas in the Anthropic Messages API shape."""
    return [{"name": t.name, "description": t.description, "input_schema": t.input_schema}
            for t in TOOLS]
