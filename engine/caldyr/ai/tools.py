"""Typed tools over the engine — the AI-native interface.

Each tool has a JSON-schema input (so it drops straight into Anthropic tool-use
or MCP) and a handler that mutates / reads an :class:`AgentSession`. Handlers
return a dict with a human-readable ``summary`` plus structured data, so the
model gets both a sentence to reason over and machine-usable fields.

`dispatch(session, name, args)` runs one tool; `anthropic_tools()` exports the
schemas in the shape the Anthropic SDK expects.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from ..core import Component
from ..core import Flowsheet as _Flowsheet
from ..economics import TEAConfig, analyze
from ..io import to_dict
from ..solver import DesignVar, optimize
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
def _list_unit_types(_s: AgentSession, _a: dict) -> dict:
    items = []
    for name, cls in sorted(REGISTRY.items()):
        unit = cls(f"_{name}")
        items.append({"type": name, "doc": (cls.__doc__ or "").strip().split("\n")[0],
                      "ports": [{"name": p.name, "direction": p.direction, "kind": p.kind}
                                for p in unit.ports]})
    return {"summary": f"{len(items)} unit types available.", "unit_types": items}


def _list_property_packages(_s: AgentSession, _a: dict) -> dict:
    pkgs = [
        {"id": "thermo:PR", "note": "Peng-Robinson cubic EOS — non-polar systems"},
        {"id": "thermo:SRK", "note": "Soave-Redlich-Kwong cubic EOS"},
        {"id": "thermo:NRTL", "note": "NRTL activity model — polar mixtures, azeotropes"},
    ]
    return {"summary": "3 property packages.", "property_packages": pkgs}


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
    unit = fs.add(REGISTRY[utype](a["id"], a.get("params", {})))
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


def _optimize(s: AgentSession, a: dict) -> dict:
    fs = s.require_flowsheet()
    design = [DesignVar(d["unit_id"], d["param"], d["lower"], d["upper"], d.get("initial"))
              for d in a["design_vars"]]

    def metric(fs_, rep, spec):
        if spec["type"] == "duty":
            return rep.duties[spec["stream"]]
        st = fs_.streams[spec["stream"]]
        if spec["type"] == "flow":
            return st.molar_flow
        return st.molar_flow * st.normalized_z().get(spec["component"], 0.0)

    obj_spec = a["objective"]
    sign = 1.0 if obj_spec.get("sense", "min") == "min" else -1.0
    cons = []
    for c in a.get("constraints", []):
        if c["op"] == ">=":
            cons.append(lambda fs_, rep, c=c: metric(fs_, rep, c["metric"]) - c["value"])
        else:
            cons.append(lambda fs_, rep, c=c: c["value"] - metric(fs_, rep, c["metric"]))

    res = optimize(fs, lambda fs_, rep: sign * metric(fs_, rep, obj_spec["metric"]),
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


def _describe_flowsheet(s: AgentSession, _a: dict) -> dict:
    from .diagnostics import describe_flowsheet
    return describe_flowsheet(s.require_flowsheet(), s.report)


def _explain_convergence(s: AgentSession, _a: dict) -> dict:
    from .diagnostics import explain_convergence
    return explain_convergence(s.require_flowsheet(), s.report)


# -- registry ---------------------------------------------------------------
TOOLS: list[Tool] = [
    Tool("list_unit_types", "List the available unit-operation types and their ports.",
         _OBJ, _list_unit_types),
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
         "Add a unit operation. params depend on the type: Heater {T_out,dP}; "
         "Mixer {dP}; Splitter {split}; Flash {T,P}; EquilibriumReactor "
         "{reaction:{stoich,key},T}; ConversionReactor {reaction,conversion,T_out}; "
         "Pump/Compressor {P_out,eta}; Valve {P_out}.",
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
