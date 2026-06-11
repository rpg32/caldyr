"""`.flow` JSON load/save. Plain text, git-friendly, exact round-trip.
Schema documented in docs/DATA_MODEL.md (schema id: "caldyr.flow/1").
"""
import json
from pathlib import Path
from typing import Any

from ..core import Component, Flowsheet, Stream
from ..unitops import REGISTRY  # importing populates the type registry

SCHEMA = "caldyr.flow/1"


def _endpoint(unit: str | None, port: str | None) -> str | None:
    return f"{unit}:{port}" if unit is not None else None


def _type_name(unit) -> str:
    for name, cls in REGISTRY.items():
        if type(unit) is cls:
            return name
    raise ValueError(f"unit {unit.id!r} of type {type(unit).__name__} is not registered")


def to_dict(flowsheet: Flowsheet, meta: dict | None = None) -> dict:
    """Serialize a Flowsheet to a schema-conformant dict."""
    components: list[dict[str, Any]] = []
    for c in flowsheet.components:
        entry: dict[str, Any] = {"id": c.id}
        if c.name:
            entry["name"] = c.name
        if c.formula is not None:
            entry["formula"] = c.formula
        if c.cas is not None:
            entry["cas"] = c.cas
        pseudo = getattr(c, "pseudo", None)
        if pseudo:                       # assay pseudo-component constants
            entry["pseudo"] = dict(pseudo)
        components.append(entry)

    units: list[dict[str, Any]] = []
    for unit in flowsheet.units.values():
        entry = {"id": unit.id, "type": _type_name(unit), "params": unit.params}
        xy = getattr(unit, "xy", None)
        if xy is not None:
            entry["xy"] = list(xy)
        units.append(entry)

    streams: list[dict[str, Any]] = []
    solved: dict[str, dict[str, Any]] = {}
    for conn in flowsheet.connections:
        sid = conn.stream_id
        s = flowsheet.streams.get(sid)
        entry = {
            "id": sid,
            "from": _endpoint(conn.from_unit, conn.from_port),
            "to": _endpoint(conn.to_unit, conn.to_port),
        }
        if conn.from_unit is None and s is not None:  # a feed carries its spec
            entry["spec"] = _spec(s)
        streams.append(entry)
        if s is not None and s.H is not None:  # cache resolved state
            solved[sid] = {
                "T": s.T, "P": s.P, "molar_flow": s.molar_flow, "z": s.z,
                "H": s.H, "phase": s.phase, "vapor_fraction": s.vapor_fraction,
            }

    doc = {
        "schema": SCHEMA,
        "meta": meta or {},
        "components": components,
        "property_package": flowsheet.property_package,
        "units": units,
        "streams": streams,
    }
    if flowsheet.logical:           # Set/Adjust logical ops (solver.logical)
        doc["logical"] = flowsheet.logical
    if flowsheet.solver_hints:      # tear guesses / tolerance override
        doc["solver_hints"] = flowsheet.solver_hints
    if solved:
        doc["solved"] = solved
    return doc


def _spec(s: Stream) -> dict:
    return {"T": s.T, "P": s.P, "molar_flow": s.molar_flow, "z": s.z}


def from_dict(data: dict) -> Flowsheet:
    """Build a Flowsheet from a parsed `.flow` dict."""
    schema = data.get("schema", "")
    if not schema.startswith("caldyr.flow/"):
        raise ValueError(f"unrecognized schema {schema!r} (expected caldyr.flow/*)")

    fs = Flowsheet(
        components=[Component(**c) for c in data.get("components", [])],
        property_package=data.get("property_package", "thermo:PR"),
        logical=list(data.get("logical", [])),
        solver_hints=dict(data.get("solver_hints", {})),
    )

    for u in data.get("units", []):
        utype = u["type"]
        if utype not in REGISTRY:
            raise ValueError(f"unknown unit type {utype!r}; registered: {sorted(REGISTRY)}")
        unit = REGISTRY[utype](u["id"], u.get("params", {}))
        if "xy" in u:
            unit.xy = tuple(u["xy"])
        fs.add(unit)

    for s in data.get("streams", []):
        fs.connect(s["id"], s.get("from"), s.get("to"))
        spec = s.get("spec")
        if spec is not None:
            fs.streams[s["id"]] = Stream(
                id=s["id"], components=list(fs.component_ids),
                T=spec.get("T"), P=spec.get("P"),
                molar_flow=spec.get("molar_flow"), z=dict(spec.get("z", {})),
            )

    for sid, st in data.get("solved", {}).items():  # restore for exact round-trips
        existing = fs.streams.get(sid)
        components = existing.components if existing else list(fs.component_ids)
        fs.streams[sid] = Stream(
            id=sid, components=components,
            T=st.get("T"), P=st.get("P"), molar_flow=st.get("molar_flow"),
            z=dict(st.get("z", {})), H=st.get("H"),
            phase=st.get("phase"), vapor_fraction=st.get("vapor_fraction"),
        )

    return fs


def load_flow(path: str | Path) -> Flowsheet:
    return from_dict(json.loads(Path(path).read_text()))


def save_flow(flowsheet: Flowsheet, path: str | Path, meta: dict | None = None) -> None:
    Path(path).write_text(json.dumps(to_dict(flowsheet, meta), indent=2))
