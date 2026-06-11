# Data model

All physics is in SI internally (K, Pa, mol/s, J). The model is deliberately small
and explicit so it serializes cleanly to `.flow` JSON and is easy for the AI layer
to manipulate.

## Component
A chemical species and the metadata needed to look up properties.
```python
@dataclass(frozen=True)
class Component:
    id: str            # canonical key, e.g. "water", "methane"
    name: str
    formula: str | None = None
    cas: str | None = None
    # property backends resolve by id/cas; no property data stored here
```

## Stream
A material stream. State is fully defined by (T or H, P, total flow, composition);
the property package resolves the rest.
```python
@dataclass
class Stream:
    id: str
    components: list[str]               # references Component.id
    T: float | None                     # K
    P: float | None                     # Pa
    molar_flow: float | None            # mol/s
    z: dict[str, float]                 # composition, mole fractions (sum→1)
    # resolved by the property package after solve:
    H: float | None = None              # J/mol (molar enthalpy)
    phase: str | None = None            # "vapor" | "liquid" | "VLE"
    vapor_fraction: float | None = None
```
An **energy stream** is a separate lightweight type carrying a duty (W) only.

## Port
A typed connection point on a unit op.
```python
@dataclass(frozen=True)
class Port:
    name: str                 # e.g. "in1", "out", "duty"
    direction: str            # "inlet" | "outlet"
    kind: str = "material"    # "material" | "energy"
```

## UnitOp
Abstract base. Holds parameters and ports; implements the solve contract.
```python
class UnitOp(ABC):
    id: str
    ports: list[Port]
    params: dict                       # model-specific (e.g. {"dP": 0, "Q": 1e5})
    @abstractmethod
    def solve(self, inlets: dict[str, Stream], pp: "PropertyPackage")
        -> dict[str, Stream]: ...
```

## Flowsheet
The directed graph.
```python
class Flowsheet:
    components: list[Component]
    property_package: str              # selected package id
    units: dict[str, UnitOp]           # id → unit op
    connections: list[Connection]      # (from_unit.port) → (to_unit.port) carrying a Stream
    def solve(self, backend="sequential") -> SolveReport: ...
```

## `.flow` file format (JSON)
Plain text, git-friendly, exact round-trip.
```json
{
  "schema": "caldyr.flow/1",
  "meta": {"name": "Mixer + Heater", "created": "..."},
  "components": [{"id": "water"}, {"id": "ethanol"}],
  "property_package": "thermo:PR",
  "units": [
    {"id": "MIX1", "type": "Mixer", "params": {"dP": 0}, "xy": [120, 200]},
    {"id": "H1",   "type": "Heater", "params": {"dP": 0, "T_out": 350}, "xy": [320, 200]}
  ],
  "streams": [
    {"id": "S1", "from": null, "to": "MIX1:in1",
     "spec": {"T": 298.15, "P": 101325, "molar_flow": 10,
              "z": {"water": 0.6, "ethanol": 0.4}}},
    {"id": "S2", "from": null, "to": "MIX1:in2",
     "spec": {"T": 320, "P": 101325, "molar_flow": 5, "z": {"water": 1.0}}},
    {"id": "S3", "from": "MIX1:out", "to": "H1:in1"},
    {"id": "S4", "from": "H1:out",  "to": null}
  ],
  "solved": { "S3": {"T": 305.1, "H": -2.81e5, "phase": "liquid"}, "...": {} }
}
```
Rules: feed streams have `from: null` and a full `spec`; internal streams carry no
spec (computed); `solved` is a cache and may be regenerated. Versioned via `schema`.
