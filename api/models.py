"""Pydantic request/response models for the Caldyr API.

Requests carry a ``.flow`` document (the same schema the engine's ``io`` layer
round-trips) so flowsheet-as-canvas and flowsheet-as-code stay equivalent.
"""
from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class FlowUnitModel(BaseModel):
    """One unit-op node in a .flow document."""
    model_config = ConfigDict(extra="allow")
    id: str
    type: str
    params: dict[str, Any] = Field(default_factory=dict)
    xy: list[float] | None = None


class FlowStreamModel(BaseModel):
    """One stream edge; `from`/`to` are 'UNIT:port' or null (boundary)."""
    model_config = ConfigDict(extra="allow", populate_by_name=True)
    id: str
    from_: str | None = Field(default=None, alias="from")
    to: str | None = None
    spec: dict[str, Any] | None = None


class FlowDocModel(BaseModel):
    """Structural validation of a .flow document at the API boundary.

    The engine's ``io.from_dict`` stays the semantic authority (unit types,
    ports, components); this model catches malformed documents with precise
    422s before they reach it. Unknown keys are allowed so the schema can
    evolve without breaking older clients.
    """
    model_config = ConfigDict(extra="allow")
    schema_id: str = Field(alias="schema", pattern=r"^caldyr\.flow/")
    components: list[dict[str, Any]] = Field(default_factory=list)
    property_package: str = "thermo:PR"
    units: list[FlowUnitModel] = Field(default_factory=list)
    streams: list[FlowStreamModel] = Field(default_factory=list)
    logical: list[dict[str, Any]] = Field(default_factory=list)
    solver_hints: dict[str, Any] = Field(default_factory=dict)


# Endpoints keep receiving the raw dict (the engine round-trips it exactly);
# FlowDocModel.model_validate(flow) is the structural gate.
FlowDoc = dict[str, Any]


class SolveRequest(BaseModel):
    flow: FlowDoc
    backend: Literal["sequential", "equation_oriented"] = "sequential"
    tol: float = 1e-6
    max_iter: int = Field(default=200, ge=1, le=10_000)
    method: Literal["wegstein", "direct"] = "wegstein"


class StreamState(BaseModel):
    id: str
    T: float | None = None
    P: float | None = None
    molar_flow: float | None = None
    z: dict[str, float] = Field(default_factory=dict)
    H: float | None = None
    phase: str | None = None
    vapor_fraction: float | None = None


class SolveReportModel(BaseModel):
    converged: bool
    iterations: int
    residual: float | None
    tol: float
    method: str
    order: list[str]
    tear_streams: list[str]
    duties: dict[str, float]
    messages: list[str]
    history: list[float] = Field(default_factory=list)  # residual per iteration


class SolveResponse(BaseModel):
    report: SolveReportModel
    streams: dict[str, StreamState]
    # per-unit design results (column profiles, FUG numbers, fuel duty, ...)
    designs: dict[str, dict[str, Any]] = Field(default_factory=dict)
    # per-component molar mass (kg/mol) so the web can render mass flow / mass
    # fractions from molar data; keyed by the flowsheet's component ids.
    molar_mass: dict[str, float] = Field(default_factory=dict)


class CostConfig(BaseModel):
    year: int = 2023
    operating_hours: float = 8000.0
    discount_rate: float = 0.10
    project_years: int = 20
    product_component: str = "ammonia"
    product_min_fraction: float = 0.5
    prices_per_kg: dict[str, float] | None = None       # raw-material/product $/kg
    utility_prices: dict[str, float] | None = None      # {utility: $/GJ}
    sizing: dict[str, float] | None = None              # SizingOptions field overrides
    factors: dict[str, float] | None = None             # CostFactors field overrides


class CostRequest(BaseModel):
    flow: FlowDoc
    backend: Literal["sequential", "equation_oriented"] = "sequential"
    config: CostConfig = Field(default_factory=CostConfig)
    tornado: bool = True
    monte_carlo: int = Field(default=0, ge=0, le=10_000)  # MC samples; 0 = skip


class EnvelopeRequest(BaseModel):
    """Phase envelope (bubble/dew T vs P) for one stream's composition."""
    flow: FlowDoc
    stream: str
    n: int = Field(default=25, ge=3, le=200)
    p_min: float | None = None  # Pa; default 0.2x stream P
    p_max: float | None = None  # Pa; default 5x stream P


class PropertyTableRequest(BaseModel):
    """Stream properties over a (T, P) grid — the HYSYS Property Table tool
    (analysis.property_table). Composition comes from a named ``stream`` in the
    .flow document, or from an explicit ``z`` (overrides the stream)."""
    flow: FlowDoc
    stream: str | None = None
    z: dict[str, float] | None = None
    T: list[float] = Field(min_length=1, max_length=200)   # K grid values
    P: list[float] = Field(min_length=1, max_length=200)   # Pa grid values
    props: list[str] | None = None              # default analysis.DEFAULT_PROPS


class ReliefRequest(BaseModel):
    """Pressure-relief valve sizing (API 520/526; analysis.relief).

    ``phase='vapor'`` needs W, T, M, Z, k, P1 (+ optional coefficients);
    ``phase='liquid'`` needs W, rho, P1, P2 (+ optional coefficients). All SI;
    P1/P2 absolute Pa. ``M`` (vapor) may be omitted and derived from a stream's
    composition when ``flow``+``stream`` are supplied."""
    phase: Literal["vapor", "liquid"]
    W: float                                    # relieving mass flow, kg/s
    P1: float                                   # upstream relieving pressure, Pa abs
    # vapor
    T: float | None = None
    M: float | None = None                      # kg/mol (or derived from stream)
    Z: float = 1.0
    k: float | None = None                      # Cp/Cv
    backpressure: float = 101_325.0
    # liquid
    rho: float | None = None                    # kg/m^3
    P2: float | None = None                     # downstream pressure, Pa abs
    # coefficients (defaults match the physics functions)
    Kd: float | None = None
    Kb: float = 1.0
    Kc: float = 1.0
    Kw: float = 1.0
    Kv: float = 1.0
    # optional composition source for M
    flow: FlowDoc | None = None
    stream: str | None = None


class PinchRequest(BaseModel):
    """Heat-integration pinch targeting on a solved flowsheet (analysis.pinch)."""
    flow: FlowDoc
    backend: Literal["sequential", "equation_oriented"] = "sequential"
    dt_min: float = Field(default=10.0, gt=0.0)   # minimum approach ΔT, K
    tol: float = 1e-6


class MetricSpec(BaseModel):
    """A scalar read off the solved flowsheet."""
    type: Literal["duty", "flow", "component_rate"]
    stream: str
    component: str | None = None


class ConstraintSpec(BaseModel):
    metric: MetricSpec
    op: Literal[">=", "<="]
    value: float


class DesignVarSpec(BaseModel):
    unit_id: str
    param: str
    lower: float
    upper: float
    initial: float | None = None


class ObjectiveSpec(BaseModel):
    sense: Literal["min", "max"] = "min"
    metric: MetricSpec


class OptimizeRequest(BaseModel):
    flow: FlowDoc
    backend: Literal["sequential", "equation_oriented"] = "sequential"
    objective: ObjectiveSpec
    design_vars: list[DesignVarSpec]
    constraints: list[ConstraintSpec] = Field(default_factory=list)
    tol: float = 1e-7
