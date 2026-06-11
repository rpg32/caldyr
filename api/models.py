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
    max_iter: int = 200
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


class CostConfig(BaseModel):
    year: int = 2023
    operating_hours: float = 8000.0
    discount_rate: float = 0.10
    project_years: int = 20
    product_component: str = "ammonia"
    prices_per_kg: dict[str, float] | None = None


class CostRequest(BaseModel):
    flow: FlowDoc
    backend: Literal["sequential", "equation_oriented"] = "sequential"
    config: CostConfig = Field(default_factory=CostConfig)
    tornado: bool = True
    monte_carlo: int = 0  # number of MC samples; 0 = skip


class EnvelopeRequest(BaseModel):
    """Phase envelope (bubble/dew T vs P) for one stream's composition."""
    flow: FlowDoc
    stream: str
    n: int = Field(default=25, ge=3, le=200)
    p_min: float | None = None  # Pa; default 0.2x stream P
    p_max: float | None = None  # Pa; default 5x stream P


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
