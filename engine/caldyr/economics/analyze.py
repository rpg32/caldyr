"""Top-level techno-economic analysis: solved flowsheet -> full TEA result.

`analyze()` runs the whole ECONOMICS.md pipeline (size -> cost -> capital -> opex
-> profitability). The financial core is factored into `evaluate_economics()` so
the Monte-Carlo and tornado tools can re-evaluate cheaply under perturbed inputs
without re-solving the flowsheet.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from ..thermo import make_package
from . import data
from .capital import CapitalEstimate, estimate_capital
from .costing import CostResult, cost_equipment
from .opex import OperatingCosts, estimate_opex
from .profitability import Profitability, profitability
from .sizing import EquipmentSize, SizingOptions, size_flowsheet

SECONDS_PER_HOUR = 3600.0


@dataclass
class TEAConfig:
    year: int = 2023
    operating_hours: float = 8000.0
    discount_rate: float = 0.10
    project_years: int = 20
    product_component: str = "ammonia"
    product_min_fraction: float = 0.5      # a "product-grade" boundary stream
    prices_per_kg: dict | None = None
    sizing: SizingOptions | None = None


@dataclass
class TEAResult:
    capital: CapitalEstimate
    opex: OperatingCosts
    profitability: Profitability
    annual_production_kg: float
    annual_revenue: float
    sizes: list[EquipmentSize]
    costs: list[CostResult]
    config: TEAConfig = field(default_factory=TEAConfig)


def _product_streams(fs, component: str, min_fraction: float):
    """Boundary product streams (no downstream unit) rich in ``component``."""
    out = []
    for conn in fs.connections:
        if conn.to_unit is None and conn.from_unit is not None:
            s = fs.streams.get(conn.stream_id)
            if s is not None and s.molar_flow and s.z.get(component, 0.0) >= min_fraction:
                out.append(s)
    return out


def annual_production_kg(fs, component: str, min_fraction: float, hours: float) -> float:
    mw = data.molar_mass(component)
    kg_s = sum(s.molar_flow * s.normalized_z().get(component, 0.0) * mw
               for s in _product_streams(fs, component, min_fraction))
    return kg_s * hours * SECONDS_PER_HOUR


def evaluate_economics(fs, sizes, cfg: TEAConfig, *, capex_multiplier: float = 1.0,
                       prices_per_kg: dict | None = None,
                       operating_hours: float | None = None,
                       discount_rate: float | None = None) -> TEAResult:
    """Cost -> capital -> opex -> profitability for fixed equipment sizes. Cheap
    (no flowsheet solve), so Monte-Carlo can call it thousands of times."""
    prices = {**data.PRICES_PER_KG, **(cfg.prices_per_kg or {}), **(prices_per_kg or {})}
    hours = operating_hours if operating_hours is not None else cfg.operating_hours
    rate = discount_rate if discount_rate is not None else cfg.discount_rate

    costs = [cost_equipment(s, cfg.year) for s in sizes]
    if capex_multiplier != 1.0:
        costs = [replace(c, purchased=c.purchased * capex_multiplier,
                         bare_module=c.bare_module * capex_multiplier,
                         bare_module_base=c.bare_module_base * capex_multiplier)
                 for c in costs]
    capital = estimate_capital(costs, cfg.year)
    opex = estimate_opex(fs, sizes, capital, operating_hours=hours, prices_per_kg=prices)

    production = annual_production_kg(fs, cfg.product_component, cfg.product_min_fraction, hours)
    if production <= 0:
        raise ValueError(
            f"no product stream rich in {cfg.product_component!r} found (need a "
            f"boundary outlet with mole fraction >= {cfg.product_min_fraction}). "
            f"Connect the product outlet to a boundary (to = null) so it is reported, "
            f"or check product_component."
        )
    if cfg.product_component not in prices:
        raise ValueError(f"no price for product {cfg.product_component!r}; "
                         f"pass prices_per_kg")
    revenue = production * prices[cfg.product_component]
    prof = profitability(capital.tci, opex.total, revenue, production,
                         discount_rate=rate, project_years=cfg.project_years)

    return TEAResult(
        capital=capital, opex=opex, profitability=prof,
        annual_production_kg=production, annual_revenue=revenue,
        sizes=sizes, costs=costs, config=cfg,
    )


def analyze(fs, report, config: TEAConfig | None = None) -> TEAResult:
    """Full techno-economic analysis of a solved flowsheet."""
    cfg = config or TEAConfig()
    pp = make_package(fs.property_package, fs.component_ids)
    opts = cfg.sizing or SizingOptions()
    sizes = size_flowsheet(fs, report, pp, opts)
    return evaluate_economics(fs, sizes, cfg)
