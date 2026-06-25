"""Operating cost: variable (raw materials, utilities) + fixed (labor, maintenance).

Variable costs come straight from the solved flowsheet — boundary feed flows
priced per kg, and the utility duties found during sizing priced per GJ. Fixed
costs and the manufacturing-cost roll-up follow Turton 4e Ch. 8:

    COM_d = 0.180·FCI + 2.73·C_OL + 1.23·(C_UT + C_WT + C_RM)

with operating labor from the Turton operator correlation
N_OL = (6.29 + 0.23·N_np)^0.5 per shift (no solids handling here).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from . import data
from .capital import CapitalEstimate
from .sizing import EquipmentSize

SECONDS_PER_HOUR = 3600.0
# Default COM/labor factors live in data.CostFactors (overridable per flowsheet).
# Unit types that count as process steps in the labor correlation (Turton: no
# pumps or vessels).
_LABOR_UNITS = {"Heater", "FiredHeater", "AirCooler", "HeatExchanger", "Compressor",
                "ConversionReactor", "EquilibriumReactor", "GibbsReactor",
                "CSTR", "PFR", "ShortcutColumn"}


@dataclass
class OperatingCosts:
    raw_materials: float           # $/yr
    utilities: float               # $/yr
    fixed: float                   # $/yr (labor, maintenance, overheads)
    total: float                   # COM_d, $/yr
    operating_hours: float
    breakdown: dict = field(default_factory=dict)


def _feed_streams(fs):
    """Boundary feed streams (no upstream unit) with a resolved composition."""
    feeds = []
    for conn in fs.connections:
        if conn.from_unit is None and conn.to_unit is not None:
            s = fs.streams.get(conn.stream_id)
            if s is not None and s.molar_flow:
                feeds.append(s)
    return feeds


def estimate_opex(fs, sizes: list[EquipmentSize], capital: CapitalEstimate, *,
                  operating_hours: float = 8000.0,
                  prices_per_kg: dict | None = None,
                  utility_prices: dict | None = None,
                  factors: data.CostFactors | None = None) -> OperatingCosts:
    prices = {**data.PRICES_PER_KG, **(prices_per_kg or {})}
    f = factors or data.CostFactors()
    seconds = operating_hours * SECONDS_PER_HOUR

    # Raw materials: each boundary feed, by component mass flow x price.
    rm = 0.0
    rm_detail: dict[str, float] = {}
    for s in _feed_streams(fs):
        z = s.normalized_z()
        for comp, frac in z.items():
            price = prices.get(comp)
            if not price:
                continue
            kg_s = s.molar_flow * frac * data.molar_mass(comp)
            cost = kg_s * price * seconds
            rm += cost
            rm_detail[comp] = rm_detail.get(comp, 0.0) + cost

    # Utilities: each sized item's duty x price per GJ.
    ut = 0.0
    ut_detail: dict[str, float] = {}
    up = utility_prices or {}
    for size in sizes:
        if not size.utility:
            continue
        price_per_gj = up.get(size.utility, data.UTILITIES[size.utility].price_per_GJ)
        cost = size.utility_duty_W / 1e9 * price_per_gj * seconds
        ut += cost
        ut_detail[size.utility] = ut_detail.get(size.utility, 0.0) + cost

    # Fixed costs (Turton COM_d decomposition; factors overridable via CostFactors).
    n_np = sum(1 for u in fs.units.values() if type(u).__name__ in _LABOR_UNITS)
    n_ol = math.sqrt(f.labor_a + f.labor_b * n_np)          # operators per shift
    labor = n_ol * f.shifts_per_operator * f.operator_salary
    fci = capital.grassroots
    com_d = f.com_fci_factor * fci + f.com_labor_factor * labor \
        + f.com_variable_factor * (ut + rm)
    fixed = com_d - (ut + rm)

    return OperatingCosts(
        raw_materials=rm,
        utilities=ut,
        fixed=fixed,
        total=com_d,
        operating_hours=operating_hours,
        breakdown={
            "raw_materials": rm_detail,
            "utilities": ut_detail,
            "operating_labor": labor,
            "operators_per_shift": n_ol,
            "maintenance_and_overhead": fixed - labor,
        },
    )
