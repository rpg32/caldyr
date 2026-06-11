"""Capital cost roll-up: bare-module items -> ISBL -> grassroots -> TCI.

Follows the Turton method (4e Ch. 7):
  * ISBL (inside battery limits) = Σ bare-module costs.
  * Total module cost  C_TM = 1.18 · ISBL   (15% contingency + 3% fee).
  * Grassroots cost    C_GR = C_TM + 0.50 · Σ Cbm(base)   (the 0.50·base term is
    the auxiliary/offsite — OSBL — facilities).
  * Total capital investment TCI = C_GR + working capital.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .costing import CostResult

CONTINGENCY_AND_FEE = 1.18      # Turton 4e Eq. 7.5 (0.15 contingency + 0.03 fee)
OFFSITE_FACTOR = 0.50           # Turton 4e grassroots: 0.50 · Σ Cbm(base)
WORKING_CAPITAL_FRACTION = 0.15  # of grassroots fixed capital (Turton ~15-20%)


@dataclass
class CapitalEstimate:
    isbl: float                # Σ bare-module (installed, this material/pressure)
    total_module: float        # 1.18 · ISBL
    osbl: float                # 0.50 · Σ bare-module(base)
    grassroots: float          # fixed capital investment
    working_capital: float
    tci: float                 # total capital investment
    year: int
    items: list[CostResult] = field(default_factory=list)


def estimate_capital(items: list[CostResult], year: int = 2023) -> CapitalEstimate:
    isbl = sum(c.bare_module for c in items)
    cbm_base = sum(c.bare_module_base for c in items)
    total_module = CONTINGENCY_AND_FEE * isbl
    osbl = OFFSITE_FACTOR * cbm_base
    grassroots = total_module + osbl
    working_capital = WORKING_CAPITAL_FRACTION * grassroots
    return CapitalEstimate(
        isbl=isbl,
        total_module=total_module,
        osbl=osbl,
        grassroots=grassroots,
        working_capital=working_capital,
        tci=grassroots + working_capital,
        year=year,
        items=list(items),
    )
