"""Discounted-cash-flow profitability: NPV, IRR, payback, and levelized cost.

The headline number for the analysis track is the **levelized cost of product**
(LCOP) — the price at which discounted revenue exactly covers annualized capital
plus operating cost:

    LCOP = (CRF · TCI + annual OPEX) / annual production
    CRF  = r(1+r)^n / ((1+r)^n − 1)     (capital recovery factor)
"""
from __future__ import annotations

from dataclasses import dataclass

from scipy.optimize import brentq


def capital_recovery_factor(rate: float, years: int) -> float:
    if rate == 0:
        return 1.0 / years
    f = (1.0 + rate) ** years
    return rate * f / (f - 1.0)


@dataclass
class Profitability:
    tci: float
    annual_revenue: float
    annual_opex: float
    net_annual_cashflow: float
    annual_production: float          # production units (e.g. kg) per year
    discount_rate: float
    project_years: int
    npv: float
    irr: float | None                 # None if no positive cash flow
    payback_years: float | None
    lcop: float                       # $ per production unit


def _npv(rate: float, tci: float, ncf: float, years: int) -> float:
    return -tci + sum(ncf / (1.0 + rate) ** t for t in range(1, years + 1))


def profitability(tci: float, annual_opex: float, annual_revenue: float,
                  annual_production: float, *, discount_rate: float = 0.10,
                  project_years: int = 20) -> Profitability:
    ncf = annual_revenue - annual_opex
    npv = _npv(discount_rate, tci, ncf, project_years)

    irr: float | None = None
    if ncf > 0:
        # NPV decreases in r; bracket a root between -0.99 and a large rate.
        lo, hi = -0.99, 10.0
        if _npv(lo, tci, ncf, project_years) > 0 > _npv(hi, tci, ncf, project_years):
            irr = float(brentq(_npv, lo, hi, args=(tci, ncf, project_years), xtol=1e-8))

    payback = tci / ncf if ncf > 0 else None
    crf = capital_recovery_factor(discount_rate, project_years)
    lcop = (crf * tci + annual_opex) / annual_production

    return Profitability(
        tci=tci, annual_revenue=annual_revenue, annual_opex=annual_opex,
        net_annual_cashflow=ncf, annual_production=annual_production,
        discount_rate=discount_rate, project_years=project_years,
        npv=npv, irr=irr, payback_years=payback, lcop=lcop,
    )
