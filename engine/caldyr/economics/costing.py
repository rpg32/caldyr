"""Equipment costing per the Turton bare-module method. See docs/ECONOMICS.md.

Pipeline for one equipment item:
  capacity attribute A  ->  purchased cost Cp0 (CEPCI 397, CS, ambient)
                        ->  pressure factor Fp, material factor Fm
                        ->  bare-module cost Cbm = Cp0 (B1 + B2 Fm Fp)
                        ->  escalate to the analysis year via CEPCI.

All correlation constants and their sources live in :mod:`caldyr.economics.data`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, log10

from . import data


# -- low-level correlation functions (pure; reusable) -----------------------
def purchased_cost(K1: float, K2: float, K3: float, A: float) -> float:
    """log10(Cp0) = K1 + K2 log10 A + K3 (log10 A)^2 ; A = capacity attribute."""
    la = log10(A)
    return 10 ** (K1 + K2 * la + K3 * la * la)


def bare_module_cost(Cp0: float, B1: float, B2: float, Fm: float, Fp: float) -> float:
    """Cbm = Cp0 (B1 + B2 Fm Fp)."""
    return Cp0 * (B1 + B2 * Fm * Fp)


def escalate_cepci(cost: float, cepci_now: float, cepci_base: float) -> float:
    return cost * (cepci_now / cepci_base)


def pressure_factor(eq_type: str, pressure_barg: float) -> float:
    """log10(Fp) = C1 + C2 log10 P + C3 (log10 P)^2, with Fp = 1 below the
    correlation's valid minimum pressure (Turton convention)."""
    c = data.PRESSURE[eq_type]
    if pressure_barg < c.pmin:
        return 1.0
    lp = log10(pressure_barg)
    return 10 ** (c.C1 + c.C2 * lp + c.C3 * lp * lp)


def vessel_pressure_factor(pressure_barg: float, diameter_m: float) -> float:
    """Pressure factor for process vessels from the membrane-stress wall-thickness
    relation (Turton 4e Eq. for vessels): based on a 940-bar allowable stress,
    0.9 weld efficiency, 3.15 mm corrosion allowance, floored at 1.0."""
    p = pressure_barg
    if p < -0.5:
        return 1.25
    fp = ((p + 1.0) * diameter_m / (2.0 * (944.0 * 0.9 - 0.6 * (p + 1.0))) + 0.00315) / 0.0063
    return max(1.0, fp)


def _clamped_fp(eq: str, pressure_barg: float, notes: list[str]) -> float:
    """Pressure factor with the correlation's max-P clamp (noted, never silent)."""
    pmax = data.PRESSURE[eq].pmax
    if pressure_barg > pmax:
        notes.append(f"pressure {pressure_barg:.0f} barg > correlation max {pmax:.0f}; "
                     f"Fp clamped at the max (refine for high-P service)")
        pressure_barg = pmax
    return pressure_factor(eq, pressure_barg)


def material_factor(eq_type: str, material: str) -> float:
    entry = data.MATERIAL.get(eq_type)
    if entry is None or material not in entry.factors:
        raise ValueError(f"no material factor for {material!r} on {eq_type!r}")
    return entry.factors[material]


def quantity_factor(eq_type: str, n: int) -> float:
    """Fq for N stacked identical items (column trays): Turton 4e Eq. A.5,
    log10(Fq) = a + b log10 N + c (log10 N)^2 for N below the full-discount
    count; Fq = 1 otherwise."""
    c = data.QUANTITY_FACTOR.get(eq_type)
    if c is None or n >= c.n_full:
        return 1.0
    ln = log10(n)
    return 10 ** (c.a + c.b * ln + c.c * ln * ln)


# -- high-level: cost a sized equipment item --------------------------------
@dataclass
class CostResult:
    """Costed equipment item, all dollar figures escalated to ``year``."""
    unit_id: str
    equipment_type: str
    attribute: float
    attribute_name: str
    n_units: int                 # parallel trains used to stay in correlation range
    purchased: float             # Cp0 (CS, ambient), escalated, all trains
    bare_module: float           # Cbm, escalated, all trains
    bare_module_base: float      # Cbm at Fm=Fp=1 (for grassroots/offsites)
    year: int
    factors: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def cost_equipment(size, year: int = 2023) -> CostResult:
    """Cost an :class:`~caldyr.economics.sizing.EquipmentSize`.

    Splits into parallel trains if the attribute exceeds the correlation's range,
    and escalates from the CEPCI=397 (2001) basis to ``year``.
    """
    eq = size.equipment_type
    notes: list[str] = []
    cepci = escalate_cepci(1.0, data.CEPCI[year], data.CEPCI[data.CEPCI_BASE_YEAR])

    pdata = data.PURCHASED[eq]
    a_total = size.attribute
    n_units = 1
    if a_total > pdata.amax:
        n_units = ceil(a_total / pdata.amax)
        notes.append(f"attribute {a_total:.3g} > {pdata.amax:.3g}; "
                     f"split into {n_units} parallel units")
    a_each = a_total / n_units
    if a_each < pdata.amin:
        notes.append(f"attribute {a_each:.3g} < {pdata.amin:.3g}; "
                     f"clamped to minimum (small-equipment floor)")
        a_each = pdata.amin

    cp0 = purchased_cost(pdata.K1, pdata.K2, pdata.K3, a_each)

    # Pressure and material factors. Direct-Fbm equipment skips the material
    # factor; if it has a pressure correlation (fired heaters: Turton 4e
    # Eq. A.6, Cbm = Cp0 Fbm Fp), the pressure factor still applies.
    if eq in data.FBM_DIRECT:
        fbm = data.FBM_DIRECT[eq].Fbm
        fm = float("nan")
        if eq in data.PRESSURE:
            fp = _clamped_fp(eq, size.pressure_barg, notes)
            cbm_each = cp0 * fbm * fp
        else:
            fp = float("nan")
            cbm_each = cp0 * fbm
        cbm_base_each = cp0 * fbm
    else:
        if eq in data.PRESSURE:
            fp = _clamped_fp(eq, size.pressure_barg, notes)
        elif eq.startswith("vessel"):
            fp = vessel_pressure_factor(size.pressure_barg, size.diameter_m)
        else:
            fp = 1.0
        fm = material_factor(eq, size.material)
        b = data.BARE_MODULE[eq]
        cbm_each = bare_module_cost(cp0, b.B1, b.B2, fm, fp)
        cbm_base_each = bare_module_cost(cp0, b.B1, b.B2, 1.0, 1.0)

    # Multiple identical stacked items (trays): Cbm = N * Cp0 * Fbm * Fq.
    qty = max(1, int(size.quantity))
    fq = quantity_factor(eq, qty)
    if qty > 1:
        notes.append(f"{qty} items; quantity factor Fq={fq:.2f}")

    return CostResult(
        unit_id=size.unit_id,
        equipment_type=eq,
        attribute=a_total,
        attribute_name=size.attribute_name,
        n_units=n_units,
        purchased=cp0 * n_units * qty * cepci,
        bare_module=cbm_each * n_units * qty * fq * cepci,
        bare_module_base=cbm_base_each * n_units * qty * fq * cepci,
        year=year,
        factors={"Fp": fp, "Fm": fm, "Fbm": cbm_each / cp0, "Fq": fq},
        warnings=notes,
    )


def six_tenths(cost_known: float, size_known: float, size_new: float, n: float = 0.6) -> float:
    """Capacity scaling: Cost_new = Cost_known * (size_new / size_known)^n."""
    return cost_known * (size_new / size_known) ** n
