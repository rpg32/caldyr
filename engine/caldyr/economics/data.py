"""Cited cost-correlation and price data for the economics layer.

Primary source for equipment correlations: **Turton, Bailie, Whiting & Shaeiwitz,
*Analysis, Synthesis, and Design of Chemical Processes*, 4th ed. (Pearson, 2012),
Appendix A** — purchased-cost constants (Table A.1), pressure-factor constants
(Table A.2), material factors (Table A.3) and bare-module factors (Table A.4).
Correlations are on the CEPCI = 397 (2001) basis.

Tables are typed dataclasses so every constant is named and a ``source`` travels
with it — nothing is an unsourced magic number. Prices are *representative*
inputs meant to be swept in sensitivity/Monte-Carlo, not fixed truths.
"""
from __future__ import annotations

from dataclasses import dataclass

# -- CEPCI plant cost index (Chemical Engineering magazine) -----------------
CEPCI = {
    2001: 397.0,   # Turton 4e correlation basis
    2018: 603.1,
    2023: 797.9,
}
CEPCI_SOURCE = "Chemical Engineering Plant Cost Index (CEPCI), annual averages"
CEPCI_BASE_YEAR = 2001


@dataclass(frozen=True)
class Purchased:
    """log10(Cp0) = K1 + K2 log10 A + K3 (log10 A)^2, A in ``attribute`` units."""
    K1: float
    K2: float
    K3: float
    amin: float
    amax: float
    attribute: str
    source: str


@dataclass(frozen=True)
class Pressure:
    """log10(Fp) = C1 + C2 log10 P + C3 (log10 P)^2, P in barg; Fp=1 below pmin."""
    C1: float
    C2: float
    C3: float
    pmin: float
    pmax: float
    source: str


@dataclass(frozen=True)
class BareModule:
    B1: float
    B2: float
    source: str


@dataclass(frozen=True)
class FbmDirect:
    Fbm: float
    source: str


@dataclass(frozen=True)
class MaterialFactors:
    factors: dict[str, float]
    source: str


@dataclass(frozen=True)
class QuantityFactor:
    """Small-quantity surcharge for stacked identical items (column trays):
    log10(Fq) = a + b log10 N + c (log10 N)^2 for N < n_full items (Fq > 1),
    and Fq = 1 for N >= n_full."""
    a: float
    b: float
    c: float
    n_full: int
    source: str


@dataclass(frozen=True)
class Utility:
    kind: str                 # "heat" | "cool" | "power"
    price_per_GJ: float
    source: str
    T_supply: float = 0.0
    T_return: float = 0.0
    U: float = 0.0            # W/m^2/K for sizing a heater/cooler as an exchanger


# -- Purchased-cost correlations (Turton 4e Table A.1) ----------------------
PURCHASED: dict[str, Purchased] = {
    "heat_exchanger": Purchased(
        4.3247, -0.3030, 0.1634, 10.0, 1000.0, "area_m2",
        "Turton 4e Table A.1, shell-and-tube floating head (10-1000 m^2)"),
    "vessel_vertical": Purchased(
        3.4974, 0.4485, 0.1074, 0.3, 520.0, "volume_m3",
        "Turton 4e Table A.1, vertical process vessel (0.3-520 m^3)"),
    "vessel_horizontal": Purchased(
        3.5565, 0.3776, 0.0905, 0.1, 628.0, "volume_m3",
        "Turton 4e Table A.1, horizontal process vessel (0.1-628 m^3)"),
    "pump_centrifugal": Purchased(
        3.3892, 0.0536, 0.1538, 1.0, 300.0, "power_kW",
        "Turton 4e Table A.1, centrifugal pump (1-300 kW shaft)"),
    "compressor_centrifugal": Purchased(
        2.2897, 1.3604, -0.1027, 450.0, 3000.0, "power_kW",
        "Turton 4e Table A.1, centrifugal compressor (450-3000 kW)"),
    # Distillation towers use the same Turton correlation as vertical process
    # vessels ("Towers: tray and packed" shares the vertical-vessel constants);
    # the trays are costed separately, per tray cross-section area.
    "tray_sieve": Purchased(
        2.9949, 0.4465, 0.3961, 0.07, 12.3, "area_m2",
        "Turton 4e Table A.1, sieve trays (0.07-12.3 m^2 per tray)"),
    "turbine_axial": Purchased(
        2.7051, 1.4398, -0.1776, 100.0, 4000.0, "power_kW",
        "Turton 4e Table A.1, axial gas turbine (100-4000 kW shaft)"),
    # Fired heater capacity is the absorbed (process) duty Q in kW.
    "fired_heater": Purchased(
        7.3488, -1.1666, 0.2028, 1000.0, 100_000.0, "duty_kW",
        "Turton 4e Table A.1, non-reactive fired heater (1,000-100,000 kW duty)"),
    # Air cooler capacity is the bare-tube heat-transfer area.
    "air_cooler": Purchased(
        4.0336, 0.2341, 0.0497, 10.0, 10_000.0, "area_m2",
        "Turton 4e Table A.1, air cooler (10-10,000 m^2 bare-tube area)"),
    # Straight pipe, costed per length with a diameter power law. The attribute
    # is A = L_m * (D_m)^0.74, so the linear correlation (K2=1, K3=0) gives
    # Cp0 = 1416 * L * D^0.74 on the CEPCI-397 basis. Source: installed cost of
    # carbon-steel pipe ~ 880 * (D, m)^0.74 GBP/m (R.K. Sinnott, Coulson &
    # Richardson's Chemical Engineering Vol. 6, 4th ed. 2005, sec. 5.5 piping
    # costs; 2004 price basis, erection/fittings included), converted at
    # 1.8 USD/GBP and deflated from CEPCI 444 (2004) to the 397 (2001) basis:
    # 880 * 1.8 * 397/444 = 1416 USD per (m * m^0.74). Order-of-magnitude data
    # meant for sensitivity sweeps, like the price tables below.
    "pipe": Purchased(
        3.1510, 1.0, 0.0, 1e-3, 1e9, "L_m_x_D_m^0.74",
        "Sinnott, C&R Vol. 6, 4e (2005), sec. 5.5: installed CS pipe "
        "880*(D m)^0.74 GBP2004/m -> 1416 USD2001 per m*D^0.74"),
    # -- Solids-handling equipment (Hameed 2025 ch. 12 units). Turton 4e has no
    # solids coverage, so these come from secondary sources, converted to the
    # CEPCI-397 (2001) USD basis. They are ORDER-OF-MAGNITUDE correlations for
    # screening/sensitivity, not vendor quotes — confidence is flagged per row.
    #
    # Cyclone, per cyclone, on the actual gas volumetric flow: Couper, Penney,
    # Fair & Walas, *Chemical Process Equipment*, 3e (2012), ch. 21 cost data
    # (after Walas 1988): heavy-duty gas cyclone, C[k$, mid-1985] =
    # 1.39 Q^0.98 with Q in kSCFM (2-40). Converted with 1 m^3/s = 2.1189
    # kSCFM (actual ~ standard flow assumed; adds <10% error below ~70 C) and
    # CEPCI 325.3 (1985) -> 397: Cp0 = 3541 * Q^0.98 USD2001, Q in m^3/s.
    # MODERATE-LOW confidence: the 0.98 exponent and 1.39 coefficient are from
    # the Walas table as recalled; verify against a physical copy before
    # relying on absolute cyclone capex.
    "cyclone": Purchased(
        3.5491, 0.98, 0.0, 0.94, 18.9, "gas_m3_s",
        "Couper et al. 3e ch. 21 (Walas): heavy-duty cyclone 1.39*Q^0.98 "
        "k$1985, Q in kSCFM (2-40) -> 3541*Q^0.98 USD2001, Q in m^3/s "
        "[order-of-magnitude]"),
    # Rotary vacuum drum filter, on filtration area: power-law anchored to the
    # rotary drum vacuum filter purchased-cost chart of Peters, Timmerhaus &
    # West, *Plant Design and Economics for Chemical Engineers*, 5e (2003),
    # ch. 14 (solids separation): ~$250k purchased at 10 m^2 in ~2002 dollars
    # (CEPCI ~396 ~ the 397 basis), with a 0.54 capacity exponent (mid-range
    # for package solids equipment; six-tenths-rule family). LOW confidence on
    # the absolute level (chart read from memory) — flagged for verification.
    "rotary_vacuum_filter": Purchased(
        4.8580, 0.54, 0.0, 1.0, 80.0, "area_m2",
        "Peters, Timmerhaus & West 5e ch. 14 chart anchor: ~$250k at 10 m^2 "
        "(USD2002~CEPCI-397), exponent 0.54 [order-of-magnitude]"),
    # Baghouse (fabric filter), on gross cloth area: EPA *Air Pollution
    # Control Cost Manual*, 6e (EPA/452/B-02-001, 2002), sec. 6 ch. 1: pulse-
    # jet flange-to-flange equipment cost ~ $8-12/ft^2 of cloth plus bags
    # ~$1-2/ft^2 (1998$, CEPCI ~390 ~ the 397 basis) -> ~ $160/m^2, linear in
    # area. MODERATE confidence (the manual's correlations are linear in cloth
    # area with a small fixed intercept this row drops).
    "baghouse": Purchased(
        2.2041, 1.0, 0.0, 10.0, 10_000.0, "area_m2",
        "EPA APC Cost Manual 6e sec. 6 ch. 1: pulse-jet baghouse ~$15/ft^2 "
        "cloth incl. bags (1998$~CEPCI-397) -> 160*A USD2001, A in m^2"),
}

# -- Pressure-factor correlations (Turton 4e Table A.2) ---------------------
PRESSURE: dict[str, Pressure] = {
    "heat_exchanger": Pressure(0.03881, -0.11272, 0.08183, 5.0, 140.0,
                               "Turton 4e Table A.2, shell-and-tube"),
    "pump_centrifugal": Pressure(-0.3935, 0.3957, -0.00226, 10.0, 100.0,
                                 "Turton 4e Table A.2, centrifugal pump"),
    # AUDIT NOTE (2026-06-11): this C-triple may belong to Table A.2's
    # *pyrolysis furnace* row; the non-reactive fired-heater row is possibly
    # (0.1347, -0.2368, 0.1021). Both give Fp ~ 1.0-1.15 below ~40 barg, so the
    # impact is small at typical process pressures; verify against a physical
    # Turton 4e before relying on >40 barg fired-heater costs.
    "fired_heater": Pressure(0.1017, -0.1957, 0.09403, 10.0, 200.0,
                             "Turton 4e Table A.2, non-reactive fired heater "
                             "(tube-side P, 10-200 barg)"),
    "air_cooler": Pressure(-0.1250, 0.15361, -0.02861, 10.0, 100.0,
                           "Turton 4e Table A.2, air cooler (tube-side P, 10-100 barg)"),
    # Vessels use a wall-stress formula (costing.vessel_pressure_factor);
    # compressors take Fbm directly (no Fp).
}

# -- Bare-module factors (Turton 4e Table A.4) ------------------------------
BARE_MODULE: dict[str, BareModule] = {
    "heat_exchanger": BareModule(1.63, 1.66, "Turton 4e Table A.4"),
    "air_cooler": BareModule(0.96, 1.21, "Turton 4e Table A.4, air cooler"),
    "vessel_vertical": BareModule(2.25, 1.82, "Turton 4e Table A.4, vertical vessel"),
    "vessel_horizontal": BareModule(1.49, 1.52, "Turton 4e Table A.4, horizontal vessel"),
    "pump_centrifugal": BareModule(1.89, 1.35, "Turton 4e Table A.4, pump"),
}

# Equipment costed with a direct bare-module factor: Cbm = Cp0 * Fbm.
FBM_DIRECT: dict[str, FbmDirect] = {
    "compressor_centrifugal": FbmDirect(2.7, "Turton 4e, centrifugal compressor + drive (CS)"),
    "tray_sieve": FbmDirect(1.0, "Turton 4e Table A.6, sieve tray, carbon steel"),
    # Turbines take a direct bare-module factor, like compressors (no Fp/Fm
    # correlation): Turton 4e Fig. A.19 gives FBM ~ 3.5 for an axial gas
    # turbine in carbon steel.
    "turbine_axial": FbmDirect(3.5, "Turton 4e Fig. A.19, axial gas turbine, CS"),
    # Fired heaters take a direct bare-module factor times the pressure factor:
    # Cbm = Cp0 * Fbm * Fp (Turton 4e Eq. A.6 for fired equipment; the steam-
    # superheat factor Ft does not apply to process fired heaters, Ft = 1).
    "fired_heater": FbmDirect(2.13, "Turton 4e Table A.7, non-reactive fired heater, CS"),
    # The pipe correlation is an *installed* cost (erection, fittings and
    # supports already included in the Sinnott figure), so Fbm = 1.
    "pipe": FbmDirect(1.0, "Sinnott C&R Vol. 6 4e sec. 5.5 — correlation is "
                           "already an installed cost"),
    # Solids equipment: simple fabricated/package items take a Hand-type
    # installation factor of ~2 on the f.o.b. purchased cost (Towler &
    # Sinnott, Chemical Engineering Design, 2e, ch. 7, installation factors
    # for miscellaneous fabricated equipment).
    "cyclone": FbmDirect(2.0, "Towler & Sinnott 2e ch. 7, Hand-type factor ~2 "
                              "for fabricated package equipment"),
    "rotary_vacuum_filter": FbmDirect(
        2.0, "Towler & Sinnott 2e ch. 7, Hand-type factor ~2 for package "
             "solids equipment"),
    # Baghouses: the EPA manual's own total-capital build-up multiplies the
    # flange-to-flange cost by ~2.17 (direct + indirect installation; EPA APC
    # Cost Manual 6e sec. 6 ch. 1 cost-factor table).
    "baghouse": FbmDirect(2.17, "EPA APC Cost Manual 6e sec. 6 ch. 1: TCI ~ "
                                "2.17 x flange-to-flange equipment cost"),
}

# Quantity factor for stacked/multiple identical items (trays): Turton 4e
# Eq. A.5 — Cbm = N * Cp0 * Fbm * Fq.
QUANTITY_FACTOR: dict[str, QuantityFactor] = {
    "tray_sieve": QuantityFactor(
        0.4771, 0.08516, -0.3473, 20,
        "Turton 4e Eq. A.5: log10 Fq = 0.4771 + 0.08516 log10 N - 0.3473 (log10 N)^2 "
        "for N < 20 trays; Fq = 1 for N >= 20"),
}

# -- Material factors Fm (Turton 4e Table A.3) ------------------------------
MATERIAL: dict[str, MaterialFactors] = {
    "heat_exchanger": MaterialFactors({"CS": 1.0, "SS": 2.73}, "Turton 4e Table A.3, S&T"),
    "vessel_vertical": MaterialFactors({"CS": 1.0, "SS": 3.12}, "Turton 4e Table A.3, vessel"),
    "vessel_horizontal": MaterialFactors({"CS": 1.0, "SS": 3.12}, "Turton 4e Table A.3, vessel"),
    "pump_centrifugal": MaterialFactors({"CS": 1.0, "SS": 2.28}, "Turton 4e Table A.3, pump"),
    # Air coolers: carbon-steel tubes only for now; add alloy factors from
    # Turton 4e Fig. A.18 (identification-number chart) if alloy service is needed.
    "air_cooler": MaterialFactors({"CS": 1.0}, "Turton 4e Fig. A.18, air cooler, CS tubes"),
}

# -- Utilities (Turton 4e Table 8.3 prices; U from Table 11.11) -------------
UTILITIES: dict[str, Utility] = {
    "cooling_water": Utility("cool", 0.378, "Turton 4e Table 8.3 / 11.11",
                             T_supply=303.15, T_return=316.15, U=850.0),
    "refrigerant_lowtemp": Utility("cool", 9.50, "Turton 4e Table 8.3, ~ -40 C refrigeration",
                                   T_supply=233.15, T_return=243.15, U=600.0),
    "low_pressure_steam": Utility("heat", 14.05, "Turton 4e Table 8.3, LPS 5 barg",
                                  T_supply=433.0, T_return=433.0, U=900.0),
    "fired_heat": Utility("heat", 11.10, "Turton 4e Table 8.3, fired-heater fuel",
                          T_supply=1100.0, T_return=1100.0, U=60.0),
    "electricity": Utility("power", 18.72, "Turton 4e Table 8.3 ($0.0674/kWh)"),
}

# -- Feed and product prices ($/kg). Representative ~2023 market values; these
# are the variables the analysis track sweeps, not fixed constants.
PRICES_PER_KG: dict[str, float] = {
    "hydrogen": 1.50,    # grey/SMR-ish merchant H2
    "nitrogen": 0.10,    # air-separation N2
    "ammonia": 0.50,     # ~ $500/tonne NH3
    "benzene": 0.95,     # ~ $950/tonne merchant benzene
    "toluene": 0.80,     # ~ $800/tonne merchant toluene
}
PRICES_SOURCE = "representative ~2023 merchant prices; intended for sensitivity sweeps"

# -- Molar masses (kg/mol), CIAAW standard atomic weights -------------------
# A small local table for the components the validated examples use; anything
# else falls back to the `chemicals` database via molar_mass() below.
MOLAR_MASS: dict[str, float] = {
    "hydrogen": 0.002016,
    "nitrogen": 0.028014,
    "ammonia": 0.017031,
    "argon": 0.039948,
    "water": 0.018015,
    "benzene": 0.078112,
    "toluene": 0.092138,
    "p-xylene": 0.106165,
}


def molar_mass(component: str) -> float:
    """Molar mass of ``component`` in kg/mol: the local table first, then the
    `chemicals` database (cached in caldyr.core.components_db). An unknown
    component raises a clear error naming it — never a bare KeyError."""
    mw = MOLAR_MASS.get(component)
    if mw is not None:
        return mw
    from ..core.components_db import UnknownComponentError
    from ..core.components_db import molar_mass as _chemicals_mw

    try:
        return _chemicals_mw(component)
    except UnknownComponentError as exc:
        raise ValueError(
            f"no molar mass available for component {component!r}: it is not in "
            f"economics.data.MOLAR_MASS and the chemicals database does not "
            f"recognize it"
        ) from exc
