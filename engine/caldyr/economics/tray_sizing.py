"""Tray hydraulic sizing: column diameter from Fair's flooding correlation.

The tower diameter is set by the flooding velocity of the governing tray.
Fair's (1961) flooding capacity chart is used in the closed form fitted by
Lygeros & Magoulas, as given in Seader, Henley & Roper, *Separation Process
Principles* 3e, sec. 6.6 (eq. 6-42)::

    C [m/s] = 0.0105 + 8.127e-4 * T_s^0.755 * exp(-1.463 * F_LV^0.842)

with ``T_s`` the tray spacing in mm and ``F_LV`` the kinetic-energy flow
ratio ``(L M_L / V M_V) * sqrt(rho_V / rho_L)`` (chart validity 0.01..1.0;
out-of-range values are clamped and noted). The capacity factor is corrected
in Fair's method by surface tension, foaming and hole-area factors
``C_F = C * F_ST * F_F * F_HA``; with no surface-tension model in the
property packages all three default to 1 (sigma = 20 dyn/cm, non-foaming,
hole/active area >= 0.10 — each documented in the notes), then::

    u_flood = C_F * sqrt((rho_L - rho_V) / rho_V)        # on the net area

The design velocity is a fraction of flooding (default **80%**, the classic
sieve-tray design point; Seader 3e sec. 6.6, Wankat 3e ch. 10) on the net
area (total cross-section minus the downcomer, default 10%). Vapor and
liquid densities come from the *per-phase* property calls
(``pp.volume_vapor`` / ``pp.volume_liquid``) at the actual stage state, so a
stage off saturation (absorbers) is handled correctly.

The governing tray is found by sizing the candidate stages a column model
exposes — for profile-bearing columns (RigorousColumn, Absorber,
ReboiledAbsorber) the top-most and bottom-most *trays*; for the
ShortcutColumn the constant-molal-overflow section loads — and taking the
largest diameter.

Sanity anchor: for the SO2/air/water absorber of Hameed (2025) sec. 9.1
(206 kmol/h gas, 1.3e5 kg/h water, 20 stages at 1 atm) this correlation
gives D ~ 1.4 m vs the book's HYSYS packed-column design of 1.285 m at 80%
of capacity — tray and packed diameters for the same service agreeing to
~10% (see tests/test_m12_tray_sizing.py).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from . import data

# Fair-chart validity range for F_LV (Seader 3e fig. 6.24).
_F_LV_MIN = 0.01
_F_LV_MAX = 1.0


@dataclass
class StageLoad:
    """Vapor/liquid loads and state of one candidate (governing) tray."""
    stage: int                  # 1-based stage label (for the notes)
    V: float                    # vapor molar flow through the tray, mol/s
    L: float                    # liquid molar flow off the tray, mol/s
    T: float                    # K
    P: float                    # Pa
    x: dict[str, float]         # liquid composition
    y: dict[str, float]         # vapor composition


@dataclass
class TrayHydraulics:
    """Result of the Fair flooding calculation at the governing tray."""
    diameter_m: float
    area_m2: float              # total tower cross-section
    stage: int                  # governing (largest-diameter) stage, 1-based
    u_flood_ms: float           # flooding velocity on the net area
    u_design_ms: float          # design velocity (flood_fraction * u_flood)
    F_LV: float
    C_F: float                  # capacity factor, m/s
    rho_V: float                # kg/m^3
    rho_L: float                # kg/m^3
    flood_fraction: float
    notes: list[str] = field(default_factory=list)


def fair_capacity(F_LV: float, tray_spacing_m: float) -> float:
    """Fair flooding capacity parameter C (m/s) — the Lygeros & Magoulas fit
    of Fair's chart quoted in Seader 3e eq. (6-42); spacing in m here."""
    ts_mm = tray_spacing_m * 1000.0
    return 0.0105 + 8.127e-4 * ts_mm ** 0.755 * math.exp(
        -1.463 * F_LV ** 0.842)


def _mw(comp: dict[str, float]) -> float:
    """Mean molar mass of a composition, kg/mol (data.molar_mass is kg/mol)."""
    return sum(frac * data.molar_mass(c)
               for c, frac in comp.items() if frac > 0.0)


def size_tray(pp: Any, load: StageLoad, *, tray_spacing_m: float = 0.61,
              flood_fraction: float = 0.8,
              downcomer_frac: float = 0.10) -> TrayHydraulics:
    """Diameter required by one tray at ``flood_fraction`` of Fair flooding."""
    notes: list[str] = []
    mw_v, mw_l = _mw(load.y), _mw(load.x)
    rho_v = mw_v / pp.volume_vapor(load.T, load.P, load.y)
    rho_l = mw_l / pp.volume_liquid(load.T, load.P, load.x)
    if rho_l <= rho_v:
        raise ValueError(
            f"tray sizing on stage {load.stage}: liquid density "
            f"{rho_l:.1f} <= vapor density {rho_v:.1f} kg/m^3 — the stage "
            f"state is unphysical (near-critical?); Fair's correlation does "
            f"not apply"
        )
    f_lv = (load.L * mw_l) / (load.V * mw_v) * math.sqrt(rho_v / rho_l)
    if not _F_LV_MIN <= f_lv <= _F_LV_MAX:
        clamped = min(max(f_lv, _F_LV_MIN), _F_LV_MAX)
        notes.append(f"F_LV={f_lv:.3g} outside the Fair chart range "
                     f"[{_F_LV_MIN}, {_F_LV_MAX}]; clamped to {clamped:.3g}")
        f_lv = clamped
    # C_F = C * F_ST * F_F * F_HA with all three factors 1 (documented):
    # sigma assumed 20 dyn/cm (F_ST=(sigma/20)^0.2=1), non-foaming (F_F=1),
    # hole/active area >= 0.10 (F_HA=1).
    c_f = fair_capacity(f_lv, tray_spacing_m)
    u_flood = c_f * math.sqrt((rho_l - rho_v) / rho_v)
    u_design = flood_fraction * u_flood
    q_v = load.V * pp.volume_vapor(load.T, load.P, load.y)   # m^3/s
    a_net = q_v / u_design
    a_total = a_net / (1.0 - downcomer_frac)
    diameter = math.sqrt(4.0 * a_total / math.pi)
    notes.insert(0, (
        f"Fair flooding (Seader 3e eq. 6-42) at stage {load.stage}: "
        f"F_LV={f_lv:.3f}, C_F={c_f:.4f} m/s, rho_V={rho_v:.2f}, "
        f"rho_L={rho_l:.1f} kg/m^3, u_flood={u_flood:.2f} m/s at "
        f"{tray_spacing_m:.2f} m spacing; designed at "
        f"{flood_fraction:.0%} of flood with {downcomer_frac:.0%} downcomer"
    ))
    return TrayHydraulics(
        diameter_m=diameter, area_m2=a_total, stage=load.stage,
        u_flood_ms=u_flood, u_design_ms=u_design, F_LV=f_lv, C_F=c_f,
        rho_V=rho_v, rho_L=rho_l, flood_fraction=flood_fraction, notes=notes,
    )


def governing_tray(pp: Any, loads: list[StageLoad], *,
                   tray_spacing_m: float = 0.61, flood_fraction: float = 0.8,
                   downcomer_frac: float = 0.10) -> TrayHydraulics:
    """Size every candidate tray and return the largest-diameter result (the
    governing tray sets the tower diameter)."""
    if not loads:
        raise ValueError("tray sizing needs at least one candidate stage load")
    results = [size_tray(pp, load, tray_spacing_m=tray_spacing_m,
                         flood_fraction=flood_fraction,
                         downcomer_frac=downcomer_frac) for load in loads]
    return max(results, key=lambda r: r.diameter_m)


def loads_from_profiles(design: dict[str, Any],
                        candidates: list[int]) -> list[StageLoad]:
    """Candidate stage loads from a column's converged ``design`` profiles
    (``T/P/L/V/x/y_profile``); ``candidates`` are 0-based stage indices."""
    return [StageLoad(
        stage=j + 1,
        V=design["V_profile"][j], L=design["L_profile"][j],
        T=design["T_profile"][j], P=design["P_profile"][j],
        x=design["x_profile"][j], y=design["y_profile"][j],
    ) for j in candidates]
