"""Packed-column hydraulic sizing: diameter from the Eckert/Strigle generalized
pressure-drop correlation (GPDC) and height from HETP rules of thumb.

A packed tower is sized differently from a tray tower (``economics.tray_sizing``):

* **Diameter** is set by the flood point of the packing. The GPDC (Eckert,
  modified by Strigle; Perry 8e Fig. 14-55) correlates a capacity parameter

      CP = C_s · F_p^0.5 · nu^0.05            (US units: C_s [ft/s], F_p [ft^-1],
                                               nu = liquid kinematic viscosity [cS])

  with the flow parameter ``F_LV = (L/G)·sqrt(rho_G/rho_L)`` (the same kinetic-
  energy ratio the tray Fair correlation uses), where ``C_s = u_s·sqrt(rho_G/
  (rho_L-rho_G))`` is the Souders-Brown capacity factor on the tower
  cross-section. The flood point is read at the flood pressure drop given by the
  Kister & Gill equation (Perry eq. 14-142)

      dP_flood = 0.12 · F_p^0.7              [in H2O / ft of packing]

  The flood gas velocity then follows from CP, and the tower is designed at a
  fraction of flood (default 70 %, the standard packed-tower design point;
  Strigle / Kister). The GPDC flood locus here is a log-log digitization of the
  published Strigle chart anchored to its worked example (Perry Example 13:
  F_LV=0.207, CP=1.01 at dP=0.38) — the same digitize-the-chart approach the
  tray sizer uses for the Fair flooding curve.

* **Height** is the theoretical-stage count times the HETP (height equivalent to
  a theoretical plate). HETP is taken from the well-tested rules of thumb in
  Perry 8e (eqs. 14-156, 14-158): ``HETP ~ 18·D_p`` (random packing, D_p the
  nominal packing size, m) or equivalently ``HETP ~ 93/a_p`` for Pall rings
  (a_p the specific surface area, m^2/m^3); structured packing uses
  ``HETP ~ 100·C/a_p + 0.10``. For water-rich systems (surface tension
  ~70 mN/m) these are **doubled** (Perry, underwetting), and for intermediate
  surface tensions multiplied by 1.5.

Packing factors ``F_p`` and specific areas ``a_p`` are from Perry 8e Tables
14-13 (random) and 14-14 (structured). Sanity anchor: the Hameed (2025)
sec. 9.1 SO2 absorber (206 kmol/h gas, 1.3e5 kg/h water, 1 atm) is HYSYS-
packed to 1.285 m at 80 % of capacity (fig. 9.10); this correlation with 50 mm
metal Pall rings at 70 % of flood gives ~1.18 m — packed and tray
(``tray_sizing`` gives ~1.4 m) diameters for the same service agreeing within
~15 % (see tests/test_m20_packed_sizing.py).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .tray_sizing import StageLoad, _mw

# Random and structured packings: (packing factor F_p [m^-1], specific surface
# area a_p [m^2/m^3], nominal size [m], structured?). Perry 8e Tables 14-13/14-14.
PACKINGS: dict[str, tuple[float, float, float, bool]] = {
    # metal Pall rings (a.k.a. Flexi-/Ballast rings) — the workhorse random packing
    "pall_metal_25mm": (183.0, 205.0, 0.025, False),
    "pall_metal_38mm": (131.0, 130.0, 0.038, False),
    "pall_metal_50mm": (89.0, 105.0, 0.050, False),
    "pall_metal_90mm": (59.0, 66.0, 0.090, False),
    # metal Intalox (IMTP) — modern high-capacity random
    "imtp_25mm": (134.0, 207.0, 0.025, False),
    "imtp_50mm": (59.0, 98.0, 0.050, False),
    # ceramic Raschig rings — classic / corrosion service
    "raschig_ceramic_25mm": (587.0, 190.0, 0.025, False),
    "raschig_ceramic_50mm": (213.0, 92.0, 0.050, False),
    # ceramic Intalox saddles
    "intalox_ceramic_25mm": (302.0, 256.0, 0.025, False),
    "intalox_ceramic_50mm": (131.0, 118.0, 0.050, False),
    # structured (corrugated metal sheet, Mellapak-style)
    "mellapak_250y": (66.0, 250.0, 0.0, True),
    "mellapak_350y": (75.0, 350.0, 0.0, True),
    "flexipac_2y": (49.0, 220.0, 0.0, True),
}

_M_PER_FT = 0.3048
_FLOOD_FRACTION = 0.70       # packed-tower design point (fraction of flood)
_NU_DEFAULT_CS = 1.0         # liquid kinematic viscosity, cS (nu^0.05 ~ 1; weak)
_F_LV_MIN = 0.005            # GPDC abscissa validity (Perry Fig. 14-55)
_F_LV_MAX = 5.0
# Surface-tension thresholds for the HETP underwetting penalty (Perry 8e):
_SIGMA_AQUEOUS = 0.045       # N/m; above this -> water-rich, double the HETP
_SIGMA_INTERMEDIATE = 0.030  # N/m; above this -> x1.5 (amines/glycols)

# Strigle GPDC, capacity parameter CP on the dP = 1.0 in-H2O/ft pressure-drop
# curve, vs the flow parameter F_LV (log-log digitization of Perry Fig. 14-55,
# anchored to its worked Example 13). Other pressure drops scale as CP ~ dP^0.43.
_GPDC_FLV = (0.01, 0.03, 0.10, 0.20, 0.50, 1.00, 2.00, 4.00)
_GPDC_CP1 = (2.30, 2.05, 1.70, 1.50, 1.15, 0.90, 0.62, 0.40)
_GPDC_DP_EXP = 0.43


@dataclass
class PackedHydraulics:
    """Result of the GPDC flood calculation for a packed bed."""
    diameter_m: float
    area_m2: float
    stage: int                  # governing (largest-diameter) candidate stage
    u_flood_ms: float           # flooding superficial gas velocity
    u_design_ms: float          # design velocity (flood_fraction * u_flood)
    F_LV: float
    capacity_param: float       # CP at flood
    dP_flood_inH2O_ft: float    # Kister-Gill flood pressure drop
    rho_V: float                # kg/m^3
    rho_L: float                # kg/m^3
    flood_fraction: float
    notes: list[str] = field(default_factory=list)


def gpdc_capacity_param(F_LV: float, dP_flood: float) -> float:
    """Capacity parameter CP = C_s·F_p^0.5·nu^0.05 at the flood pressure drop,
    from the digitized Strigle GPDC. ``dP_flood`` in in-H2O per ft of packing."""
    x = min(max(F_LV, _F_LV_MIN), _F_LV_MAX)
    lx = math.log(x)
    # piecewise log-log interpolation on the dP = 1.0 curve
    cp1 = _GPDC_CP1[0]
    for i in range(1, len(_GPDC_FLV)):
        if x <= _GPDC_FLV[i] or i == len(_GPDC_FLV) - 1:
            x0, x1 = _GPDC_FLV[i - 1], _GPDC_FLV[i]
            y0, y1 = _GPDC_CP1[i - 1], _GPDC_CP1[i]
            frac = (lx - math.log(x0)) / (math.log(x1) - math.log(x0))
            cp1 = math.exp(math.log(y0) + frac * (math.log(y1) - math.log(y0)))
            break
    return cp1 * dP_flood ** _GPDC_DP_EXP


def get_packing(name: str) -> tuple[float, float, float, bool]:
    try:
        return PACKINGS[name]
    except KeyError as exc:
        raise ValueError(
            f"unknown packing {name!r}; choose one of {sorted(PACKINGS)}"
        ) from exc


def size_packed(pp: Any, load: StageLoad, packing: str, *,
                flood_fraction: float = _FLOOD_FRACTION,
                nu_cS: float = _NU_DEFAULT_CS) -> PackedHydraulics:
    """Tower diameter required by one candidate stage at ``flood_fraction`` of
    the GPDC flood point for ``packing``."""
    Fp_m, _ap, _size, _struct = get_packing(packing)
    Fp_ft = Fp_m * _M_PER_FT                      # m^-1 -> ft^-1
    mw_v, mw_l = _mw(load.y), _mw(load.x)
    rho_v = mw_v / pp.volume_vapor(load.T, load.P, load.y)
    rho_l = mw_l / pp.volume_liquid(load.T, load.P, load.x)
    notes: list[str] = []
    if rho_l <= rho_v:
        raise ValueError(
            f"packed sizing on stage {load.stage}: liquid density {rho_l:.1f} "
            f"<= vapor density {rho_v:.1f} kg/m^3 — the stage state is "
            f"unphysical (near-critical?); the GPDC does not apply"
        )
    f_lv = (load.L * mw_l) / (load.V * mw_v) * math.sqrt(rho_v / rho_l)
    if not _F_LV_MIN <= f_lv <= _F_LV_MAX:
        clamped = min(max(f_lv, _F_LV_MIN), _F_LV_MAX)
        notes.append(f"F_LV={f_lv:.3g} outside the GPDC range "
                     f"[{_F_LV_MIN}, {_F_LV_MAX}]; clamped to {clamped:.3g}")
        f_lv = clamped
    dP_flood = 0.12 * Fp_ft ** 0.7                # Kister-Gill, in-H2O/ft
    cp = gpdc_capacity_param(f_lv, dP_flood)
    # CP = C_s·F_p^0.5·nu^0.05 -> C_s [ft/s] -> SI; u_flood from Souders-Brown.
    c_s_ft = cp / (Fp_ft ** 0.5 * nu_cS ** 0.05)
    c_s = c_s_ft * _M_PER_FT                      # ft/s -> m/s
    u_flood = c_s * math.sqrt((rho_l - rho_v) / rho_v)
    u_design = flood_fraction * u_flood
    q_v = load.V * pp.volume_vapor(load.T, load.P, load.y)    # m^3/s
    area = q_v / u_design
    diameter = math.sqrt(4.0 * area / math.pi)
    notes.insert(0, (
        f"GPDC flood (Perry 8e Fig. 14-55) at stage {load.stage}, {packing}: "
        f"F_LV={f_lv:.3f}, F_p={Fp_m:.0f} m^-1, dP_flood={dP_flood:.2f} in-H2O/ft "
        f"(Kister-Gill), CP={cp:.3f}, u_flood={u_flood:.2f} m/s; designed at "
        f"{flood_fraction:.0%} of flood"))
    return PackedHydraulics(
        diameter_m=diameter, area_m2=area, stage=load.stage,
        u_flood_ms=u_flood, u_design_ms=u_design, F_LV=f_lv,
        capacity_param=cp, dP_flood_inH2O_ft=dP_flood,
        rho_V=rho_v, rho_L=rho_l, flood_fraction=flood_fraction, notes=notes)


def governing_packed(pp: Any, loads: list[StageLoad], packing: str, *,
                     flood_fraction: float = _FLOOD_FRACTION,
                     nu_cS: float = _NU_DEFAULT_CS) -> PackedHydraulics:
    """Largest-diameter (governing) result over the candidate stage loads."""
    if not loads:
        raise ValueError("packed sizing needs at least one candidate stage load")
    results = [size_packed(pp, ld, packing, flood_fraction=flood_fraction,
                           nu_cS=nu_cS) for ld in loads]
    return max(results, key=lambda r: r.diameter_m)


def hetp_m(packing: str, *, aqueous: bool = False,
           sigma_N_m: float | None = None) -> tuple[float, list[str]]:
    """Height equivalent to a theoretical plate (m), from the Perry 8e rules of
    thumb. ``HETP = 18·D_p`` (random, D_p the nominal size) or ``93/a_p`` for
    Pall rings; structured ``= 100/a_p + 0.10``. Doubled for water-rich systems
    (surface tension > ~45 mN/m) and x1.5 for intermediate (> ~30 mN/m), per the
    underwetting penalty. Pass ``aqueous=True`` (a convenience for water systems)
    or an explicit ``sigma_N_m``."""
    Fp_m, ap, size, struct = get_packing(packing)
    notes: list[str] = []
    if struct:
        base = 100.0 / ap + 0.10
        notes.append(f"structured HETP = 100/a_p + 0.10 (a_p={ap:.0f} m^2/m^3)")
    elif size > 0.0:
        base = 18.0 * size
        notes.append(f"random HETP = 18·D_p (D_p={size * 1000:.0f} mm)")
    else:
        base = 93.0 / ap
        notes.append(f"random HETP = 93/a_p (a_p={ap:.0f} m^2/m^3)")
    factor = 1.0
    sigma = sigma_N_m if sigma_N_m is not None else (
        _SIGMA_AQUEOUS + 0.01 if aqueous else 0.0)
    if sigma >= _SIGMA_AQUEOUS:
        factor, why = 2.0, "water-rich (sigma>45 mN/m): HETP x2 (underwetting)"
    elif sigma >= _SIGMA_INTERMEDIATE:
        factor, why = 1.5, "intermediate sigma (30-45 mN/m): HETP x1.5"
    else:
        why = "low-sigma organic: no underwetting penalty"
    notes.append(why)
    return base * factor, notes


def packing_volume_m3(diameter_m: float, height_m: float) -> float:
    """Volume of packing in the bed (tower cross-section x packed height)."""
    return math.pi / 4.0 * diameter_m ** 2 * height_m
