"""Pressure-relief-valve orifice sizing per API Standard 520 Part I / API 526.

The HYSYS "Safety Analysis -> PSV sizing" equivalent (Hameed, *Chemical Process
Simulations using Aspen HYSYS*, Wiley 2025, sec. 7.3): given a relieving rate
and the relieving conditions, compute the required effective discharge area
and pick the standard API 526 orifice letter.

Vapor, critical (choked) flow — API 520 Part I (9th ed. 2014 / 10th ed. 2020),
sec. 5.6.3, sizing equation in SI form::

    A = W / (C Kd P1 Kb Kc) sqrt(T Z / M)   [A mm^2, W kg/h, P1 kPaa, M g/mol]
    C = 0.03948 sqrt( k (2/(k+1))^((k+1)/(k-1)) )

:func:`relief_vapor` implements the identical relation in pure SI units (kg/s,
Pa, m^2) via the choked ideal-(Z-corrected)-gas mass flux

    G = P1 sqrt( k M / (Z R T) * (2/(k+1))^((k+1)/(k-1)) ),   A = W / (Kd Kb Kc G)

(the two forms are algebraically identical — API's 0.03948 is exactly
3.6 / sqrt(1000 R) = 0.039481 in those mixed units — so they agree to the
rounding of the published 0.03948 constant, ~0.02%). Only
critical flow is implemented: a backpressure above the critical-flow pressure
``P1 (2/(k+1))^(k/(k-1))`` raises :class:`ReliefSizingError` (API 520's
subcritical equation, sec. 5.6.4, is future work).

Liquid (capacity-certified valves) — API 520 Part I sec. 5.8, SI form
``A = 11.78 Q / (Kd Kw Kc Kv) sqrt(G / (p1 - p2))`` [A mm^2, Q L/min, dp kPa];
:func:`relief_liquid` implements the same relation in pure SI via the orifice
equation ``A = W / (Kd Kw Kc Kv sqrt(2 rho (P1 - P2)))``.

Standard orifices — API Standard 526 (7th ed. 2017), letters D through T with
the effective discharge areas (in^2) used with the API 520 effective
coefficients: D 0.110, E 0.196, F 0.307, G 0.503, H 0.785, J 1.287, K 1.838,
L 2.853, M 3.60, N 4.34, P 6.38, Q 11.05, R 16.0, T 26.0.

This is a *rating/analysis* utility (a relief valve is sized against a
scenario, not solved in the flowsheet graph), so it lives in
``caldyr.analysis`` rather than ``caldyr.unitops``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

_R = 8.314462618              # J/mol/K
_IN2_TO_M2 = 6.4516e-4

#: API 526 standard effective orifices, (letter, effective area in^2), D -> T.
API_526_ORIFICES: tuple[tuple[str, float], ...] = (
    ("D", 0.110), ("E", 0.196), ("F", 0.307), ("G", 0.503), ("H", 0.785),
    ("J", 1.287), ("K", 1.838), ("L", 2.853), ("M", 3.60), ("N", 4.34),
    ("P", 6.38), ("Q", 11.05), ("R", 16.0), ("T", 26.0),
)


class ReliefSizingError(ValueError):
    """Invalid or out-of-scope relief-sizing inputs (e.g. subcritical vapor
    flow, which API 520 sec. 5.6.4 covers but this v1 does not)."""


@dataclass
class ReliefResult:
    """Required orifice area and the selected API 526 standard orifice.

    ``orifice`` is ``None`` when the requirement exceeds the largest standard
    orifice (T, 26 in^2) — multiple valves are then required (noted).
    ``capacity_used`` = required / selected area (the HYSYS "Capacity Used %"
    divided by 100).
    """
    area_m2: float                       # required effective discharge area
    orifice: str | None                  # API 526 letter, D..T
    orifice_area_m2: float | None        # selected standard effective area
    capacity_used: float | None          # required / selected
    phase: str                           # "vapor" | "liquid"
    critical: bool | None                # vapor flow regime (None for liquid)
    details: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _select_orifice(area_m2: float, result: ReliefResult) -> ReliefResult:
    for letter, area_in2 in API_526_ORIFICES:
        std = area_in2 * _IN2_TO_M2
        if std >= area_m2:
            result.orifice = letter
            result.orifice_area_m2 = std
            result.capacity_used = area_m2 / std
            return result
    result.notes.append(
        f"required area {area_m2 * 1e4:.1f} cm^2 exceeds the largest API 526 "
        f"orifice (T, {26.0 * _IN2_TO_M2 * 1e4:.1f} cm^2); use multiple valves "
        f"in parallel"
    )
    return result


def relief_vapor(W: float, T: float, M: float, Z: float, k: float, P1: float, *,
                 backpressure: float = 101_325.0, Kd: float = 0.975,
                 Kb: float = 1.0, Kc: float = 1.0) -> ReliefResult:
    """Required orifice area for vapor relief in critical (choked) flow,
    API 520 Part I sec. 5.6.3 (see module docstring for the equation).

    Parameters (SI):
      W : relieving mass flow, kg/s
      T : relieving temperature, K
      M : molar mass, kg/mol
      Z : compressibility factor at relieving conditions
      k : ideal-gas specific-heat ratio Cp/Cv
      P1 : upstream relieving pressure, Pa **absolute** (set pressure +
           allowable overpressure + atmospheric)
      backpressure : downstream total backpressure, Pa absolute (default 1 atm)
      Kd : effective discharge coefficient (API 520 default 0.975 for a relief
           valve without a rupture disk)
      Kb : backpressure correction (1.0 for conventional valves in critical
           flow; balanced-bellows values from API 520 Fig. 30)
      Kc : combination factor (1.0 alone; 0.9 with an upstream rupture disk)
    """
    for name, value in (("W", W), ("T", T), ("M", M), ("Z", Z), ("P1", P1),
                        ("Kd", Kd), ("Kb", Kb), ("Kc", Kc)):
        if value <= 0.0 or not math.isfinite(value):
            raise ReliefSizingError(f"relief_vapor: {name} must be positive (got {value})")
    if k <= 1.0:
        raise ReliefSizingError(f"relief_vapor: k = Cp/Cv must exceed 1 (got {k})")
    if backpressure < 0.0:
        raise ReliefSizingError(
            f"relief_vapor: backpressure must be >= 0 Pa absolute (got {backpressure})"
        )

    p_crit_ratio = (2.0 / (k + 1.0)) ** (k / (k - 1.0))
    p_crit = p_crit_ratio * P1
    if backpressure > p_crit:
        raise ReliefSizingError(
            f"relief_vapor: backpressure {backpressure:.4g} Pa exceeds the "
            f"critical-flow pressure {p_crit:.4g} Pa (= {p_crit_ratio:.3f} P1) "
            f"— the flow is subcritical. API 520 sec. 5.6.4 subcritical sizing "
            f"is not implemented in this v1."
        )

    # Choked mass flux through the orifice (API 520 critical vapor equation,
    # pure-SI form; identical to the 13160/0.03948 kg/h-kPa-mm^2 form).
    flux = P1 * math.sqrt(
        k * M / (Z * _R * T) * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0))
    )
    area = W / (Kd * Kb * Kc * flux)
    result = ReliefResult(
        area_m2=area, orifice=None, orifice_area_m2=None, capacity_used=None,
        phase="vapor", critical=True,
        details={"mass_flux_kg_m2_s": flux, "P_critical_Pa": p_crit,
                 "critical_pressure_ratio": p_crit_ratio, "k": k, "Z": Z},
    )
    return _select_orifice(area, result)


def relief_liquid(W: float, rho: float, P1: float, P2: float, *,
                  Kd: float = 0.65, Kw: float = 1.0, Kc: float = 1.0,
                  Kv: float = 1.0) -> ReliefResult:
    """Required orifice area for liquid relief (capacity-certified valves),
    API 520 Part I sec. 5.8 (see module docstring for the equation).

    Parameters (SI):
      W : relieving mass flow, kg/s
      rho : liquid density at relieving conditions, kg/m^3
      P1 : upstream relieving pressure, Pa absolute
      P2 : total backpressure, Pa absolute
      Kd : effective discharge coefficient (API 520 default 0.65 for liquids)
      Kw : backpressure correction (1.0 for conventional valves /
           atmospheric discharge; API 520 Fig. 31 for balanced valves)
      Kc : combination factor (0.9 with an upstream rupture disk)
      Kv : viscosity correction (1.0 for nonviscous; API 520 Fig. 36 /
           Eq. 5.10 once a trial area fixes the orifice Reynolds number)
    """
    for name, value in (("W", W), ("rho", rho), ("P1", P1),
                        ("Kd", Kd), ("Kw", Kw), ("Kc", Kc), ("Kv", Kv)):
        if value <= 0.0 or not math.isfinite(value):
            raise ReliefSizingError(f"relief_liquid: {name} must be positive (got {value})")
    if P2 < 0.0 or P2 >= P1:
        raise ReliefSizingError(
            f"relief_liquid: need 0 <= P2 < P1 (got P1={P1}, P2={P2})"
        )

    area = W / (Kd * Kw * Kc * Kv * math.sqrt(2.0 * rho * (P1 - P2)))
    result = ReliefResult(
        area_m2=area, orifice=None, orifice_area_m2=None, capacity_used=None,
        phase="liquid", critical=None,
        details={"dP_Pa": P1 - P2, "rho_kg_m3": rho},
    )
    return _select_orifice(area, result)
