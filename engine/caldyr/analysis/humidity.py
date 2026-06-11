"""Psychrometrics: moist-air state from one humidity specification.

The HYSYS "Saturator / humidity calculations" workflow of Hameed, *Chemical
Process Simulations using Aspen HYSYS* (Wiley 2025), §2.4 — humidity is the
mass of water vapor per unit mass of *dry* air (book Eq. 2.5), relative
humidity the ratio of the water partial pressure to its saturation value
(Eq. 2.6).

Backed by CoolProp's ``HAPropsSI`` (the ASHRAE RP-1485 real-gas humid-air
model — psychrometric-chart quality), not by the flowsheet property package:
psychrometrics is a property *analysis*, so it lives here next to the other
analysis tools.
"""
from __future__ import annotations

P_ATM = 101_325.0


def humidity(
    T: float,
    P: float = P_ATM,
    *,
    rh: float | None = None,
    w: float | None = None,
    t_wb: float | None = None,
    t_dp: float | None = None,
) -> dict[str, float]:
    """Resolve the full moist-air state at dry-bulb ``T`` (K) and ``P`` (Pa)
    from exactly one humidity spec.

    Parameters
    ----------
    rh : relative humidity as a *fraction* (0–1, not percent).
    w : humidity ratio, kg water / kg dry air (the book's H, Eq. 2.5).
    t_wb : wet-bulb temperature, K.
    t_dp : dew-point temperature, K.

    Returns
    -------
    dict with keys
      * ``T``, ``P`` — the inputs, K and Pa
      * ``rh`` — relative humidity (0–1)
      * ``w`` — humidity ratio, kg water / kg dry air
      * ``t_dp`` — dew point, K
      * ``t_wb`` — wet-bulb temperature, K
      * ``h`` — moist-air specific enthalpy, J / kg *dry* air (ASHRAE
        reference: h = 0 for dry air at 0 °C)
      * ``v`` — moist volume, m^3 / kg dry air
      * ``p_w`` — water vapor partial pressure, Pa
    """
    from CoolProp.HumidAirProp import HAPropsSI

    specs = {"R": rh, "W": w, "Twb": t_wb, "Tdp": t_dp}
    given = {k: v for k, v in specs.items() if v is not None}
    if len(given) != 1:
        raise ValueError(
            f"humidity() needs exactly one of rh / w / t_wb / t_dp; got "
            f"{len(given)}: {sorted(given) or 'none'}"
        )
    ((key, value),) = given.items()
    if key == "R" and not 0.0 <= float(value) <= 1.0:
        raise ValueError(
            f"rh={value} must be a fraction in [0, 1] (60% relative humidity "
            f"is rh=0.6, not 60)"
        )
    if key == "W" and float(value) < 0.0:
        raise ValueError(f"w={value} (kg water / kg dry air) must be >= 0")

    args = ("T", float(T), "P", float(P), key, float(value))
    return {
        "T": float(T),
        "P": float(P),
        "rh": float(HAPropsSI("R", *args)),
        "w": float(HAPropsSI("W", *args)),
        "t_dp": float(HAPropsSI("Tdp", *args)),
        "t_wb": float(HAPropsSI("Twb", *args)),
        "h": float(HAPropsSI("Hda", *args)),
        "v": float(HAPropsSI("Vda", *args)),
        "p_w": float(HAPropsSI("P_w", *args)),
    }
