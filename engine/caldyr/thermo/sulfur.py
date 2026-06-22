"""Liquid-sulfur condensed-phase properties for Claus sulfur-recovery condensers.

The gas-phase thermo of a Claus plant (H2S/SO2/S2/S8/H2O/N2/...) lives in the
NASA ideal-gas package (:mod:`caldyr.thermo.nasa_pkg`), which carries the sulfur
*vapour* allotropes. What that package cannot do is condense elemental sulfur:
``nasa_gas.yaml`` has no liquid phase, and caldyr's cubic-EOS package cannot even
build (``chemicals`` has no critical constants for the S2 dimer). A
:class:`~caldyr.unitops.sulfur_condenser.SulfurCondenser` therefore needs an
independent liquid-sulfur model — that is what this module provides.

Both correlations come from the *same* ``thermo``/``chemicals`` stack the rest of
the engine uses (Caleb Bell, MIT-licensed), so they sit on a consistent data
basis:

* **Vapour pressure of liquid sulfur** — ``thermo.VaporPressure`` for elemental
  sulfur (CAS 7704-34-9; Poling Antoine constants). It reproduces the normal
  boiling point to <1 % (1.02e5 Pa at 717.8 K) and gives the very low vapour
  pressures (~10-60 Pa) at Claus condenser temperatures (~400-435 K) that set the
  residual elemental-sulfur loss to the tail gas.
* **Heat of vaporisation of sulfur** — ``thermo.EnthalpyVaporization`` for the
  same species (~10.4 kJ per mol of S atoms near the condensers), the latent heat
  released as sulfur condenses.

Both are reported by ``thermo`` on a **per-mole-of-S-atom** basis (e.g. the
liquid heat capacity is ~32 J/mol-S/K ≈ 1 J/g/K), which is how this module
exposes them. Callers working in S8 multiply by 8.
"""
from __future__ import annotations

from functools import lru_cache

_SULFUR_CAS = "7704-34-9"        # elemental sulfur (the liquid/condensed phase)


@lru_cache(maxsize=1)
def _vapor_pressure_model():
    from thermo import VaporPressure

    return VaporPressure(CASRN=_SULFUR_CAS)


@lru_cache(maxsize=1)
def _hvap_model():
    from thermo import EnthalpyVaporization

    return EnthalpyVaporization(CASRN=_SULFUR_CAS)


def liquid_sulfur_psat(T: float) -> float:
    """Vapour pressure of liquid sulfur (Pa) at temperature ``T`` (K).

    This is the *total* equilibrium pressure of sulfur vapour over the liquid;
    at Claus condenser temperatures the vapour is overwhelmingly S8, so a
    condenser treats this as the S8 partial pressure that saturates the exit gas
    (the residual elemental-sulfur loss). Validated at the normal boiling point:
    ``liquid_sulfur_psat(717.8) ≈ 1.0e5 Pa``.
    """
    p = _vapor_pressure_model().T_dependent_property(T)
    if p is None:
        raise ValueError(
            f"liquid-sulfur vapour pressure is undefined at T={T:.2f} K "
            f"(thermo's Antoine fit is out of range)"
        )
    return float(p)


def sulfur_hvap_per_atom(T: float) -> float:
    """Heat of vaporisation of sulfur (J per mol of S atoms) at ``T`` (K).

    Multiply by 8 for a per-mol-S8 latent heat. Near the condensers (~430 K)
    this is ~10.4 kJ/mol-S (≈ 83 kJ/mol-S8)."""
    h = _hvap_model().T_dependent_property(T)
    if h is None:
        raise ValueError(
            f"liquid-sulfur heat of vaporisation is undefined at T={T:.2f} K"
        )
    return float(h)
