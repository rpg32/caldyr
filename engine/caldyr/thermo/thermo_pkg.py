"""Cubic-EOS property package wrapping `thermo` (Caleb Bell, MIT).

Backs streams with a cubic-EOS VLE flash (Peng-Robinson by default, SRK
optional). Best for non-polar / lightly-polar systems (hydrocarbons, light
gases). For strongly non-ideal polar mixtures (alcohols + water, with
azeotropes) prefer the activity-coefficient package — see
:mod:`caldyr.thermo.activity_pkg`.

Enthalpies share `thermo`'s internally consistent reference state, so energy
balances close to machine precision across the flowsheet.
"""
from __future__ import annotations

from functools import lru_cache

from ._flasher import FlasherPackage, formation_props


@lru_cache(maxsize=64)
def _build_flasher(components: tuple[str, ...], method: str):
    """Build (and cache) a cubic-EOS thermo flasher for an ordered component
    tuple. Returns ``(flasher, formation_enthalpies, formation_gibbs)``.

    `ChemicalConstantsPackage.from_IDs` hits a database and is slow, so the
    result is memoized per (components, method) — flashers are pure and reusable.
    A single-component system uses `FlashPureVLS` (the multicomponent `FlashVL`
    divides by N-1 internally and breaks for pure fluids); both share the same
    cubic-EOS enthalpy reference, so mixed/pure results stay consistent.
    """
    from thermo import (
        CEOSGas,
        CEOSLiquid,
        ChemicalConstantsPackage,
        FlashPureVLS,
        FlashVL,
        PRMIX,
        SRKMIX,
    )

    eos = {"PR": PRMIX, "SRK": SRKMIX}[method]
    constants, props = ChemicalConstantsPackage.from_IDs(list(components))
    eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas = CEOSGas(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    liq = CEOSLiquid(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    hf, gf = formation_props(components, constants.Hfgs, constants.Gfgs)
    if len(components) == 1:
        return FlashPureVLS(constants, props, gas=gas, liquids=[liq], solids=[]), hf, gf
    return FlashVL(constants, props, liquid=liq, gas=gas), hf, gf


@lru_cache(maxsize=64)
def _build_flasher_3p(components: tuple[str, ...], method: str):
    """Build (and cache) a *three-phase* cubic-EOS flasher (`thermo.FlashVLN`,
    a gas plus two trial liquids of the same EOS — the standard VLLE setup from
    the thermo docs). Returns ``(flasher, molar_masses_kg_per_mol)``; the molar
    masses let the caller order the liquids by mass density (light vs heavy).
    Requires two or more components (a pure fluid cannot split into two liquids).
    """
    from thermo import CEOSGas, CEOSLiquid, ChemicalConstantsPackage, FlashVLN, PRMIX, SRKMIX

    if len(components) < 2:
        raise ValueError(
            "a three-phase (VLLE) flash needs at least two components; "
            f"got {list(components)}"
        )
    eos = {"PR": PRMIX, "SRK": SRKMIX}[method]
    constants, props = ChemicalConstantsPackage.from_IDs(list(components))
    eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas = CEOSGas(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    liq = CEOSLiquid(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    flasher = FlashVLN(constants, props, liquids=[liq, liq], gas=gas)
    mws = [float(mw) / 1000.0 for mw in constants.MWs]      # g/mol -> kg/mol
    return flasher, mws


class ThermoPackage(FlasherPackage):
    """PR/SRK cubic-EOS property package over a fixed, ordered component list."""

    SUPPORTED = ("PR", "SRK")

    def __init__(self, components: list[str], method: str = "PR") -> None:
        method = method.upper()
        if method not in self.SUPPORTED:
            raise ValueError(
                f"unsupported EOS method {method!r}; expected one of {self.SUPPORTED}"
            )
        self.method = method
        flasher, hf, gf = _build_flasher(tuple(components), method)
        self._init(components, flasher, hf, gf)

    def _build_3p(self):
        """Three-phase (VLLE) flasher: FlashVLN with two trial cubic-EOS liquids.
        Built lazily on the first flash_pt_3p/flash_ph_3p call and cached."""
        return _build_flasher_3p(tuple(self.components), self.method)

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "ThermoPackage":
        """Build from a flowsheet `property_package` string like ``"thermo:PR"``."""
        backend, _, method = spec.partition(":")
        if backend != "thermo":
            raise ValueError(
                f"ThermoPackage cannot build backend {backend!r} (got {spec!r})"
            )
        return cls(components, method or "PR")
