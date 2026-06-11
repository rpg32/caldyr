"""Cubic-EOS property package wrapping `thermo` (Caleb Bell, MIT).

Backs streams with a cubic-EOS VLE flash (Peng-Robinson by default, SRK
optional). Best for non-polar / lightly-polar systems (hydrocarbons, light
gases). For strongly non-ideal polar mixtures (alcohols + water, with
azeotropes) prefer the activity-coefficient package — see
:mod:`caldyr.thermo.activity_pkg`.

Enthalpies share `thermo`'s internally consistent reference state, so energy
balances close to machine precision across the flowsheet.

Petroleum pseudo-components (assay boiling-point cuts — see
:mod:`caldyr.assay`) are supported: any component id registered in the
pseudo-component registry (:mod:`caldyr.core.components_db`) gets its
Tc/Pc/omega/MW/Tb and ideal-gas Cp from the registry instead of the databank,
merged into one `ChemicalConstantsPackage` with the databank species. Their
formation enthalpy/Gibbs default to 0 (documented on
:class:`caldyr.core.component.Component`): energy balances stay exact because
the formation offset cancels at conserved composition, but reactions on
pseudo-components are unsupported.
"""
from __future__ import annotations

from functools import lru_cache

from ..core.components_db import pseudo_signature
from ._flasher import FlasherPackage, formation_props

# Ideal-gas Cp polynomial validity span given to thermo for pseudo-components
# (thermo extrapolates linearly outside). The Kesler-Lee Cp correlation behind
# the coefficients is fitted to ~255-922 K; the wider span avoids hard edges
# mid-flash while staying a documented extrapolation.
_PSEUDO_CP_T_RANGE = (200.0, 1500.0)


def _build_constants(components: tuple[str, ...], pseudo_sig: tuple):
    """Build ``(ChemicalConstantsPackage, PropertyCorrelationsPackage)`` for an
    ordered component list that may mix databank species with registered
    pseudo-components (``pseudo_sig`` as built by
    :func:`caldyr.core.components_db.pseudo_signature`).

    Databank species keep their full `from_IDs` constants and Cp correlations;
    pseudo entries contribute user/assay constants, a polynomial ideal-gas Cp,
    and Hf = Gf = Sf = 0 (no formation data exists for a lumped cut — and
    setting 0 instead of None deliberately suppresses the "no formation
    enthalpy" warning that real species get: for pseudos this is by design,
    since reactions on them are unsupported, not missing data).
    """
    from thermo import (
        ChemicalConstantsPackage,
        HeatCapacityGas,
        PropertyCorrelationsPackage,
    )

    pseudo = {cid: dict(items) for cid, items in pseudo_sig}
    if not pseudo:
        return ChemicalConstantsPackage.from_IDs(list(components))

    real = [c for c in components if c not in pseudo]
    rcon, rprops = (None, None)
    if real:
        rcon, rprops = ChemicalConstantsPackage.from_IDs(real)
    ridx = {c: i for i, c in enumerate(real)}

    def real_attr(obj, attr: str, cid: str):
        assert obj is not None     # cid not in pseudo => real list was non-empty
        return getattr(obj, attr)[ridx[cid]]

    def pick(attr: str, cid: str, pkey: str, pdefault=None):
        if cid in pseudo:
            value = pseudo[cid].get(pkey, pdefault)
            return None if value is None else float(value)
        return real_attr(rcon, attr, cid)

    names = list(components)
    cass = [None if c in pseudo else real_attr(rcon, "CASs", c) for c in components]
    mws = [float(pseudo[c]["MW"]) * 1000.0 if c in pseudo       # registry kg/mol
           else real_attr(rcon, "MWs", c) for c in components]  # thermo wants g/mol
    tbs = [pick("Tbs", c, "Tb") for c in components]
    tcs = [pick("Tcs", c, "Tc") for c in components]
    pcs = [pick("Pcs", c, "Pc") for c in components]
    omegas = [pick("omegas", c, "omega") for c in components]
    hfgs = [pick("Hfgs", c, "Hf", 0.0) for c in components]
    gfgs = [pick("Gfgs", c, "Gf", 0.0) for c in components]
    # Sf = (Hf - Gf)/T298 keeps the three thermodynamically consistent.
    sfgs = [
        (h - g) / 298.15 if (c in pseudo and h is not None and g is not None)
        else (None if c in pseudo else real_attr(rcon, "Sfgs", c))
        for c, h, g in zip(components, hfgs, gfgs)
    ]

    cps = []
    tmin, tmax = _PSEUDO_CP_T_RANGE
    for c in components:
        if c not in pseudo:
            cps.append(real_attr(rprops, "HeatCapacityGases", c))
            continue
        coeffs = pseudo[c].get("Cp_ig")
        if coeffs is None:
            raise ValueError(
                f"pseudo-component {c!r} has no ideal-gas Cp coefficients "
                f"('Cp_ig'); characterize it via caldyr.assay (which estimates "
                f"them with the Kesler-Lee correlation) or supply "
                f"Cp_ig=[a0, a1, a2] (J/mol/K, ascending powers of T)"
            )
        descending = [float(a) for a in coeffs][::-1]   # thermo wants highest first
        cps.append(HeatCapacityGas(poly_fit=(tmin, tmax, descending)))

    constants = ChemicalConstantsPackage(
        names=names, CASs=cass, MWs=mws, Tbs=tbs, Tcs=tcs, Pcs=pcs,
        omegas=omegas, Hfgs=hfgs, Gfgs=gfgs, Sfgs=sfgs,
    )
    props = PropertyCorrelationsPackage(
        constants=constants, HeatCapacityGases=cps, skip_missing=True,
    )
    return constants, props


@lru_cache(maxsize=64)
def _build_flasher(components: tuple[str, ...], method: str, pseudo_sig: tuple = ()):
    """Build (and cache) a cubic-EOS thermo flasher for an ordered component
    tuple. Returns ``(flasher, formation_enthalpies, formation_gibbs)``.

    `ChemicalConstantsPackage.from_IDs` hits a database and is slow, so the
    result is memoized per (components, method, pseudo constants) — flashers
    are pure and reusable, and keying on the pseudo constants themselves means
    a re-characterized assay (same ids, new constants) never reuses a stale
    flasher. A single-component system uses `FlashPureVLS` (the multicomponent
    `FlashVL` divides by N-1 internally and breaks for pure fluids); both share
    the same cubic-EOS enthalpy reference, so mixed/pure results stay
    consistent.
    """
    from thermo import CEOSGas, CEOSLiquid, FlashPureVLS, FlashVL, PRMIX, SRKMIX

    eos = {"PR": PRMIX, "SRK": SRKMIX}[method]
    constants, props = _build_constants(components, pseudo_sig)
    eos_kwargs = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas = CEOSGas(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    liq = CEOSLiquid(eos, eos_kwargs, HeatCapacityGases=props.HeatCapacityGases)
    hf, gf = formation_props(components, constants.Hfgs, constants.Gfgs)
    if len(components) == 1:
        return FlashPureVLS(constants, props, gas=gas, liquids=[liq], solids=[]), hf, gf
    return FlashVL(constants, props, liquid=liq, gas=gas), hf, gf


@lru_cache(maxsize=64)
def _build_flasher_3p(components: tuple[str, ...], method: str, pseudo_sig: tuple = ()):
    """Build (and cache) a *three-phase* cubic-EOS flasher (`thermo.FlashVLN`,
    a gas plus two trial liquids of the same EOS — the standard VLLE setup from
    the thermo docs). Returns ``(flasher, molar_masses_kg_per_mol)``; the molar
    masses let the caller order the liquids by mass density (light vs heavy).
    Requires two or more components (a pure fluid cannot split into two liquids).
    Pseudo-components are supported exactly as in :func:`_build_flasher`.
    """
    from thermo import CEOSGas, CEOSLiquid, FlashVLN, PRMIX, SRKMIX

    if len(components) < 2:
        raise ValueError(
            "a three-phase (VLLE) flash needs at least two components; "
            f"got {list(components)}"
        )
    eos = {"PR": PRMIX, "SRK": SRKMIX}[method]
    constants, props = _build_constants(components, pseudo_sig)
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
        flasher, hf, gf = _build_flasher(
            tuple(components), method, pseudo_signature(components)
        )
        self._init(components, flasher, hf, gf)

    def _build_3p(self):
        """Three-phase (VLLE) flasher: FlashVLN with two trial cubic-EOS liquids.
        Built lazily on the first flash_pt_3p/flash_ph_3p call and cached."""
        return _build_flasher_3p(
            tuple(self.components), self.method, pseudo_signature(self.components)
        )

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "ThermoPackage":
        """Build from a flowsheet `property_package` string like ``"thermo:PR"``."""
        backend, _, method = spec.partition(":")
        if backend != "thermo":
            raise ValueError(
                f"ThermoPackage cannot build backend {backend!r} (got {spec!r})"
            )
        return cls(components, method or "PR")
