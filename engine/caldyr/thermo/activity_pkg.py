"""Activity-coefficient (gamma-phi) property package wrapping `thermo`.

The liquid is modelled with an excess-Gibbs activity-coefficient model (NRTL)
and the vapor as an ideal gas — the standard low-pressure gamma-phi approach.
Unlike a cubic EOS, this captures strongly non-ideal polar liquids and the
azeotropes they form (e.g. ethanol + water). Binary interaction parameters come
from the ChemSep dataset bundled with `thermo` (via its IPDB); pairs without
data fall back to an ideal-solution liquid, and that fallback is warned about
rather than applied silently.

Validated against the ethanol/water minimum-boiling azeotrope (~89 mol% ethanol,
78.2 C at 1 atm — DePriester / Gmehling DECHEMA), which a cubic EOS cannot
represent.
"""
from __future__ import annotations

import warnings
from functools import lru_cache

from ._flasher import FlasherPackage, formation_props

# ChemSep NRTL stores tau as b_ij / T and the non-randomness as alpha_ij.
_NRTL_TABLE = "ChemSep NRTL"


def _has_offdiagonal(matrix) -> bool:
    return any(abs(matrix[i][j]) > 0.0
               for i in range(len(matrix)) for j in range(len(matrix)) if i != j)


@lru_cache(maxsize=64)
def _build_flasher(components: tuple[str, ...], model: str):
    """Build (and cache) a gamma-phi thermo flasher for an ordered component
    tuple. Pure-component systems use `FlashPureVLS` (the activity model is
    irrelevant for one component: gamma = 1)."""
    from thermo import (
        ChemicalConstantsPackage,
        FlashPureVLS,
        FlashVL,
        GibbsExcessLiquid,
        IdealGas,
        NRTL,
    )
    from thermo.interaction_parameters import IPDB

    if model != "NRTL":
        raise ValueError(f"unsupported activity model {model!r}; expected 'NRTL'")

    constants, props = ChemicalConstantsPackage.from_IDs(list(components))
    gas = IdealGas(HeatCapacityGases=props.HeatCapacityGases)
    hf, gf = formation_props(components, constants.Hfgs, constants.Gfgs)
    liquid_kwargs = dict(
        VaporPressures=props.VaporPressures,
        HeatCapacityGases=props.HeatCapacityGases,
        VolumeLiquids=props.VolumeLiquids,
        equilibrium_basis="Psat",
        caloric_basis="Psat",
    )

    if len(components) == 1:
        liquid = GibbsExcessLiquid(**liquid_kwargs)
        return FlashPureVLS(constants, props, gas=gas, liquids=[liquid], solids=[]), hf, gf

    cas = constants.CASs
    tau_bs = IPDB.get_ip_asymmetric_matrix(_NRTL_TABLE, cas, "bij")
    alpha_cs = IPDB.get_ip_asymmetric_matrix(_NRTL_TABLE, cas, "alphaij")
    if not _has_offdiagonal(tau_bs):
        warnings.warn(
            f"no ChemSep NRTL parameters for {list(components)}; the liquid "
            f"falls back to an ideal solution (no activity correction). Results "
            f"for polar mixtures may be inaccurate.",
            stacklevel=2,
        )
    n = len(components)
    ge_model = NRTL(T=298.15, xs=[1.0 / n] * n, tau_bs=tau_bs, alpha_cs=alpha_cs)
    liquid = GibbsExcessLiquid(GibbsExcessModel=ge_model, **liquid_kwargs)
    return FlashVL(constants, props, liquid=liquid, gas=gas), hf, gf


class ActivityPackage(FlasherPackage):
    """NRTL (gamma) + ideal-gas (phi) property package over a fixed component
    list. Selected by ``property_package`` strings like ``"thermo:NRTL"``."""

    SUPPORTED = ("NRTL",)

    def __init__(self, components: list[str], model: str = "NRTL") -> None:
        model = model.upper()
        if model not in self.SUPPORTED:
            raise ValueError(
                f"unsupported activity model {model!r}; expected one of {self.SUPPORTED}"
            )
        self.model = model
        flasher, hf, gf = _build_flasher(tuple(components), model)
        self._init(components, flasher, hf, gf)

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "ActivityPackage":
        """Build from a flowsheet `property_package` string like ``"thermo:NRTL"``."""
        backend, _, model = spec.partition(":")
        if backend != "thermo":
            raise ValueError(
                f"ActivityPackage cannot build backend {backend!r} (got {spec!r})"
            )
        return cls(components, model or "NRTL")
