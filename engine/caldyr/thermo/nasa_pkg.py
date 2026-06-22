"""Ideal-gas property package over Cantera's NASA polynomials (``nasa_gas.yaml``).

Selected as ``"nasa:gas"`` (alias ``"nasa:claus"``). This is the enabling thermo
for Claus sulfur recovery and other high-temperature combustion/equilibrium
flowsheets whose species the cubic-EOS packages cannot carry. In particular the
sulfur dimer **S2** — the dominant elemental-sulfur allotrope above ~800 K, and
unavoidable in a Claus thermal reactor — has *no* critical constants in
``chemicals``, so ``make_package("thermo:PR", [... "S2" ...])`` crashes building
the EOS. The NASA polynomials need only the species' standard thermochemistry,
which Cantera bundles for the whole Claus slate (H2S, SO2, SO3, S2, S8, S, SO,
COS, CS2, plus H2O/O2/N2/CO2/CO/H2/CH4).

**Scope — ideal gas, single phase.** Every state this package resolves is a
vapour: it provides ideal-gas enthalpy/entropy/volume and PT/PH/PS flashes, which
is all a combustion/equilibrium reactor train needs. It deliberately does *not*
do VLE — there is no liquid phase in ``nasa_gas.yaml``. Sulfur condensation (the
one condensed phase a Claus plant forms) is handled explicitly by the
:class:`~caldyr.unitops.sulfur_condenser.SulfurCondenser` using
:mod:`caldyr.thermo.sulfur`, exactly as the :class:`GibbsReactor` keeps its
equilibrium composition in Cantera and its state in the package. Bubble/dew and
three-phase flashes raise a typed ``NotImplementedError`` pointing at the cubic
packages.

**Enthalpy basis.** Cantera's NASA polynomials are formation-inclusive (referred
to the elements in their standard states at 298.15 K), the *same* convention as
the engine's ``thermo:*`` packages — verified: H2S comes out at −20.5 kJ/mol and
H2O at −241.8 kJ/mol, their formation enthalpies. So reactor heats of reaction
and adiabatic flame temperatures fall out of plain enthalpy balances, and every
stream in a NASA-package flowsheet shares one absolute basis. (A NASA-package
flowsheet should not be mixed with ``thermo:*`` streams; within it everything is
self-consistent.)
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

from .base import PhaseResult, ThreePhaseResult

_R = 8.314462618          # J/mol/K
_KMOL = 1000.0            # mol per kmol (Cantera molar properties are per kmol)

# CAS -> nasa_gas.yaml species name, for the stable species a combustion / Claus
# flowsheet would carry. Caldyr ids resolve to CAS via `chemicals`, so any
# synonym the database knows lands on the right Cantera species. Reactive
# radicals (S, SO, OH, O, H) are deliberately not exposed as flowsheet
# components; a reactor's equilibrium runs over exactly the mapped components, so
# fidelity is the modeller's choice of component list (as with GibbsReactor).
_CAS_TO_NASA: dict[str, str] = {
    "1333-74-0": "H2",        # hydrogen
    "7727-37-9": "N2",        # nitrogen
    "7782-44-7": "O2",        # oxygen
    "7732-18-5": "H2O",       # water
    "124-38-9": "CO2",        # carbon dioxide
    "630-08-0": "CO",         # carbon monoxide
    "74-82-8": "CH4",         # methane
    "7783-06-4": "H2S",       # hydrogen sulfide
    "7446-09-5": "SO2",       # sulfur dioxide
    "7446-11-9": "SO3",       # sulfur trioxide
    "23550-45-0": "S2",       # disulfur (high-T allotrope)
    "10544-50-0": "S8",       # cyclo-octasulfur (low-T allotrope)
    "463-58-1": "COS",        # carbonyl sulfide
    "75-15-0": "CS2",         # carbon disulfide
}


class NasaSpeciesError(ValueError):
    """A flowsheet component could not be mapped onto the bundled Cantera
    ``nasa_gas.yaml`` species set supported by the NASA ideal-gas package."""


@lru_cache(maxsize=512)
def _nasa_name(component_id: str) -> str | None:
    """nasa_gas.yaml species name for a caldyr component id, or None if the id
    does not resolve / is outside the mapped species set."""
    from chemicals.identifiers import CAS_from_any

    try:
        cas = CAS_from_any(component_id)
    except ValueError:
        return None
    return _CAS_TO_NASA.get(cas)


@lru_cache(maxsize=32)
def _build_solution(species_names: tuple[str, ...]) -> Any:
    """Ideal-gas Cantera Solution restricted to ``species_names`` (a subset of
    nasa_gas.yaml). Cached and reused; callers always set the full state before
    reading a property."""
    import cantera as ct

    by_name = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    return ct.Solution(thermo="ideal-gas",
                       species=[by_name[n] for n in species_names])


class NasaGasPackage:
    """Ideal-gas :class:`~caldyr.thermo.base.PropertyPackage` over Cantera's
    NASA polynomials. Built by ``make_package("nasa:gas", components)``."""

    def __init__(self, components: list[str], method: str = "gas") -> None:
        self.components = list(components)
        self.method = method
        self._name: dict[str, str] = {}          # component id -> nasa species
        unmapped: list[str] = []
        for cid in components:
            name = _nasa_name(cid)
            if name is None:
                unmapped.append(cid)
            else:
                self._name[cid] = name
        if unmapped:
            raise NasaSpeciesError(
                f"the NASA ideal-gas package ('nasa:{method}') cannot carry "
                f"component(s) {unmapped!r}; supported species are "
                f"{sorted(set(_CAS_TO_NASA.values()))}. (Liquid elemental sulfur "
                f"is represented as S8 and handled by the SulfurCondenser, not as "
                f"a separate component.)"
            )
        # Stable Cantera species order for this component set.
        self._species = tuple(self._name[c] for c in self.components)

    # -- internals -----------------------------------------------------------
    def _solution(self) -> Any:
        return _build_solution(self._species)

    def _X(self, z: dict[str, float]) -> dict[str, float]:
        """Cantera mole-fraction dict (by species name) from a caldyr
        composition (by component id). Normalised; missing -> 0."""
        x = {self._name[c]: max(z.get(c, 0.0), 0.0) for c in self.components}
        total = sum(x.values())
        if total <= 0.0:
            raise ValueError(f"composition sums to {total}; expected > 0")
        return {k: v / total for k, v in x.items()}

    def _set_TPX(self, T: float, P: float, z: dict[str, float]) -> Any:
        gas = self._solution()
        gas.TPX = T, P, self._X(z)
        return gas

    def _vapor_result(self, T: float, P: float, H: float,
                      z: dict[str, float]) -> PhaseResult:
        comp = {c: z.get(c, 0.0) for c in self.components}
        return PhaseResult(T=T, P=P, H=H, phase="vapor", vapor_fraction=1.0,
                           x=None, y=comp, H_liquid=None, H_vapor=H)

    # -- PropertyPackage protocol --------------------------------------------
    def enthalpy(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar enthalpy (J/mol), formation-inclusive NASA basis."""
        return self._set_TPX(T, P, z).enthalpy_mole / _KMOL

    def entropy(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar entropy (J/mol/K), NASA basis."""
        return self._set_TPX(T, P, z).entropy_mole / _KMOL

    def volume(self, T: float, P: float, z: dict[str, float]) -> float:
        """Molar volume (m^3/mol). Ideal gas: v = RT/P."""
        # _set_TPX validates the composition; ideal-gas volume is exact.
        self._X(z)
        return _R * T / P

    def flash_pt(self, T: float, P: float, z: dict[str, float]) -> PhaseResult:
        H = self._set_TPX(T, P, z).enthalpy_mole / _KMOL
        return self._vapor_result(T, P, H, z)

    def flash_ph(self, P: float, H: float, z: dict[str, float]) -> PhaseResult:
        """Find T such that the ideal-gas enthalpy equals ``H`` (J/mol) at fixed
        composition, via Cantera's HP state setter."""
        gas = self._solution()
        X = self._X(z)
        gas.TPX = 298.15, P, X                       # seed
        gas.HPX = H * _KMOL / gas.mean_molecular_weight, P, X
        return self._vapor_result(gas.T, P, H, z)

    def flash_ps(self, P: float, S: float, z: dict[str, float]) -> PhaseResult:
        gas = self._solution()
        X = self._X(z)
        gas.TPX = 298.15, P, X
        gas.SPX = S * _KMOL / gas.mean_molecular_weight, P, X
        H = gas.enthalpy_mole / _KMOL
        return self._vapor_result(gas.T, P, H, z)

    # -- forced-phase properties (vapour only; this is an ideal-gas package) --
    def k_values(self, T: float, P: float, x: dict[str, float],
                 y: dict[str, float]) -> dict[str, float]:
        raise NotImplementedError(
            "nasa:gas is an ideal-gas (single-phase) package; it has no VLE "
            "K-values. Use a 'thermo:*' package for distillation/absorption."
        )

    def enthalpy_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        return self.enthalpy(T, P, y)

    def enthalpy_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        raise NotImplementedError(
            "nasa:gas has no liquid phase (ideal-gas, NASA polynomials). Liquid "
            "elemental sulfur is modelled by the SulfurCondenser "
            "(caldyr.thermo.sulfur)."
        )

    def volume_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        return self.volume(T, P, y)

    def volume_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        raise NotImplementedError("nasa:gas has no liquid phase.")

    # -- unsupported (no condensed phase / no VLE) ----------------------------
    def bubble_dew(self, P: float, z: dict[str, float]) -> tuple[float, float]:
        raise NotImplementedError(
            "nasa:gas is an ideal-gas package with no condensed phase; "
            "bubble/dew points are undefined. Use a 'thermo:*' package."
        )

    def bubble_point(self, P: float, x: dict[str, float]) -> PhaseResult:
        raise NotImplementedError(
            "nasa:gas is an ideal-gas package with no condensed phase; "
            "bubble_point is undefined. Use a 'thermo:*' package."
        )

    def flash_pt_3p(self, T: float, P: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "nasa:gas does not support three-phase (VLLE) flashes; "
            "select 'thermo:PR' or 'thermo:SRK'."
        )

    def flash_ph_3p(self, P: float, H: float, z: dict[str, float]) -> ThreePhaseResult:
        raise NotImplementedError(
            "nasa:gas does not support three-phase (VLLE) flashes; "
            "select 'thermo:PR' or 'thermo:SRK'."
        )

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "NasaGasPackage":
        backend, _, method = spec.partition(":")
        if backend != "nasa":
            raise ValueError(f"NasaGasPackage builds only 'nasa:*' (got {spec!r})")
        return cls(components, method or "gas")
