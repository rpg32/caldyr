from .activity_pkg import ActivityPackage
from .amine_pkg import AmineAcidGasPackage
from .base import PhaseResult, PropertyPackage, ThreePhaseResult
from .coolprop_pkg import CoolPropWaterPackage
from .nasa_pkg import NasaGasPackage
from .thermo_pkg import ThermoPackage

# Which backend builds each `thermo:<method>` selector.
_CUBIC = {"PR", "SRK"}
_ACTIVITY = {"NRTL", "UNIFAC", "UNIFAC-LLE"}

#: Every selector :func:`make_package` accepts, with a one-line "when to use".
#: The single source of truth for what property packages exist — the AI tool
#: layer (and anything else listing packages) serves from here.
AVAILABLE_PACKAGES: list[dict[str, str]] = [
    {"id": "thermo:PR", "name": "Peng-Robinson (cubic EOS)",
     "use": "non-polar gases/hydrocarbons; supports petroleum assay pseudo-components"},
    {"id": "thermo:SRK", "name": "Soave-Redlich-Kwong (cubic EOS)",
     "use": "non-polar systems; supports petroleum assay pseudo-components"},
    {"id": "thermo:NRTL", "name": "NRTL gamma-phi",
     "use": "polar mixtures / azeotropes"},
    {"id": "thermo:UNIFAC", "name": "Modified UNIFAC (Dortmund) gamma-phi",
     "use": "predictive VLE+LLE / heteroazeotropes (3-phase distillation)"},
    {"id": "thermo:UNIFAC-LLE", "name": "UNIFAC-LLE gamma-phi",
     "use": "liquid-liquid extraction (Magnussen LLE table)"},
    {"id": "coolprop:Water", "name": "Steam tables (IAPWS-95 via CoolProp)",
     "use": "pure-water / steam flowsheets only"},
    {"id": "amine:DEA", "name": "Reactive acid-gas, DEA (modified Kent-Eisenberg)",
     "use": "CO2/H2S amine sweetening — absorbers, strippers, regenerators"},
    {"id": "amine:MDEA", "name": "Reactive acid-gas, MDEA (modified Kent-Eisenberg)",
     "use": "CO2/H2S amine sweetening — absorbers, strippers, regenerators"},
    {"id": "amine:MEA", "name": "Reactive acid-gas, MEA (modified Kent-Eisenberg)",
     "use": "CO2/H2S amine sweetening — absorbers, strippers, regenerators"},
    {"id": "nasa:gas", "name": "Ideal gas over NASA polynomials (Cantera data)",
     "use": "combustion/Claus species incl. sulfur allotropes S2/S8"},
]


def _validate_component_ids(components: list[str]) -> None:
    """Check every component id resolves in the chemicals database *before* the
    thermo backend hits a cryptic lookup failure deep inside its constants
    loader. A bad id raises ``UnknownComponentError`` naming it. Lookups are
    cached, so this costs nothing on repeat solves."""
    from ..core.components_db import UnknownComponentError, resolve_component

    bad: list[str] = []
    for cid in components:
        try:
            resolve_component(cid)
        except UnknownComponentError:
            bad.append(cid)
    if bad:
        raise UnknownComponentError(
            f"flowsheet component id(s) {bad!r} could not be resolved to chemical "
            f"species; use a name/CAS/formula the chemicals database recognizes "
            f"(see caldyr.core.components_db.COMMON_COMPONENTS for a curated list)"
        )


def _reject_pseudo_components(spec: str, components: list[str]) -> None:
    """Pseudo-components (assay cuts) carry only Tc/Pc/omega/MW/Cp constants —
    enough for a cubic EOS, but NRTL needs binary interaction parameters and
    the CoolProp backend needs a reference EOS, neither of which can exist for
    a lumped cut. Fail loudly with the supported alternative."""
    from ..core.components_db import is_pseudo_component

    pseudos = [c for c in components if is_pseudo_component(c)]
    if pseudos:
        raise ValueError(
            f"property package {spec!r} does not support petroleum "
            f"pseudo-components (got {pseudos}); use a cubic EOS package — "
            f"'thermo:PR' or 'thermo:SRK'"
        )


def make_package(spec: str, components: list[str]) -> PropertyPackage:
    """Build the property package selected by a flowsheet's ``property_package``
    string for a given ordered component list.

    Supported selectors:
      * ``thermo:PR`` / ``thermo:SRK`` — cubic EOS (non-polar systems). Also the
        only packages supporting petroleum pseudo-components (assay cuts; see
        :mod:`caldyr.assay`).
      * ``thermo:NRTL`` — activity-coefficient liquid (polar systems, azeotropes).
      * ``thermo:UNIFAC`` — predictive Modified UNIFAC (Dortmund) gamma-phi liquid
        with automatic group assignment; the simultaneous VLE+LLE model for
        heteroazeotropic (3-phase) distillation. ``thermo:UNIFAC-LLE`` is the
        Magnussen LLE-fit variant.
      * ``coolprop:Water`` — pure-water steam tables (CoolProp IAPWS-95);
        single-component water flowsheets only.
      * ``amine:DEA`` / ``amine:MDEA`` / ``amine:MEA`` — reactive acid-gas
        (modified Kent-Eisenberg) package for amine gas sweetening (CO2/H2S
        absorbers, strippers and regenerators); see
        :mod:`caldyr.thermo.amine_pkg`.
      * ``nasa:gas`` (alias ``nasa:claus``) — ideal-gas package over Cantera's
        NASA polynomials, carrying combustion/Claus species (incl. the sulfur
        allotropes S2/S8 the cubic EOS cannot); see
        :mod:`caldyr.thermo.nasa_pkg`.
    """
    backend, _, method = spec.partition(":")
    if backend == "nasa":
        _reject_pseudo_components(spec, components)
        _validate_component_ids(components)
        return NasaGasPackage(components, method or "gas")
    if backend == "amine":
        _reject_pseudo_components(spec, components)
        _validate_component_ids(components)
        return AmineAcidGasPackage(components, method or "DEA")
    if backend == "coolprop":
        if (method or "").lower() != "water":
            raise ValueError(
                f"unknown coolprop method {method!r} in {spec!r}; the only "
                f"supported coolprop selector is 'coolprop:Water'"
            )
        _reject_pseudo_components(spec, components)
        _validate_component_ids(components)
        return CoolPropWaterPackage(components)
    if backend != "thermo":
        raise ValueError(f"unknown property package backend in {spec!r}")
    _validate_component_ids(components)
    method = (method or "PR").upper()
    if method in _CUBIC:
        return ThermoPackage(components, method)
    if method in _ACTIVITY:
        _reject_pseudo_components(spec, components)
        return ActivityPackage(components, method)
    raise ValueError(
        f"unknown property method {method!r} in {spec!r}; "
        f"expected one of {sorted(_CUBIC | _ACTIVITY)}"
    )


__all__ = [
    "AVAILABLE_PACKAGES",
    "PhaseResult",
    "ThreePhaseResult",
    "PropertyPackage",
    "ThermoPackage",
    "ActivityPackage",
    "CoolPropWaterPackage",
    "AmineAcidGasPackage",
    "NasaGasPackage",
    "make_package",
]
