from .activity_pkg import ActivityPackage
from .base import PhaseResult, PropertyPackage, ThreePhaseResult
from .thermo_pkg import ThermoPackage

# Which backend builds each `thermo:<method>` selector.
_CUBIC = {"PR", "SRK"}
_ACTIVITY = {"NRTL"}


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


def make_package(spec: str, components: list[str]) -> PropertyPackage:
    """Build the property package selected by a flowsheet's ``property_package``
    string for a given ordered component list.

    Supported selectors:
      * ``thermo:PR`` / ``thermo:SRK`` — cubic EOS (non-polar systems).
      * ``thermo:NRTL`` — activity-coefficient liquid (polar systems, azeotropes).
    """
    backend, _, method = spec.partition(":")
    if backend != "thermo":
        raise ValueError(f"unknown property package backend in {spec!r}")
    _validate_component_ids(components)
    method = (method or "PR").upper()
    if method in _CUBIC:
        return ThermoPackage(components, method)
    if method in _ACTIVITY:
        return ActivityPackage(components, method)
    raise ValueError(
        f"unknown property method {method!r} in {spec!r}; "
        f"expected one of {sorted(_CUBIC | _ACTIVITY)}"
    )


__all__ = [
    "PhaseResult",
    "ThreePhaseResult",
    "PropertyPackage",
    "ThermoPackage",
    "ActivityPackage",
    "make_package",
]
