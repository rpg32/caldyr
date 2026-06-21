"""Component identity resolution + a curated catalog of common species.

:func:`resolve_component` turns any identifier the `chemicals` library
(Caleb Bell, MIT) recognizes — a common name ("water"), CAS number
("7732-18-5"), formula ("H2O"), InChI, or SMILES — into a typed
:class:`~caldyr.core.component.Component` with its name/formula/CAS filled in.
An unknown identifier raises :class:`UnknownComponentError` naming the bad id
(no silent failures), which is also how flowsheet component validation surfaces
typos before the thermo layer hits a cryptic database miss.

:data:`COMMON_COMPONENTS` is a curated catalog of ~120 industrially common
species (light gases, C1-C12 alkanes/alkenes/aromatics, alcohols, ketones,
acids, esters, amines, chlorinated solvents, refrigerants). The id strings are
what a flowsheet uses; name/formula/CAS were resolved through `chemicals` 1.5.2
at author time and hardcoded so importing the catalog costs nothing at startup.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

from .component import Component


class UnknownComponentError(ValueError):
    """A component identifier could not be resolved to a chemical species."""


# -- petroleum pseudo-component registry --------------------------------------
# Pseudo-components (assay boiling-point cuts; see caldyr.assay) have no entry
# in the chemicals databank — their constants are supplied by the user or the
# assay characterization. The registry maps component id -> constants dict
# (keys documented on caldyr.core.component.Component.pseudo) and is consulted
# by resolve_component / molar_mass and by the thermo layer when building
# flashers. It is populated as a side effect of constructing a Component with
# a `pseudo` payload (including `.flow` round-trips through Component(**c)).
_PSEUDO_REGISTRY: dict[str, dict] = {}

_PSEUDO_REQUIRED_KEYS = ("MW", "Tb", "SG", "Tc", "Pc", "omega")


def register_pseudo_component(component_id: str, constants: dict) -> None:
    """Register (or overwrite) the constants of a pseudo-component.

    Validates the required keys (MW kg/mol, Tb K, SG, Tc K, Pc Pa, omega) are
    present and positive where physical, raising ``ValueError`` naming the
    offender — a typo here would otherwise surface as a cryptic flash failure.
    """
    missing = [k for k in _PSEUDO_REQUIRED_KEYS if k not in constants]
    if missing:
        raise ValueError(
            f"pseudo-component {component_id!r} is missing required constant(s) "
            f"{missing}; required: {list(_PSEUDO_REQUIRED_KEYS)} "
            f"(MW kg/mol, Tb K, SG, Tc K, Pc Pa, omega)"
        )
    for key in ("MW", "Tb", "SG", "Tc", "Pc"):
        if not float(constants[key]) > 0.0:
            raise ValueError(
                f"pseudo-component {component_id!r} has non-positive {key} = "
                f"{constants[key]!r}"
            )
    _PSEUDO_REGISTRY[component_id] = dict(constants)


def pseudo_constants(component_id: str) -> dict | None:
    """The registered constants dict of a pseudo-component, or None if
    ``component_id`` is not a registered pseudo-component."""
    found = _PSEUDO_REGISTRY.get(component_id)
    return dict(found) if found is not None else None


def is_pseudo_component(component_id: str) -> bool:
    return component_id in _PSEUDO_REGISTRY


def _hashable(value):
    return tuple(value) if isinstance(value, (list, tuple)) else float(value)


def pseudo_signature(components: list[str] | tuple[str, ...]) -> tuple:
    """A hashable signature of the registered pseudo constants among
    ``components`` — () when none are pseudo. Used as part of thermo flasher
    cache keys so re-characterizing an assay (same ids, new constants) can
    never reuse a stale flasher."""
    return tuple(
        (cid, tuple(sorted((k, _hashable(v)) for k, v in _PSEUDO_REGISTRY[cid].items())))
        for cid in components
        if cid in _PSEUDO_REGISTRY
    )


# Process abbreviations the chemicals index does not resolve on its own
# (notably "MEA"); mapped to CAS so the common gas-sweetening amine ids work.
_ABBREV_CAS = {
    "MEA": "141-43-5", "DEA": "111-42-2", "MDEA": "105-59-9",
}


@lru_cache(maxsize=1024)
def _search(identifier: str) -> Any:
    """Cached `chemicals` metadata lookup (the chemicals index search is not
    free, and the same ids are resolved repeatedly across solves)."""
    from chemicals.identifiers import search_chemical

    return search_chemical(_ABBREV_CAS.get(identifier.strip().upper(), identifier))


def resolve_component(identifier: str) -> Component:
    """Resolve ``identifier`` (name / CAS / formula / InChI / SMILES) to a
    :class:`Component` via the `chemicals` database.

    The returned component keeps ``identifier`` as its ``id`` (streams and
    property packages key on the id the flowsheet was built with) and fills
    ``name``, ``formula`` and ``cas`` from the database.

    Raises :class:`UnknownComponentError` with the offending id if the lookup
    fails.
    """
    if not identifier or not identifier.strip():
        raise UnknownComponentError("component identifier is empty")
    pseudo = pseudo_constants(identifier)
    if pseudo is not None:        # registered pseudo-component: no databank hit
        return Component(id=identifier, name=identifier, pseudo=pseudo)
    try:
        meta = _search(identifier)
    except ValueError as exc:
        raise UnknownComponentError(
            f"unknown component {identifier!r}: not recognized by the chemicals "
            f"database (try a CAS number, a formula, or a common name — see "
            f"caldyr.core.components_db.COMMON_COMPONENTS for a curated list)"
        ) from exc
    return Component(
        id=identifier,
        name=str(meta.common_name),
        formula=str(meta.formula),
        cas=str(meta.CASs),
    )


def molar_mass(identifier: str) -> float:
    """Molar mass of ``identifier`` in **kg/mol**, from the `chemicals`
    database (or the pseudo-component registry for registered assay cuts).
    Raises :class:`UnknownComponentError` for unknown ids."""
    pseudo = pseudo_constants(identifier)
    if pseudo is not None:
        return float(pseudo["MW"])
    try:
        meta = _search(identifier)
    except ValueError as exc:
        raise UnknownComponentError(
            f"unknown component {identifier!r}: cannot look up its molar mass "
            f"in the chemicals database"
        ) from exc
    return float(meta.MW) / 1000.0          # chemicals reports g/mol


# -- curated catalog ----------------------------------------------------------
# name/formula/cas resolved through `chemicals` 1.5.2 (search_chemical) at
# author time, 2026-06; ids chosen so that resolve_component(id) returns the
# same species (verified by test_m9_components_db). Formulas are in the Hill
# order chemicals reports (e.g. ammonia = "H3N").
COMMON_COMPONENTS: list[dict[str, str]] = [
    # light / inorganic gases
    {"id": "hydrogen", "name": "hydrogen", "formula": "H2", "cas": "1333-74-0"},
    {"id": "nitrogen", "name": "nitrogen", "formula": "N2", "cas": "7727-37-9"},
    {"id": "oxygen", "name": "oxygen", "formula": "O2", "cas": "7782-44-7"},
    {"id": "argon", "name": "argon", "formula": "Ar", "cas": "7440-37-1"},
    {"id": "helium", "name": "helium", "formula": "He", "cas": "7440-59-7"},
    {"id": "neon", "name": "neon", "formula": "Ne", "cas": "7440-01-9"},
    {"id": "carbon monoxide", "name": "carbon monoxide", "formula": "CO", "cas": "630-08-0"},
    {"id": "carbon dioxide", "name": "carbon dioxide", "formula": "CO2", "cas": "124-38-9"},
    {"id": "water", "name": "water", "formula": "H2O", "cas": "7732-18-5"},
    {"id": "ammonia", "name": "ammonia", "formula": "H3N", "cas": "7664-41-7"},
    {"id": "hydrogen sulfide", "name": "hydrogen sulfide", "formula": "H2S", "cas": "7783-06-4"},
    {"id": "sulfur dioxide", "name": "sulfur dioxide", "formula": "O2S", "cas": "7446-09-5"},
    {"id": "sulfur trioxide", "name": "sulfur trioxide", "formula": "O3S", "cas": "7446-11-9"},
    {"id": "carbonyl sulfide", "name": "carbonyl sulfide", "formula": "COS", "cas": "463-58-1"},
    {"id": "nitric oxide", "name": "nitric oxide", "formula": "NO", "cas": "10102-43-9"},
    {"id": "nitrogen dioxide", "name": "nitrogen dioxide", "formula": "NO2", "cas": "10102-44-0"},
    {"id": "nitrous oxide", "name": "nitrous oxide", "formula": "N2O", "cas": "10024-97-2"},
    {"id": "chlorine", "name": "chlorine", "formula": "Cl2", "cas": "7782-50-5"},
    {"id": "hydrogen chloride", "name": "hydrochloric acid", "formula": "ClH", "cas": "7647-01-0"},
    {"id": "hydrogen cyanide", "name": "hydrogen cyanide", "formula": "CHN", "cas": "74-90-8"},
    # alkanes (C1-C12, plus naphthenes)
    {"id": "methane", "name": "methane", "formula": "CH4", "cas": "74-82-8"},
    {"id": "ethane", "name": "ethane", "formula": "C2H6", "cas": "74-84-0"},
    {"id": "propane", "name": "propane", "formula": "C3H8", "cas": "74-98-6"},
    {"id": "n-butane", "name": "butane", "formula": "C4H10", "cas": "106-97-8"},
    {"id": "isobutane", "name": "isobutane", "formula": "C4H10", "cas": "75-28-5"},
    {"id": "n-pentane", "name": "pentane", "formula": "C5H12", "cas": "109-66-0"},
    {"id": "isopentane", "name": "isopentane", "formula": "C5H12", "cas": "78-78-4"},
    {"id": "neopentane", "name": "neopentane", "formula": "C5H12", "cas": "463-82-1"},
    {"id": "n-hexane", "name": "hexane", "formula": "C6H14", "cas": "110-54-3"},
    {"id": "n-heptane", "name": "heptane", "formula": "C7H16", "cas": "142-82-5"},
    {"id": "n-octane", "name": "octane", "formula": "C8H18", "cas": "111-65-9"},
    {"id": "isooctane", "name": "2,2,4-trimethylpentane", "formula": "C8H18", "cas": "540-84-1"},
    {"id": "n-nonane", "name": "nonane", "formula": "C9H20", "cas": "111-84-2"},
    {"id": "n-decane", "name": "decane", "formula": "C10H22", "cas": "124-18-5"},
    {"id": "n-undecane", "name": "undecane", "formula": "C11H24", "cas": "1120-21-4"},
    {"id": "n-dodecane", "name": "dodecane", "formula": "C12H26", "cas": "112-40-3"},
    {"id": "cyclopentane", "name": "cyclopentane", "formula": "C5H10", "cas": "287-92-3"},
    {"id": "cyclohexane", "name": "cyclohexane", "formula": "C6H12", "cas": "110-82-7"},
    {"id": "methylcyclohexane", "name": "methylcyclohexane", "formula": "C7H14",
     "cas": "108-87-2"},
    # alkenes / alkynes / dienes
    {"id": "ethylene", "name": "ethene", "formula": "C2H4", "cas": "74-85-1"},
    {"id": "propylene", "name": "propene", "formula": "C3H6", "cas": "115-07-1"},
    {"id": "1-butene", "name": "1-butene", "formula": "C4H8", "cas": "106-98-9"},
    {"id": "isobutylene", "name": "isobutene", "formula": "C4H8", "cas": "115-11-7"},
    {"id": "cis-2-butene", "name": "cis-2-butene", "formula": "C4H8", "cas": "590-18-1"},
    {"id": "trans-2-butene", "name": "trans-2-butene", "formula": "C4H8", "cas": "624-64-6"},
    {"id": "1,3-butadiene", "name": "1,3-butadiene", "formula": "C4H6", "cas": "106-99-0"},
    {"id": "1-pentene", "name": "1-pentene", "formula": "C5H10", "cas": "109-67-1"},
    {"id": "1-hexene", "name": "1-hexene", "formula": "C6H12", "cas": "592-41-6"},
    {"id": "1-octene", "name": "1-octene", "formula": "C8H16", "cas": "111-66-0"},
    {"id": "acetylene", "name": "acetylene", "formula": "C2H2", "cas": "74-86-2"},
    {"id": "styrene", "name": "styrene", "formula": "C8H8", "cas": "100-42-5"},
    # aromatics
    {"id": "benzene", "name": "benzene", "formula": "C6H6", "cas": "71-43-2"},
    {"id": "toluene", "name": "toluene", "formula": "C7H8", "cas": "108-88-3"},
    {"id": "ethylbenzene", "name": "ethylbenzene", "formula": "C8H10", "cas": "100-41-4"},
    {"id": "o-xylene", "name": "o-xylene", "formula": "C8H10", "cas": "95-47-6"},
    {"id": "m-xylene", "name": "m-xylene", "formula": "C8H10", "cas": "108-38-3"},
    {"id": "p-xylene", "name": "p-xylene", "formula": "C8H10", "cas": "106-42-3"},
    {"id": "cumene", "name": "cumene", "formula": "C9H12", "cas": "98-82-8"},
    {"id": "naphthalene", "name": "naphthalene", "formula": "C10H8", "cas": "91-20-3"},
    {"id": "phenol", "name": "phenol", "formula": "C6H6O", "cas": "108-95-2"},
    {"id": "aniline", "name": "aniline", "formula": "C6H7N", "cas": "62-53-3"},
    {"id": "pyridine", "name": "pyridine", "formula": "C5H5N", "cas": "110-86-1"},
    # alcohols / glycols
    {"id": "methanol", "name": "methanol", "formula": "CH4O", "cas": "67-56-1"},
    {"id": "ethanol", "name": "ethanol", "formula": "C2H6O", "cas": "64-17-5"},
    {"id": "1-propanol", "name": "1-propanol", "formula": "C3H8O", "cas": "71-23-8"},
    {"id": "isopropanol", "name": "isopropanol", "formula": "C3H8O", "cas": "67-63-0"},
    {"id": "1-butanol", "name": "1-butanol", "formula": "C4H10O", "cas": "71-36-3"},
    {"id": "isobutanol", "name": "2-methyl-1-propanol", "formula": "C4H10O", "cas": "78-83-1"},
    {"id": "2-butanol", "name": "2-butanol", "formula": "C4H10O", "cas": "78-92-2"},
    {"id": "ethylene glycol", "name": "ethylene glycol", "formula": "C2H6O2", "cas": "107-21-1"},
    {"id": "propylene glycol", "name": "1,2-propanediol", "formula": "C3H8O2", "cas": "57-55-6"},
    {"id": "glycerol", "name": "glycerol", "formula": "C3H8O3", "cas": "56-81-5"},
    # aldehydes / ketones
    {"id": "formaldehyde", "name": "formaldehyde", "formula": "CH2O", "cas": "50-00-0"},
    {"id": "acetaldehyde", "name": "acetaldehyde", "formula": "C2H4O", "cas": "75-07-0"},
    {"id": "acetone", "name": "acetone", "formula": "C3H6O", "cas": "67-64-1"},
    {"id": "methyl ethyl ketone", "name": "2-butanone", "formula": "C4H8O", "cas": "78-93-3"},
    {"id": "methyl isobutyl ketone", "name": "4-methyl-2-pentanone", "formula": "C6H12O",
     "cas": "108-10-1"},
    {"id": "furfural", "name": "2-furaldehyde", "formula": "C5H4O2", "cas": "98-01-1"},
    # acids
    {"id": "formic acid", "name": "formic acid", "formula": "CH2O2", "cas": "64-18-6"},
    {"id": "acetic acid", "name": "acetic acid", "formula": "C2H4O2", "cas": "64-19-7"},
    {"id": "propionic acid", "name": "propionic acid", "formula": "C3H6O2", "cas": "79-09-4"},
    {"id": "acrylic acid", "name": "acrylic acid", "formula": "C3H4O2", "cas": "79-10-7"},
    {"id": "sulfuric acid", "name": "sulfuric acid", "formula": "H2O4S", "cas": "7664-93-9"},
    {"id": "nitric acid", "name": "nitric acid", "formula": "HNO3", "cas": "7697-37-2"},
    # esters / ethers / oxides
    {"id": "methyl acetate", "name": "methyl acetate", "formula": "C3H6O2", "cas": "79-20-9"},
    {"id": "ethyl acetate", "name": "ethyl acetate", "formula": "C4H8O2", "cas": "141-78-6"},
    {"id": "n-butyl acetate", "name": "butyl acetate", "formula": "C6H12O2", "cas": "123-86-4"},
    {"id": "vinyl acetate", "name": "vinyl acetate", "formula": "C4H6O2", "cas": "108-05-4"},
    {"id": "methyl methacrylate", "name": "methyl methacrylate", "formula": "C5H8O2",
     "cas": "80-62-6"},
    {"id": "dimethyl ether", "name": "dimethyl ether", "formula": "C2H6O", "cas": "115-10-6"},
    {"id": "diethyl ether", "name": "diethyl ether", "formula": "C4H10O", "cas": "60-29-7"},
    {"id": "methyl tert-butyl ether", "name": "tert-butyl methyl ether", "formula": "C5H12O",
     "cas": "1634-04-4"},
    {"id": "tetrahydrofuran", "name": "tetrahydrofuran", "formula": "C4H8O", "cas": "109-99-9"},
    {"id": "1,4-dioxane", "name": "1,4-dioxane", "formula": "C4H8O2", "cas": "123-91-1"},
    {"id": "ethylene oxide", "name": "oxirane", "formula": "C2H4O", "cas": "75-21-8"},
    {"id": "propylene oxide", "name": "propylene oxide", "formula": "C3H6O", "cas": "75-56-9"},
    # nitrogen chemicals / common dipolar solvents
    {"id": "acetonitrile", "name": "acetonitrile", "formula": "C2H3N", "cas": "75-05-8"},
    {"id": "acrylonitrile", "name": "acrylonitrile", "formula": "C3H3N", "cas": "107-13-1"},
    {"id": "monoethanolamine", "name": "ethanolamine", "formula": "C2H7NO", "cas": "141-43-5"},
    {"id": "diethanolamine", "name": "diethanolamine", "formula": "C4H11NO2", "cas": "111-42-2"},
    {"id": "methylamine", "name": "methylamine", "formula": "CH5N", "cas": "74-89-5"},
    {"id": "dimethylamine", "name": "dimethylamine", "formula": "C2H7N", "cas": "124-40-3"},
    {"id": "trimethylamine", "name": "trimethylamine", "formula": "C3H9N", "cas": "75-50-3"},
    {"id": "urea", "name": "urea", "formula": "CH4N2O", "cas": "57-13-6"},
    {"id": "dimethylformamide", "name": "n,n-dimethylformamide", "formula": "C3H7NO",
     "cas": "68-12-2"},
    {"id": "dimethyl sulfoxide", "name": "dimethyl sulfoxide", "formula": "C2H6OS",
     "cas": "67-68-5"},
    # chlorinated
    {"id": "dichloromethane", "name": "dichloromethane", "formula": "CH2Cl2", "cas": "75-09-2"},
    {"id": "chloroform", "name": "chloroform", "formula": "CHCl3", "cas": "67-66-3"},
    {"id": "carbon tetrachloride", "name": "carbon tetrachloride", "formula": "CCl4",
     "cas": "56-23-5"},
    {"id": "vinyl chloride", "name": "ethene, chloro-", "formula": "C2H3Cl", "cas": "75-01-4"},
    {"id": "1,2-dichloroethane", "name": "1,2-dichloroethane", "formula": "C2H4Cl2",
     "cas": "107-06-2"},
    {"id": "chlorobenzene", "name": "chlorobenzene", "formula": "C6H5Cl", "cas": "108-90-7"},
    # refrigerants (ASHRAE R-numbers resolve in chemicals except R-32/R-125,
    # whose ids are therefore the chemical names — "R-32" is unrecognized and
    # "R-125" mis-resolves, so the names are the reliable identifiers)
    {"id": "R-22", "name": "difluorochloromethane", "formula": "CHClF2", "cas": "75-45-6"},
    {"id": "difluoromethane", "name": "difluoromethane", "formula": "CH2F2", "cas": "75-10-5"},
    {"id": "pentafluoroethane", "name": "pentafluoroethane", "formula": "C2HF5",
     "cas": "354-33-6"},
    {"id": "R-134a", "name": "norflurane", "formula": "C2H2F4", "cas": "811-97-2"},
    {"id": "R-143a", "name": "1,1,1-Trifluoroethane", "formula": "C2H3F3", "cas": "420-46-2"},
    {"id": "R-152a", "name": "1,1-difluoroethane", "formula": "C2H4F2", "cas": "75-37-6"},
    {"id": "R-1234yf", "name": "2,3,3,3-Tetrafluoropropene", "formula": "C3H2F4",
     "cas": "754-12-1"},
    {"id": "R-1234ze(E)", "name": "(1E)-1,3,3,3-Tetrafluoro-1-propene", "formula": "C3H2F4",
     "cas": "29118-24-9"},
]
