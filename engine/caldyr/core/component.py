from dataclasses import dataclass, field


@dataclass(frozen=True)
class Component:
    """A chemical species. Property backends resolve by id/cas; no property
    data is stored on the component itself — EXCEPT for petroleum
    pseudo-components, which by definition have no databank entry.

    ``pseudo`` carries the characterized constants of a petroleum
    pseudo-component (a boiling-point cut from an assay — see
    :mod:`caldyr.assay`), keyed:

    ====== ======================================== =========
    key    meaning                                  unit
    ====== ======================================== =========
    MW     molar mass                               kg/mol
    Tb     normal boiling point                     K
    SG     specific gravity (60F/60F)               --
    Tc     critical temperature                     K
    Pc     critical pressure                        Pa
    omega  acentric factor                          --
    Cp_ig  ideal-gas Cp poly coeffs [a0, a1, a2]    J/mol/K (ascending powers of T[K])
    Hf     ideal-gas formation enthalpy (optional)  J/mol (default 0 — see note)
    Gf     ideal-gas formation Gibbs (optional)     J/mol (default 0 — see note)
    ====== ======================================== =========

    Constructing a Component with ``pseudo`` registers the constants in the
    process-wide pseudo-component registry
    (:func:`caldyr.core.components_db.register_pseudo_component`), so that
    property packages — which receive only component *id* strings — can build
    flashers for it. Re-registering the same id with different constants
    overwrites (last writer wins); thermo flasher caches key on the constants
    themselves, so stale flashers are never reused.

    Note on formation properties: pseudo-components default to Hf = Gf = 0
    (their true values are unknowable for a lumped cut). Energy balances stay
    exact because the formation offset cancels wherever composition is
    conserved — but REACTIONS involving pseudo-components are unsupported
    (heats of reaction / equilibrium constants would be wrong) and the thermo
    layer makes no attempt to support them.
    """
    id: str                 # canonical key, e.g. "water" or "NBP_123C"
    name: str = ""
    formula: str | None = None
    cas: str | None = None
    pseudo: dict | None = field(default=None, compare=True)

    def __post_init__(self) -> None:
        if self.pseudo:
            from .components_db import register_pseudo_component

            register_pseudo_component(self.id, self.pseudo)
