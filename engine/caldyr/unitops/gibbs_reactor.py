"""Multi-reaction chemical-equilibrium reactor (Gibbs minimization via Cantera).

Where the single-reaction :class:`EquilibriumReactor` solves one mass-action
extent, the GibbsReactor minimizes the total Gibbs energy of an ideal-gas
mixture over *all* species simultaneously (Cantera's ``equilibrate("TP")``), so
coupled reaction networks — steam-methane reforming + water-gas shift, ammonia
+ cracking side reactions — come out without enumerating the reactions at all.
Species that cannot react (no shared elements, e.g. argon) are simply carried
through by the minimization; no explicit ``inerts`` list is needed.

**Species mapping.** Caldyr component ids are resolved to CAS numbers through
the `chemicals` database and mapped onto the GRI-Mech 3.0 species set bundled
with Cantera (``gri30.yaml`` — H2/N2/NH3/CH4/CO/CO2/H2O/O2/Ar and the C1-C3 /
N-chemistry around them). The reactor builds an ideal-gas
:class:`cantera.Solution` over *only the mapped flowsheet components*, so the
equilibrium atom balance closes exactly over the streams the engine sees (no
trace radicals leak out of the unit). A component the table cannot map raises a
typed :class:`CanteraSpeciesError` naming it — unless its inlet flow is zero,
in which case it passes through untouched (it cannot participate anyway).

**Enthalpy/duty basis.** Only the outlet *composition* comes from Cantera; the
outlet state (enthalpy, phase) is computed with the flowsheet's own property
package, exactly like the other reactors. The reported duty therefore closes on
caldyr's formation-inclusive enthalpy basis — consistent with every other unit
in the flowsheet — and not on Cantera's NASA-polynomial basis (the two agree to
within the small differences between the underlying thermochemical data).
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .reaction import reactor_outlet

_FLOW_EPS = 1e-12     # mol/s below which an unmappable component is "absent"

# CAS -> GRI-Mech 3.0 species name, for the stable species a flowsheet would
# carry. Caldyr ids resolve to CAS via `chemicals`, so any synonym the database
# knows ("water", "H2O", "7732-18-5") lands on the right Cantera species.
# Radical/intermediate gri30 species are deliberately not exposed.
_CAS_TO_GRI30: dict[str, str] = {
    "1333-74-0": "H2",        # hydrogen
    "7727-37-9": "N2",        # nitrogen
    "7782-44-7": "O2",        # oxygen
    "7440-37-1": "AR",        # argon
    "7732-18-5": "H2O",       # water
    "7664-41-7": "NH3",       # ammonia
    "74-82-8": "CH4",         # methane
    "630-08-0": "CO",         # carbon monoxide
    "124-38-9": "CO2",        # carbon dioxide
    "67-56-1": "CH3OH",       # methanol
    "50-00-0": "CH2O",        # formaldehyde
    "75-07-0": "CH3CHO",      # acetaldehyde
    "74-86-2": "C2H2",        # acetylene
    "74-85-1": "C2H4",        # ethylene
    "74-84-0": "C2H6",        # ethane
    "74-98-6": "C3H8",        # propane
    "7722-84-1": "H2O2",      # hydrogen peroxide
    "74-90-8": "HCN",         # hydrogen cyanide
    "10102-43-9": "NO",       # nitric oxide
    "10102-44-0": "NO2",      # nitrogen dioxide
    "10024-97-2": "N2O",      # nitrous oxide
    "463-51-4": "CH2CO",      # ketene
}


class CanteraSpeciesError(ValueError):
    """A flowing component could not be mapped onto the bundled Cantera
    (gri30.yaml) species set."""


@lru_cache(maxsize=512)
def _gri30_name(component_id: str) -> str | None:
    """GRI-Mech 3.0 species name for a caldyr component id, or None if the id
    does not resolve / is not in the mapped species set."""
    from chemicals.identifiers import CAS_from_any

    try:
        cas = CAS_from_any(component_id)
    except ValueError:
        return None
    return _CAS_TO_GRI30.get(cas)


@lru_cache(maxsize=32)
def _build_solution(species_names: tuple[str, ...]) -> Any:
    """Ideal-gas Cantera Solution restricted to ``species_names`` (a subset of
    gri30.yaml). Cached: Solutions are reusable, and we always set TPX before
    use. Restricting the species set keeps the equilibrium atom balance exact
    over the flowsheet components (no trace radicals)."""
    import cantera as ct

    by_name = {s.name: s for s in ct.Species.list_from_file("gri30.yaml")}
    return ct.Solution(thermo="ideal-gas",
                       species=[by_name[n] for n in species_names])


@register("GibbsReactor")
class GibbsReactor(UnitOp):
    """Isothermal multi-reaction equilibrium reactor (Gibbs minimization).

    Params:
      * ``T`` (required) — outlet temperature, K (the equilibrium is evaluated
        at outlet T and P; the duty to hold it is reported on ``duty``).
      * ``P`` (optional) — outlet pressure, Pa; defaults to the inlet pressure.
      * ``dP`` (optional) — pressure drop, Pa, subtracted from the above.

    See the module docstring for the species mapping and the enthalpy basis.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"GibbsReactor {self.id!r}: missing/empty inlet on 'in1'")
        if self.params.get("T") is None:
            raise ValueError(
                f"GibbsReactor {self.id!r}: params['T'] is required — the Gibbs "
                f"equilibrium is evaluated at a specified outlet temperature"
            )
        T = float(self.params["T"])
        _, P_in, n_in = inlet.require_state()
        P_out = float(self.params.get("P", P_in)) - float(self.params.get("dP", 0.0))
        if P_out <= 0.0:
            raise ValueError(f"GibbsReactor {self.id!r}: outlet pressure {P_out} Pa <= 0")
        z_in = inlet.normalized_z()

        # Map flowsheet components onto gri30 species. Unmapped components with
        # flow are a hard, typed error; unmapped zero-flow components are inert
        # bystanders of the equilibrium (they cannot participate or form).
        mapped: dict[str, str] = {}            # component id -> gri30 name
        unmapped_flowing: list[str] = []
        for comp in inlet.components:
            name = _gri30_name(comp)
            if name is not None:
                mapped[comp] = name
            elif n_in * z_in.get(comp, 0.0) > _FLOW_EPS:
                unmapped_flowing.append(comp)
        if unmapped_flowing:
            raise CanteraSpeciesError(
                f"GibbsReactor {self.id!r}: component(s) {unmapped_flowing!r} have "
                f"no mapping onto the bundled Cantera gri30.yaml species set; "
                f"supported species are "
                f"{sorted(set(_CAS_TO_GRI30.values()))}"
            )
        if not mapped or all(n_in * z_in.get(c, 0.0) <= _FLOW_EPS for c in mapped):
            raise CanteraSpeciesError(
                f"GibbsReactor {self.id!r}: no mappable components are flowing"
            )

        gas = _build_solution(tuple(mapped.values()))
        gas.TPX = T, P_out, {mapped[c]: z_in.get(c, 0.0) for c in mapped}
        mw_in = gas.mean_molecular_weight           # kg/kmol
        gas.equilibrate("TP")
        # Mass is conserved by the equilibrium; total moles follow from the
        # mean molar mass change (n_out * MW_out = n_in * MW_in).
        n_out = n_in * mw_in / gas.mean_molecular_weight

        x = gas.mole_fraction_dict(threshold=0.0)
        moles_out = {c: n_out * float(x.get(mapped[c], 0.0)) for c in mapped}
        for comp in inlet.components:                # unmapped zero-flow comps
            moles_out.setdefault(comp, 0.0)

        out, duty = reactor_outlet(self.id, inlet, pp, moles_out, P_out, T)
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}
