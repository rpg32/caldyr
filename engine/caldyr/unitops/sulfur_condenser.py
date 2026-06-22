"""Sulfur condenser for a Claus train: cool the process gas and knock out liquid
elemental sulfur.

After each Claus reaction stage the gas is cooled (the thermal stage's cooler is
the waste-heat boiler) so that the elemental sulfur formed condenses to a liquid
and is drained, shifting the next stage's equilibrium toward yet more sulfur. The
gas leaving is *saturated* in sulfur vapour at the condenser temperature — the
small residual that escapes is the per-stage sulfur "loss", and the loss from the
final condenser sets the plant's tail-gas emission.

Model
-----
* **Phase split.** Only *elemental* sulfur (the vapour allotropes S2/S8 in the
  feed) can condense; the compound sulfur species (H2S, SO2, COS, CS2) and the
  inerts stay in the gas. The exit gas holds elemental sulfur at its saturation
  partial pressure ``p_S = liquid_sulfur_psat(T)`` (:mod:`caldyr.thermo.sulfur`),
  modelled as S8 (the dominant allotrope at condenser temperatures): the residual
  S8 vapour satisfies ``y_S8 = p_S / P``. Everything above that condenses to
  liquid. Elemental sulfur is tracked as **S8** throughout — the liquid product
  is S8 and any residual elemental vapour is re-expressed as S8 — so the sulfur
  atom balance is exact (S2 is re-speciated to S8; the next catalytic converter
  re-equilibrates the allotropes regardless).
* **Energy.** The duty is ``H_gas_out + H_liquid_out − H_gas_in``: the sensible
  cooling of the gas (gas enthalpies from the property package) plus the latent
  heat released by the condensing sulfur. Liquid-S8 enthalpy is the gas-phase S8
  enthalpy at the condenser temperature minus the heat of vaporisation
  (``8 × sulfur_hvap_per_atom(T)``), keeping the liquid on the same NASA basis as
  the gas so the balance closes.

Requires the ``nasa:gas`` property package (it must carry S8); ``params['T']`` is
the condenser outlet temperature (K), ``params['P']`` / ``params['dP']`` the
outlet pressure (defaults to the inlet pressure).
"""
from __future__ import annotations

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from ..thermo.sulfur import liquid_sulfur_psat, sulfur_hvap_per_atom
from .base import register

# Caldyr ids accepted as elemental sulfur vapour, with their S-atom count.
_ELEMENTAL_S = {"S2": 2, "S8": 8, "sulfur dimer": 2, "cyclooctasulfur": 8}
_S8_ATOMS = 8


class SulfurCondenserError(ValueError):
    """A SulfurCondenser could not be set up or solved."""


def _elemental_id(components: list[str]) -> tuple[str, dict[str, int]]:
    """Pick the component id used to carry liquid/residual sulfur (S8 preferred)
    and the map of feed component id -> S-atom count for the elemental vapours
    present. Raises if no elemental-sulfur component is in the list."""
    from chemicals.identifiers import CAS_from_any

    # CAS -> atoms, robust to id synonyms.
    cas_atoms = {"23550-45-0": 2, "10544-50-0": 8}
    present: dict[str, int] = {}
    s8_id = None
    for c in components:
        try:
            cas = CAS_from_any(c)
        except ValueError:
            continue
        if cas in cas_atoms:
            present[c] = cas_atoms[cas]
            if cas == "10544-50-0":
                s8_id = c
    if s8_id is None:
        raise SulfurCondenserError(
            "SulfurCondenser requires S8 in the component list to carry liquid "
            f"and residual elemental sulfur; got {components!r}")
    return s8_id, present


@register("SulfurCondenser")
class SulfurCondenser(UnitOp):
    """Cool Claus process gas to ``params['T']`` and drain liquid sulfur.

    Ports ``in1`` (process gas), ``gas`` (cooled, sulfur-saturated gas),
    ``liquid`` (liquid sulfur product, as S8) and ``duty``.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("gas", "outlet"),
                Port("liquid", "outlet"), Port("duty", "outlet", "energy")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise SulfurCondenserError(
                f"SulfurCondenser {self.id!r}: missing/empty inlet on 'in1'")
        if self.params.get("T") is None:
            raise SulfurCondenserError(
                f"SulfurCondenser {self.id!r}: params['T'] (condenser outlet "
                f"temperature, K) is required")

        T_in, P_in, n_in = inlet.require_state()
        T = float(self.params["T"])
        P = float(self.params.get("P", P_in)) - float(self.params.get("dP", 0.0))
        if P <= 0.0:
            raise SulfurCondenserError(
                f"SulfurCondenser {self.id!r}: outlet pressure {P} Pa <= 0")
        components = list(inlet.components)
        s8_id, elemental = _elemental_id(components)
        z_in = inlet.normalized_z()

        moles_in = {c: n_in * z_in.get(c, 0.0) for c in components}
        # Total elemental-sulfur atoms entering as vapour (S2 + S8 + ...).
        s_atoms_in = sum(moles_in.get(c, 0.0) * k for c, k in elemental.items())
        # Non-elemental gas (everything that stays vapour): keep unchanged.
        non_elem = {c: m for c, m in moles_in.items() if c not in elemental}
        n_non_elem = sum(non_elem.values())

        # Saturation: residual elemental sulfur leaves as S8 at p_S8 = Psat(T).
        p_sat = liquid_sulfur_psat(T)
        if p_sat >= P:
            # Above the sulfur boiling point at this P nothing condenses.
            n_s8_vap = s_atoms_in / _S8_ATOMS
        else:
            # y_S8 = n_s8_vap / (n_non_elem + n_s8_vap) = p_sat / P.
            n_s8_vap = n_non_elem * p_sat / (P - p_sat)
        s_atoms_vap = min(_S8_ATOMS * n_s8_vap, s_atoms_in)   # cannot exceed feed
        n_s8_vap = s_atoms_vap / _S8_ATOMS
        s_atoms_liq = s_atoms_in - s_atoms_vap
        n_s8_liq = s_atoms_liq / _S8_ATOMS

        # -- gas product: inerts + residual S8 vapour --------------------------
        gas_moles = dict(non_elem)
        gas_moles[s8_id] = gas_moles.get(s8_id, 0.0) + n_s8_vap
        for c in elemental:                              # S2 etc fully re-speciated
            if c != s8_id:
                gas_moles.setdefault(c, 0.0)
        n_gas = sum(gas_moles.values())
        z_gas = {c: gas_moles.get(c, 0.0) / n_gas for c in components}
        h_gas = pp.enthalpy(T, P, z_gas)

        gas = Stream(id=f"{self.id}.gas", components=components,
                     T=T, P=P, molar_flow=n_gas, z=z_gas,
                     H=h_gas, phase="vapor", vapor_fraction=1.0)

        # -- liquid sulfur product (as S8) ------------------------------------
        # Liquid-S8 enthalpy on the NASA basis: gas-phase S8 at T minus latent.
        z_s8 = {c: (1.0 if c == s8_id else 0.0) for c in components}
        h_s8_gas = pp.enthalpy_vapor(T, P, z_s8)
        h_s8_liq = h_s8_gas - _S8_ATOMS * sulfur_hvap_per_atom(T)
        liquid = Stream(id=f"{self.id}.liquid", components=components,
                        T=T, P=P, molar_flow=n_s8_liq, z=dict(z_s8),
                        H=h_s8_liq, phase="liquid", vapor_fraction=0.0)

        # -- duty (heat removed: negative) ------------------------------------
        h_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z_in)
        duty = n_gas * h_gas + n_s8_liq * h_s8_liq - n_in * h_in

        return {"gas": gas, "liquid": liquid,
                "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}
