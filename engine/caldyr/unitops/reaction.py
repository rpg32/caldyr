"""Reaction stoichiometry + shared reactor outlet/energy logic.

Reactions are stored as plain dicts in unit-op ``params`` (so ``.flow`` JSON
round-trips without custom encoders); :class:`Reaction` is the typed view a
reactor builds from that dict. The enthalpy basis is formation-inclusive (see
:mod:`caldyr.thermo._flasher`), so reactor energy balances carry the heat of
reaction automatically — an adiabatic reactor's temperature change falls out of
a plain PH flash.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from ..core import Stream

_R = 8.314462618          # J/mol/K


@dataclass(frozen=True)
class Reaction:
    """A single stoichiometric reaction. ``stoich`` maps component id to its
    coefficient (negative for reactants, positive for products); ``key`` is the
    limiting reactant a conversion is measured against."""
    stoich: dict[str, float]
    key: str | None = None

    @classmethod
    def from_param(cls, d: dict) -> "Reaction":
        return cls(stoich={k: float(v) for k, v in d["stoich"].items()}, key=d.get("key"))

    @property
    def dn(self) -> float:
        """Change in total moles per unit extent, Σ ν_i."""
        return sum(self.stoich.values())

    def extent_for_conversion(self, moles: dict[str, float], conversion: float) -> float:
        """Reaction extent ξ that consumes ``conversion`` of the key reactant
        present in ``moles``."""
        if self.key is None:
            raise ValueError("reaction needs a 'key' reactant for a conversion spec")
        nu_key = self.stoich[self.key]
        if nu_key >= 0:
            raise ValueError(f"key {self.key!r} must be a reactant (negative coefficient)")
        return conversion * moles.get(self.key, 0.0) / -nu_key


@dataclass(frozen=True)
class KineticReaction:
    """A power-law rate expression on top of a stoichiometric reaction:

        r = k0 · exp(-Ea / RT) · Π C_i^order_i          [mol/(m^3·s)]

    Param-dict form (JSON-friendly, like :class:`Reaction`)::

        {"stoich": {...}, "key": <reactant>, "k0": ..., "Ea": ...,
         "orders": {component: order}}        # orders default to {key: 1}

    ``k0``'s units must make r come out in mol/(m^3·s) given concentrations in
    mol/m^3 (e.g. 1/s for first order, m^3/(mol·s) for second). ``Ea`` is in
    J/mol. Concentrations are clamped at zero so a component driven to
    exhaustion stops the reaction instead of producing NaNs.

    **Reversible reactions.** Supplying ``k0_rev``/``Ea_rev`` (and optional
    ``orders_rev`` over the products, defaulting to the product stoichiometry)
    adds a reverse term so the *net* rate is

        r = k_f · Π C_i^order_i  −  k_r · Π C_j^order_rev,j

    which vanishes — and pins the conversion — at the equilibrium
    ``Π C_prod^ν / Π C_react^ν = k_f/k_r``. This is what an esterification or
    etherification in a reactive-distillation column needs: a forward-only fast
    reaction would drive to complete (unphysical) conversion. Without the
    reverse parameters the rate is forward-only (unchanged).
    """
    stoich: dict[str, float]
    key: str
    k0: float
    Ea: float
    orders: dict[str, float] = field(default_factory=dict)
    k0_rev: float = 0.0
    Ea_rev: float = 0.0
    orders_rev: dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_param(cls, d: dict) -> "KineticReaction":
        key = d.get("key")
        if key is None:
            raise ValueError("kinetic reaction needs a 'key' reactant")
        stoich = {k: float(v) for k, v in d["stoich"].items()}
        if stoich.get(key, 0.0) >= 0:
            raise ValueError(
                f"kinetic reaction key {key!r} must be a reactant "
                f"(negative coefficient)")
        missing = [f for f in ("k0", "Ea") if f not in d]
        if missing:
            raise ValueError(f"kinetic reaction is missing field(s) {missing} "
                             f"(power law needs k0 and Ea)")
        orders = {k: float(v) for k, v in d.get("orders", {key: 1.0}).items()}
        # Reverse term: default the orders to the product stoichiometry.
        k0_rev = float(d.get("k0_rev", 0.0))
        Ea_rev = float(d.get("Ea_rev", 0.0))
        default_rev = {c: nu for c, nu in stoich.items() if nu > 0}
        orders_rev = {k: float(v)
                      for k, v in d.get("orders_rev", default_rev).items()}
        return cls(stoich=stoich, key=key, k0=float(d["k0"]), Ea=float(d["Ea"]),
                   orders=orders, k0_rev=k0_rev, Ea_rev=Ea_rev,
                   orders_rev=orders_rev)

    def rate(self, conc: dict[str, float], T: float) -> float:
        """Net rate in mol/(m^3·s) at concentrations ``conc`` (mol/m^3) and
        T (K): forward minus reverse (reverse is zero unless ``k0_rev`` is set)."""
        r = self.k0 * math.exp(-self.Ea / (_R * T))
        for comp, order in self.orders.items():
            r *= max(conc.get(comp, 0.0), 0.0) ** order
        if self.k0_rev > 0.0:
            rr = self.k0_rev * math.exp(-self.Ea_rev / (_R * T))
            for comp, order in self.orders_rev.items():
                rr *= max(conc.get(comp, 0.0), 0.0) ** order
            r -= rr
        return r


def concentrations(pp, T: float, P: float,
                   moles: dict[str, float]) -> dict[str, float]:
    """Molar concentrations C_i = z_i / v(T,P,z) in mol/m^3, with the bulk
    molar volume from the property package.

    This is the *single* concentration basis the kinetic reactors use, for
    vapor and liquid alike — a real-gas/compressed-liquid volume, not the
    ideal-gas P/RT shortcut, so kinetics stay consistent with the rest of the
    flowsheet's thermodynamics. (Negative round-off moles are clamped to 0.)
    """
    moles = {c: max(m, 0.0) for c, m in moles.items()}
    n_tot = sum(moles.values())
    if n_tot <= 0.0:
        raise ValueError("cannot evaluate concentrations: total moles <= 0")
    z = {c: m / n_tot for c, m in moles.items()}
    v = pp.volume(T, P, z)                       # m^3/mol, bulk
    return {c: zi / v for c, zi in z.items()}


def apply_extents(moles: dict[str, float], reactions, extents) -> dict[str, float]:
    """Return moles after applying each reaction at its extent (in order)."""
    out = dict(moles)
    for rxn, xi in zip(reactions, extents):
        for comp, nu in rxn.stoich.items():
            out[comp] = out.get(comp, 0.0) + nu * xi
    # numerical floor: tiny negatives from round-off become zero
    return {c: (v if abs(v) > 1e-12 else 0.0) for c, v in out.items()}


def reactor_outlet(unit_id: str, inlet: Stream, pp, moles_out: dict[str, float],
                   P_out: float, T_spec: float | None):
    """Build the outlet stream + duty from a post-reaction mole inventory.

    Isothermal (``T_spec`` given): outlet at that temperature; duty is the heat
    that must be added/removed. Adiabatic (``T_spec`` is None): conserve total
    enthalpy flow via a PH flash; duty is zero and the temperature moves on its
    own (exothermic reactions heat the stream).
    """
    n_out = sum(moles_out.values())
    if n_out <= 0:
        raise ValueError(f"reactor {unit_id!r} produced non-positive total moles")
    z_out = {c: m / n_out for c, m in moles_out.items()}

    T_in, P_in, n_in = inlet.require_state()
    z_in = inlet.normalized_z()
    h_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z_in)
    H_in_total = n_in * h_in

    if T_spec is not None:
        res = pp.flash_pt(T_spec, P_out, z_out)
        duty = n_out * res.H - H_in_total          # heat added to hold T_spec
    else:
        res = pp.flash_ph(P_out, H_in_total / n_out, z_out)
        duty = 0.0

    out = Stream(
        id=f"{unit_id}.out", components=list(inlet.components),
        T=res.T, P=res.P, molar_flow=n_out, z=z_out,
        H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
    )
    return out, duty
