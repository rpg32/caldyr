"""Claus sulfur-recovery reactor: thermal furnace and catalytic converter.

The Claus process turns the H2S in an acid-gas stream (the product of an amine
sweetening unit — see :mod:`caldyr.thermo.amine_pkg`) into elemental sulfur in
two reaction regimes, both of which this one unit op covers by chemical
equilibrium (Cantera ``equilibrate`` over the mapped flowsheet species, exactly
the :class:`GibbsReactor` pattern but over the sulfur slate of
``nasa_gas.yaml``):

* **Thermal stage (reaction furnace) — adiabatic.** A third of the H2S is burnt
  with air,  ``H2S + 3/2 O2 -> SO2 + H2O``, and the SO2 then reacts with the
  remaining H2S, ``2 H2S + SO2 -> 3/x S_x + 2 H2O``. At the ~1300-1600 K flame
  temperature the sulfur is the dimer S2. Run with ``params['T']`` omitted: the
  outlet temperature and composition come out together from a constant-enthalpy
  (``equilibrate('HP')``) flame calculation, so the adiabatic flame temperature
  is predicted, not assumed. The duty is zero.

* **Catalytic converter — isothermal.** Over alumina at ~470-570 K the Claus
  reaction ``2 H2S + SO2 -> 3/8 S8 + 2 H2O`` runs toward more sulfur (it is
  exothermic, so a lower temperature gives higher conversion); the sulfur is now
  the ring S8. Run with ``params['T']`` set: ``equilibrate('TP')`` gives the
  converter outlet and the duty to hold the temperature.

**Thermo.** The reactor builds its own ideal-gas Cantera solution over the
mapped components and takes the *composition* (and, when adiabatic, the
temperature) from it. The outlet stream's enthalpy is then evaluated through the
flowsheet's property package, which for a Claus plant is ``nasa:gas`` — the same
NASA basis Cantera uses — so the duty closes consistently (see
:mod:`caldyr.thermo.nasa_pkg`). Components that cannot be mapped onto
``nasa_gas.yaml`` are a typed error if they flow, and pass through untouched if
their inlet flow is zero.
"""
from __future__ import annotations

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from ..thermo.nasa_pkg import _build_solution, _nasa_name
from .base import register

_FLOW_EPS = 1e-12         # mol/s below which an unmappable component is "absent"
_KMOL = 1000.0


class ClausReactorError(ValueError):
    """A ClausReactor could not be set up or solved (bad spec / unmappable
    flowing component)."""


@register("ClausReactor")
class ClausReactor(UnitOp):
    """Claus reaction-furnace (adiabatic) or catalytic converter (isothermal),
    by Cantera chemical equilibrium over the sulfur slate.

    Params:
      * ``T`` (optional) — outlet temperature, K. Omit for the adiabatic thermal
        furnace (flame temperature predicted); set it for a catalytic converter.
      * ``P`` (optional) — outlet pressure, Pa; defaults to the ``in1`` pressure.
      * ``dP`` (optional) — pressure drop, Pa, subtracted from the above.

    Ports ``in1`` (acid gas), ``in2`` (optional second feed, e.g. combustion
    air), ``out`` and ``duty``.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("in2", "inlet"),
                Port("out", "outlet"), Port("duty", "outlet", "energy")]

    # -- feed handling -------------------------------------------------------
    def _combine_inlets(self, inlets: dict[str, Stream], pp
                        ) -> tuple[list[str], dict[str, float], float, float]:
        """Sum the wired material inlets into (component order, moles per
        component, total inlet enthalpy [W], reference pressure)."""
        streams = [s for k, s in inlets.items()
                   if isinstance(s, Stream) and s is not None and s.molar_flow]
        if not streams:
            raise ClausReactorError(
                f"ClausReactor {self.id!r}: no non-empty inlet on in1/in2")
        components = list(streams[0].components)
        moles: dict[str, float] = {c: 0.0 for c in components}
        H_in = 0.0
        pressures: list[float] = []
        for s in streams:
            T, P, n = s.require_state()
            z = s.normalized_z()
            h = s.H if s.H is not None else pp.enthalpy(T, P, z)
            H_in += n * h
            for c in s.components:
                moles[c] = moles.get(c, 0.0) + n * z.get(c, 0.0)
            pressures.append(P)
        return components, moles, H_in, min(pressures)

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        components, moles_in, H_in, P_in = self._combine_inlets(inlets, pp)
        n_in = sum(moles_in.values())
        P_out = float(self.params.get("P", P_in)) - float(self.params.get("dP", 0.0))
        if P_out <= 0.0:
            raise ClausReactorError(
                f"ClausReactor {self.id!r}: outlet pressure {P_out} Pa <= 0")

        # Map components -> nasa_gas species; flowing & unmappable is an error.
        mapped: dict[str, str] = {}
        unmapped_flowing: list[str] = []
        for c in components:
            name = _nasa_name(c)
            if name is not None:
                mapped[c] = name
            elif moles_in.get(c, 0.0) > _FLOW_EPS:
                unmapped_flowing.append(c)
        if unmapped_flowing:
            raise ClausReactorError(
                f"ClausReactor {self.id!r}: component(s) {unmapped_flowing!r} "
                f"have no mapping onto nasa_gas.yaml; supported sulfur/combustion "
                f"species only (use the 'nasa:gas' property package)")
        if all(moles_in.get(c, 0.0) <= _FLOW_EPS for c in mapped):
            raise ClausReactorError(
                f"ClausReactor {self.id!r}: no mappable components are flowing")

        gas = _build_solution(tuple(mapped.values()))
        X_in = {mapped[c]: moles_in.get(c, 0.0) / n_in for c in mapped}

        # Mean MW of the feed (kg/kmol) before reaction, for the n_out balance.
        gas.TPX = 298.15, P_out, X_in
        mw_in = float(gas.mean_molecular_weight)

        T_spec = self.params.get("T")
        if T_spec is None:                                   # adiabatic furnace
            gas.HPX = (H_in / n_in) * _KMOL / mw_in, P_out, X_in
            gas.equilibrate("HP")
            T_out = float(gas.T)
        else:                                                # isothermal converter
            T_out = float(T_spec)
            gas.TPX = T_out, P_out, X_in
            gas.equilibrate("TP")

        # Mass is conserved by the equilibrium; total moles follow the mean
        # molar-mass change: n_out * MW_out = n_in * MW_in.
        n_out = n_in * mw_in / float(gas.mean_molecular_weight)

        x_out = gas.mole_fraction_dict(threshold=0.0)
        moles_out = {c: n_out * float(x_out.get(mapped[c], 0.0)) for c in mapped}
        for c in components:                                 # zero-flow passthrough
            moles_out.setdefault(c, 0.0)
        n_out = sum(moles_out.values())
        z_out = {c: m / n_out for c, m in moles_out.items()}

        H_out_molar = pp.enthalpy(T_out, P_out, z_out)
        duty = 0.0 if T_spec is None else n_out * H_out_molar - H_in

        out = Stream(id=f"{self.id}.out", components=list(components),
                     T=T_out, P=P_out, molar_flow=n_out, z=z_out,
                     H=H_out_molar, phase="vapor", vapor_fraction=1.0)
        return {"out": out,
                "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}
