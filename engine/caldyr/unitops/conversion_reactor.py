from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .reaction import Reaction, apply_extents, reactor_outlet


@register("ConversionReactor")
class ConversionReactor(UnitOp):
    """Stoichiometric reactor with specified fractional conversion(s).

    Two param forms:
      * single — ``params['reaction'] = {"stoich": {...}, "key": <reactant>}``
        with ``params['conversion'] = X``;
      * multiple — ``params['reactions'] = [{"stoich": ..., "key": ...,
        "conversion": X}, ...]`` (each reaction carries its own conversion).
    Each reaction's extent consumes fraction ``X`` of its key reactant present
    when it is applied (reactions apply in order). Isothermal if ``params['T_out']``
    is given (duty reported on the energy port), otherwise adiabatic.
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet"),
                Port("duty", "outlet", "energy")]

    def _reactions(self) -> list[tuple[Reaction, float]]:
        if "reactions" in self.params:
            return [(Reaction.from_param(d), float(d["conversion"]))
                    for d in self.params["reactions"]]
        return [(Reaction.from_param(self.params["reaction"]),
                 float(self.params["conversion"]))]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"ConversionReactor {self.id!r}: missing/empty inlet on 'in1'")

        _, P_in, n_in = inlet.require_state()
        z_in = inlet.normalized_z()
        moles = {c: n_in * z_in.get(c, 0.0) for c in inlet.components}

        reactions, extents = [], []
        running = dict(moles)
        for rxn, conv in self._reactions():
            xi = rxn.extent_for_conversion(running, conv)
            reactions.append(rxn)
            extents.append(xi)
            running = apply_extents(running, [rxn], [xi])  # sequential basis

        moles_out = apply_extents(moles, reactions, extents)
        P_out = P_in - float(self.params.get("dP", 0.0))
        T_spec = self.params.get("T_out")
        out, duty = reactor_outlet(
            self.id, inlet, pp, moles_out, P_out,
            None if T_spec is None else float(T_spec),
        )
        return {"out": out, "duty": EnergyStream(id=f"{self.id}.duty", duty=duty)}
