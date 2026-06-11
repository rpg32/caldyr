from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("ComponentSplitter")
class ComponentSplitter(UnitOp):
    """Black-box component separator (the pragmatic "perfect membrane"): each
    component is split between the ``overhead`` and ``bottoms`` outlets by a
    specified fraction, with no physics implied. Useful as a placeholder for a
    separation that will later be a real column/membrane/absorber.

    Params (JSON-friendly; ``.flow`` round-trips):
      * ``splits`` — ``{component_id: fraction to overhead}``; each in [0, 1].
      * ``default_split`` — fraction to overhead for components not listed in
        ``splits`` (default 0.0, i.e. unlisted components go to bottoms).
      * ``T_overhead`` / ``T_bottoms`` — outlet temperatures, K (default: keep
        the inlet temperature).
      * ``P_overhead`` / ``P_bottoms`` — outlet pressures, Pa (default: inlet
        pressure).

    Because the split is unphysical, the energy books are kept honest with a
    ``duty`` energy outlet reporting the net enthalpy change
    ``Q = sum(n_out * h_out) - n_in * h_in`` (positive = heat added), so
    flowsheet energy balances still close exactly.
    """

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("overhead", "outlet"),
            Port("bottoms", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(
                f"ComponentSplitter {self.id!r}: missing or empty inlet on 'in1'"
            )

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)

        splits = dict(self.params.get("splits", {}))
        unknown = set(splits) - set(comps)
        if unknown:
            raise ValueError(
                f"ComponentSplitter {self.id!r}: 'splits' has components not in "
                f"the flowsheet: {sorted(unknown)} (known: {comps})"
            )
        default = float(self.params.get("default_split", 0.0))
        frac = {c: float(splits.get(c, default)) for c in comps}
        bad = {c: v for c, v in frac.items() if not 0.0 <= v <= 1.0}
        if bad:
            raise ValueError(
                f"ComponentSplitter {self.id!r}: split fractions must lie in "
                f"[0, 1]; got {bad}"
            )

        f = {c: n * z.get(c, 0.0) for c in comps}            # feed rates, mol/s
        ov = {c: f[c] * frac[c] for c in comps}
        bt = {c: f[c] - ov[c] for c in comps}
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        def outlet(name: str, flows: dict[str, float]) -> tuple[Stream, float]:
            """Build one outlet; also return its total enthalpy flow n*h (W)."""
            n_out = sum(flows.values())
            t = float(self.params.get(f"T_{name}", T_in))
            p = float(self.params.get(f"P_{name}", P_in))
            # A zero-flow outlet keeps the inlet composition (its state is
            # irrelevant to the balances but should still be well-defined).
            z_out = ({c: v / n_out for c, v in flows.items()} if n_out > 0.0
                     else dict(z))
            res = pp.flash_pt(t, p, z_out)
            stream = Stream(
                id=f"{self.id}.{name}", components=comps,
                T=res.T, P=p, molar_flow=n_out, z=z_out,
                H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
            )
            return stream, n_out * res.H

        overhead, h_ov = outlet("overhead", ov)
        bottoms, h_bt = outlet("bottoms", bt)
        duty = h_ov + h_bt - n * H_in
        return {
            "overhead": overhead,
            "bottoms": bottoms,
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }
