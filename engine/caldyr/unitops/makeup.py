from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Makeup")
class Makeup(UnitOp):
    """A **make-up controller**: top up one component of a process stream to a
    target molar flow by injecting a pure make-up of that component, then pass
    the combined stream on. The injected rate is computed analytically each
    solve (``design['makeup_flow']``), so it acts as an *in-loop* inventory
    controller — the robust way to close a recycle whose solvent/water inventory
    would otherwise drift (an open-steam amine regenerator loses water overhead;
    holding the circulating water with a make-up is what makes the lean-amine
    recycle converge — see `examples/24`).

    This is the analytic equivalent of a make-up *Source* on an Adjust (vary the
    make-up rate until the circulating component flow hits its target), but it is
    solved exactly per sweep rather than by an outer root find — so it does not
    destabilize the recycle the way an outer controller on an unstable inventory
    loop does.

    Ports: ``in1`` (the process stream), ``out`` (topped-up stream).

    Params:
      * ``component`` — the make-up species id (required; must be in the stream).
      * ``target`` — the desired outlet molar flow of ``component``, mol/s
        (required, >= 0). If the inlet already carries more than ``target``,
        nothing is added (a *surplus* needs a purge, not a make-up — use a
        Splitter) and ``design['makeup_flow']`` is 0.
      * ``T`` — make-up temperature, K (default: the inlet temperature).
      * ``P`` — make-up pressure, Pa (default: the inlet pressure).
    """

    design: dict | None = None

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out", "outlet")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Makeup {self.id!r}: missing or empty inlet on 'in1'")
        comp = self.params.get("component")
        if comp is None or comp not in inlet.components:
            raise ValueError(
                f"Makeup {self.id!r}: 'component' must be one of the stream "
                f"components {list(inlet.components)}; got {comp!r}"
            )
        try:
            target = float(self.params["target"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Makeup {self.id!r}: numeric 'target' (mol/s) is required "
                f"(got {self.params.get('target')!r})"
            ) from exc
        if target < 0.0:
            raise ValueError(f"Makeup {self.id!r}: target={target} must be >= 0")

        T_in, P_in, n_in = inlet.require_state()
        z_in = inlet.normalized_z()
        components = list(inlet.components)
        moles = {c: n_in * z_in.get(c, 0.0) for c in components}
        add = max(0.0, target - moles[comp])

        T_mk = float(self.params.get("T", T_in))
        P_mk = float(self.params.get("P", P_in))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z_in)
        h_mk = pp.enthalpy(T_mk, P_mk, {comp: 1.0})

        moles[comp] += add
        n_out = n_in + add
        z_out = {c: moles[c] / n_out for c in components}
        H_out = (n_in * H_in + add * h_mk) / n_out          # J/mol
        res = pp.flash_ph(P_in, H_out, z_out)
        out = Stream(
            id=f"{self.id}.out", components=components,
            T=res.T, P=res.P, molar_flow=n_out, z=z_out,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        self.design = {"makeup_flow": add, "component": comp, "target": target}
        return {"out": out}
