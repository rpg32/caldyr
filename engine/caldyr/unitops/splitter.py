from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Splitter")
class Splitter(UnitOp):
    """Flow splitter: divide one inlet into two outlets with identical intensive
    state (T, P, composition, phase). Only the molar flow is partitioned.

    ``params['split']`` is the fraction of the inlet flow sent to ``out1`` (the
    rest goes to ``out2``); it must lie in [0, 1].
    """

    def define_ports(self) -> list[Port]:
        return [Port("in1", "inlet"), Port("out1", "outlet"), Port("out2", "outlet")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or inlet.molar_flow is None:
            raise ValueError(f"Splitter {self.id!r}: missing inlet on 'in1'")

        split = float(self.params.get("split", 0.5))
        if not 0.0 <= split <= 1.0:
            raise ValueError(f"Splitter {self.id!r}: split={split} must be in [0, 1]")

        _, _, n = inlet.require_state()

        def branch(suffix: str, frac: float) -> Stream:
            return Stream(
                id=f"{self.id}.{suffix}",
                components=list(inlet.components),
                T=inlet.T, P=inlet.P, molar_flow=n * frac, z=dict(inlet.z),
                H=inlet.H, phase=inlet.phase, vapor_fraction=inlet.vapor_fraction,
            )

        return {"out1": branch("out1", split), "out2": branch("out2", 1.0 - split)}
