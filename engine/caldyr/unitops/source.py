from ..core import Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Source")
class Source(UnitOp):
    """A material **source**: a boundary feed expressed as a unit op, so its
    rate/state is a *parameter* the logical ops (Set/Adjust) and the optimizer
    can drive — unlike a plain `Flowsheet.feed`, whose spec is fixed.

    The canonical use is an **adjustable makeup/injection stream**: e.g. an
    Adjust that varies ``molar_flow`` to hold a downstream spec (the amine plant's
    circulating-water makeup in `examples/24`). With no inlets, a Source is a root
    of the solve order; it contributes no equipment cost.

    Params (JSON-friendly; ``.flow`` round-trips):
      * ``molar_flow`` — total molar flow, mol/s (required, >= 0).
      * ``T`` — temperature, K (required).
      * ``P`` — pressure, Pa (required).
      * ``z`` — composition ``{component: mole_fraction}`` (required; normalized).
    """

    def define_ports(self) -> list[Port]:
        return [Port("out", "outlet")]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        try:
            n = float(self.params["molar_flow"])
            T = float(self.params["T"])
            P = float(self.params["P"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Source {self.id!r}: 'molar_flow', 'T' and 'P' are required "
                f"numeric params (got molar_flow={self.params.get('molar_flow')!r}, "
                f"T={self.params.get('T')!r}, P={self.params.get('P')!r})"
            ) from exc
        if n < 0.0 or T <= 0.0 or P <= 0.0:
            raise ValueError(
                f"Source {self.id!r}: need molar_flow >= 0 and T, P > 0 "
                f"(got molar_flow={n}, T={T}, P={P})"
            )
        z = {c: float(v) for c, v in (self.params.get("z") or {}).items()}
        total = sum(z.values())
        if total <= 0.0:
            raise ValueError(
                f"Source {self.id!r}: composition 'z' must have positive total "
                f"(got {self.params.get('z')!r})"
            )
        z = {c: v / total for c, v in z.items()}
        res = pp.flash_pt(T, P, z)
        out = Stream(
            id=f"{self.id}.out", components=list(z.keys()),
            T=res.T, P=res.P, molar_flow=n, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out": out}
