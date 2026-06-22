"""Liquid-liquid decanter — the settling drum of a heterogeneous (3-phase)
distillation.

A decanter receives a (condensed) two-liquid stream and lets it settle into a
**light** (organic) and a **heavy** (aqueous) liquid layer by density, the
operation that lets a heterogeneous-azeotrope distillation cross the
distillation boundary: the column overhead is the heteroazeotrope, the decanter
splits it, one layer is refluxed (kept in the column) and the other is drawn as
product. The phase split is the property package's VLLE flash
(``flash_pt_3p``), so the unit op stays free of phase-split numerics.

This is the close cousin of :class:`ThreePhaseSeparator` (a VLLE drum with a
vapor space): both flash at a specified temperature. The Decanter adds the
**reflux-drum** operation heterogeneous-azeotrope distillation needs — an
optional ``reflux_fraction`` of one settled layer is returned on a ``reflux``
port and the remainder drawn on ``product`` — and treats any flashed vapor as a
vented stream rather than a primary product (a decanter is liquid-liquid).

Three-phase (VLLE) flashes are implemented by the cubic-EOS packages
(``thermo:PR`` / ``thermo:SRK``) and the activity packages (``thermo:NRTL`` /
``thermo:UNIFAC``, via an isoactivity liquid-liquid split — UNIFAC is the
predictive simultaneous VLE+LLE model for heteroazeotropes). A fully miscible
feed degrades gracefully: everything leaves on ``liquid_light`` with an empty
``liquid_heavy`` (and ``reflux``/``product``) — not an error.
"""
from typing import Any

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("Decanter")
class Decanter(UnitOp):
    """Liquid-liquid (VLLE) decanter / settling drum.

    Params:
      * ``T`` — decanter temperature, K (**required**; a PT three-phase flash,
        like :class:`ThreePhaseSeparator` — an adiabatic PH split is not
        supported because upstream enthalpies live on the two-phase flash
        surface). The duty to hold T is reported on the ``duty`` port.
      * ``P`` — pressure, Pa (default: inlet pressure).
      * ``reflux_fraction`` — optional, in [0, 1]: the fraction of the
        ``reflux_layer`` returned on the ``reflux`` port; the rest of that layer
        leaves on ``product``. The OTHER layer always leaves on its
        ``liquid_light`` / ``liquid_heavy`` port. Omit (or 0) to take both
        layers straight off the light/heavy ports (``reflux``/``product`` empty).
      * ``reflux_layer`` — ``"light"`` (organic, default) or ``"heavy"``
        (aqueous): which settled layer feeds the reflux split.

    The split diagnostics (betas, layer densities, T) are stored on ``result``.
    """

    result: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("liquid_light", "outlet"),
            Port("liquid_heavy", "outlet"),
            Port("reflux", "outlet"),
            Port("product", "outlet"),
            Port("vapor", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(
                f"Decanter {self.id!r}: missing or empty inlet on 'in1'")
        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)
        P = float(self.params.get("P", P_in))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        if self.params.get("T") is None:
            raise ValueError(
                f"Decanter {self.id!r}: params['T'] is required — an adiabatic "
                f"(PH) three-phase flash is not supported (upstream enthalpies "
                f"live on the two-phase flash surface); specify the decanter "
                f"temperature")
        if not hasattr(pp, "flash_pt_3p"):
            raise ValueError(
                f"Decanter {self.id!r}: the property package does not implement "
                f"the three-phase (VLLE) flash — use thermo:PR / thermo:SRK / "
                f"thermo:NRTL / thermo:UNIFAC")
        layer = str(self.params.get("reflux_layer", "light"))
        if layer not in ("light", "heavy"):
            raise ValueError(
                f"Decanter {self.id!r}: reflux_layer={layer!r} must be "
                f"'light' or 'heavy'")
        rf_raw = self.params.get("reflux_fraction")
        rf = 0.0 if rf_raw is None else float(rf_raw)
        if not 0.0 <= rf <= 1.0:
            raise ValueError(
                f"Decanter {self.id!r}: reflux_fraction={rf} must be in [0, 1]")

        res = pp.flash_pt_3p(float(self.params["T"]), P, z)
        duty = n * (res.H - H_in)
        self.result = {
            "T": res.T, "P": P, "beta_vapor": res.beta_vapor,
            "beta_light": res.beta_light, "beta_heavy": res.beta_heavy,
            "rho_light": res.rho_light, "rho_heavy": res.rho_heavy,
        }

        def stream(name: str, beta: float, comp: dict[str, float] | None,
                   h: float | None, phase: str, vf: float) -> Stream:
            return Stream(
                id=f"{self.id}.{name}", components=comps,
                T=res.T, P=P, molar_flow=n * beta,
                z=dict(comp) if comp is not None else dict(z),
                H=h if h is not None else res.H, phase=phase, vapor_fraction=vf)

        out: dict[str, PortStream] = {
            "vapor": stream("vapor", res.beta_vapor, res.y, res.H_vapor,
                            "vapor", 1.0),
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }
        # The two settled liquids; the reflux_layer optionally feeds the
        # reflux/product split, the other layer leaves on its own port.
        light = ("liquid_light", res.beta_light, res.x_light, res.H_light)
        heavy = ("liquid_heavy", res.beta_heavy, res.x_heavy, res.H_heavy)
        keep, other = (light, heavy) if layer == "light" else (heavy, light)

        if rf > 0.0:
            kname, kbeta, kcomp, kh = keep
            out["reflux"] = stream("reflux", kbeta * rf, kcomp, kh, "liquid", 0.0)
            out["product"] = stream("product", kbeta * (1.0 - rf), kcomp, kh,
                                    "liquid", 0.0)
            out[kname] = stream(kname, 0.0, kcomp, kh, "liquid", 0.0)
        else:
            kname, kbeta, kcomp, kh = keep
            out[kname] = stream(kname, kbeta, kcomp, kh, "liquid", 0.0)
            out["reflux"] = stream("reflux", 0.0, kcomp, kh, "liquid", 0.0)
            out["product"] = stream("product", 0.0, kcomp, kh, "liquid", 0.0)
        oname, obeta, ocomp, oh = other
        out[oname] = stream(oname, obeta, ocomp, oh, "liquid", 0.0)
        return out
