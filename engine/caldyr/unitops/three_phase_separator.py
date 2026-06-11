"""Three-phase (vapor / light liquid / heavy liquid) separator.

A VLLE flash drum — the classic horizontal "decanter with a vapor space" used
wherever two partially-miscible liquids coexist with a gas (e.g. a
water/hydrocarbon condensate separator). The flash itself comes from the
property package's three-phase methods (``flash_pt_3p`` / ``flash_ph_3p``,
built on `thermo`'s ``FlashVLN``), so the unit op stays free of phase-split
numerics.

Only the cubic-EOS packages (``thermo:PR`` / ``thermo:SRK``) implement the
three-phase flash for now; the NRTL activity package raises a clear
``NotImplementedError``.

The two liquids are identified by **mass density**: ``liquid_light`` is the
less-dense liquid (e.g. the organic phase), ``liquid_heavy`` the denser (e.g.
the aqueous phase). A system with fewer phases degrades gracefully: a fully
miscible liquid leaves on ``liquid_light`` with an *empty* (zero-flow)
``liquid_heavy`` stream — not an error.
"""
from typing import Any

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


@register("ThreePhaseSeparator")
class ThreePhaseSeparator(UnitOp):
    """Three-phase (VLLE) separator.

    Operating spec: ``params['T']`` is **required** (a PT three-phase flash;
    ``params['P']`` defaults to the inlet pressure, like :class:`FlashDrum`),
    with the heat needed to hold T reported on the ``duty`` energy port.

    Unlike the two-phase FlashDrum there is deliberately **no adiabatic (PH)
    mode**: stream enthalpies elsewhere in the engine are computed on the
    two-phase (VL) flash surface, and feeding that enthalpy to a PH *three*-
    phase flash mixes two slightly different enthalpy surfaces — for a
    liquid-liquid feed `thermo`'s FlashVLN can then converge to unphysical
    states (negative phase fractions were observed). A missing ``T`` raises a
    clear error instead. (For the same reason, the small enthalpy-of-demixing
    difference between the VL and VLN surfaces of an LL feed shows up honestly
    in the reported duty rather than being hidden.)

    Phase-split diagnostics (betas, liquid densities, T) are stored on the
    unit's ``result`` attribute after each solve.
    """

    result: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("vapor", "outlet"),
            Port("liquid_light", "outlet"),
            Port("liquid_heavy", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(
                f"ThreePhaseSeparator {self.id!r}: missing or empty inlet on 'in1'"
            )

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)
        P = float(self.params.get("P", P_in))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        if self.params.get("T") is None:
            raise ValueError(
                f"ThreePhaseSeparator {self.id!r}: params['T'] is required — "
                f"an adiabatic (PH) three-phase flash is not supported because "
                f"upstream enthalpies live on the two-phase flash surface (see "
                f"the class docstring); specify the separator temperature"
            )
        res = pp.flash_pt_3p(float(self.params["T"]), P, z)
        duty = n * (res.H - H_in)

        def stream(name: str, beta: float, comp: dict[str, float] | None,
                   h: float | None, phase: str, vf: float) -> Stream:
            return Stream(
                id=f"{self.id}.{name}", components=comps,
                T=res.T, P=P, molar_flow=n * beta,
                z=dict(comp) if comp is not None else dict(z),
                H=h if h is not None else res.H,
                phase=phase, vapor_fraction=vf,
            )

        self.result = {
            "T": res.T, "P": P,
            "beta_vapor": res.beta_vapor, "beta_light": res.beta_light,
            "beta_heavy": res.beta_heavy,
            "rho_light": res.rho_light, "rho_heavy": res.rho_heavy,
        }
        return {
            "vapor": stream("vapor", res.beta_vapor, res.y, res.H_vapor,
                            "vapor", 1.0),
            "liquid_light": stream("liquid_light", res.beta_light, res.x_light,
                                   res.H_light, "liquid", 0.0),
            "liquid_heavy": stream("liquid_heavy", res.beta_heavy, res.x_heavy,
                                   res.H_heavy, "liquid", 0.0),
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }
