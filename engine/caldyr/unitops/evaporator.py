"""Evaporator: a heated flash — vaporize part of a liquid feed at pressure P.

Reference: Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley
2025), §5.2 "Simulation of Evaporator" — there modeled as a heat exchanger
(heated by condensing saturated steam) followed by a flash separator; here the
two are folded into one unit: a liquid inlet, an operating pressure, ONE
thermal specification, and vapor/liquid outlets plus the heating duty on an
energy port.

Specifications (``params``) — give ``P`` (Pa; defaults to the inlet pressure)
plus exactly one of:

* ``T`` — outlet temperature, K: a PT flash; duty follows from the energy
  balance.
* ``duty`` — heat input Q, W: a PH flash at H_in + Q/n.
* ``vapor_fraction`` — molar vapor fraction of the combined outlet (0–1): the
  flash temperature/enthalpy is found by root-finding. For a multicomponent
  feed the unit brackets T between the bubble and dew points and solves
  ``VF(T) = spec`` with ``scipy.optimize.brentq``; for a (near-)pure fluid the
  isobaric two-phase region collapses to a single temperature (T_bubble =
  T_dew = Tsat), so T is not a usable unknown — the target enthalpy is taken
  directly from the saturated-phase enthalpies, H = (1-VF)·h_liq + VF·h_vap,
  and resolved with a PH flash.

The duty reported on the ``duty`` port is the net heat input
Q = n·(H_out - H_in) (W, positive when evaporating). Economics: sized as a
vertical vessel whose heating duty draws a hot utility (see
:mod:`caldyr.economics.sizing`).
"""
from __future__ import annotations

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register

# Bubble/dew spread below which the feed is treated as a pure fluid (K).
_PURE_DT = 1e-6


@register("Evaporator")
class Evaporator(UnitOp):
    """Heated flash drum (book §5.2). See module docstring for the specs."""

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("vapor", "outlet"),
            Port("liquid", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def _spec(self) -> tuple[str, float]:
        given: dict[str, float] = {}
        for k in ("T", "duty", "vapor_fraction"):
            v = self.params.get(k)
            if v is not None:
                given[k] = float(v)
        if len(given) != 1:
            raise ValueError(
                f"Evaporator {self.id!r}: specify exactly one of 'T', 'duty' or "
                f"'vapor_fraction' (got {sorted(given) or 'none'})"
            )
        ((key, value),) = given.items()
        if key == "vapor_fraction" and not 0.0 <= value <= 1.0:
            raise ValueError(
                f"Evaporator {self.id!r}: vapor_fraction={value} must be in [0, 1]"
            )
        return key, value

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"Evaporator {self.id!r}: missing or empty inlet on 'in1'")

        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        P = float(self.params.get("P", P_in))
        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)
        key, value = self._spec()

        if key == "T":
            res = pp.flash_pt(value, P, z)
        elif key == "duty":
            res = pp.flash_ph(P, H_in + value / n, z)
        else:
            res = self._flash_vf(pp, P, z, value)
        duty = n * (res.H - H_in)

        vf = res.vapor_fraction
        y = res.y if res.y is not None else z
        x = res.x if res.x is not None else z
        h_vap = res.H_vapor if res.H_vapor is not None else res.H
        h_liq = res.H_liquid if res.H_liquid is not None else res.H

        vapor = Stream(
            id=f"{self.id}.vapor", components=list(inlet.components),
            T=res.T, P=P, molar_flow=n * vf, z=dict(y),
            H=h_vap, phase="vapor", vapor_fraction=1.0,
        )
        liquid = Stream(
            id=f"{self.id}.liquid", components=list(inlet.components),
            T=res.T, P=P, molar_flow=n * (1.0 - vf), z=dict(x),
            H=h_liq, phase="liquid", vapor_fraction=0.0,
        )
        return {
            "vapor": vapor,
            "liquid": liquid,
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }

    def _flash_vf(self, pp, P: float, z: dict[str, float], vf: float):
        """Resolve the vapor-fraction spec (see module docstring)."""
        t_bub, t_dew = pp.bubble_dew(P, z)

        if t_dew - t_bub < _PURE_DT:
            # Pure fluid: T is pinned at Tsat for any 0 < VF < 1, so solve in
            # enthalpy instead — exact, no iteration needed.
            sat = pp.bubble_point(P, z)
            if sat.H_liquid is None or sat.H_vapor is None:
                raise ValueError(
                    f"Evaporator {self.id!r}: property package returned no "
                    f"saturated phase enthalpies at P={P:.4g} Pa"
                )
            h_target = (1.0 - vf) * sat.H_liquid + vf * sat.H_vapor
            return pp.flash_ph(P, h_target, z)

        from scipy.optimize import brentq

        def f(t: float) -> float:
            return pp.flash_pt(float(t), P, z).vapor_fraction - vf

        # VF rises monotonically from 0 at the bubble point to 1 at the dew
        # point, so [t_bub, t_dew] brackets every spec in [0, 1].
        t_flash = float(brentq(f, t_bub, t_dew, xtol=1e-8))
        return pp.flash_pt(t_flash, P, z)
