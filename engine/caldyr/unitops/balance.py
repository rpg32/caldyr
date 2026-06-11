"""Balance: the HYSYS general-purpose material/heat balance logical operation.

Reference: Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley
2025), §6.3 — "The Balance operation provides a general-purpose heat and
material balance facility"; the five balance types of its Parameters tab map
to ``params['mode']`` here. HYSYS's Balance is non-directional (it computes
whichever attached stream is unknown); in Caldyr's directional solve contract
the inlets ``in1..inN`` are the known streams and the single outlet ``out1``
is the computed unknown.

Modes — what is conserved, and what is left free:

* ``"mole"`` (book §6.3.1, *Component Mole Flow*): per-component molar flows
  are conserved (``out1`` carries the summed component moles of all inlets).
  No energy balance is made, and T/P are **not** passed — the outlet's
  T/P/H/phase are left unset, exactly as HYSYS does ("this operation does not
  pass pressure or temperature"); a downstream unit consuming it raises the
  engine's underspecified-stream error unless the state is supplied.
* ``"mass"`` (§6.3.2, *Mass Flow*): only total **mass** is conserved. The
  outlet composition must be supplied as ``params['z_out']`` (the book: "the
  composition must be specified for all streams") and the outlet molar flow is
  computed from the mass balance. Energy, moles and chemical species are NOT
  conserved — the book's use case is a reactor of unknown stoichiometry
  (alkylation units, hydrotreaters). T/P/H are not passed.
* ``"heat"`` (§6.3.3, *Heat Flow*): an overall **heat** balance only. The
  outlet ``out1`` is an *energy* port carrying the total enthalpy flow of the
  inlets, Σ nᵢ·Hᵢ (W) — "transfer the enthalpy of a process stream into a
  second energy stream". All inlets must be fully specified material streams;
  nothing material is passed. Note the engine's enthalpy basis is
  formation-inclusive, so the number is an absolute enthalpy flow.
* ``"mole_heat"`` (§6.3.4, *Component Mole and Heat Flow*): material (per
  component, molar) and energy are balanced **independently** — the outlet
  carries the summed component moles AND the summed enthalpy flow. Because a
  usable Caldyr stream needs a resolved state, the outlet is PH-flashed at
  ``params['P']`` (default: the lowest inlet pressure) to recover T/phase —
  the one place we go beyond HYSYS, which leaves T/P unset. Should not be used
  around reactors (the book's warning; with conserved moles *and* the
  formation-inclusive enthalpy basis both balances are honest here, but the
  HYSYS semantics are "not for changing species").
* ``"mass_heat"`` (§6.3.5, *Mass and Heat Flow*): total mass AND energy are
  conserved; moles and species are not. Outlet composition from
  ``params['z_out']``, molar flow from the mass balance, molar enthalpy from
  the energy balance (n_out·H_out = Σ nᵢ·Hᵢ), then PH-flashed at ``params['P']``
  (default lowest inlet P) to recover T/phase.

Ports are param-driven: ``params['n_inlets']`` (default 1) declares inlets
``in1..inN``; ``out1`` is a material outlet, or an energy outlet in ``"heat"``
mode.
"""
from __future__ import annotations

import math

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.components_db import molar_mass
from ..core.unitop import PortStream
from .base import register

MODES = ("mole", "mass", "heat", "mole_heat", "mass_heat")


@register("Balance")
class Balance(UnitOp):
    """HYSYS-style Balance logical op (see module docstring for the modes)."""

    def define_ports(self) -> list[Port]:
        n = int(self.params.get("n_inlets", 1))
        if n < 1:
            raise ValueError(f"Balance {self.id!r}: n_inlets={n} must be >= 1")
        mode = self._mode()
        out_kind = "energy" if mode == "heat" else "material"
        ports = [Port(f"in{i}", "inlet") for i in range(1, n + 1)]
        ports.append(Port("out1", "outlet", out_kind))
        return ports

    # -- helpers ---------------------------------------------------------------
    def _mode(self) -> str:
        mode = str(self.params.get("mode", "mole"))
        if mode not in MODES:
            raise ValueError(
                f"Balance {self.id!r}: unknown mode {mode!r}; expected one of {MODES}"
            )
        return mode

    def _inlet_streams(self, inlets: dict[str, Stream]) -> list[Stream]:
        streams: list[Stream] = []
        for port in self.ports:
            if port.direction != "inlet":
                continue
            s = inlets.get(port.name)
            if s is None:
                raise ValueError(
                    f"Balance {self.id!r}: missing inlet stream on {port.name!r} "
                    f"(declared n_inlets={self.params.get('n_inlets', 1)})"
                )
            streams.append(s)
        return streams

    def _z_out(self) -> dict[str, float]:
        z_out = self.params.get("z_out")
        if not isinstance(z_out, dict) or not z_out:
            raise ValueError(
                f"Balance {self.id!r}: mode {self._mode()!r} conserves only "
                f"mass/heat, so the outlet composition must be given as "
                f"params['z_out'] = {{component: mole fraction}} (the book: "
                f"'the composition must be specified for all streams')"
            )
        total = sum(float(v) for v in z_out.values())
        if total <= 0.0:
            raise ValueError(f"Balance {self.id!r}: z_out sums to {total}; expected > 0")
        return {c: float(v) / total for c, v in z_out.items()}

    @staticmethod
    def _mole_sum(streams: list[Stream]) -> tuple[float, dict[str, float]]:
        """Total molar flow and per-component molar flows over ``streams``.
        Mole-balance modes need only flow + composition (HYSYS: "the direction
        of flow ... is of no consequence"; T/P play no part)."""
        moles: dict[str, float] = {}
        n_total = 0.0
        for s in streams:
            if s.molar_flow is None:
                raise ValueError(f"Balance: inlet {s.id!r} has no molar_flow")
            zi = s.normalized_z()
            n_total += s.molar_flow
            for c, frac in zi.items():
                moles[c] = moles.get(c, 0.0) + s.molar_flow * frac
        return n_total, moles

    @staticmethod
    def _mass_flow(streams: list[Stream]) -> float:
        """Total mass flow (kg/s) over ``streams``."""
        total = 0.0
        for s in streams:
            if s.molar_flow is None:
                raise ValueError(f"Balance: inlet {s.id!r} has no molar_flow")
            total += s.molar_flow * sum(
                f * molar_mass(c) for c, f in s.normalized_z().items()
            )
        return total

    @staticmethod
    def _enthalpy_flow(streams: list[Stream], pp) -> float:
        """Total enthalpy flow Σ nᵢ·Hᵢ (W) over fully-specified inlets."""
        total = 0.0
        for s in streams:
            T, P, n = s.require_state()
            h = s.H if s.H is not None else pp.enthalpy(T, P, s.normalized_z())
            total += n * h
        return total

    # -- solve -----------------------------------------------------------------
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        mode = self._mode()
        streams = self._inlet_streams(inlets)
        components = list(streams[0].components)

        if mode == "heat":
            duty = self._enthalpy_flow(streams, pp)
            return {"out1": EnergyStream(id=f"{self.id}.out1", duty=duty)}

        if mode == "mole":
            n_total, moles = self._mole_sum(streams)
            z = {c: moles.get(c, 0.0) / n_total for c in components}
            # Component moles conserved; T/P/H deliberately NOT passed (HYSYS
            # §6.3.1) — the outlet state is free until something specifies it.
            return {"out1": Stream(id=f"{self.id}.out1", components=components,
                                   T=None, P=None, molar_flow=n_total, z=z)}

        if mode == "mass":
            z_out = self._z_out()
            mw_out = sum(f * molar_mass(c) for c, f in z_out.items())
            n_out = self._mass_flow(streams) / mw_out
            # Only total mass survives; moles/species/energy are not conserved
            # and T/P/H are not passed (HYSYS §6.3.2).
            return {"out1": Stream(id=f"{self.id}.out1", components=components,
                                   T=None, P=None, molar_flow=n_out, z=z_out)}

        # The two combined modes: material balance + independent energy balance,
        # then a PH flash to give the outlet a usable resolved state.
        h_flow = self._enthalpy_flow(streams, pp)
        p_out = float(self.params.get("P") or min(s.P for s in streams if s.P is not None))
        if mode == "mole_heat":
            n_out, moles = self._mole_sum(streams)
            z = {c: moles.get(c, 0.0) / n_out for c in components}
        else:                                                       # mass_heat
            z = self._z_out()
            n_out = self._mass_flow(streams) / sum(f * molar_mass(c) for c, f in z.items())

        h_out = h_flow / n_out                                      # J/mol
        if not math.isfinite(h_out):
            raise ValueError(f"Balance {self.id!r}: non-finite outlet enthalpy {h_out}")
        res = pp.flash_ph(p_out, h_out, z)
        out = Stream(
            id=f"{self.id}.out1", components=components,
            T=res.T, P=res.P, molar_flow=n_out, z=z,
            H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
        )
        return {"out1": out}
