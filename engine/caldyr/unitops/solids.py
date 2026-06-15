"""Solids operations: gas cyclone, rotary vacuum filter, baghouse (fabric) filter.

The HYSYS "solid operations" palette of Hameed, *Chemical Process Simulations
using Aspen HYSYS* (Wiley 2025), ch. 12: the cyclone (sec. 12.1), the rotary
vacuum filter (sec. 12.2) and the baghouse filter (sec. 12.3).

The v1 particle model — solids without a stream-model change
------------------------------------------------------------
Caldyr streams carry molar component flows only (no particle attributes), and
``core/`` is deliberately untouched. Solids are therefore modeled pragmatically:

* the **solid is an ordinary component** of the stream (carbon dust in air is
  simply the component ``"carbon"``);
* the **particle size distribution lives on the unit**, as the ``psd`` param —
  a list of ``{"d_microns": ..., "mass_frac": ...}`` bins, exactly how the
  book's worked example enters it (Hameed sec. 12.1.2 step 9);
* the unit computes a per-bin (grade) efficiency, collapses it to an overall
  collection efficiency with the bin mass fractions, and splits the designated
  solid component(s) between its outlets by that overall efficiency.

Honest limitations of this v1: the PSD does **not** propagate between units
(an outlet stream has no particle data), so chained solids operations each
re-specify their own PSD — physically the downstream PSD is finer, since the
cyclone preferentially removes the coarse bins (the per-bin splits are
reported on ``unit.design`` so a user can hand-chain them); all designated
solid components share one PSD; and outlet enthalpies come from flashing the
solid component on the flowsheet EOS — adequate for closing the energy books
across an isothermal split, not a real solid-phase thermodynamic model. Like
the ComponentSplitter, every unit carries a ``duty`` energy outlet reporting
the (small) net enthalpy change of the split so flowsheet energy balances
still close exactly.

Physics sources
---------------
* **Cyclone** — Lapple's cut-diameter model: ``d50 = sqrt(9 mu W / (2 pi N_e
  v_i (rho_p - rho_g)))`` with grade efficiency ``eta_j = 1 / (1 +
  (d50/d_j)^2)`` and effective turns ``N_e = (L_b + L_c/2) / H`` (Cooper &
  Alley, *Air Pollution Control: A Design Approach*, 4e, ch. 4, Eqs. 4.4-4.9;
  Perry's 8e sec. 17). Pressure drop by Shepherd-Lapple velocity heads
  ``dP = 0.5 rho_g v_i^2 N_H`` with ``N_H = 16 H W / D_e^2`` (Cooper & Alley
  Eqs. 4.10-4.11; the book's Eqs. 12.2-12.3 are the same form). Standard
  geometry families from Cooper & Alley Table 4.2 (Lapple 1951, Stairmand
  1951); the book's HYSYS "High Efficiency" / "High Output" sizing ratios
  (Figs. 12.5-12.6) match Stairmand's two designs and are included.
* **Rotary vacuum filter** — constant-pressure cake filtration with
  negligible medium resistance, applied to a continuous rotary drum: the
  cycle-averaged filtrate flux is ``Q/A = sqrt(2 c dP f / (alpha mu t_c))``
  where ``f`` is the submerged fraction of the cycle (McCabe, Smith &
  Harriott, *Unit Operations of Chemical Engineering*, 7e, ch. 29, continuous
  filtration; Perry's 8e sec. 18; the book's Eqs. 12.8-12.10 are the same
  equation).
* **Baghouse filter** — cloth area from the gross air-to-cloth ratio
  ``A = Q / v_f`` (Cooper & Alley ch. 6, Table 6.1; EPA *Air Pollution
  Control Cost Manual*, 6e, EPA/452/B-02-001, sec. 6 ch. 1) and pressure drop
  from the classic filter-drag model ``dP = S_E v_f + K2 c_i v_f^2 t``
  (Cooper & Alley ch. 6; Billings & Wilder via EPA manual), giving the
  filtration time to a maximum allowed dP.

Gas/liquid viscosity comes from the same thermo/chemicals databank path the
PipeSegment uses (:func:`caldyr.unitops.pipe.mixture_viscosity`), evaluated on
the carrier (non-solid) components only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.components_db import molar_mass
from ..core.unitop import PortStream
from .base import register
from .pipe import mixture_viscosity


class SolidsOperationError(ValueError):
    """A solids unit op cannot model the requested separation (bad PSD or
    solids spec, impossible pressure drop, no carrier phase, ...)."""


# -- shared helpers -----------------------------------------------------------
def _solid_ids(label: str, params: dict, components: list[str]) -> list[str]:
    """The designated solid component id(s): ``params['solids']`` as a string
    or list of strings, all of which must exist in the flowsheet."""
    raw = params.get("solids")
    ids = [raw] if isinstance(raw, str) else list(raw or [])
    if not ids:
        raise SolidsOperationError(
            f"{label}: parameter 'solids' is required — the component id (or "
            f"list of ids) of the solid(s) to separate, e.g. 'carbon'"
        )
    unknown = [s for s in ids if s not in components]
    if unknown:
        raise SolidsOperationError(
            f"{label}: 'solids' lists components not in the flowsheet: "
            f"{unknown} (known: {components})"
        )
    return ids


def _parse_psd(label: str, raw: Any) -> list[tuple[float, float]]:
    """Validate and normalize a PSD param into ``[(d_m, mass_frac), ...]``
    sorted by size. Bins are ``{"d_microns": d, "mass_frac": w}`` dicts; the
    mass fractions must sum to 1 within 2% (renormalized exactly)."""
    if not isinstance(raw, list) or not raw:
        raise SolidsOperationError(
            f"{label}: parameter 'psd' is required — a non-empty list of "
            f"{{'d_microns': ..., 'mass_frac': ...}} bins (a single-particle "
            f"diameter is a one-bin psd)"
        )
    bins: list[tuple[float, float]] = []
    for i, b in enumerate(raw):
        try:
            d = float(b["d_microns"])
            w = float(b["mass_frac"])
        except (TypeError, KeyError) as exc:
            raise SolidsOperationError(
                f"{label}: psd bin {i} must be a dict with 'd_microns' and "
                f"'mass_frac' keys; got {b!r}"
            ) from exc
        if d <= 0.0 or not math.isfinite(d):
            raise SolidsOperationError(
                f"{label}: psd bin {i} has non-positive d_microns ({d})")
        if w < 0.0:
            raise SolidsOperationError(
                f"{label}: psd bin {i} has negative mass_frac ({w})")
        bins.append((d * 1e-6, w))
    total = sum(w for _, w in bins)
    if not 0.98 <= total <= 1.02:
        raise SolidsOperationError(
            f"{label}: psd mass fractions sum to {total:.4f}; they must sum "
            f"to 1 (within 2%)"
        )
    return sorted(((d, w / total) for d, w in bins), key=lambda b: b[0])


def _carrier_state(label: str, inlet: Stream, solid_ids: list[str], pp,
                   phase: str) -> dict[str, Any]:
    """Bulk properties of the carrier (non-solid) phase at inlet conditions:
    normalized composition, molar/volumetric flow, density, viscosity. The
    solid's volume contribution is neglected (dilute-dust/slurry carrier)."""
    T, P, n = inlet.require_state()
    z = inlet.normalized_z()
    z_raw = {c: z.get(c, 0.0) for c in inlet.components
             if c not in solid_ids and z.get(c, 0.0) > 0.0}
    frac = sum(z_raw.values())
    if frac <= 0.0:
        kind = "gas" if phase == "vapor" else "liquid"
        raise SolidsOperationError(
            f"{label}: the inlet carries no {kind} — every component with "
            f"positive flow is a designated solid; a carrier {kind} phase is "
            f"required"
        )
    z_c = {c: v / frac for c, v in z_raw.items()}
    comps = tuple(z_c)
    zs = [z_c[c] for c in comps]
    vm = pp.volume(T, P, z_c)                      # m^3/mol at T, P
    if phase == "liquid":
        vol_liq = getattr(pp, "volume_liquid", None)
        if vol_liq is not None:
            vm = vol_liq(T, P, z_c)
    mw = sum(zi * molar_mass(c) for c, zi in zip(comps, zs))   # kg/mol
    n_carrier = n * frac
    return {
        "z": z_c, "n": n_carrier, "Vdot": n_carrier * vm,
        "rho": mw / vm, "mu": mixture_viscosity(comps, zs, T, phase),
    }


def _solids_mass_rate(inlet: Stream, solid_ids: list[str]) -> float:
    """Total mass rate (kg/s) of the designated solid components in a stream."""
    z = inlet.normalized_z()
    _, _, n = inlet.require_state()
    return sum(n * z.get(s, 0.0) * molar_mass(s) for s in solid_ids)


def _outlet(unit_id: str, name: str, comps: list[str], flows: dict[str, float],
            t: float, p: float, z_fallback: dict[str, float],
            pp) -> tuple[Stream, float]:
    """Build one outlet stream from per-component molar flows (mol/s), flashed
    at (t, p); also return its total enthalpy flow n*h (W) for the duty books.
    A zero-flow outlet keeps the fallback composition (state well-defined,
    irrelevant to the balances)."""
    n_out = sum(flows.values())
    z_out = ({c: v / n_out for c, v in flows.items()} if n_out > 0.0
             else dict(z_fallback))
    res = pp.flash_pt(t, p, z_out)
    stream = Stream(
        id=f"{unit_id}.{name}", components=comps,
        T=res.T, P=p, molar_flow=n_out, z=z_out,
        H=res.H, phase=res.phase, vapor_fraction=res.vapor_fraction,
    )
    return stream, n_out * res.H


# -- cyclone ------------------------------------------------------------------
@dataclass(frozen=True)
class CycloneGeometry:
    """Standard cyclone proportions as fractions of the body diameter D
    (Cooper & Alley 4e Table 4.2). ``Ne`` is Lapple's number of effective
    turns, ``NH`` the Shepherd-Lapple inlet velocity heads (K = 16, standard
    tangential inlet)."""
    name: str
    H: float        # inlet height / D
    W: float        # inlet width / D
    De: float       # gas outlet diameter / D
    S: float        # gas outlet (vortex finder) length / D
    Lb: float       # body (cylinder) length / D
    Lc: float       # cone length / D
    Dd: float       # solids outlet diameter / D

    @property
    def Ne(self) -> float:
        """Effective turns N_e = (L_b + L_c/2) / H (Cooper & Alley Eq. 4.4)."""
        return (self.Lb + self.Lc / 2.0) / self.H

    @property
    def NH(self) -> float:
        """Inlet velocity heads N_H = 16 H W / De^2 (Cooper & Alley Eq. 4.11,
        Shepherd-Lapple K = 16 for a standard tangential inlet)."""
        return 16.0 * self.H * self.W / (self.De * self.De)


#: Standard geometry families (Cooper & Alley 4e Table 4.2). "Stairmand_HE" /
#: "Stairmand_HT" are the book's HYSYS "High Efficiency" / "High Output"
#: sizing ratios (Hameed Figs. 12.5-12.6: inlet width 0.2/0.375, inlet height
#: 0.5/0.75, gas outlet 0.5/0.75, solids outlet 0.375).
CYCLONE_GEOMETRIES: dict[str, CycloneGeometry] = {
    "Lapple": CycloneGeometry("Lapple", H=0.5, W=0.25, De=0.5, S=0.625,
                              Lb=2.0, Lc=2.0, Dd=0.25),
    "Stairmand_HE": CycloneGeometry("Stairmand_HE", H=0.5, W=0.2, De=0.5,
                                    S=0.5, Lb=1.5, Lc=2.5, Dd=0.375),
    "Stairmand_HT": CycloneGeometry("Stairmand_HT", H=0.75, W=0.375, De=0.75,
                                    S=0.875, Lb=1.5, Lc=2.5, Dd=0.375),
}


def lapple_d50(mu: float, w_inlet: float, n_e: float, v_i: float,
               rho_p: float, rho_g: float) -> float:
    """Lapple cut diameter (m): the particle size collected with 50%
    efficiency, ``d50 = sqrt(9 mu W / (2 pi N_e v_i (rho_p - rho_g)))``
    (Cooper & Alley 4e Eq. 4.7; Perry's 8e sec. 17)."""
    if rho_p <= rho_g:
        raise SolidsOperationError(
            f"particle density {rho_p:.1f} kg/m^3 must exceed the gas density "
            f"{rho_g:.2f} kg/m^3 for centrifugal collection"
        )
    if v_i <= 0.0:
        raise SolidsOperationError(f"inlet velocity must be positive, got {v_i}")
    return math.sqrt(9.0 * mu * w_inlet / (2.0 * math.pi * n_e * v_i * (rho_p - rho_g)))


def lapple_grade_efficiency(d50: float, d: float) -> float:
    """Lapple's grade-efficiency curve ``eta = 1 / (1 + (d50/d)^2)``
    (Cooper & Alley 4e Eq. 4.9, Theodore-DePaola algebraic fit)."""
    return 1.0 / (1.0 + (d50 / d) ** 2)


@register("Cyclone")
class Cyclone(UnitOp):
    """Gas cyclone: split solid dust components from a carrier gas with
    Lapple cut-diameter physics (see module docstring for the model and the
    v1 PSD-as-param particle representation).

    Ports: ``gas_in`` -> ``gas_out`` + ``solids_out`` (+ ``duty`` energy
    outlet closing the enthalpy books of the split).

    Parameters
    ----------
    solids : component id (str) or list — the dust component(s). Required.
    psd : list of ``{"d_microns", "mass_frac"}`` bins. Required.
    particle_density : kg/m^3 of the solid. Required (Lapple needs rho_p).
    geometry : "Lapple" (default) | "Stairmand_HE" | "Stairmand_HT" —
        standard proportion family (Cooper & Alley Table 4.2; the Stairmand
        pair matches the book's HYSYS High Efficiency / High Output ratios).
    body_diameter : m — size the cyclone(s) directly, with ``n_cyclones``
        identical units in parallel (default 1); **or**
    inlet_velocity : m/s — size the body diameter from a design inlet
        velocity (classic design point 15-27 m/s, ~18 m/s optimum; Cooper &
        Alley ch. 4) given ``n_cyclones``: D = sqrt(Q_each / (v H_r W_r)).
        Exactly one of ``body_diameter`` / ``inlet_velocity`` is required.
    n_cyclones : parallel identical cyclones sharing the flow (default 1).

    After a solve, ``unit.design`` carries d50 (m and microns), the per-bin
    grade efficiencies, overall_efficiency, dP_Pa (Shepherd-Lapple),
    inlet_velocity_m_s, body_diameter_m, n_cyclones, Ne/NH, gas density &
    viscosity, volumetric flow and the solids capture/emission mass rates.
    """

    design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("gas_in", "inlet"),
            Port("gas_out", "outlet"),
            Port("solids_out", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        label = f"Cyclone {self.id!r}"
        inlet = inlets.get("gas_in")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"{label}: missing or empty inlet on 'gas_in'")
        t_in, p_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)

        solids = _solid_ids(label, self.params, comps)
        psd = _parse_psd(label, self.params.get("psd"))
        rho_p = self.params.get("particle_density")
        if rho_p is None or float(rho_p) <= 0.0:
            raise SolidsOperationError(
                f"{label}: parameter 'particle_density' (kg/m^3, > 0) is "
                f"required — the Lapple cut diameter scales with rho_p"
            )
        rho_p = float(rho_p)
        geo_name = str(self.params.get("geometry", "Lapple"))
        geo = CYCLONE_GEOMETRIES.get(geo_name)
        if geo is None:
            raise SolidsOperationError(
                f"{label}: unknown geometry {geo_name!r}; available: "
                f"{sorted(CYCLONE_GEOMETRIES)}"
            )
        n_cyc = int(self.params.get("n_cyclones", 1))
        if n_cyc < 1:
            raise SolidsOperationError(f"{label}: n_cyclones must be >= 1, got {n_cyc}")

        gas = _carrier_state(label, inlet, solids, pp, "vapor")
        q_each = gas["Vdot"] / n_cyc

        # Sizing: explicit body diameter, or from a target inlet velocity.
        d_body = self.params.get("body_diameter")
        v_target = self.params.get("inlet_velocity")
        if d_body is not None:
            d_body = float(d_body)
            if d_body <= 0.0:
                raise SolidsOperationError(
                    f"{label}: body_diameter must be positive, got {d_body}")
        elif v_target is not None:
            v_target = float(v_target)
            if v_target <= 0.0:
                raise SolidsOperationError(
                    f"{label}: inlet_velocity must be positive, got {v_target}")
            d_body = math.sqrt(q_each / (v_target * geo.H * geo.W))
        else:
            raise SolidsOperationError(
                f"{label}: specify either 'body_diameter' (m, with "
                f"'n_cyclones') or a design 'inlet_velocity' (m/s) to size "
                f"the cyclone"
            )
        w_in, h_in = geo.W * d_body, geo.H * d_body
        v_i = q_each / (h_in * w_in)

        d50 = lapple_d50(gas["mu"], w_in, geo.Ne, v_i, rho_p, gas["rho"])
        grade = [{"d_microns": d * 1e6, "mass_frac": w,
                  "efficiency": lapple_grade_efficiency(d50, d)}
                 for d, w in psd]
        eta = sum(g["mass_frac"] * g["efficiency"] for g in grade)

        dp = 0.5 * gas["rho"] * v_i * v_i * geo.NH      # Shepherd-Lapple
        p_out = p_in - dp
        if p_out <= 0.0:
            raise SolidsOperationError(
                f"{label}: Shepherd-Lapple pressure drop {dp:.4g} Pa exceeds "
                f"the inlet pressure {p_in:.4g} Pa (inlet velocity "
                f"{v_i:.1f} m/s, N_H={geo.NH:.1f}) — slow the cyclone down "
                f"(larger body_diameter or more n_cyclones)"
            )

        feed = {c: n * z.get(c, 0.0) for c in comps}
        captured = {c: (feed[c] * eta if c in solids else 0.0) for c in comps}
        to_gas = {c: feed[c] - captured[c] for c in comps}
        m_in = _solids_mass_rate(inlet, solids)

        h_in_total = n * (inlet.H if inlet.H is not None
                          else pp.enthalpy(t_in, p_in, z))
        gas_out, h_gas = _outlet(self.id, "gas_out", comps, to_gas,
                                 t_in, p_out, z, pp)
        sol_out, h_sol = _outlet(self.id, "solids_out", comps, captured,
                                 t_in, p_out, z, pp)

        self.design = {
            "geometry": geo.name,
            "body_diameter_m": d_body,
            "n_cyclones": n_cyc,
            "inlet_width_m": w_in,
            "inlet_height_m": h_in,
            "inlet_velocity_m_s": v_i,
            "Q_m3_s": gas["Vdot"],
            "Q_per_cyclone_m3_s": q_each,
            "gas_density_kg_m3": gas["rho"],
            "gas_viscosity_Pa_s": gas["mu"],
            "particle_density_kg_m3": rho_p,
            "Ne_turns": geo.Ne,
            "NH_velocity_heads": geo.NH,
            "d50_m": d50,
            "d50_microns": d50 * 1e6,
            "grade": grade,
            "overall_efficiency": eta,
            "dP_Pa": dp,
            "dust_loading_kg_m3": m_in / gas["Vdot"],
            "solids_in_kg_s": m_in,
            "solids_captured_kg_s": m_in * eta,
            "solids_emitted_kg_s": m_in * (1.0 - eta),
        }
        return {
            "gas_out": gas_out,
            "solids_out": sol_out,
            "duty": EnergyStream(id=f"{self.id}.duty",
                                 duty=h_gas + h_sol - h_in_total),
        }


# -- rotary vacuum filter -----------------------------------------------------
@register("RotaryVacuumFilter")
class RotaryVacuumFilter(UnitOp):
    """Continuous rotary-drum vacuum filter: split a slurry into a wet cake
    and a filtrate, and size the drum from constant-pressure cake-filtration
    theory (see module docstring; McCabe-Smith-Harriott 7e ch. 29).

    Ports: ``slurry_in`` -> ``filtrate_out`` + ``cake_out`` (+ ``duty``).

    Split model (the book's HYSYS unit assumes 100% solids removal and
    computes solvent retention in the cake — Hameed sec. 12.2):

    * ``solids_capture`` (default 1.0, the book's assumption) of each solid
      component reports to the cake;
    * ``cake_moisture`` (default 0.5 — the 50%-moisture wet-basis cake of the
      book's Exercise 12.2 / McCabe 7e ch. 29 worked example) is the liquid
      mass fraction of the wet cake; the corresponding liquid is taken from
      the carrier at the filtrate composition.

    Sizing parameters (all required unless noted): ``pressure_drop`` (Pa
    across the cake), ``cycle_time_s`` (one drum revolution),
    ``submergence`` (fraction of the cycle the drum is submerged, default
    0.20 — the book's sec. 12.2 example; McCabe's example uses 0.30),
    ``alpha`` (specific cake resistance, m/kg — no honest default exists;
    typical incompressible CaCO3 cakes run ~1e10-1e12 m/kg, Perry's 8e
    sec. 18). ``filtrate_viscosity`` (Pa s) overrides the databank value.
    ``drum_radius_m`` (optional): report the drum width for a given radius,
    as HYSYS does (area = 2 pi R W); otherwise a square drum (W = 2R) is
    suggested.

    ``unit.design``: filtration area required (m^2), cycle-averaged filtrate
    flux, solids concentration c (kg solids per m^3 filtrate), drum
    radius/width suggestion, cake and filtrate rates.
    """

    design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("slurry_in", "inlet"),
            Port("filtrate_out", "outlet"),
            Port("cake_out", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def _required(self, label: str, name: str, why: str) -> float:
        value = self.params.get(name)
        if value is None or float(value) <= 0.0:
            raise SolidsOperationError(
                f"{label}: parameter {name!r} (> 0) is required — {why}")
        return float(value)

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        label = f"RotaryVacuumFilter {self.id!r}"
        inlet = inlets.get("slurry_in")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"{label}: missing or empty inlet on 'slurry_in'")
        t_in, p_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)

        solids = _solid_ids(label, self.params, comps)
        capture = float(self.params.get("solids_capture", 1.0))
        if not 0.0 < capture <= 1.0:
            raise SolidsOperationError(
                f"{label}: solids_capture must lie in (0, 1], got {capture}")
        moisture = float(self.params.get("cake_moisture", 0.5))
        if not 0.0 <= moisture < 1.0:
            raise SolidsOperationError(
                f"{label}: cake_moisture (wet-basis liquid mass fraction) "
                f"must lie in [0, 1), got {moisture}"
            )
        dp = self._required(label, "pressure_drop",
                            "the constant-pressure driving force across the cake (Pa)")
        t_cycle = self._required(label, "cycle_time_s", "one drum revolution (s)")
        alpha = self._required(
            label, "alpha",
            "the specific cake resistance (m/kg); measure it or take a "
            "literature value (Perry's 8e sec. 18)")
        f_sub = float(self.params.get("submergence", 0.20))
        if not 0.0 < f_sub < 1.0:
            raise SolidsOperationError(
                f"{label}: submergence must lie in (0, 1), got {f_sub}")
        if dp >= p_in:
            raise SolidsOperationError(
                f"{label}: pressure_drop {dp:.4g} Pa >= inlet pressure "
                f"{p_in:.4g} Pa — the vacuum side cannot go below absolute zero"
            )

        liq = _carrier_state(label, inlet, solids, pp, "liquid")
        m_solids_in = _solids_mass_rate(inlet, solids)
        if m_solids_in <= 0.0:
            raise SolidsOperationError(
                f"{label}: the inlet carries none of the designated solid(s) "
                f"{solids} — nothing to filter"
            )

        # Split: solids to cake by capture; cake liquid from the moisture spec.
        m_cake_solids = capture * m_solids_in
        m_cake_liquid = m_cake_solids * moisture / (1.0 - moisture)
        m_liq_total = sum(n * z.get(c, 0.0) * molar_mass(c)
                          for c in comps if c not in solids)
        if m_cake_liquid >= m_liq_total:
            raise SolidsOperationError(
                f"{label}: the cake moisture spec needs {m_cake_liquid:.4g} "
                f"kg/s of liquid but the slurry only carries "
                f"{m_liq_total:.4g} kg/s — reduce cake_moisture ({moisture})"
            )
        phi = m_cake_liquid / m_liq_total      # mole fraction of each liquid taken
        feed = {c: n * z.get(c, 0.0) for c in comps}
        cake = {c: feed[c] * (capture if c in solids else phi) for c in comps}
        filtrate = {c: feed[c] - cake[c] for c in comps}

        # Continuous-filtration sizing (McCabe 7e ch. 29, negligible medium
        # resistance): cycle-averaged flux = sqrt(2 c dP f / (alpha mu t_c)).
        n_liq_filtrate = sum(v for c, v in filtrate.items() if c not in solids)
        q_filtrate = n_liq_filtrate * liq["Vdot"] / liq["n"]   # m^3/s
        mu = float(self.params.get("filtrate_viscosity", liq["mu"]))
        c_solids = m_cake_solids / q_filtrate                  # kg/m^3 filtrate
        flux = math.sqrt(2.0 * c_solids * dp * f_sub / (alpha * mu * t_cycle))
        area = q_filtrate / flux

        # Drum suggestion: HYSYS asks for the radius and reports the width
        # (area = 2 pi R W, the full drum face); default to a "square" drum.
        radius = self.params.get("drum_radius_m")
        radius = float(radius) if radius is not None else math.sqrt(area / (4.0 * math.pi))
        width = area / (2.0 * math.pi * radius)

        h_in_total = n * (inlet.H if inlet.H is not None
                          else pp.enthalpy(t_in, p_in, z))
        filt_out, h_f = _outlet(self.id, "filtrate_out", comps, filtrate,
                                t_in, p_in - dp, z, pp)
        cake_out, h_c = _outlet(self.id, "cake_out", comps, cake,
                                t_in, p_in, z, pp)

        self.design = {
            "area_m2": area,
            "flux_m3_m2_s": flux,
            "c_solids_kg_m3_filtrate": c_solids,
            "alpha_m_kg": alpha,
            "filtrate_viscosity_Pa_s": mu,
            "pressure_drop_Pa": dp,
            "cycle_time_s": t_cycle,
            "submergence": f_sub,
            "drum_radius_m": radius,
            "drum_width_m": width,
            "filtrate_m3_s": q_filtrate,
            "cake_solids_kg_s": m_cake_solids,
            "cake_liquid_kg_s": m_cake_liquid,
            "cake_moisture": moisture,
            "solids_capture": capture,
        }
        return {
            "filtrate_out": filt_out,
            "cake_out": cake_out,
            "duty": EnergyStream(id=f"{self.id}.duty",
                                 duty=h_f + h_c - h_in_total),
        }


# -- baghouse filter ----------------------------------------------------------
@register("BaghouseFilter")
class BaghouseFilter(UnitOp):
    """Baghouse (fabric) filter: near-total dust capture, sized on the gross
    air-to-cloth ratio with the classic filter-drag pressure-drop model (see
    module docstring; Cooper & Alley ch. 6, EPA Cost Manual 6e sec. 6 ch. 1).

    Ports: ``gas_in`` -> ``gas_out`` + ``solids_out`` (+ ``duty``).

    Parameters
    ----------
    solids : component id (str) or list — the dust component(s). Required.
    efficiency : overall collection efficiency, default 0.999 (well-run
        baghouses routinely exceed 99.9%; Cooper & Alley ch. 6). HYSYS uses
        an internal grade curve; this v1 takes the bulk value directly.
    face_velocity : superficial gas velocity through the cloth, m/s
        (= air-to-cloth ratio, m^3/s per m^2). Default 0.01 m/s (~2 ft/min:
        woven fabric, shaker/reverse-air service; pulse-jet felts run
        0.025-0.06 m/s. Cooper & Alley Table 6.1; EPA manual sec. 6 ch. 1).
    S_E : conditioned-fabric drag, Pa per (m/s). Default 2.5e4 (~0.5 in
        H2O/(ft/min), mid-range; Cooper & Alley ch. 6).
    K2 : specific dust-cake resistance coefficient, Pa·s·m/kg. Default 5.0e4
        (~5 in H2O/(ft/min)/(lb/ft^2), mid-range of the 1-30 tabulated span;
        Cooper & Alley ch. 6).
    dP_max : bag-cleaning trigger pressure drop, Pa. Default 2000 Pa
        (~8 in H2O) — also the gas-side dP applied to ``gas_out``
        (conservative: the dirty-bag ceiling, not the cycle average).
    bag_diameter_m / bag_length_m : per-bag geometry for the bag count
        suggestion (defaults 0.15 m x 3.0 m, typical pulse-jet bags).

    ``unit.design``: cloth area (A = Q/v_f), bag count, inlet dust loading,
    clean/dirty pressure drops, and ``filtration_time_s`` — the time after
    cleaning at which the drag model reaches dP_max (the book's sec. 12.3
    case-study variable; infinite when the inlet carries no dust).
    """

    design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("gas_in", "inlet"),
            Port("gas_out", "outlet"),
            Port("solids_out", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        label = f"BaghouseFilter {self.id!r}"
        inlet = inlets.get("gas_in")
        if inlet is None or not inlet.molar_flow:
            raise ValueError(f"{label}: missing or empty inlet on 'gas_in'")
        t_in, p_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)

        solids = _solid_ids(label, self.params, comps)
        eta = float(self.params.get("efficiency", 0.999))
        if not 0.0 < eta <= 1.0:
            raise SolidsOperationError(
                f"{label}: efficiency must lie in (0, 1], got {eta}")
        v_f = float(self.params.get("face_velocity", 0.01))
        if v_f <= 0.0:
            raise SolidsOperationError(
                f"{label}: face_velocity must be positive, got {v_f}")
        s_e = float(self.params.get("S_E", 2.5e4))
        k2 = float(self.params.get("K2", 5.0e4))
        if s_e < 0.0 or k2 <= 0.0:
            raise SolidsOperationError(
                f"{label}: S_E must be >= 0 and K2 > 0 (got S_E={s_e}, K2={k2})")
        dp_max = float(self.params.get("dP_max", 2000.0))
        dp_clean = s_e * v_f
        if dp_max <= dp_clean:
            raise SolidsOperationError(
                f"{label}: dP_max {dp_max:.4g} Pa must exceed the clean-bag "
                f"drop S_E*v = {dp_clean:.4g} Pa — raise dP_max or lower the "
                f"face_velocity"
            )
        if dp_max >= p_in:
            raise SolidsOperationError(
                f"{label}: dP_max {dp_max:.4g} Pa >= inlet pressure {p_in:.4g} Pa")
        bag_d = float(self.params.get("bag_diameter_m", 0.15))
        bag_l = float(self.params.get("bag_length_m", 3.0))
        if bag_d <= 0.0 or bag_l <= 0.0:
            raise SolidsOperationError(
                f"{label}: bag_diameter_m and bag_length_m must be positive "
                f"(got {bag_d}, {bag_l})"
            )

        gas = _carrier_state(label, inlet, solids, pp, "vapor")
        area = gas["Vdot"] / v_f
        n_bags = math.ceil(area / (math.pi * bag_d * bag_l))
        m_in = _solids_mass_rate(inlet, solids)
        loading = m_in / gas["Vdot"]                       # kg dust / m^3 gas
        # Filter drag: dP(t) = S_E v + K2 c v^2 t -> time to reach dP_max.
        t_filter = ((dp_max - dp_clean) / (k2 * loading * v_f * v_f)
                    if loading > 0.0 else math.inf)

        feed = {c: n * z.get(c, 0.0) for c in comps}
        captured = {c: (feed[c] * eta if c in solids else 0.0) for c in comps}
        to_gas = {c: feed[c] - captured[c] for c in comps}

        h_in_total = n * (inlet.H if inlet.H is not None
                          else pp.enthalpy(t_in, p_in, z))
        gas_out, h_gas = _outlet(self.id, "gas_out", comps, to_gas,
                                 t_in, p_in - dp_max, z, pp)
        sol_out, h_sol = _outlet(self.id, "solids_out", comps, captured,
                                 t_in, p_in - dp_max, z, pp)

        self.design = {
            "cloth_area_m2": area,
            "face_velocity_m_s": v_f,
            "n_bags": n_bags,
            "bag_diameter_m": bag_d,
            "bag_length_m": bag_l,
            "Q_m3_s": gas["Vdot"],
            "gas_density_kg_m3": gas["rho"],
            "dust_loading_kg_m3": loading,
            "overall_efficiency": eta,
            "S_E_Pa_s_m": s_e,
            "K2_Pa_s_m_kg": k2,
            "dP_clean_Pa": dp_clean,
            "dP_max_Pa": dp_max,
            "filtration_time_s": t_filter,
            "solids_in_kg_s": m_in,
            "solids_captured_kg_s": m_in * eta,
            "solids_emitted_kg_s": m_in * (1.0 - eta),
        }
        return {
            "gas_out": gas_out,
            "solids_out": sol_out,
            "duty": EnergyStream(id=f"{self.id}.duty",
                                 duty=h_gas + h_sol - h_in_total),
        }
