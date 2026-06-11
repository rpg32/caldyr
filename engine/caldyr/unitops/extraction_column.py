"""Counter-current liquid-liquid extraction column (stage-wise LLE cascade).

A tray (or mixer-settler cascade) extractor: two partially-miscible liquids
contact counter-currently over ``n_stages`` theoretical stages; the solute
transfers from the feed phase to the solvent phase (Hameed, *Chemical Process
Simulations using Aspen Hysys*, Wiley 2025, sec. 9.4; Seader, Henley & Roper,
*Separation Process Principles* 3e, ch. 8).

**Port convention (the book's / HYSYS's orientation, Hameed 2025 fig. 9.39):**
the **heavier** liquid enters the **top** (``feed_in``, stage 1) and flows
down; the **lighter** liquid enters the **bottom** (``solvent_in``, stage
n_stages) and rises. Whatever enters at the bottom leaves at the top:
``extract_out`` is the top product (the phase that travelled up with the
bottom-fed solvent — the book's "Ovhd Light Liquid" / Rich-Sol), and
``raffinate_out`` is the bottom product (the descending feed phase, the
book's "Bottoms Heavy Liquid"). The book's own rule applies: *"the higher
liquid density must enter from the top of the column, and lighter liquid
density must enter from the bottom"* — when the *solvent* is the heavier
liquid (e.g. water washing a solute out of a hydrocarbon), feed the solvent
on ``feed_in`` and the solute-rich light feed on ``solvent_in``; the
solute-rich solvent then leaves on ``raffinate_out``. The ports are
positional (top/bottom), not semantic, exactly like HYSYS's. The model does
not enforce gravity — the feed liquid densities (as computed by the property
package) are reported in ``design['rho_feed_in']/['rho_solvent_in']`` so a
mis-oriented column can be spotted, but they are not policed: cubic-EOS
liquid densities for water are poor enough (PR puts water near toluene's
density) that a hard check would reject physically correct setups.

**Isothermal v1**: the stage temperature is ``params['T']`` (required, like
the ThreePhaseSeparator) at pressure ``params['P']`` (default: feed
pressure). LLE is only weakly temperature-dependent stage to stage and the
heat effects of liquid-liquid mixing are small, so a single operating
temperature is the standard first model; the heat needed to hold it is
reported honestly on the ``duty`` energy port (like a ComponentSplitter),
including the small enthalpy-of-demixing offset between the VL and VLL
enthalpy surfaces (see ThreePhaseSeparator for the same note).

**Per-stage equilibrium** comes from the property package's three-phase
flash ``flash_pt_3p`` (`thermo`'s FlashVLN through the existing public API) —
so only the cubic-EOS packages (``thermo:PR`` / ``thermo:SRK``) are
supported; the NRTL activity package raises a clear error. **Phase tracking
is by composition continuity, not density**: FlashVLN labels its two liquids
light/heavy by mass density, but PR's liquid densities can be nearly
degenerate (water vs toluene differ by < 0.1% under PR), so each stage's two
phases are assigned to the descending/ascending streams by closeness (L1
distance) to the previous compositions of those streams — the descending
family is seeded by the top feed and the ascending family by the bottom
solvent.

**Algorithm** (isothermal sum-rates-like substitution): with descending
streams ``d_j`` (leaving stage j downward) and ascending streams ``a_j``
(leaving stage j upward), each Gauss-Seidel sweep recomputes every stage
top-to-bottom: combine the stage inputs ``d_{j-1} + a_{j+1}`` (the column
feeds at the end stages), flash at (T, P), split into the two liquid phases,
and damp the update. Iterate until the stage component flows stop moving
(max change <= 1e-9 of the total feed). Convergence is geometric with ratio
~ min(E, 1/E) in the extraction factor E, so a well-posed column converges
in tens of sweeps; an oscillating profile triggers damping, and a
non-convergent or degenerate (single-phase / vaporizing) column raises
:class:`ExtractionColumnError` with diagnostics — never a silent wrong
answer.

Component mass balances close to machine precision (the raffinate is the
feeds minus the extract); the duty port closes the energy balance exactly.
The full converged stage profiles (ascending/descending flows and
compositions) are stored on ``unit.design``.

**Validation honesty (PR vs NRTL/experiment)**: see
tests/test_m13_extraction.py. For the book's polar acetone-water /
3-methylhexane example (Hameed 2025 sec. 9.4), PR expels essentially *all*
non-water species from the aqueous phase (K_acetone ~ 1e5 organic/aqueous vs
a realistic O(1) distribution), so the acetone split is structural, not
NRTL-quantitative — though the book's own headline result (raffinate water
mole fraction = 1.0) is reproduced, since HYSYS-PRSV behaves the same way
there. The quantitative burden is carried by a system PR distributes
finitely (methanol between toluene and water) checked against the Kremser
closed form.
"""
from __future__ import annotations

import math
from typing import Any

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register

_MAX_ITER = 300          # default sweep cap ('max_iter' overrides)
_TOL_FLOW = 1e-9         # stage component-flow convergence, fraction of feed
_FLOW_FLOOR = 0.0        # an absent phase is honestly zero (checked at the end)
_VAPOR_TOL = 1e-9        # any vapor fraction above this is an error


class ExtractionColumnError(ValueError):
    """Specification or convergence error in an ExtractionColumn (bad stage
    counts, vaporizing or single-phase stages, sweep failure, ...)."""


def _l1(a: dict[str, float], b: dict[str, float], comps: list[str]) -> float:
    return sum(abs(a.get(c, 0.0) - b.get(c, 0.0)) for c in comps)


@register("ExtractionColumn")
class ExtractionColumn(UnitOp):
    """Counter-current liquid-liquid extraction column. See the module
    docstring for the stage convention (heavy liquid in at the top on
    ``feed_in``, light liquid in at the bottom on ``solvent_in``; the
    ascending phase leaves on ``extract_out`` at the top, the descending one
    on ``raffinate_out`` at the bottom) and the algorithm.

    Ports: ``feed_in`` (top-stage inlet), ``solvent_in`` (bottom-stage
    inlet), ``extract_out`` (top outlet), ``raffinate_out`` (bottom outlet),
    ``duty`` (energy: heat to hold the operating temperature).

    Params (JSON-friendly scalars; ``.flow`` round-trips):
      * ``n_stages`` — theoretical stages. Required, >= 1.
      * ``T`` — operating temperature, K. **Required** (isothermal v1).
      * ``P`` — operating pressure, Pa (default: feed_in pressure).
      * ``max_iter`` — sweep cap (default 300).
    """

    design: dict[str, Any] | None = None

    def __init__(self, id: str, params: dict | None = None) -> None:
        super().__init__(id, params)
        self._warm: dict[str, Any] | None = None
        self._cache_key: tuple | None = None
        self._cache_out: dict[str, PortStream] | None = None
        self._cache_design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("feed_in", "inlet"),
            Port("solvent_in", "inlet"),
            Port("extract_out", "outlet"),
            Port("raffinate_out", "outlet"),
            Port("duty", "outlet", "energy"),
        ]

    def _read_params(self, P_in: float) -> tuple[int, float, float, int]:
        try:
            n_stages = int(self.params["n_stages"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: integer 'n_stages' is "
                f"required (got {self.params.get('n_stages')!r})"
            ) from exc
        if n_stages < 1:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: n_stages={n_stages} must be "
                f">= 1"
            )
        if self.params.get("T") is None:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: params['T'] is required — "
                f"the column is isothermal in v1 (LLE is weakly "
                f"T-dependent); specify the operating temperature"
            )
        T = float(self.params["T"])
        if T <= 0.0:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: T={T} K must be > 0")
        P = float(self.params.get("P") or P_in)
        if P <= 0.0:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: P={P} Pa must be > 0")
        max_iter = int(self.params.get("max_iter", _MAX_ITER))
        return n_stages, T, P, max_iter

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        feed = inlets.get("feed_in")
        solvent = inlets.get("solvent_in")
        if feed is None or not feed.molar_flow:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: missing or empty inlet on "
                f"'feed_in'"
            )
        if solvent is None or not solvent.molar_flow:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: missing or empty inlet on "
                f"'solvent_in'"
            )
        T_f, P_f, F = feed.require_state()
        T_s, P_s, S = solvent.require_state()
        z_f, z_s = feed.normalized_z(), solvent.normalized_z()
        comps = list(feed.components)
        N, T, P, max_iter = self._read_params(P_f)

        key = (repr(sorted(self.params.items())),
               T_f, P_f, F, tuple(sorted(z_f.items())), feed.H,
               T_s, P_s, S, tuple(sorted(z_s.items())), solvent.H)
        if key == self._cache_key and self._cache_out is not None:
            assert self._cache_design is not None
            self.design = _copy_design(self._cache_design)
            return _copy_out(self._cache_out)

        f = {c: F * z_f.get(c, 0.0) for c in comps}      # top (descending) feed
        s = {c: S * z_s.get(c, 0.0) for c in comps}      # bottom (ascending) feed
        f_tot = {c: f[c] + s[c] for c in comps}
        F_total = F + S

        res = self._solve_cascade(pp, comps, N, T, P, f, s, F_total, max_iter)
        a, d = res["a"], res["d"]

        # -- products: extract by the converged top ascending stream, raffinate
        # by difference — every component balance closes to machine precision.
        e_flows = {c: a[0].get(c, 0.0) for c in comps}
        r_flows = {c: f_tot[c] - e_flows[c] for c in comps}
        neg = min(r_flows.values())
        if neg < -1e-7 * F_total:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: the converged extract would "
                f"carry more of a component than the feeds supply (raffinate "
                f"flow {neg:.3g} mol/s) — tighten tolerances or check the "
                f"specifications"
            )
        if neg < 0.0:
            scale_to = sum(f_tot.values()) - sum(e_flows.values())
            r_flows = {c: max(v, 0.0) for c, v in r_flows.items()}
            scale = scale_to / sum(r_flows.values())
            r_flows = {c: v * scale for c, v in r_flows.items()}
        E_out = sum(e_flows.values())
        R_out = sum(r_flows.values())
        if E_out <= 0.0 or R_out <= 0.0:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: a product phase vanished at "
                f"the converged point (extract {E_out:.3g}, raffinate "
                f"{R_out:.3g} mol/s) — the two feeds are fully miscible at "
                f"T={T:.1f} K, P={P:.4g} Pa, so the column cannot operate"
            )
        x_e = {c: v / E_out for c, v in e_flows.items()}
        x_r = {c: v / R_out for c, v in r_flows.items()}

        # Honest isothermal duty: outlet enthalpies on the package's liquid
        # surface at the product compositions; the duty closes the overall
        # energy balance exactly.
        h_f = feed.H if feed.H is not None else pp.enthalpy(T_f, P_f, z_f)
        h_s = (solvent.H if solvent.H is not None
               else pp.enthalpy(T_s, P_s, z_s))
        h_e = pp.enthalpy_liquid(T, P, x_e)
        h_r = pp.enthalpy_liquid(T, P, x_r)
        duty = E_out * h_e + R_out * h_r - (F * h_f + S * h_s)

        extract = Stream(
            id=f"{self.id}.extract_out", components=comps,
            T=T, P=P, molar_flow=E_out, z=x_e,
            H=h_e, phase="liquid", vapor_fraction=0.0,
        )
        raffinate = Stream(
            id=f"{self.id}.raffinate_out", components=comps,
            T=T, P=P, molar_flow=R_out, z=x_r,
            H=h_r, phase="liquid", vapor_fraction=0.0,
        )

        rho_f, rho_s = self._feed_densities(pp, T, P, z_f, z_s, comps)
        self.design = {
            "P": P, "T": T, "n_stages": N,
            # tray-count key for the economics sizer (no condenser/reboiler).
            "N": float(N),
            "extract_flows": dict(e_flows), "raffinate_flows": dict(r_flows),
            "x_extract": dict(x_e), "x_raffinate": dict(x_r),
            "extract_total": E_out, "raffinate_total": R_out,
            # Fraction of the total feed of each component leaving in the
            # extract (the recovery, for components fed only at the top).
            "recovery": {c: (e_flows[c] / f_tot[c]) if f_tot[c] > 0.0 else 0.0
                         for c in comps},
            "T_profile": [T] * N, "P_profile": [P] * N,
            # Ascending (extract-side) and descending (raffinate-side) stage
            # streams leaving each stage, top -> bottom.
            "E_profile": [sum(row.values()) for row in a],
            "R_profile": [sum(row.values()) for row in d],
            "x_extract_profile": [_norm_or_zero(row, comps) for row in a],
            "x_raffinate_profile": [_norm_or_zero(row, comps) for row in d],
            "rho_feed_in": rho_f, "rho_solvent_in": rho_s,
            "iterations": res["iterations"],
            "max_dF_rel": res["max_dF_rel"], "damping": res["damping"],
            "flash_calls": res["flash_calls"],
            "duty": duty,
        }
        out: dict[str, PortStream] = {
            "extract_out": extract,
            "raffinate_out": raffinate,
            "duty": EnergyStream(id=f"{self.id}.duty", duty=duty),
        }
        self._warm = {"n_stages": N, "comps": list(comps),
                      "z": {c: f_tot[c] / F_total for c in comps},
                      "a": [dict(row) for row in a],
                      "d": [dict(row) for row in d]}
        self._cache_key = key
        self._cache_out = _copy_out(out)
        self._cache_design = _copy_design(self.design)
        return out

    # -- the counter-current cascade ----------------------------------------
    def _solve_cascade(self, pp, comps: list[str], N: int, T: float, P: float,
                       f: dict[str, float], s: dict[str, float],
                       F_total: float, max_iter: int) -> dict[str, Any]:
        """Gauss-Seidel substitution over the stage LLE flashes (see the
        module docstring). Returns the converged ascending/descending stage
        streams and the iteration diagnostics."""
        a, d = self._initial_profiles(N, comps, f, s, F_total)
        tol = _TOL_FLOW * F_total
        omega = 1.0
        dF_prev = math.inf
        worsened = 0
        converged = False
        n_flash = 0
        it = 0
        dF = math.inf
        for it in range(1, max_iter + 1):
            dF = 0.0
            n_two_phase = 0
            for j in range(N):
                m = {c: (d[j - 1][c] if j > 0 else f[c])
                     + (a[j + 1][c] if j < N - 1 else s[c]) for c in comps}
                m_tot = sum(m.values())
                if m_tot <= 0.0:
                    raise ExtractionColumnError(
                        f"ExtractionColumn {self.id!r}: stage {j + 1} ran "
                        f"dry during iteration {it}"
                    )
                z_j = {c: v / m_tot for c, v in m.items()}
                try:
                    res = pp.flash_pt_3p(T, P, z_j)
                except NotImplementedError as exc:
                    raise ExtractionColumnError(
                        f"ExtractionColumn {self.id!r}: the property package "
                        f"does not implement the three-phase (LLE) flash — "
                        f"select 'thermo:PR' or 'thermo:SRK' as the "
                        f"flowsheet property package ({exc})"
                    ) from exc
                n_flash += 1
                liq_phases = self._liquid_phases(pp, res, T, P, j)
                if len(liq_phases) >= 2:
                    n_two_phase += 1
                a_new, d_new = self._split_phases(liq_phases, m, m_tot,
                                                  comps, a[j], d[j], f, s)
                for c in comps:
                    da, dd = a_new[c] - a[j][c], d_new[c] - d[j][c]
                    dF = max(dF, abs(da), abs(dd))
                    a[j][c] += omega * da
                    d[j][c] += omega * dd
            # Every stage flashing to a single liquid means there is no
            # two-liquid region anywhere along the cascade — the feed and
            # solvent are mutually soluble, and continuing the sweep can only
            # flip-flop the one phase between the two streams. Refuse early
            # with the real diagnosis instead of a convergence failure.
            if n_two_phase == 0:
                raise ExtractionColumnError(
                    f"ExtractionColumn {self.id!r}: every stage flashed to a "
                    f"single liquid phase (iteration {it}) — the feed and "
                    f"solvent are fully miscible at T={T:.1f} K, "
                    f"P={P:.4g} Pa, so liquid-liquid extraction cannot "
                    f"operate; change the solvent or the temperature"
                )
            # Adaptive damping if the sweep residual oscillates upward.
            if dF > dF_prev:
                worsened += 1
                if worsened >= 2:
                    omega = max(0.25, 0.5 * omega)
                    worsened = 0
            else:
                worsened = 0
            dF_prev = dF
            if dF <= tol:
                converged = True
                break

        if not converged:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: the counter-current sweep "
                f"did not converge in {max_iter} iterations (max stage-flow "
                f"change {dF:.3g} mol/s vs {tol:.3g}, damping={omega}). "
                f"Check that both liquid phases exist at T={T:.1f} K and "
                f"that the solvent rate is sensible, or raise 'max_iter'"
            )
        # A stage that lost one of its phases at convergence means the column
        # is operating outside the two-liquid region — refuse to mislabel it.
        for j in range(N):
            if sum(a[j].values()) <= 0.0 or sum(d[j].values()) <= 0.0:
                raise ExtractionColumnError(
                    f"ExtractionColumn {self.id!r}: stage {j + 1} collapsed "
                    f"to a single liquid phase at the converged point "
                    f"(ascending {sum(a[j].values()):.3g}, descending "
                    f"{sum(d[j].values()):.3g} mol/s) — the feed/solvent "
                    f"pair is miscible at this stage's composition; change "
                    f"T, the solvent, or the solvent rate"
                )
        return {"a": a, "d": d, "iterations": it,
                "max_dF_rel": dF / F_total, "damping": omega,
                "flash_calls": n_flash}

    def _liquid_phases(self, pp, res, T: float, P: float,
                       j: int) -> list[tuple[float, dict[str, float]]]:
        """The stage's liquid phases as ``(beta, x)`` pairs, after auditing a
        claimed vapor phase. `thermo`'s FlashVLN sometimes mislabels a
        hydrocarbon-rich *liquid* as the gas phase (e.g. the toluene-rich
        liquid of water/toluene/methanol at 25 C, 1 atm comes back as a
        "vapor" of 91% toluene whose dew point is ~110 C). The audit is by
        dew point: a *real* equilibrium vapor at (T, P) sits exactly at its
        dew temperature, so a "vapor" whose dew temperature lies well above
        the stage temperature cannot be a stable gas and is reclassified as
        the liquid it is. A genuine vapor (stage actually boiling) raises —
        a liquid-liquid extractor must run all-liquid."""
        phases: list[tuple[float, dict[str, float]]] = []
        for beta, x in ((res.beta_light, res.x_light),
                        (res.beta_heavy, res.x_heavy)):
            if beta > 0.0 and x is not None:
                phases.append((beta, x))
        if res.beta_vapor > _VAPOR_TOL and res.y is not None:
            _, dew_t = pp.bubble_dew(P, res.y)
            if T < dew_t - 1.0:
                phases.append((res.beta_vapor, res.y))   # mislabeled liquid
            else:
                raise ExtractionColumnError(
                    f"ExtractionColumn {self.id!r}: stage {j + 1} vaporizes "
                    f"at T={T:.1f} K, P={P:.4g} Pa (vapor fraction "
                    f"{res.beta_vapor:.3g}, vapor dew point {dew_t:.1f} K) "
                    f"— a liquid-liquid extractor must run all-liquid; "
                    f"lower T or raise P"
                )
        if len(phases) > 2:
            raise ExtractionColumnError(
                f"ExtractionColumn {self.id!r}: stage {j + 1} split into "
                f"{len(phases)} liquid-like phases at T={T:.1f} K, "
                f"P={P:.4g} Pa — a two-liquid counter-current cascade "
                f"cannot route a third liquid phase"
            )
        return phases

    @staticmethod
    def _split_phases(liq_phases: list[tuple[float, dict[str, float]]],
                      m: dict[str, float], m_tot: float,
                      comps: list[str], a_prev: dict[str, float],
                      d_prev: dict[str, float], f: dict[str, float],
                      s: dict[str, float],
                      ) -> tuple[dict[str, float], dict[str, float]]:
        """Assign the stage's liquid phases to the ascending/descending
        streams by composition continuity (L1 distance to the streams'
        previous compositions; falling back to the column feeds when a stream
        is still empty) — NOT by the flasher's density labels, which can be
        nearly degenerate or inverted under PR (see the module docstring)."""
        phases = [({c: m_tot * beta * x.get(c, 0.0) for c in comps}, x)
                  for beta, x in liq_phases]
        zero = {c: 0.0 for c in comps}
        if not phases:
            return zero, dict(m)            # degenerate: caught at convergence
        ref_a = _norm_or_zero(a_prev, comps)
        if sum(ref_a.values()) <= 0.0:
            ref_a = _norm_or_zero(s, comps)
        ref_d = _norm_or_zero(d_prev, comps)
        if sum(ref_d.values()) <= 0.0:
            ref_d = _norm_or_zero(f, comps)
        if len(phases) == 1:
            flows, x = phases[0]
            if _l1(x, ref_a, comps) < _l1(x, ref_d, comps):
                return flows, zero
            return zero, flows
        (fl0, x0), (fl1, x1) = phases
        straight = _l1(x0, ref_a, comps) + _l1(x1, ref_d, comps)
        swapped = _l1(x1, ref_a, comps) + _l1(x0, ref_d, comps)
        if straight <= swapped:
            return fl0, fl1
        return fl1, fl0

    def _initial_profiles(self, N: int, comps: list[str],
                          f: dict[str, float], s: dict[str, float],
                          F_total: float,
                          ) -> tuple[list[dict[str, float]],
                                     list[dict[str, float]]]:
        """Stage starting streams: the last converged profiles when the
        layout and combined feed are close (recycle sweeps, EO re-solves),
        else the feeds replicated down the column."""
        w = self._warm
        z_tot = {c: (f[c] + s[c]) / F_total for c in comps}
        if (w is not None and w["n_stages"] == N and w["comps"] == comps
                and max(abs(w["z"].get(c, 0.0) - z_tot[c]) for c in comps)
                < 0.1):
            return ([dict(row) for row in w["a"]],
                    [dict(row) for row in w["d"]])
        return ([dict(s) for _ in range(N)], [dict(f) for _ in range(N)])

    @staticmethod
    def _feed_densities(pp, T: float, P: float, z_f: dict[str, float],
                        z_s: dict[str, float], comps: list[str],
                        ) -> tuple[float | None, float | None]:
        """Mass densities (kg/m^3) of the two feeds as liquids at (T, P) —
        diagnostics for the heavy-on-top orientation rule (reported, not
        policed; see the module docstring)."""
        try:
            from ..core.components_db import molar_mass
            rho = []
            for z in (z_f, z_s):
                mw = sum(z.get(c, 0.0) * molar_mass(c) for c in comps)
                rho.append(mw / pp.volume_liquid(T, P, z))
            return rho[0], rho[1]
        except Exception:                    # diagnostics must never fail a solve
            return None, None


def _norm_or_zero(row: dict[str, float], comps: list[str]) -> dict[str, float]:
    tot = sum(row.get(c, 0.0) for c in comps)
    if tot <= 0.0:
        return {c: 0.0 for c in comps}
    return {c: row.get(c, 0.0) / tot for c in comps}


def _copy_out(out: dict[str, PortStream]) -> dict[str, PortStream]:
    return {
        name: (st.with_() if isinstance(st, Stream)
               else EnergyStream(id=st.id, duty=st.duty))
        for name, st in out.items()
    }


def _copy_design(design: dict[str, Any]) -> dict[str, Any]:
    return {k: (list(v) if isinstance(v, list) else
                dict(v) if isinstance(v, dict) else v)
            for k, v in design.items()}
