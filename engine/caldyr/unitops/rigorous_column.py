"""Rigorous tray-by-tray distillation column: bubble-point (Wang-Henke) MESH.

Solves the full MESH equations (Material balance, Equilibrium, Summation,
Heat balance) for a simple column — one feed, a condenser, a reboiler, two
products — with the bubble-point method of Wang & Henke (1966), the standard
workhorse for narrow- and medium-boiling distillation (Seader, Henley & Roper,
*Separation Process Principles*, 3e, ch. 10.3):

1. **Initialize** from shortcut (FUG-style) results: a Fenske split of the
   feed sized to the distillate rate estimates the product compositions, the
   stage compositions are interpolated linearly between them (which puts the
   temperature profile on a near-linear ramp between the products' bubble
   points), and the internal traffic is constant-molal-overflow from R, D and
   the feed quality q.
2. **Component balances**: with the current ``K_ij``, ``L_j``, ``V_j``, each
   component's stage balances form a tridiagonal system, solved by the Thomas
   algorithm.
3. **Equilibrium + summation**: normalize each stage's liquid; one
   saturated-liquid (VF=0) flash per stage then yields the new stage
   temperature (the bubble point of x_j), the stage K-values (incipient vapor
   over liquid) and both saturated phase molar enthalpies.
4. **Heat balances**: the stage energy balances reduce (after eliminating L_j
   with the column-section total balance ``L_j = V_{j+1} + S_j - D``) to a
   forward recurrence for the vapor traffic ``V_j``; updates are damped if the
   temperature profile starts oscillating.
5. Repeat 2-4 until the stage temperatures and vapor traffic stop moving
   (max |dT| <= 1e-4 K and max |dV|/V <= 1e-6). A non-convergent or
   degenerate iteration raises :class:`RigorousColumnError` with diagnostics —
   never a silent wrong answer.

**Stage convention**: stages are numbered from the top, ``1 .. n_stages``,
**including the condenser as stage 1 and the reboiler as stage n_stages**.
``feed_stage`` uses the same numbering, so it must lie in
``2 .. n_stages - 1``. A ShortcutColumn's ``N`` (which counts the reboiler but
not a total condenser) therefore corresponds to ``n_stages = N + 1``, and its
``feed_stage`` (counted from the top tray) to ``feed_stage + 1`` here.

Products leave at stage pressure: the distillate as liquid at the stage-1
bubble point (total condenser) or vapor at its dew point
(``partial_condenser=True``); the bottoms as liquid at the stage-N bubble
point. The condenser duty comes from the stage-1 energy balance (Heater sign
convention: heat added to the process is positive, so Q_cond <= 0); the
reboiler duty closes the overall energy balance exactly::

    F h_F + Q_cond + Q_reb = D h_D + B h_B

and the residual against the stage-N (reboiler) balance is reported in
``design['energy_residual_rel']`` as a convergence diagnostic.

The full converged profile (per-stage T, P, L, V, x, y) is stored on
``unit.design`` so column profiles can be plotted; the FUG-compatible keys
(``N``, ``V_top``, ``x_D`` ...) are stored too, so the unit is sized and
costed exactly like a ShortcutColumn (tower + trays + condenser + reboiler).

Out of scope (this round): absorber/stripper mode (no condenser/reboiler),
multiple feeds, side draws, and pumparounds.
"""
from __future__ import annotations

import math
from typing import Any

from scipy.optimize import brentq

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


class RigorousColumnError(ValueError):
    """Specification or convergence error in a RigorousColumn (bad stage
    counts, infeasible distillate rate, MESH iteration failure, ...)."""


_MAX_ITER = 200          # default MESH iteration cap (param 'max_iter' overrides)
_TOL_T = 1e-4            # K, stage-temperature convergence
_TOL_V = 1e-6            # relative vapor-traffic convergence
_X_FLOOR = 1e-15         # floor for stage mole fractions (Thomas can go negative)
_V_FLOOR_FRAC = 1e-8     # vapor/liquid traffic floor, as a fraction of the feed
_WARM_Z_TOL = 0.1        # max |dz| of the feed for warm-starting from last solve
_FENSKE_CLAMP = 50.0     # caps the logistic exponent in the Fenske initializer


def _thomas(a: list[float], b: list[float], c: list[float],
            d: list[float]) -> list[float]:
    """Solve the tridiagonal system (a=sub, b=diag, c=super) by the Thomas
    algorithm. O(n); the standard solver for MESH component balances."""
    n = len(b)
    cp = [0.0] * n
    dp = [0.0] * n
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for j in range(1, n):
        denom = b[j] - a[j] * cp[j - 1]
        cp[j] = c[j] / denom if j < n - 1 else 0.0
        dp[j] = (d[j] - a[j] * dp[j - 1]) / denom
    x = [0.0] * n
    x[-1] = dp[-1]
    for j in range(n - 2, -1, -1):
        x[j] = dp[j] - cp[j] * x[j + 1]
    return x


@register("RigorousColumn")
class RigorousColumn(UnitOp):
    """Rigorous (MESH, bubble-point) distillation column. See the module
    docstring for the algorithm and the stage-numbering convention.

    Params (JSON-friendly scalars; ``.flow`` round-trips):
      * ``n_stages`` — theoretical stages, **including** the condenser
        (stage 1) and the reboiler (stage n_stages). Required, >= 3.
      * ``feed_stage`` — feed stage in the same numbering; 2..n_stages-1.
      * ``reflux_ratio`` — molar reflux ratio L_1/D, > 0. Required.
      * exactly one of ``distillate_rate`` (mol/s) or ``distillate_to_feed``
        (molar D/F fraction).
      * ``P`` — top-stage pressure, Pa (default: feed pressure).
      * ``dP_stage`` — linear pressure rise per stage going down, Pa
        (default 0: uniform column pressure).
      * ``partial_condenser`` — vapor distillate at its dew point (default
        False: total condenser, liquid distillate at its bubble point).
      * ``max_iter`` — MESH iteration cap (default 200).
    """

    design: dict[str, Any] | None = None

    def __init__(self, id: str, params: dict | None = None) -> None:
        super().__init__(id, params)
        # Warm-start memory (last converged liquid profile) and an exact-repeat
        # cache: the equation-oriented solver re-calls solve() with identical
        # inlets dozens of times per Newton step, and a recycle sweep with
        # near-identical ones — both must not pay full MESH price every time.
        self._warm: dict[str, Any] | None = None
        self._cache_key: tuple | None = None
        self._cache_out: dict[str, PortStream] | None = None
        self._cache_design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("distillate", "outlet"),
            Port("bottoms", "outlet"),
            Port("condenser_duty", "outlet", "energy"),
            Port("reboiler_duty", "outlet", "energy"),
        ]

    # -- parameter validation --------------------------------------------------
    def _read_params(self, F: float, P_in: float):
        try:
            n_stages = int(self.params["n_stages"])
            feed_stage = int(self.params["feed_stage"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: integer 'n_stages' and 'feed_stage' "
                f"are required (got n_stages={self.params.get('n_stages')!r}, "
                f"feed_stage={self.params.get('feed_stage')!r})"
            ) from exc
        if n_stages < 3:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: n_stages={n_stages} must be >= 3 — "
                f"the count includes the condenser (stage 1) and the reboiler "
                f"(stage n_stages), so 3 is one tray plus condenser and reboiler"
            )
        if not 2 <= feed_stage <= n_stages - 1:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: feed_stage={feed_stage} is out of "
                f"range — it must lie between stage 2 (below the condenser) and "
                f"stage {n_stages - 1} (above the reboiler) for n_stages={n_stages}"
            )

        try:
            R = float(self.params["reflux_ratio"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: 'reflux_ratio' is required "
                f"(got {self.params.get('reflux_ratio')!r})"
            ) from exc
        if R <= 0.0:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: reflux_ratio={R} must be > 0 — "
                f"a column with no reflux cannot rectify"
            )

        dr = self.params.get("distillate_rate")
        dtf = self.params.get("distillate_to_feed")
        if (dr is None) == (dtf is None):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: specify exactly one of "
                f"'distillate_rate' (mol/s) or 'distillate_to_feed' (fraction); "
                f"got distillate_rate={dr!r}, distillate_to_feed={dtf!r}"
            )
        D = float(dr) if dr is not None else float(dtf or 0.0) * F
        if not 0.0 < D < F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: distillate rate D={D:.6g} mol/s must "
                f"lie strictly between 0 and the feed rate F={F:.6g} mol/s "
                f"(both products must exist)"
            )

        P = float(self.params.get("P") or P_in)
        if P <= 0.0:
            raise RigorousColumnError(f"RigorousColumn {self.id!r}: P={P} Pa must be > 0")
        dP = float(self.params.get("dP_stage", 0.0))
        if dP < 0.0:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: dP_stage={dP} Pa must be >= 0 "
                f"(pressure increases going down the column)"
            )
        partial = bool(self.params.get("partial_condenser", False))
        max_iter = int(self.params.get("max_iter", _MAX_ITER))
        return n_stages, feed_stage, R, D, P, dP, partial, max_iter

    # -- solve -------------------------------------------------------------------
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: missing or empty inlet on 'in1'"
            )
        T_in, P_in, F = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)
        n_stages, feed_stage, R, D, P, dP, partial, max_iter = \
            self._read_params(F, P_in)

        # Exact-repeat cache (see __init__): same params + same inlet state
        # -> return copies of the previous result without re-iterating MESH.
        key = (repr(sorted(self.params.items())), T_in, P_in, F,
               tuple(sorted(z.items())), inlet.H)
        if key == self._cache_key and self._cache_out is not None:
            assert self._cache_design is not None
            self.design = {k: (list(v) if isinstance(v, list) else
                               dict(v) if isinstance(v, dict) else v)
                           for k, v in self._cache_design.items()}
            return {
                name: (s.with_() if isinstance(s, Stream)
                       else EnergyStream(id=s.id, duty=s.duty))
                for name, s in self._cache_out.items()
            }

        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)
        active = [c for c in comps if z.get(c, 0.0) > 0.0]
        if len(active) < 2:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: the feed carries "
                f"{len(active)} component(s) ({active}); distillation needs at "
                f"least two"
            )

        N = n_stages
        f0 = feed_stage - 1                              # 0-based feed stage
        P_j = [P + j * dP for j in range(N)]             # stage pressures, top->down
        B = F - D
        f = {c: F * z[c] for c in active}                # component feed flows

        # Feed quality q at the feed stage (isenthalpic flash to stage pressure).
        feed_res = pp.flash_ph(P_j[f0], H_in, {c: z[c] for c in active})
        q = 1.0 - feed_res.vapor_fraction

        x = self._initial_x(pp, P, f, D, N, z)
        L, V = self._initial_traffic(N, f0, R, D, F, q, partial)

        # Stage-by-stage bubble points of the initial profile seed T, K, h.
        T, K, hL, hV, n_flash = self._stage_bubble_points(pp, P_j, x, active, 0)

        # -- MESH loop (Wang-Henke) ---------------------------------------------
        S = [F if j >= f0 else 0.0 for j in range(N)]    # cumulative feed above j
        v_floor = _V_FLOOR_FRAC * F
        omega = 1.0                                      # damping on the V update
        dT_prev = math.inf
        worsened = 0
        converged = False
        it = 0
        dT = dV = math.inf
        for it in range(1, max_iter + 1):
            x = self._component_balances(N, f0, active, z, F, D, K, L, V, partial)
            T_new, K, hL, hV, n_flash = \
                self._stage_bubble_points(pp, P_j, x, active, n_flash)
            dT = max(abs(tn - to) for tn, to in zip(T_new, T))
            T = T_new

            # Damp the traffic update if the temperature profile oscillates.
            if dT > dT_prev:
                worsened += 1
                if worsened >= 2:
                    omega = max(0.25, 0.5 * omega)
                    worsened = 0
            else:
                worsened = 0
            dT_prev = dT

            V_new = self._energy_balances(N, f0, F, D, H_in, S, hL, hV, V, v_floor)
            dV = max(abs(vn - vo) for vn, vo in zip(V_new[1:], V[1:])) / max(V_new[1:])
            V = [V[0]] + [max(vo + omega * (vn - vo), v_floor)
                          for vo, vn in zip(V[1:], V_new[1:])]
            L = [V[j + 1] + S[j] - D for j in range(N - 1)] + [B]

            if dT <= _TOL_T and dV <= _TOL_V:
                converged = True
                break

        if not converged:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: MESH (bubble-point) iteration did "
                f"not converge in {max_iter} iterations (max |dT|={dT:.3g} K vs "
                f"{_TOL_T} K, max |dV|/V={dV:.3g} vs {_TOL_V}, damping={omega}). "
                f"Check the specs (R={R}, D={D:.4g} mol/s, n_stages={N}) — a "
                f"distillate rate far from what the reflux can enrich, or a "
                f"wide-boiling feed, may need more iterations ('max_iter') or a "
                f"different method"
            )
        if any(lj <= v_floor for lj in L) or any(vj <= v_floor for vj in V[1:]):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: a column section dried up at the "
                f"converged point (min L={min(L):.3g}, min V={min(V[1:]):.3g} "
                f"mol/s) — the specified R={R} and D={D:.4g} mol/s are "
                f"infeasible for this feed"
            )

        # -- products (mass balance closed exactly by difference) ----------------
        y = [{c: K[j][c] * x[j][c] for c in active} for j in range(N)]
        for j in range(N):
            tot = sum(y[j].values())
            y[j] = {c: v / tot for c, v in y[j].items()}

        x_d = dict(y[0]) if partial else dict(x[0])
        d_flows = {c: D * x_d[c] for c in active}
        b_flows = {c: f[c] - d_flows[c] for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: converged distillate would carry "
                f"more of a component than the feed supplies "
                f"(bottoms flow {neg:.3g} mol/s) — tighten tolerances or check "
                f"the distillate-rate spec"
            )
        if neg < 0.0:           # trace negatives from finite tolerance: clip
            b_flows = {c: max(v, 0.0) for c, v in b_flows.items()}
            scale = B / sum(b_flows.values())
            b_flows = {c: v * scale for c, v in b_flows.items()}
        x_b = {c: v / B for c, v in b_flows.items()}

        z_dist = {c: x_d.get(c, 0.0) for c in comps}
        z_bot = {c: x_b.get(c, 0.0) for c in comps}
        bub_d, dew_d = pp.bubble_dew(P_j[0], z_dist)
        res_d = pp.flash_pt(dew_d if partial else bub_d, P_j[0], z_dist)
        bub_b, _ = pp.bubble_dew(P_j[-1], z_bot)
        res_b = pp.flash_pt(bub_b, P_j[-1], z_bot)

        distillate = Stream(
            id=f"{self.id}.distillate", components=comps,
            T=res_d.T, P=P_j[0], molar_flow=D, z=z_dist,
            H=res_d.H, phase=res_d.phase, vapor_fraction=res_d.vapor_fraction,
        )
        bottoms = Stream(
            id=f"{self.id}.bottoms", components=comps,
            T=res_b.T, P=P_j[-1], molar_flow=B, z=z_bot,
            H=res_b.H, phase=res_b.phase, vapor_fraction=res_b.vapor_fraction,
        )

        # -- duties (Heater sign convention: positive heats the process) ----------
        # Condenser from the stage-1 energy balance; reboiler closes the overall
        # balance exactly (enthalpies are absolute / formation-inclusive):
        #   F h_F + Q_cond + Q_reb = D h_D + B h_B
        if partial:
            q_cond = L[0] * hL[0] + D * res_d.H - V[1] * hV[1]
        else:
            q_cond = (L[0] + D) * res_d.H - V[1] * hV[1]
        q_reb = D * res_d.H + B * res_b.H - F * H_in - q_cond
        # Diagnostic: the independently computed stage-N (reboiler) balance.
        q_reb_stage = V[N - 1] * hV[N - 1] + B * hL[N - 1] - L[N - 2] * hL[N - 2]
        e_scale = max(abs(q_reb), abs(q_cond), 1.0)
        energy_residual = abs(q_reb - q_reb_stage) / e_scale

        self._warm = {"n_stages": N, "partial": partial, "active": list(active),
                      "z": dict(z), "x": [dict(row) for row in x]}

        self.design = {
            # FUG-compatible keys (the economics sizer reads these — a
            # RigorousColumn is costed exactly like a ShortcutColumn). "N" is
            # the equilibrium-stage count in the ShortcutColumn convention
            # (counts the reboiler, not the condenser): n_stages - 1.
            "N": float(N - 1), "P": P_j[0], "x_D": dict(x_d), "x_B": dict(x_b),
            "T_top": res_d.T, "T_top_dew": dew_d, "T_bottom": res_b.T,
            "V_top": V[1], "Q_condenser": q_cond, "Q_reboiler": q_reb,
            "D": D, "B": B, "R": R, "q": q, "feed_stage": feed_stage,
            "partial_condenser": partial,
            "distillate_flows": dict(d_flows), "bottoms_flows": dict(b_flows),
            # Rigorous extras: the full converged stage profiles (top -> bottom,
            # length n_stages; stage 1 = condenser, stage n_stages = reboiler).
            "n_stages": N,
            "T_profile": list(T),
            "P_profile": list(P_j),
            "L_profile": list(L),
            "V_profile": list(V),
            "x_profile": [{c: row.get(c, 0.0) for c in comps} for row in x],
            "y_profile": [{c: row.get(c, 0.0) for c in comps} for row in y],
            # Convergence diagnostics.
            "iterations": it, "max_dT": dT, "max_dV_rel": dV,
            "damping": omega, "flash_calls": n_flash,
            "energy_residual_rel": energy_residual,
        }
        out: dict[str, PortStream] = {
            "distillate": distillate,
            "bottoms": bottoms,
            "condenser_duty": EnergyStream(id=f"{self.id}.condenser_duty", duty=q_cond),
            "reboiler_duty": EnergyStream(id=f"{self.id}.reboiler_duty", duty=q_reb),
        }
        self._cache_key = key
        self._cache_out = {
            name: (s.with_() if isinstance(s, Stream)
                   else EnergyStream(id=s.id, duty=s.duty))
            for name, s in out.items()
        }
        self._cache_design = {k: (list(v) if isinstance(v, list) else
                                  dict(v) if isinstance(v, dict) else v)
                              for k, v in self.design.items()}
        return out

    # -- initialization ------------------------------------------------------------
    def _initial_x(self, pp, P: float, f: dict[str, float], D: float, N: int,
                   z: dict[str, float]) -> list[dict[str, float]]:
        """Initial liquid-composition profile.

        Warm path: reuse the last converged profile when the column layout is
        unchanged and the feed composition is close (the equation-oriented
        solver and recycle sweeps re-solve with near-identical feeds).
        Cold path: a Fenske-style split — d_i/b_i proportional to
        alpha_i^N_min with the split factor solved so sum(d) = D (Fenske 1932;
        Seader 3e eq. 9-12 rearranged) — then linear interpolation between the
        estimated product compositions, the shortcut initialization the
        bubble-point method calls for.
        """
        active = list(f)
        w = self._warm
        if (w is not None and w["n_stages"] == N and w["active"] == active
                and max(abs(w["z"].get(c, 0.0) - z.get(c, 0.0)) for c in active)
                < _WARM_Z_TOL):
            return [dict(row) for row in w["x"]]

        res = pp.bubble_point(P, {c: v / sum(f.values()) for c, v in f.items()})
        if res.y is None or res.x is None:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: no vapor-liquid equilibrium for the "
                f"feed at P={P:.4g} Pa (bubble flash returned a single phase); "
                f"check the pressure"
            )
        K_feed = {c: res.y[c] / res.x[c] for c in active}
        k_min = min(K_feed.values())
        ln_alpha = {c: math.log(K_feed[c] / k_min) for c in active}
        n_min = max(2.0, 0.5 * (N - 1))      # effective Fenske stages (N ~ 2 N_min)

        def overhead(ln_s: float) -> float:
            tot = 0.0
            for c in active:
                t = ln_s + n_min * ln_alpha[c]
                t = max(min(t, _FENSKE_CLAMP), -_FENSKE_CLAMP)
                tot += f[c] / (1.0 + math.exp(-t))
            return tot - D

        ln_s = float(brentq(overhead, -300.0, 300.0, xtol=1e-10))
        d = {}
        for c in active:
            t = max(min(ln_s + n_min * ln_alpha[c], _FENSKE_CLAMP), -_FENSKE_CLAMP)
            d[c] = f[c] / (1.0 + math.exp(-t))
        x_d = {c: max(d[c] / D, _X_FLOOR) for c in active}
        b_tot = sum(f.values()) - D
        x_b = {c: max((f[c] - d[c]) / b_tot, _X_FLOOR) for c in active}

        x: list[dict[str, float]] = []
        for j in range(N):
            wgt = j / (N - 1)
            row = {c: (1.0 - wgt) * x_d[c] + wgt * x_b[c] for c in active}
            tot = sum(row.values())
            x.append({c: v / tot for c, v in row.items()})
        return x

    @staticmethod
    def _initial_traffic(N: int, f0: int, R: float, D: float, F: float,
                         q: float, partial: bool) -> tuple[list[float], list[float]]:
        """Constant-molal-overflow traffic from R, D and the feed quality q.
        V[0] is the condenser vapor product (0 for a total condenser); V[1] is
        fixed by the condenser balance at D(1+R) throughout the iteration."""
        v_rect = D * (1.0 + R)
        v_strip = max(v_rect - (1.0 - q) * F, 0.2 * F)
        V = [D if partial else 0.0] + \
            [v_rect if j <= f0 else v_strip for j in range(1, N)]
        S = [F if j >= f0 else 0.0 for j in range(N)]
        L = [V[j + 1] + S[j] - D for j in range(N - 1)] + [F - D]
        return [max(lj, _V_FLOOR_FRAC * F) for lj in L], V

    # -- MESH building blocks ---------------------------------------------------
    @staticmethod
    def _component_balances(N: int, f0: int, active: list[str],
                            z: dict[str, float], F: float, D: float,
                            K: list[dict[str, float]], L: list[float],
                            V: list[float], partial: bool) -> list[dict[str, float]]:
        """One tridiagonal solve per component (Thomas algorithm) for the stage
        liquid compositions, then floor and renormalize each stage (the
        summation equations). Coefficients per Seader 3e eqs. 10-8..10-12 with
        no side draws: the only liquid draw is the distillate off stage 1 of a
        total condenser (U_1 = D); a partial condenser's vapor product is
        V_1 = D."""
        u1 = 0.0 if partial else D
        cols: dict[str, list[float]] = {}
        for c in active:
            a = [0.0] + [L[j - 1] for j in range(1, N)]
            b = [-(L[0] + u1 + V[0] * K[0][c])] + \
                [-(L[j] + V[j] * K[j][c]) for j in range(1, N)]
            cc = [V[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
            d = [-F * z[c] if j == f0 else 0.0 for j in range(N)]
            cols[c] = _thomas(a, b, cc, d)
        x: list[dict[str, float]] = []
        for j in range(N):
            row = {c: max(cols[c][j], _X_FLOOR) for c in active}
            tot = sum(row.values())
            x.append({c: v / tot for c, v in row.items()})
        return x

    def _stage_bubble_points(self, pp, P_j: list[float], x: list[dict[str, float]],
                             active: list[str], n_flash: int):
        """One saturated-liquid flash per stage: returns the new temperature
        profile, K-values, and saturated liquid/vapor molar enthalpies."""
        T: list[float] = []
        K: list[dict[str, float]] = []
        hL: list[float] = []
        hV: list[float] = []
        for j, row in enumerate(x):
            res = pp.bubble_point(P_j[j], row)
            n_flash += 1
            if (res.y is None or res.H_liquid is None or res.H_vapor is None
                    or not math.isfinite(res.T)):
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: bubble-point flash failed on "
                    f"stage {j + 1} (T={res.T!r}, composition {row}) at "
                    f"P={P_j[j]:.4g} Pa"
                )
            T.append(res.T)
            K.append({c: res.y[c] / row[c] for c in active})
            hL.append(res.H_liquid)
            hV.append(res.H_vapor)
        return T, K, hL, hV, n_flash

    def _energy_balances(self, N: int, f0: int, F: float, D: float, h_F: float,
                         S: list[float], hL: list[float], hV: list[float],
                         V: list[float], v_floor: float) -> list[float]:
        """Forward recurrence for the vapor traffic from the stage energy
        balances (Seader 3e eqs. 10-28..10-30), with L eliminated via the
        section total balance L_j = V_{j+1} + S_j - D. V_2 stays fixed at
        D(1+R) by the condenser balance."""
        V_new = list(V)
        for j in range(1, N - 1):           # 0-based stage j (= stage j+1, 1-based)
            feed_j = F if j == f0 else 0.0
            denom = hV[j + 1] - hL[j]
            if denom <= 0.0:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: degenerate energy balance on "
                    f"stage {j + 2} (saturated vapor enthalpy below the liquid's, "
                    f"dH={denom:.4g} J/mol) — the property package returned an "
                    f"unphysical state"
                )
            num = ((S[j] - D) * hL[j] - (S[j - 1] - D) * hL[j - 1]
                   - feed_j * h_F - V_new[j] * (hL[j - 1] - hV[j]))
            V_new[j + 1] = max(num / denom, v_floor)
        return V_new
