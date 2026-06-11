"""Shortcut distillation column: Fenske-Underwood-Gilliland (FUG) design.

The classic shortcut method for a simple column (one feed, two products,
constant relative volatility, constant molal overflow):

* **Fenske (1932)** — minimum stages at total reflux from the key-component
  recoveries::

      N_min = ln[(d_LK/b_LK)(b_HK/d_HK)] / ln(alpha_LK)

  Non-key components are then distributed with the same equation at N_min:
  ``d_i/b_i = (d_HK/b_HK) * alpha_i^N_min``.

* **Underwood (1948)** — minimum reflux. First the root(s) theta of::

      sum_i alpha_i z_Fi / (alpha_i - theta) = 1 - q

  with alpha_HK < theta < alpha_LK (Brent's method between the poles), then::

      R_min + 1 = sum_i alpha_i x_Di / (alpha_i - theta)

  The feed quality q comes from the feed's thermal state at column pressure
  (an isenthalpic flash to P): q = 1 - vapor fraction.

* **Gilliland (1940)**, in the analytic Molokanov et al. (1972) form — actual
  stages at R = rr_factor * R_min::

      X = (R - R_min)/(R + 1)
      Y = 1 - exp[ (1 + 54.4 X)/(11 + 117.2 X) * (X - 1)/sqrt(X) ]
      N = (N_min + Y)/(1 - Y)

  N counts equilibrium stages including the (partial) reboiler but not a total
  condenser.

* **Kirkbride (1944)** — feed-stage location::

      N_R/N_S = [ (B/D) (z_HK,F/z_LK,F) (x_B,LK/x_D,HK)^2 ]^0.206

Relative volatilities alpha_i = K_i/K_HK come from the flowsheet's property
package, evaluated at the top (distillate composition) and bottom (bottoms
composition) of the column at P, combined as a geometric mean, and iterated to
self-consistency with the Fenske component distribution.

Products leave at column pressure ``P``: the distillate at its bubble point
(total condenser) or dew point (``partial_condenser=True``), the bottoms at its
bubble point. The condenser duty condenses the overhead vapor V = D(1+R)
(Heater sign convention: heat added is positive, so Q_cond <= 0); the reboiler
duty closes the overall energy balance exactly::

    F h_F + Q_cond + Q_reb = D h_D + B h_B

A ``dP`` param is accepted for `.flow` compatibility but ignored — the shortcut
treats the column at uniform pressure.

All design results (N_min, R_min, R, N, feed stage, component splits, duties)
are stored on the unit's ``design`` attribute after each solve.
"""
import math
from typing import Any

from scipy.optimize import brentq

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register


class ShortcutColumnError(ValueError):
    """Specification or feasibility error in a ShortcutColumn (bad keys,
    inverted volatility, infeasible recoveries, no Underwood root, ...)."""


# Fractions of the bubble->dew span tried when hunting a two-phase flash for
# K-value evaluation (the midpoint almost always works).
_TWO_PHASE_FRACS = (0.5, 0.25, 0.75, 0.1, 0.9)
_EXP_CLAMP = 35.0          # caps Fenske split factors at ~1e15 (avoids overflow)


def _k_values(pp, P: float, z: dict[str, float]) -> dict[str, float]:
    """K_i = y_i/x_i for the positive components of ``z``, from a two-phase
    flash of that composition at pressure ``P`` (between bubble and dew)."""
    bub, dew = pp.bubble_dew(P, z)
    span = max(dew - bub, 1e-6)
    for frac in _TWO_PHASE_FRACS:
        res = pp.flash_pt(bub + frac * span, P, z)
        if res.x and res.y and 1e-9 < res.vapor_fraction < 1.0 - 1e-9:
            return {c: res.y[c] / res.x[c]
                    for c, v in z.items() if v > 0.0 and res.x.get(c, 0.0) > 0.0}
    raise ShortcutColumnError(
        f"could not find a two-phase state at P={P:.4g} Pa for composition {z} "
        f"(bubble={bub:.2f} K, dew={dew:.2f} K); K-values unavailable"
    )


def _normalized(flows: dict[str, float]) -> dict[str, float]:
    total = sum(flows.values())
    return {c: v / total for c, v in flows.items()} if total > 0 else {}


@register("ShortcutColumn")
class ShortcutColumn(UnitOp):
    """Shortcut (FUG) distillation column. See the module docstring for the
    governing equations.

    Params (JSON-friendly scalars; ``.flow`` round-trips):
      * ``light_key`` / ``heavy_key`` — component ids (required).
      * ``recovery_light`` — fraction of LK recovered to distillate (0.99).
      * ``recovery_heavy`` — fraction of HK recovered to bottoms (0.99).
      * ``rr_factor`` — R/R_min, must be > 1 (default 1.3).
      * ``P`` — column pressure, Pa (default: feed pressure).
      * ``partial_condenser`` — vapor distillate at its dew point (default
        False: total condenser, liquid distillate at its bubble point).
    """

    design: dict[str, Any] | None = None

    def define_ports(self) -> list[Port]:
        return [
            Port("in1", "inlet"),
            Port("distillate", "outlet"),
            Port("bottoms", "outlet"),
            Port("condenser_duty", "outlet", "energy"),
            Port("reboiler_duty", "outlet", "energy"),
        ]

    # -- parameter validation ------------------------------------------------
    def _read_params(self, comps: list[str], z: dict[str, float], P_in: float):
        lk = self.params.get("light_key")
        hk = self.params.get("heavy_key")
        if not lk or not hk:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: both 'light_key' and 'heavy_key' are "
                f"required (got light_key={lk!r}, heavy_key={hk!r})"
            )
        for name, key in (("light_key", lk), ("heavy_key", hk)):
            if key not in comps:
                raise ShortcutColumnError(
                    f"ShortcutColumn {self.id!r}: {name}={key!r} is not in the "
                    f"flowsheet component list {comps}"
                )
            if z.get(key, 0.0) <= 0.0:
                raise ShortcutColumnError(
                    f"ShortcutColumn {self.id!r}: {name}={key!r} has zero flow in "
                    f"the feed; pick keys present in the feed"
                )
        if lk == hk:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: light and heavy key are both {lk!r}; "
                f"they must be different components"
            )

        r_lk = float(self.params.get("recovery_light", 0.99))
        r_hk = float(self.params.get("recovery_heavy", 0.99))
        for name, r in (("recovery_light", r_lk), ("recovery_heavy", r_hk)):
            if not 0.0 < r < 1.0:
                raise ShortcutColumnError(
                    f"ShortcutColumn {self.id!r}: {name}={r} must lie strictly in "
                    f"(0, 1) — a recovery of exactly 0 or 1 needs infinite stages"
                )
        if r_lk * r_hk / ((1.0 - r_lk) * (1.0 - r_hk)) <= 1.0:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: recoveries (LK {r_lk}, HK {r_hk}) give "
                f"a Fenske separation factor <= 1 — no enrichment is asked of the "
                f"column; raise the recoveries"
            )

        rr = float(self.params.get("rr_factor", 1.3))
        if rr <= 1.0:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: rr_factor={rr} must exceed 1 "
                f"(R = rr_factor * R_min; at or below R_min the column needs "
                f"infinite stages)"
            )

        P = float(self.params.get("P") or P_in)
        if P <= 0.0:
            raise ShortcutColumnError(f"ShortcutColumn {self.id!r}: P={P} Pa must be > 0")
        return lk, hk, r_lk, r_hk, rr, P, bool(self.params.get("partial_condenser", False))

    # -- solve ----------------------------------------------------------------
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        inlet = inlets.get("in1")
        if inlet is None or not inlet.molar_flow:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: missing or empty inlet on 'in1'"
            )
        T_in, P_in, n = inlet.require_state()
        z = inlet.normalized_z()
        comps = list(inlet.components)
        lk, hk, r_lk, r_hk, rr, P, partial = self._read_params(comps, z, P_in)

        H_in = inlet.H if inlet.H is not None else pp.enthalpy(T_in, P_in, z)

        # Feed quality q from the feed's thermal state at column pressure
        # (isenthalpic flash to P): q = 1 - vapor fraction.
        feed_res = pp.flash_ph(P, H_in, z)
        q = 1.0 - feed_res.vapor_fraction

        f = {c: n * z.get(c, 0.0) for c in comps}      # component feed flows, mol/s
        active = [c for c in comps if f[c] > 0.0]

        # K-values at feed conditions: volatility-orders the components and
        # checks the keys make sense before any design math.
        K_feed = _k_values(pp, P, {c: z.get(c, 0.0) for c in comps})
        if K_feed[lk] <= K_feed[hk]:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: light key {lk!r} (K={K_feed[lk]:.4g}) is "
                f"not more volatile than heavy key {hk!r} (K={K_feed[hk]:.4g}) at "
                f"P={P:.4g} Pa; swap the keys or check the pressure"
            )
        a_feed = {c: K_feed[c] / K_feed[hk] for c in active}

        # Key splits are fixed by the recovery specs.
        d = {c: 0.0 for c in active}
        d[lk] = r_lk * f[lk]
        d[hk] = (1.0 - r_hk) * f[hk]
        # Initial non-key guess from feed volatilities.
        for c in active:
            if c in (lk, hk):
                continue
            if a_feed[c] >= a_feed[lk]:
                d[c] = f[c]                  # lighter than LK -> distillate
            elif a_feed[c] <= 1.0:
                d[c] = 0.0                   # heavier than HK -> bottoms
            else:
                d[c] = 0.5 * f[c]            # sandwich component

        # Fenske separation factor and d_HK/b_HK from the key recoveries.
        sep = (r_lk / (1.0 - r_lk)) * (r_hk / (1.0 - r_hk))
        db_hk = (1.0 - r_hk) / r_hk

        # Iterate alpha <-> non-key distribution to self-consistency (alphas
        # move little, so a few passes suffice).
        alpha = dict(a_feed)
        n_min = 0.0
        for _ in range(3):
            b = {c: f[c] - d[c] for c in active}
            a_top = self._alphas(pp, P, _normalized(d), hk)
            a_bot = self._alphas(pp, P, _normalized(b), hk)
            for c in active:
                at, ab = a_top.get(c), a_bot.get(c)
                if at and ab:
                    alpha[c] = math.sqrt(at * ab)      # geometric mean top/bottom
                else:
                    alpha[c] = at or ab or a_feed[c]
            if alpha[lk] <= 1.0:
                raise ShortcutColumnError(
                    f"ShortcutColumn {self.id!r}: relative volatility "
                    f"alpha_LK={alpha[lk]:.4g} <= 1 between {lk!r} and {hk!r} at "
                    f"P={P:.4g} Pa (azeotrope or inverted keys); this separation is "
                    f"infeasible by simple distillation"
                )
            n_min = math.log(sep) / math.log(alpha[lk])
            for c in active:                            # Fenske non-key distribution
                if c in (lk, hk):
                    continue
                t = n_min * math.log(alpha[c]) + math.log(db_hk)
                t = max(min(t, _EXP_CLAMP), -_EXP_CLAMP)
                d[c] = f[c] / (1.0 + math.exp(-t))

        b = {c: f[c] - d[c] for c in active}
        D, B = sum(d.values()), sum(b.values())
        x_d, x_b = _normalized(d), _normalized(b)
        zf = {c: z[c] for c in active}

        # -- Underwood ------------------------------------------------------
        theta, r_min = self._underwood(zf, x_d, alpha, q, lk)
        if r_min <= 0.0:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: Underwood gave R_min={r_min:.4g} <= 0 — "
                f"the requested split happens without reflux; check keys/recoveries"
            )
        R = rr * r_min

        # -- Gilliland (Molokanov) and Kirkbride ------------------------------
        X = (R - r_min) / (R + 1.0)
        Y = 1.0 - math.exp((1.0 + 54.4 * X) / (11.0 + 117.2 * X)
                           * (X - 1.0) / math.sqrt(X))
        N = (n_min + Y) / (1.0 - Y)

        ratio = ((B / D) * (zf[hk] / zf[lk]) * (x_b[lk] / x_d[hk]) ** 2) ** 0.206
        n_rect = N * ratio / (1.0 + ratio)             # stages above the feed
        feed_stage = max(1, min(round(n_rect), round(N) - 1))

        # -- product streams ---------------------------------------------------
        z_dist = {c: x_d.get(c, 0.0) for c in comps}
        z_bot = {c: x_b.get(c, 0.0) for c in comps}
        bub_d, dew_d = pp.bubble_dew(P, z_dist)
        res_d = pp.flash_pt(dew_d if partial else bub_d, P, z_dist)
        bub_b, _ = pp.bubble_dew(P, z_bot)
        res_b = pp.flash_pt(bub_b, P, z_bot)

        distillate = Stream(
            id=f"{self.id}.distillate", components=comps,
            T=res_d.T, P=P, molar_flow=D, z=z_dist,
            H=res_d.H, phase=res_d.phase, vapor_fraction=res_d.vapor_fraction,
        )
        bottoms = Stream(
            id=f"{self.id}.bottoms", components=comps,
            T=res_b.T, P=P, molar_flow=B, z=z_bot,
            H=res_b.H, phase=res_b.phase, vapor_fraction=res_b.vapor_fraction,
        )

        # -- duties (Heater sign convention: positive heats the process) -------
        V_top = D * (1.0 + R)                          # overhead vapor to condenser
        h_vap_top = pp.enthalpy(dew_d, P, z_dist)      # saturated overhead vapor
        if partial:
            # Only the reflux L = R*D is condensed; distillate leaves as vapor.
            h_liq_top = pp.enthalpy(bub_d, P, z_dist)
            q_cond = R * D * (h_liq_top - h_vap_top)
        else:
            q_cond = V_top * (res_d.H - h_vap_top)     # condense V to bubble point
        # Reboiler closes the overall energy balance exactly (enthalpies are on
        # the engine's absolute, formation-inclusive basis):
        #   F h_F + Q_cond + Q_reb = D h_D + B h_B
        q_reb = D * res_d.H + B * res_b.H - n * H_in - q_cond

        self.design = {
            "N_min": n_min, "R_min": r_min, "R": R, "N": N,
            "feed_stage": feed_stage, "rectifying_stages": n_rect,
            "theta": theta, "q": q, "alpha": dict(alpha), "P": P,
            "D": D, "B": B,
            "distillate_flows": dict(d), "bottoms_flows": dict(b),
            "x_D": dict(x_d), "x_B": dict(x_b),
            "T_top": res_d.T, "T_top_dew": dew_d, "T_bottom": res_b.T,
            "V_top": V_top, "Q_condenser": q_cond, "Q_reboiler": q_reb,
            "partial_condenser": partial,
        }
        return {
            "distillate": distillate,
            "bottoms": bottoms,
            "condenser_duty": EnergyStream(id=f"{self.id}.condenser_duty", duty=q_cond),
            "reboiler_duty": EnergyStream(id=f"{self.id}.reboiler_duty", duty=q_reb),
        }

    # -- helpers ---------------------------------------------------------------
    @staticmethod
    def _alphas(pp, P: float, comp: dict[str, float], hk: str) -> dict[str, float]:
        """alpha_i = K_i/K_HK at the two-phase state of ``comp`` at P. The heavy
        key is always present in both products (recoveries < 1), so K_HK exists."""
        if not comp:
            return {}
        K = _k_values(pp, P, comp)
        return {c: v / K[hk] for c, v in K.items()}

    def _underwood(self, zf: dict[str, float], x_d: dict[str, float],
                   alpha: dict[str, float], q: float, lk: str) -> tuple[float, float]:
        """Solve sum(alpha_i z_i/(alpha_i - theta)) = 1 - q for theta between
        alpha_HK (=1) and alpha_LK, then R_min from the distillate composition.

        With sandwich components (alphas between the keys) the interval has
        interior poles; each pole-bounded subinterval is searched and the most
        conservative (largest) R_min is used — a documented simplification of the
        full multi-theta Underwood system.
        """
        a_lk = alpha[lk]

        def phi(theta: float) -> float:
            return sum(alpha[c] * zf[c] / (alpha[c] - theta) for c in zf) - (1.0 - q)

        poles = sorted({alpha[c] for c in zf if 1.0 < alpha[c] < a_lk})
        bounds = [1.0, *poles, a_lk]
        eps = 1e-9 * (a_lk - 1.0)
        thetas: list[float] = []
        for lo, hi in zip(bounds, bounds[1:]):
            a, bnd = lo + eps, hi - eps
            if a >= bnd:
                continue
            fa, fb = phi(a), phi(bnd)
            if fa == 0.0:
                thetas.append(a)
            elif fa * fb < 0.0:
                thetas.append(float(brentq(phi, a, bnd, xtol=1e-12)))
        if not thetas:
            raise ShortcutColumnError(
                f"ShortcutColumn {self.id!r}: no Underwood root found between "
                f"alpha_HK=1 and alpha_LK={a_lk:.4g} (q={q:.4g}); check the feed "
                f"thermal state and key selection"
            )

        best_theta, best_rmin = thetas[0], -math.inf
        for theta in thetas:
            rmin = sum(alpha[c] * x_d.get(c, 0.0) / (alpha[c] - theta) for c in zf) - 1.0
            if rmin > best_rmin:
                best_theta, best_rmin = theta, rmin
        return best_theta, best_rmin
