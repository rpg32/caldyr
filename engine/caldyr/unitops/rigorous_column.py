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

**Multiple feeds**: pass ``feeds=[{"stage": j1}, {"stage": j2}, ...]`` and the
column exposes one inlet port per entry (``in1``, ``in2``, ...), each feed
landing on its stage; without the param the unit keeps its classic single
``in1`` + ``feed_stage`` form, byte-for-byte compatible with existing
flowsheets. **Side draws**: ``side_draws=[{"stage": j, "phase": "liquid" |
"vapor", "rate": mol/s}, ...]`` adds outlet ports ``side1``, ``side2``, ...
drawing saturated liquid (or vapor) off the stage at exactly the given molar
rate. Feed terms enter the tridiagonal right-hand side on their stages and
draw terms the diagonal / total balances (Seader 3e eqs. 10-8..10-12 with
U_j/W_j nonzero), so the same Thomas solve covers every configuration.

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

Absorber/stripper mode (no condenser/reboiler) lives in
:mod:`caldyr.unitops.absorber` (sum-rates method — the bubble-point method
here is for narrow/medium-boiling distillation). Out of scope (this round):
pumparounds.

**Steam-stripped main fractionator** (``reboiled=False``): a crude atmospheric
tower (after Hameed 2025 sec. 10.2.1) has a condenser but NO reboiler — the
bottom vapor is open stripping steam fed on the bottom stage, and the tower has
liquid side draws and pumparounds (modeled as ``stage_duties``). The single
specification is the distillate rate; the reflux ratio is an OUTPUT
(``reflux_ratio`` only seeds the initial traffic), reported in ``design['R']``.

The recommended method for this form is ``method="bubble_point"`` (the default).
Each stage's temperature is the K-value bubble point of its liquid
(:func:`_bubble_T_kvalues` — a secant on ln(sum K x) with phi-phi K-values, no
flash object, so it is fast and immune to the PVF-flash failure modes on
wide-boiling pseudo-component + water liquids). The vapor traffic comes from
*envelope* energy balances (each V_j from its own around-stage-j-down-to-the-
bottom balance) on a sensible basis: each molar enthalpy has its
composition-weighted formation-enthalpy offset subtracted, which keeps the
envelope pivot hLs_{j-1} - hVs_j ~ -(latent heat) < 0 even across the steam
stage (where the raw formation-inclusive offsets — water -242 kJ/mol vs the
pseudo-components 0 — would make it sign-indefinite). Both the temperature and
the traffic updates are under-relaxed, with the damping reduced whenever the
combined temperature+traffic change stalls, which is what lets the iteration
escape the period-2 limit cycle the open-steam bottom stage otherwise sits in.

A second method, ``method="sum_rates"`` (Burningham-Otto: liquid traffic from
the component-flow sums, temperatures from a simultaneous Newton solve of the
stage energy balances), is also provided but is **not recommended** — it
limit-cycles on equilibrium-dominated sections (the documented Seader 3e ch.
10.4 division of labor). Use ``bubble_point``.

The condenser stage is pinned to saturation (the bubble temperature of the
distillate liquid), its energy balance absorbed by the condenser duty. The
bottoms leave as liquid at the converged bottom-stage temperature (a
steam-stripped liquid is below its dry-EOS bubble point, so no saturation flash
is imposed). Mass balances close machine-exact by difference; the condenser
duty closes the overall energy balance exactly, and the independently computed
condenser-stage balance is reported in ``design['energy_residual_rel']``.

**Property-package range.** Because the energy balances need the saturated
vapor enthalpy to exceed the liquid's (positive latent heat), this mode is only
as good as the property package's pseudo-component enthalpies: with the cubic-
EOS ``thermo`` backends, heavy petroleum pseudo-components (NBP above roughly
the light-gas-oil range) can return an unphysical hV < hL from a bubble-point
flash, which no MESH method can resolve. The crude-tower example and test
therefore use a *light* crude (naphtha-through-light-gas-oil), whose cuts stay
in the physically valid range; a full resid-bearing crude needs a property
backend with corrected pseudo-component enthalpies (out of scope here).
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import brentq

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .reaction import KineticReaction


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
_RXN_RELAX = 0.5         # under-relaxation on the lagged reaction source (reactive
#                          distillation) — eases the generation in as the profile
#                          settles; at convergence the used == newly-computed rate
_RXN_MAX_CONV = 0.95     # per-stage extent cap: fraction of a reactant's liquid
#                          traffic a single stage may consume (flash conditioning)

# sum-rates mode (steam-stripped / non-reboiled main fractionator)
_SR_TOL_FLOW = 1e-7      # relative traffic convergence
_SR_ENERGY_TOL = 1e-9    # relative stage-energy residual required at convergence
_SR_NEWTON_MAX = 30      # inner Newton iterations on the energy balances
_SR_NEWTON_TOL = 1e-10   # inner Newton relative residual target
_SR_STEP_MAX = 15.0      # K, per-stage Newton step clamp
_SR_DH_DT = 0.05         # K, finite-difference step for dh/dT
_SR_T_OUTER_MAX = 30.0   # K, per-stage temperature move per OUTER iteration
_SR_L_RATIO = 4.0        # max factor the liquid traffic may move per iteration
_SR_DUTY_RAMP = 12       # iterations over which stage duties ramp to full
_BUBBLE_T_LO = 150.0     # K, admissible window for the K-value bubble secant
_BUBBLE_T_HI = 1500.0

# inside-out mode (Boston-Sullivan; the robust wide-boiling / resid method)
_IO_OUTER_MAX = 60       # outer (rigorous model-refit) iterations cap
_IO_OUTER_TOL_T = 1e-3   # K, outer stage-temperature convergence
_IO_OUTER_TOL_F = 1e-5   # relative outer traffic convergence
_IO_B_LO = 300.0         # K, clamp on the base-K Clausius slope ln Kb = a - b/T
_IO_B_HI = 60000.0
_IO_B_DEFAULT = 6000.0   # K, fallback slope when the finite difference is noisy
_IO_DH_DT = 0.1          # K, finite-difference step for the model slopes
_IO_DUTY_RAMP = 8        # outer iterations over which stage duties ramp to full
# inner damped-Newton (on the 2N-1 tear variables, frozen inside models)
_IO_NEWTON_MAX = 40      # inner Newton iterations per outer pass
_IO_NEWTON_TOL = 1e-9    # inner residual infinity-norm convergence
_IO_LS_MAX = 16          # backtracking line-search halvings
_IO_DT_PERT = 0.02       # K, Jacobian finite-difference step on temperatures
_IO_DV_FRAC = 1e-4       # relative Jacobian step on the vapour flows
_IO_DV_MIN = 1e-3        # mol/s, floor on the vapour-flow Jacobian step
_IO_NEWTON_ACCEPT = 1e-5  # inner residual required before outer convergence
_IO_TR_DT = 60.0         # K, trust-region half-width on the stage temperatures
_IO_TR_VR = 4.0          # trust-region ratio bound on the vapour flows
_IO_RCOND = 1e-8         # truncated-SVD cutoff for the regularized Newton step

# Naphtali-Sandholm simultaneous-correction method (per-component flow variables)
_NS_TOL = 1e-9           # scaled residual infinity-norm convergence
_NS_MAX = 60             # Newton iterations cap (param 'max_iter' overrides)
_NS_FLOOR = 1e-12        # mol/s, floor on the component flows
_NS_DL_FRAC = 1e-6       # relative finite-difference step on the flows
_NS_DL_MIN = 1e-10       # mol/s, floor on the flow finite-difference step
_NS_DT = 0.01            # K, finite-difference step on the temperatures
_NS_LS_MAX = 24          # backtracking line-search halvings
_NS_TRUST = 0.5          # max fractional reduction a flow may take in one step
_NS_STEP_REL = 1e-4      # only flows above this fraction of the stage total
#                          constrain the fractional Newton step
_NS_TOL_REBOILED = 1e-8  # reboiled/total-condenser NS scaled-residual tol: the
#                          extra global reboiler-duty unknown + the appended
#                          bubble-point row floor the float64 endgame near ~2e-9
#                          (the same marginality as the resid tower); 1e-8 is
#                          still very tight -- products match the bubble-point
#                          method to ~5 digits and the balances close ~1e-8
_NS_WARM_ITERS = 40      # inside-out homotopy warm-start iteration cap
_NS_CONT_STEP0 = 0.3     # draw-rate continuation (broader-basin fallback for
#                          draw-heavy towers): initial fraction step, ramping the
#                          side draws to 100% of target, warm-starting each NS
#                          solve from the previous; the step halves on a failed
#                          target and grows on success
_NS_CONT_GROW = 1.5      # ... step growth factor after a converged target
_NS_CONT_STEP_MIN = 0.02  # ... give up if the step must shrink below this
_NS_CONT_IT = 40         # ... per-step iteration cap for intermediate targets
#                          (a warm-started step converges in ~7-11 Newton iters
#                          or not at all; only the final target gets max_iter)
_NS_LM_INIT = 1e-3       # Levenberg-Marquardt damping: initial value
_NS_LM_MIN = 1e-10       # ... lower bound (-> Gauss-Newton, quadratic endgame)
_NS_LM_MAX = 1e10        # ... upper bound (-> short scaled gradient descent)
_NS_LM_UP = 4.0          # ... growth factor on a rejected trial
_NS_LM_DOWN = 3.0        # ... shrink factor on an accepted step
_NS_REFINE = 2           # iterative-refinement passes on the Newton linear
#                          solve (recovers digits lost to the ill-conditioned
#                          MESH Jacobian, so the endgame plunge is reproducible)
_DECANT_BETA_MIN = 1e-4  # a layer thinner than this is treated as "no split" by
#                          the decanting condenser (degenerate single-liquid
#                          overhead) -> fall back to a fixed-fraction reflux
_DECANT_FALLBACK_R = 2.0  # reflux ratio used while the decant has not yet formed
#                          (early Newton iterates), so the solver stays in a sane
#                          internal-traffic regime until the two-liquid split sets


class _NoVLEError(ValueError):
    """Internal: the feed has no two-phase state at the given pressure."""


def fenske_profile(pp, P: float, f: dict[str, float], D: float,
                   N: int, K_feed: dict[str, float] | None = None,
                   ) -> list[dict[str, float]]:
    """Initial stage liquid-composition profile for an N-stage column from a
    Fenske-style split: d_i/b_i proportional to alpha_i^N_min with the split
    factor solved so sum(d) = D (Fenske 1932; Seader 3e eq. 9-12 rearranged),
    then linear interpolation between the estimated product compositions —
    the shortcut initialization the tray-by-tray methods call for.

    ``f`` carries the (total) component feed flows in mol/s; ``D`` the
    overhead product rate. ``K_feed`` optionally supplies the feed K-values
    (e.g. from an already-computed feed flash) so no bubble-point flash is
    needed here. Shared by the RigorousColumn and the ReboiledAbsorber
    initializers.
    """
    active = list(f)
    if K_feed is None:
        res = pp.bubble_point(P, {c: v / sum(f.values()) for c, v in f.items()})
        if res.y is None or res.x is None:
            raise _NoVLEError(f"no VLE for the feed at P={P:.4g} Pa")
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


def _bubble_T_kvalues(pp, P: float, x: dict[str, float], T0: float,
                      ) -> tuple[float, dict[str, float], dict[str, float]]:
    """Bubble temperature of liquid ``x`` at ``P`` by a secant iteration on
    ``ln(sum_i K_i x_i) = 0`` with phi-phi K-values (``pp.k_values``) and the
    incipient vapor refreshed as ``y_i = K_i x_i / sum`` each evaluation — the
    classical K-value bubble point (Seader 3e sec. 4.4). No flash object is
    built, so it is both fast (two lnphi evaluations per step) and immune to
    the PVF-flash failure modes on wide-boiling pseudo-component + water
    liquids. Returns ``(T, K, y)``.
    """
    T = min(max(float(T0), _BUBBLE_T_LO + 1.0), _BUBBLE_T_HI - 1.0)
    y = {c: max(v, _X_FLOOR) for c, v in x.items()}
    tot = sum(y.values())
    y = {c: v / tot for c, v in y.items()}

    def eval_at(t: float):
        # settle the incipient-vapor composition so f(T) is single-valued
        nonlocal y
        s, K = 1.0, {}
        for _ in range(4):
            K = pp.k_values(t, P, x, y)
            s = sum(K[c] * x[c] for c in x)
            y_new = {c: max(K[c] * x[c], _X_FLOOR) / s for c in x}
            dy = max(abs(y_new[c] - y[c]) for c in x)
            y = y_new
            if dy < 1e-12:
                break
        return math.log(s), K

    # f = ln(sum K x) rises with T (more volatility); bracket the zero by
    # marching, then close in with Illinois regula falsi.
    f, K = eval_at(T)
    if abs(f) <= 1e-9:
        return T, K, y
    t_a, f_a = T, f
    step = 25.0 if f < 0.0 else -25.0
    bracketed = False
    for _ in range(80):
        t_b = min(max(t_a + step, _BUBBLE_T_LO), _BUBBLE_T_HI)
        f_b, K = eval_at(t_b)
        if abs(f_b) <= 1e-9:
            return t_b, K, y
        if (f_a < 0.0) != (f_b < 0.0):
            bracketed = True
            break
        if t_b in (_BUBBLE_T_LO, _BUBBLE_T_HI):
            raise RigorousColumnError(
                f"no K-value bubble point for the stage liquid at "
                f"P={P:.4g} Pa within {_BUBBLE_T_LO:.0f}-{_BUBBLE_T_HI:.0f} K "
                f"(ln sum Kx = {f_b:.3g} at the bound)"
            )
        t_a, f_a = t_b, f_b
    if not bracketed:
        raise RigorousColumnError(
            f"K-value bubble-point bracketing failed at P={P:.4g} Pa "
            f"(last ln sum Kx = {f_b:.3g} at T={t_b:.1f} K)"
        )
    if t_a > t_b:
        t_a, f_a, t_b, f_b = t_b, f_b, t_a, f_a
    side = 0
    t_m, K_m = t_a, K
    for _ in range(80):
        t_m = t_a + (t_b - t_a) * (-f_a) / (f_b - f_a)
        if not (t_a < t_m < t_b):
            t_m = 0.5 * (t_a + t_b)
        f_m, K_m = eval_at(t_m)
        if abs(f_m) <= 1e-9 or (t_b - t_a) < 1e-8:
            return t_m, K_m, y
        if (f_m < 0.0) == (f_a < 0.0):
            t_a, f_a = t_m, f_m
            if side == -1:
                f_b *= 0.5            # Illinois: halve the stagnant end
            side = -1
        else:
            t_b, f_b = t_m, f_m
            if side == 1:
                f_a *= 0.5
            side = 1
    return t_m, K_m, y


def _newton_T_draws(
    pp, N: int, P_j: list[float],
    x: list[dict[str, float]], y: list[dict[str, float]],
    L: list[float], V: list[float],
    fh_stage: list[float], u: list[float], w: list[float],
    T0: list[float], t_bounds: tuple[float, float],
) -> tuple[list[float], list[float], list[float], float, int]:
    """Solve the stage energy balances of stages 2..N (0-based 1..N-1) for
    their temperatures by Newton's method with a tridiagonal Jacobian — the
    Burningham-Otto temperature step (Seader 3e eqs. 10-65..10-67) extended
    with liquid/vapor side draws (``u``/``w``) and a *fixed* stage-1
    temperature (the condenser, pinned to saturation; its balance is absorbed
    by the condenser duty). Best-effort: returns the lowest-residual iterate —
    the outer sum-rates loop requires a tight residual before declaring
    convergence. Returns ``(T, hL, hV, residual_rel, n_prop_calls)``.
    """
    T = list(T0)
    t_lo, t_hi = t_bounds
    n_prop = 0
    best: tuple[float, list[float], list[float], list[float]] | None = None
    for _ in range(_SR_NEWTON_MAX):
        hL = [pp.enthalpy_liquid(T[j], P_j[j], x[j]) for j in range(N)]
        hV = [pp.enthalpy_vapor(T[j], P_j[j], y[j]) for j in range(N)]
        cpL = [0.0] + [
            (pp.enthalpy_liquid(T[j] + _SR_DH_DT, P_j[j], x[j]) - hL[j]) / _SR_DH_DT
            for j in range(1, N)]
        cpV = [0.0] + [
            (pp.enthalpy_vapor(T[j] + _SR_DH_DT, P_j[j], y[j]) - hV[j]) / _SR_DH_DT
            for j in range(1, N)]
        n_prop += 4 * N - 2

        H = [0.0] * N
        scale = [1.0] * N
        for j in range(1, N):
            inflow = (L[j - 1] * hL[j - 1]
                      + (V[j + 1] * hV[j + 1] if j < N - 1 else 0.0)
                      + fh_stage[j])
            outflow = (L[j] + u[j]) * hL[j] + (V[j] + w[j]) * hV[j]
            H[j] = inflow - outflow
            scale[j] = max(abs((L[j] + u[j]) * hL[j]),
                           abs((V[j] + w[j]) * hV[j]),
                           abs(fh_stage[j]), 1.0)
        resid_rel = max(abs(H[j]) / scale[j] for j in range(1, N))
        if best is None or resid_rel < best[0]:
            best = (resid_rel, list(T), hL, hV)
        if resid_rel <= _SR_NEWTON_TOL:
            return T, hL, hV, resid_rel, n_prop

        # Tridiagonal Newton system over stages 1..N-1 (0-based); the
        # coupling of stage 1 to the fixed condenser temperature is dropped.
        a = [0.0] + [L[j - 1] * cpL[j - 1] for j in range(2, N)]
        b = [-((L[j] + u[j]) * cpL[j] + (V[j] + w[j]) * cpV[j])
             for j in range(1, N)]
        c = [V[j + 1] * cpV[j + 1] for j in range(1, N - 1)] + [0.0]
        dT = _thomas(a, b, c, [-H[j] for j in range(1, N)])
        T = [T[0]] + [
            min(max(tj + max(min(d, _SR_STEP_MAX), -_SR_STEP_MAX), t_lo), t_hi)
            for tj, d in zip(T[1:], dT)]
    assert best is not None
    resid_rel, T, hL, hV = best
    return T, hL, hV, resid_rel, n_prop


@register("RigorousColumn")
class RigorousColumn(UnitOp):
    """Rigorous (MESH, bubble-point) distillation column. See the module
    docstring for the algorithm and the stage-numbering convention.

    Params (JSON-friendly scalars; ``.flow`` round-trips):
      * ``n_stages`` — theoretical stages, **including** the condenser
        (stage 1) and the reboiler (stage n_stages). Required, >= 3.
      * ``feed_stage`` — feed stage in the same numbering; 2..n_stages-1
        (single-feed form; the feed arrives on port ``in1``).
      * ``feeds`` — ``[{"stage": j}, ...]`` for multiple feeds: one inlet
        port per entry (``in1``, ``in2``, ...), each on its stage
        (2..n_stages-1). Mutually exclusive with ``feed_stage``.
      * ``side_draws`` — ``[{"stage": j, "phase": "liquid"|"vapor",
        "rate": mol/s}, ...]``: outlet ports ``side1``, ``side2``, ...
        drawing off stage j (2..n_stages-1) at exactly the given rate.
      * ``reflux_ratio`` — molar reflux ratio L_1/D, > 0. Required.
      * exactly one of ``distillate_rate`` (mol/s) or ``distillate_to_feed``
        (molar D/F fraction of the *total* feed).
      * ``P`` — top-stage pressure, Pa (default: feed pressure).
      * ``dP_stage`` — linear pressure rise per stage going down, Pa
        (default 0: uniform column pressure).
      * ``partial_condenser`` — vapor distillate at its dew point (default
        False: total condenser, liquid distillate at its bubble point).
      * ``max_iter`` — MESH iteration cap (default 200).
      * ``reactions`` — reactive distillation: a list of kinetic-reaction dicts
        (see :class:`caldyr.unitops.reaction.KineticReaction`, plus a
        ``"stages": [first, last]`` 1-based inclusive tray range) run on their
        stages, with ``tray_holdup`` the per-stage liquid holdup volume (m^3).
        The per-stage generation enters the component balances and the heat of
        reaction rides the formation-inclusive stage enthalpies. Reversible
        kinetics (``k0_rev``/``Ea_rev``) bound the conversion at equilibrium.
        Implemented for the default reboiled bubble-point column; see
        :meth:`_read_reactions`.
      * ``reboiled`` — default True. ``False`` removes the reboiler: stage
        ``n_stages`` becomes an ordinary (adiabatic) stage that may receive a
        feed — open stripping steam, the crude-tower configuration. The
        column then has a single specification (the distillate rate);
        ``reflux_ratio`` only seeds the initial traffic and the converged
        reflux ratio is reported in ``design['R']``. Requires a non-reboiled
        method (``"bubble_point"``, ``"inside_out"``).
      * ``method`` — the MESH temperature/traffic update:

        - ``"bubble_point"`` + ``reboiled=True`` (the **default**): the
          Wang-Henke bubble-point method — a saturated-liquid flash per stage
          sets the temperature and the vapor traffic comes from a forward
          energy recurrence. The narrow/medium-boiling distillation workhorse.
        - ``"bubble_point"`` + ``reboiled=False``: the same per-stage K-value
          bubble temperatures, but the vapor traffic comes from *envelope*
          energy balances (each V_j from its own around-stage-j-to-bottom
          balance) run on a formation-offset-conditioned sensible basis. The
          **robust method for a steam-stripped main fractionator** (open
          stripping steam, side draws, pumparound stage duties); it converges
          the crude atmospheric tower (see ``examples/20_crude_tower`` and
          ``tests/test_m15_crude_column``).
        - ``"naphtali_sandholm"`` + ``reboiled=False`` (partial condenser): the
          **Naphtali-Sandholm simultaneous-correction** method (see
          :meth:`_solve_ns`) — all MESH equations solved at once by a damped
          Newton with per-component flow variables. The **robust method for a
          full resid-bearing crude atmospheric tower** (the case bubble_point
          and inside_out cannot solve); warm-started by inside_out. Reproduces
          the bubble-point solution on the light tower and closes the resid
          tower to machine precision (see ``tests/test_m15_resid_column``).
        - ``"inside_out"`` + ``reboiled=False``: an inside-out method with a
          damped-Newton inner loop on the (T, V) tear variables (see
          :meth:`_inside_out_loops`). Converges narrow-to-medium wide-boiling
          steam-stripped towers and reproduces the bubble-point solution to
          machine precision on the light crude tower; on a full resid tower it
          recovers the physics but stalls (a draw-stage degeneracy that
          ``naphtali_sandholm`` resolves) — so it is used as the NS warm start.
        - ``"sum_rates"`` + ``reboiled=False``: a Burningham-Otto variant.
          **Provided but not recommended**: on equilibrium-dominated sections
          it limit-cycles (Seader 3e ch. 10.4). Prefer ``"bubble_point"``.

        ``reboiled=True`` with a non-default method is rejected.
      * ``stage_duties`` — ``[{"stage": j, "duty": W}, ...]``: heat added
        (negative = removed) directly on stage j (2..n_stages-1). This is
        the standard reduced model of a **pumparound**: the circulating
        liquid's net effect on the column is heat removal at the draw/return
        stages, which condenses the internal liquid the side draws need
        (Kister, *Distillation Design*, ch. 13 treats pumparounds as
        intermediate condensers). Requires ``reboiled=False`` (the reboiled
        bubble-point vapor recurrence does not carry stage duties); it is
        honored by both non-reboiled methods, which fold the (ramped) duty
        into each stage's feed-enthalpy flow.
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

    def warm_start_from(self, other: "RigorousColumn") -> None:
        """**Stage-count continuation** for the decanting-condenser column: seed
        this column's Naphtali-Sandholm warm start from another, *coarser*
        converged decant column by interpolating its stage profile onto this
        column's stage count.

        A cold solve of a long (50-62-stage) heteroazeotropic entrainer column is
        intractable — the entrainer-rich seed is too far from the basin and the
        FD-Jacobian Newton grinds for many minutes. But a short (~30-stage) column
        cold-solves in ~1 min, and its converged profile, linearly interpolated
        onto more stages, lands the long column squarely in its basin (it then
        re-solves in a handful of warm Newton iterations). This is the stage-count
        analogue of the distillate-rate continuation the entrainer examples
        already use; do ``col62.warm_start_from(col30)`` once, then drive the rest
        with the usual distillate-rate / inventory continuation.

        ``other`` must be a converged decant column (its ``_warm`` carries the
        decant profile). The feed composition signature is copied across, so the
        two columns must dehydrate the same mixture (the seed only needs to be in
        the basin, not identical — feed *stage* locations may differ)."""
        w = getattr(other, "_warm", None)
        if not w or w.get("method") != "decant":
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: warm_start_from needs a CONVERGED "
                f"decanting-condenser column to interpolate from (the source has "
                f"no decant warm profile — solve it first)"
            )
        try:
            n2 = int(self.params["n_stages"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: integer 'n_stages' is required "
                f"to interpolate a warm start onto"
            ) from exc
        active = list(w["active"])

        def _interp(src: list[float], j: int, n1: int) -> float:
            p = j * (n1 - 1) / (n2 - 1) if n2 > 1 else 0.0
            lo = int(p)
            hi = min(lo + 1, n1 - 1)
            fr = p - lo
            return src[lo] * (1.0 - fr) + src[hi] * fr

        n1 = int(w["n_stages"])
        rows = w["x"]
        x2: list[dict[str, float]] = []
        for j in range(n2):
            p = j * (n1 - 1) / (n2 - 1) if n2 > 1 else 0.0
            lo = int(p)
            hi = min(lo + 1, n1 - 1)
            fr = p - lo
            row = {c: rows[lo].get(c, 0.0) * (1.0 - fr) + rows[hi].get(c, 0.0) * fr
                   for c in active}
            tot = sum(row.values()) or 1.0
            x2.append({c: v / tot for c, v in row.items()})
        self._warm = {
            "n_stages": n2, "method": "decant", "active": active,
            "z": dict(w["z"]), "x": x2,
            "T": [_interp(w["T"], j, n1) for j in range(n2)],
            "L": [_interp(w["L"], j, n1) for j in range(n2)],
            "V": [_interp(w["V"], j, n1) for j in range(n2)],
        }

    def define_ports(self) -> list[Port]:
        feeds = self.params.get("feeds")
        draws = self.params.get("side_draws")
        if feeds is not None and not isinstance(feeds, list):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: 'feeds' must be a list of "
                f"{{'stage': j}} dicts (got {feeds!r})"
            )
        if draws is not None and not isinstance(draws, list):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: 'side_draws' must be a list of "
                f"{{'stage': j, 'phase': ..., 'rate': ...}} dicts "
                f"(got {draws!r})"
            )
        n_feeds = max(1, len(feeds or []))
        ports = [Port(f"in{i + 1}", "inlet") for i in range(n_feeds)]
        ports += [Port("distillate", "outlet"), Port("bottoms", "outlet")]
        ports += [Port(f"side{i + 1}", "outlet")
                  for i in range(len(draws or []))]
        ports += [Port("condenser_duty", "outlet", "energy"),
                  Port("reboiler_duty", "outlet", "energy")]
        return ports

    # -- parameter validation --------------------------------------------------
    def _read_params(self, F: float, P_in: float):
        try:
            n_stages = int(self.params["n_stages"])
        except (KeyError, TypeError, ValueError) as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: integer 'n_stages' is required "
                f"(got n_stages={self.params.get('n_stages')!r})"
            ) from exc
        if n_stages < 3:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: n_stages={n_stages} must be >= 3 — "
                f"the count includes the condenser (stage 1) and the reboiler "
                f"(stage n_stages), so 3 is one tray plus condenser and reboiler"
            )

        method = str(self.params.get("method", "bubble_point"))
        if method not in ("bubble_point", "sum_rates", "inside_out",
                          "naphtali_sandholm"):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: method={method!r} must be "
                f"'bubble_point', 'inside_out', 'naphtali_sandholm', or "
                f"'sum_rates'"
            )
        reboiled = bool(self.params.get("reboiled", True))
        if reboiled and method in ("sum_rates", "inside_out"):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: method={method!r} is "
                f"implemented for the non-reboiled (steam-stripped main "
                f"fractionator) form only — set reboiled=False, or use the "
                f"default bubble_point method for a reboiled column"
            )
        # method="naphtali_sandholm" + reboiled=True is the total-condenser +
        # reboiler simultaneous-correction form (specs D + reflux_ratio, same as
        # the bubble-point method; the reboiler/condenser duties are recovered).
        # a non-reboiled column's bottom stage is an ordinary stage and may
        # receive a feed (the stripping steam)
        max_feed_stage = n_stages if not reboiled else n_stages - 1

        feeds_param = self.params.get("feeds")
        if feeds_param is not None:
            if self.params.get("feed_stage") is not None:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: give either 'feed_stage' "
                    f"(single feed on in1) or 'feeds' (stages for in1, in2, "
                    f"...), not both"
                )
            if not feeds_param:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: 'feeds' must name at least "
                    f"one feed stage"
                )
            feed_stages = []
            for k, entry in enumerate(feeds_param):
                try:
                    stage = int(entry["stage"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise RigorousColumnError(
                        f"RigorousColumn {self.id!r}: feeds[{k}] needs an "
                        f"integer 'stage' (got {entry!r})"
                    ) from exc
                feed_stages.append(stage)
        else:
            try:
                feed_stages = [int(self.params["feed_stage"])]
            except (KeyError, TypeError, ValueError) as exc:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: integer 'feed_stage' is "
                    f"required (got "
                    f"feed_stage={self.params.get('feed_stage')!r})"
                ) from exc
        for stage in feed_stages:
            if not 2 <= stage <= max_feed_stage:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: feed_stage={stage} is out of "
                    f"range — it must lie between stage 2 (below the condenser) "
                    f"and stage {max_feed_stage} "
                    f"({'the bottom stage' if not reboiled else 'above the reboiler'}) "
                    f"for n_stages={n_stages}"
                )

        draws: list[tuple[int, str, float]] = []
        for k, entry in enumerate(self.params.get("side_draws") or []):
            try:
                stage = int(entry["stage"])
                rate = float(entry["rate"])
            except (KeyError, TypeError, ValueError) as exc:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: side_draws[{k}] needs an "
                    f"integer 'stage' and a numeric 'rate' (got {entry!r})"
                ) from exc
            phase = str(entry.get("phase", "liquid"))
            if phase not in ("liquid", "vapor"):
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: side_draws[{k}] phase="
                    f"{phase!r} must be 'liquid' or 'vapor'"
                )
            if not 2 <= stage <= n_stages - 1:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: side_draws[{k}] "
                    f"stage={stage} is out of range — draws must come off a "
                    f"tray (2..{n_stages - 1}); use 'distillate'/'bottoms' "
                    f"for the end stages"
                )
            if rate <= 0.0:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: side_draws[{k}] "
                    f"rate={rate} mol/s must be > 0"
                )
            draws.append((stage, phase, rate))

        duties: list[tuple[int, float]] = []
        for k, entry in enumerate(self.params.get("stage_duties") or []):
            try:
                stage = int(entry["stage"])
                duty = float(entry["duty"])
            except (KeyError, TypeError, ValueError) as exc:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: stage_duties[{k}] needs an "
                    f"integer 'stage' and a numeric 'duty' in W "
                    f"(got {entry!r})"
                ) from exc
            if reboiled:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: 'stage_duties' (pumparound "
                    f"heat) is supported on the non-reboiled form only — set "
                    f"reboiled=False (the classic reboiled vapor recurrence "
                    f"does not carry stage duties)"
                )
            if not 2 <= stage <= n_stages - 1:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: stage_duties[{k}] "
                    f"stage={stage} is out of range — duties go on a tray "
                    f"(2..{n_stages - 1}); the condenser duty is computed and "
                    f"a bottom-stage duty would be a reboiler"
                )
            duties.append((stage, duty))

        if not reboiled and self.params.get("reflux_ratio") is None:
            R = 1.0          # initial-traffic seed only; the converged R is an output
        else:
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
        draw_total = sum(rate for _, _, rate in draws)
        if draw_total > 0.0 and not D + draw_total < F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: distillate D={D:.6g} plus side "
                f"draws {draw_total:.6g} mol/s must stay below the total feed "
                f"F={F:.6g} mol/s (the bottoms must exist)"
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

        # Integrated decanting condenser (heteroazeotropic entrainer column).
        decant = bool(self.params.get("decant_condenser", False))
        condenser_T = 0.0
        reflux_organic = True
        if decant:
            if method != "naphtali_sandholm" or not reboiled:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: 'decant_condenser' requires "
                    f"method='naphtali_sandholm' and reboiled=True (the "
                    f"integrated decant is formulated on the simultaneous-"
                    f"correction reboiled column); got method={method!r}, "
                    f"reboiled={reboiled}"
                )
            if partial:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: 'decant_condenser' and "
                    f"'partial_condenser' are mutually exclusive (the decanting "
                    f"condenser is a total condenser that settles two liquids)"
                )
            try:
                condenser_T = float(self.params["condenser_T"])
            except (KeyError, TypeError, ValueError) as exc:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: 'decant_condenser' needs a "
                    f"numeric 'condenser_T' (K) — the temperature the overhead "
                    f"is condensed and decanted at (got "
                    f"{self.params.get('condenser_T')!r})"
                ) from exc
            if not _BUBBLE_T_LO < condenser_T < _BUBBLE_T_HI:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: condenser_T={condenser_T} K is "
                    f"out of range ({_BUBBLE_T_LO:.0f}-{_BUBBLE_T_HI:.0f} K)"
                )
            layer = str(self.params.get("reflux_layer", "organic")).lower()
            if layer in ("organic", "light"):
                reflux_organic = True
            elif layer in ("aqueous", "heavy"):
                reflux_organic = False
            else:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: reflux_layer={layer!r} must be "
                    f"'organic'/'light' or 'aqueous'/'heavy'"
                )
        return (n_stages, feed_stages, draws, R, D, P, dP, partial, max_iter,
                method, reboiled, duties, decant, condenser_T, reflux_organic)

    # -- solve -------------------------------------------------------------------
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        n_feeds = max(1, len(self.params.get("feeds") or []))
        feed_streams: list[Stream] = []
        for i in range(n_feeds):
            inlet = inlets.get(f"in{i + 1}")
            if inlet is None or not inlet.molar_flow:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: missing or empty inlet on "
                    f"'in{i + 1}'"
                )
            feed_streams.append(inlet)
        comps = list(feed_streams[0].components)
        F_k: list[float] = []
        z_k: list[dict[str, float]] = []
        for s in feed_streams:
            _, _, flow = s.require_state()
            F_k.append(flow)
            z_k.append(s.normalized_z())
        F = sum(F_k)
        z = {c: sum(Fi * zi.get(c, 0.0) for Fi, zi in zip(F_k, z_k)) / F
             for c in comps}                              # combined feed
        P_in = feed_streams[0].require_state()[1]
        (n_stages, feed_stages, draws, R, D, P, dP, partial, max_iter,
         method, reboiled, duties, decant, condenser_T,
         reflux_organic) = self._read_params(F, P_in)
        reactions, holdup = self._read_reactions(n_stages, reboiled, method)

        # Exact-repeat cache (see __init__): same params + same inlet states
        # -> return copies of the previous result without re-iterating MESH.
        key = (repr(sorted(self.params.items())),) + tuple(
            (s.T, s.P, s.molar_flow, tuple(sorted(s.z.items())), s.H)
            for s in feed_streams)
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

        H_k = [s.H if s.H is not None
               else pp.enthalpy(*s.require_state()[:2], s.normalized_z())
               for s in feed_streams]
        active = [c for c in comps if z.get(c, 0.0) > 0.0]
        if len(active) < 2:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: the feed carries "
                f"{len(active)} component(s) ({active}); distillation needs at "
                f"least two"
            )

        N = n_stages
        P_j = [P + j * dP for j in range(N)]             # stage pressures, top->down
        draw_total = sum(rate for _, _, rate in draws)
        B = F - D - draw_total
        f = {c: F * z[c] for c in active}                # total component feeds

        # Per-stage feed component flows / enthalpy flows, and the per-feed
        # quality q (isenthalpic flash to the feed-stage pressure).
        fz_stage: list[dict[str, float]] = [{} for _ in range(N)]
        fh_stage = [0.0] * N
        feed_info: list[tuple[int, float, float]] = []   # (stage0, F_k, q_k)
        feed_flashes = []                                # PhaseResults, same order
        for stage, Fi, zi, Hi in zip(feed_stages, F_k, z_k, H_k):
            j = stage - 1
            for c in active:
                if zi.get(c, 0.0) > 0.0:
                    fz_stage[j][c] = fz_stage[j].get(c, 0.0) + Fi * zi[c]
            fh_stage[j] += Fi * Hi
            res_k = pp.flash_ph(P_j[j], Hi,
                                {c: zi[c] for c in active if zi.get(c, 0.0) > 0.0})
            feed_info.append((j, Fi, 1.0 - res_k.vapor_fraction))
            feed_flashes.append(res_k)
        q = feed_info[0][2]                              # first feed's quality

        # Per-stage liquid/vapor side draws (mol/s; never on stage 1 or N).
        u = [0.0] * N
        w = [0.0] * N
        for stage, phase, rate in draws:
            if phase == "liquid":
                u[stage - 1] += rate
            else:
                w[stage - 1] += rate

        if not reboiled:
            feed_h = sum(Fi * Hi for Fi, Hi in zip(F_k, H_k))
            out: dict[str, PortStream] = self._solve_nonreboiled(
                pp, comps, active, N, P_j, F, f, fz_stage, fh_stage,
                feed_info, feed_flashes, u, w, draws, duties, R, D, B, q,
                partial, max_iter, feed_h, method)
            self._store_cache(key, out)
            return out

        if method == "naphtali_sandholm":
            feed_h = sum(Fi * Hi for Fi, Hi in zip(F_k, H_k))
            if decant:
                out = self._solve_decant_ns(
                    pp, comps, active, N, P_j, F, f, z, fz_stage, fh_stage,
                    feed_info, u, w, draws, R, D, B, q, condenser_T,
                    reflux_organic, max_iter, feed_h)
            else:
                out = self._solve_reboiled_ns(
                    pp, comps, active, N, P_j, F, f, z, fz_stage, fh_stage,
                    feed_info, feed_flashes, u, w, draws, R, D, B, q, max_iter,
                    feed_h)
            self._store_cache(key, out)
            return out

        # Cumulative (feeds - draws - distillate) above-and-including stage j:
        # the total balance reads L_j = V_{j+1} + A_j (A_{N-1} = B). A reactive
        # stage adds its net mole change (zero for a mole-conserving reaction);
        # `gen_stage[j]` is recomputed each MESH iteration from the liquid profile.
        def _cumulative_A(gen_stage: list[dict[str, float]]) -> list[float]:
            A_local = [0.0] * N
            run_local = -D
            for j in range(N):
                run_local += (sum(fz_stage[j].values())
                              + sum(gen_stage[j].values()) - u[j] - w[j])
                A_local[j] = run_local
            return A_local

        no_gen: list[dict[str, float]] = [{} for _ in range(N)]
        A = _cumulative_A(no_gen)

        # Energy-balance basis. A non-reactive column conditions its stage energy
        # balances far better on a SENSIBLE basis — subtract each stream's
        # composition-weighted formation-enthalpy offset (admissible: a
        # per-component constant cancels against the conserved component balances
        # at convergence; the same reconditioning the non-reboiled envelope form
        # uses). This keeps the recurrence pivot hV_{j+1} - hL_j ~ the true latent
        # heat > 0 even where stage compositions swing between components with very
        # different formation enthalpies (water -242 vs cyclohexane -123 kJ/mol in
        # a heteroazeotrope); on the bare formation-inclusive basis that pivot goes
        # sign-indefinite and the recurrence reports a spurious "degenerate energy
        # balance". A REACTIVE column is LEFT on the formation-inclusive basis —
        # that is exactly how it carries the heat of reaction — so reactive results
        # stay byte-for-byte unchanged.
        use_sensible = not reactions
        _form_fn = getattr(pp, "formation_enthalpies", None)
        if use_sensible and callable(_form_fn):
            _hf_all = _form_fn()
            hf_c = {c: float(_hf_all.get(c, 0.0)) for c in active}
        else:
            hf_c = {c: 0.0 for c in active}

        def _hf_mix(row: dict[str, float]) -> float:
            return sum(hf_c[c] * row.get(c, 0.0) for c in active)

        fh_stage_s = [
            fh_stage[j] - sum(Fk * _hf_mix(zk)
                              for st, Fk, zk in zip(feed_stages, F_k, z_k)
                              if st - 1 == j)
            for j in range(N)]

        x = self._initial_x(pp, P, f, D, N, z)
        L, V = self._initial_traffic(N, feed_info, R, D, F, partial, A, w)

        # Stage-by-stage bubble points of the initial profile seed T, K, h.
        T, K, hL, hV, n_flash = self._stage_bubble_points(pp, P_j, x, active, 0)

        # -- MESH loop (Wang-Henke) ---------------------------------------------
        v_floor = _V_FLOOR_FRAC * F
        omega = 1.0                                      # damping on the V update
        dT_prev = math.inf
        worsened = 0
        converged = False
        it = 0
        dT = dV = math.inf
        gen_stage = no_gen
        B_eff = B
        for it in range(1, max_iter + 1):
            # Reaction generation from the current (x, T) profile, injected as an
            # extra per-stage source into the component balances (the heat of
            # reaction is carried automatically by the formation-inclusive stage
            # enthalpies). A re-derives from the net mole change each iteration.
            # The lagged source is under-relaxed for stability (at convergence the
            # used generation equals the freshly computed rate).
            if reactions:
                gen_new = self._stage_generation(pp, reactions, holdup, x, T,
                                                 P_j, L, active)
                gen_stage = [
                    {c: gen_stage[j].get(c, 0.0)
                     + _RXN_RELAX * (gen_new[j].get(c, 0.0) - gen_stage[j].get(c, 0.0))
                     for c in set(gen_stage[j]) | set(gen_new[j])}
                    for j in range(N)]
                src_stage = [{c: fz_stage[j].get(c, 0.0) + gen_stage[j].get(c, 0.0)
                              for c in set(fz_stage[j]) | set(gen_stage[j])}
                             for j in range(N)]
                A = _cumulative_A(gen_stage)
                B_eff = A[N - 1]
            else:
                src_stage = fz_stage
            x = self._component_balances(N, active, src_stage, D, K, L, V,
                                         partial, u, w)
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

            if use_sensible:
                y_eb = []
                for j in range(N):
                    yj = {c: K[j][c] * x[j][c] for c in active}
                    ytot = sum(yj.values()) or 1.0
                    y_eb.append({c: v / ytot for c, v in yj.items()})
                hL_eb = [hL[j] - _hf_mix(x[j]) for j in range(N)]
                hV_eb = [hV[j] - _hf_mix(y_eb[j]) for j in range(N)]
                fh_eb = fh_stage_s
            else:
                hL_eb, hV_eb, fh_eb = hL, hV, fh_stage
            V_new = self._energy_balances(N, fh_eb, A, u, w, hL_eb, hV_eb, V,
                                          v_floor)
            dV = max(abs(vn - vo) for vn, vo in zip(V_new[1:], V[1:])) / max(V_new[1:])
            V = [V[0]] + [max(vo + omega * (vn - vo), v_floor)
                          for vo, vn in zip(V[1:], V_new[1:])]
            L = [max(V[j + 1] + A[j], v_floor) for j in range(N - 1)] + [B_eff]

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
        # Side-draw flows: the draw leaves at exactly its specified rate with
        # the converged stage composition (saturated liquid or vapor).
        draw_flows: list[dict[str, float]] = []
        for stage, phase, rate in draws:
            comp = x[stage - 1] if phase == "liquid" else y[stage - 1]
            draw_flows.append({c: rate * comp[c] for c in active})
        # Net component feed including reaction generation across all stages:
        # the overall balance is f + Σ_j gen_j = D x_D + B x_B + Σ draws.
        gen_tot = {c: sum(gen_stage[j].get(c, 0.0) for j in range(N))
                   for c in active}
        f_eff = {c: f[c] + gen_tot[c] for c in active}
        b_flows = {c: f_eff[c] - d_flows[c]
                   - sum(df[c] for df in draw_flows) for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: converged distillate/side draws "
                f"would carry more of a component than the feed (plus reaction) "
                f"supplies (bottoms flow {neg:.3g} mol/s) — tighten tolerances or "
                f"check the rate specs"
            )
        if neg < 0.0:
            # Trace negatives from finite MESH tolerance (a light component
            # essentially absent from the bottoms): move the excess out of
            # the distillate's entry for that component and park it on the
            # distillate's share of the bottoms' dominant component — every
            # per-component balance stays machine-exact and D is untouched.
            donor = max(b_flows, key=lambda c: b_flows[c])
            for c, v in b_flows.items():
                if v < 0.0:
                    d_flows[c] += v
                    d_flows[donor] -= v
                    b_flows[c] = 0.0
            b_flows[donor] = (f_eff[donor] - d_flows[donor]
                              - sum(df[donor] for df in draw_flows))
            x_d = {c: v / D for c, v in d_flows.items()}
        x_b = {c: v / B_eff for c, v in b_flows.items()}

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
            T=res_b.T, P=P_j[-1], molar_flow=B_eff, z=z_bot,
            H=res_b.H, phase=res_b.phase, vapor_fraction=res_b.vapor_fraction,
        )

        # Side-draw product streams: saturated at the converged stage state.
        side_streams: list[Stream] = []
        side_h = 0.0
        for i, (stage, phase, rate) in enumerate(draws):
            j = stage - 1
            comp = x[j] if phase == "liquid" else y[j]
            h_draw = hL[j] if phase == "liquid" else hV[j]
            side_h += rate * h_draw
            side_streams.append(Stream(
                id=f"{self.id}.side{i + 1}", components=comps,
                T=T[j], P=P_j[j], molar_flow=rate,
                z={c: comp.get(c, 0.0) for c in comps},
                H=h_draw, phase=phase,
                vapor_fraction=0.0 if phase == "liquid" else 1.0,
            ))

        # -- duties (Heater sign convention: positive heats the process) ----------
        # Condenser from the stage-1 energy balance; reboiler closes the overall
        # balance exactly (enthalpies are absolute / formation-inclusive):
        #   sum_k F_k h_Fk + Q_cond + Q_reb = D h_D + B h_B + sum draws
        if partial:
            q_cond = L[0] * hL[0] + D * res_d.H - V[1] * hV[1]
        else:
            q_cond = (L[0] + D) * res_d.H - V[1] * hV[1]
        feed_h = sum(Fi * Hi for Fi, Hi in zip(F_k, H_k))
        q_reb = D * res_d.H + B_eff * res_b.H + side_h - feed_h - q_cond
        # Diagnostic: the independently computed stage-N (reboiler) balance.
        q_reb_stage = V[N - 1] * hV[N - 1] + B_eff * hL[N - 1] - L[N - 2] * hL[N - 2]
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
            "D": D, "B": B_eff, "R": R, "q": q, "feed_stage": feed_stages[0],
            "partial_condenser": partial,
            "feeds": [{"stage": j + 1, "F": Fi, "q": qi}
                      for j, Fi, qi in feed_info],
            "side_draws": [{"stage": stage, "phase": phase, "rate": rate}
                           for stage, phase, rate in draws],
            "distillate_flows": dict(d_flows), "bottoms_flows": dict(b_flows),
            "side_draw_flows": [dict(df) for df in draw_flows],
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
        out = {
            "distillate": distillate,
            "bottoms": bottoms,
            "condenser_duty": EnergyStream(id=f"{self.id}.condenser_duty", duty=q_cond),
            "reboiler_duty": EnergyStream(id=f"{self.id}.reboiler_duty", duty=q_reb),
        }
        for i, s in enumerate(side_streams):
            out[f"side{i + 1}"] = s
        self._store_cache(key, out)
        return out

    # -- reboiled column on Naphtali-Sandholm (total condenser + reboiler) --------
    def _solve_reboiled_ns(self, pp, comps: list[str], active: list[str],
                           N: int, P_j: list[float], F: float,
                           f: dict[str, float], z: dict[str, float],
                           fz_stage: list[dict[str, float]],
                           fh_stage: list[float],
                           feed_info: list[tuple[int, float, float]],
                           feed_flashes: list,
                           u: list[float], w: list[float],
                           draws: list[tuple[int, str, float]],
                           R: float, D: float, B: float, q: float,
                           max_iter: int, feed_h: float,
                           ) -> dict[str, PortStream]:
        """Reboiled distillation column on the Naphtali-Sandholm simultaneous-
        correction method (Seader 3e §10.4) — the total-condenser + reboiler
        form of :meth:`_solve_ns`. Specs are the distillate rate ``D`` and the
        reflux ratio ``R`` (the same pair as the bubble-point method); the
        reboiler and condenser duties are recovered outputs.

        Two end-stage changes vs the steam-stripped main fractionator (handled
        inside ``_solve_ns`` under ``total_condenser=True``):

        * **stage 0 = total condenser** — the distillate is the liquid draw
          ``U_0 = D``, the stage-0 vapour is pinned to zero, the reflux
          ``L_0 = R*D`` is pinned (its 3rd row), and an appended global
          equation sets ``T_0`` to the condenser-liquid bubble point.
        * **stage N-1 = reboiler** — its energy balance carries one extra global
          unknown, the reboiler duty, boiling up the internal vapour (no open
          stripping steam). The duty is the recovered output.

        Wide-boiling reboiled columns (vacuum towers, reboiled fractionators)
        are where the simultaneous-correction method earns its keep over the
        bubble-point recurrence; this routes them to NS. The mass balance closes
        by difference (machine-exact) and the recovered duties close the overall
        energy balance.
        """
        # Total condenser: the distillate leaves stage 0 as the liquid draw D.
        u = list(u)
        u[0] = D

        # Seed: a Fenske composition ramp + constant-molal-overflow traffic with
        # a total condenser (V_0 = 0) and the reflux R.
        x = self._initial_x(pp, P_j[0], f, D, N, z)
        L, V = self._initial_traffic(N, feed_info, R, D, F, False,
                                     [B] * N, w)
        T, K, _hL, _hV, _nf = self._stage_bubble_points(pp, P_j, x, active, 0)
        y = [{c: max(K[j][c] * x[j].get(c, 0.0), _X_FLOOR) for c in active}
             for j in range(N)]
        for j in range(N):
            tot = sum(y[j].values())
            y[j] = {c: v / tot for c, v in y[j].items()}
        t_bounds = (max(min(T) - 120.0, _BUBBLE_T_LO),
                    min(max(T) + 300.0, _BUBBLE_T_HI))

        (T, x, y, K, L, V, hL, hV, it, dF, n_prop, converged) = \
            self._solve_ns(
                pp, N, P_j, active, fz_stage, fh_stage, u, w, [0.0] * N, D,
                partial=False, T=T, x=x, y=y, L=L, V=V, t_bounds=t_bounds,
                max_iter=max_iter, total_condenser=True, reflux=R,
                tol=_NS_TOL_REBOILED)
        if not converged:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: reboiled naphtali_sandholm did not "
                f"converge in {it} iterations (scaled residual {dF:.3g} vs "
                f"{_NS_TOL}). Check the specs (D={D:.4g} mol/s, reflux_ratio="
                f"{R:.4g}, n_stages={N}) or raise 'max_iter'"
            )

        # -- products (liquid distillate; mass balance closed by difference) -----
        x_d = dict(x[0])
        d_flows = {c: D * x_d.get(c, 0.0) for c in active}
        draw_flows: list[dict[str, float]] = []
        for stage, phase, rate in draws:
            comp = x[stage - 1] if phase == "liquid" else y[stage - 1]
            draw_flows.append({c: rate * comp.get(c, 0.0) for c in active})
        b_flows = {c: f[c] - d_flows[c]
                   - sum(df[c] for df in draw_flows) for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: converged distillate/side draws "
                f"would carry more of a component than the feed supplies "
                f"(bottoms {neg:.3g} mol/s) — check D and the reflux ratio"
            )
        if neg < 0.0:
            donor = max(b_flows, key=lambda c: b_flows[c])
            for c, v in b_flows.items():
                if v < 0.0:
                    d_flows[c] += v
                    d_flows[donor] -= v
                    b_flows[c] = 0.0
            b_flows[donor] = (f[donor] - d_flows[donor]
                              - sum(df[donor] for df in draw_flows))
            x_d = {c: v / D for c, v in d_flows.items()}
        x_b = {c: v / B for c, v in b_flows.items()}
        z_dist = {c: x_d.get(c, 0.0) for c in comps}
        z_bot = {c: x_b.get(c, 0.0) for c in comps}

        side_streams: list[Stream] = []
        side_h = 0.0
        for i, (stage, phase, rate) in enumerate(draws):
            j = stage - 1
            comp = x[j] if phase == "liquid" else y[j]
            h_draw = hL[j] if phase == "liquid" else hV[j]
            side_h += rate * h_draw
            side_streams.append(Stream(
                id=f"{self.id}.side{i + 1}", components=comps,
                T=T[j], P=P_j[j], molar_flow=rate,
                z={c: comp.get(c, 0.0) for c in comps},
                H=h_draw, phase=phase,
                vapor_fraction=0.0 if phase == "liquid" else 1.0))

        # Duties (both recovered). Reboiler from the bottom-stage energy balance
        # (a bubble-point stage with the duty as the NS global free var):
        #   Q_reb = B h_B + V_{N-1} h_V,N-1 - L_{N-2} h_L,N-2.
        # Condenser closes the overall balance (absolute / formation-inclusive
        # enthalpies):  feed_h + Q_cond + Q_reb = D h_D + B h_B + sum draws.
        h_D, h_B = hL[0], hL[N - 1]
        q_reb = B * h_B + V[N - 1] * hV[N - 1] - L[N - 2] * hL[N - 2]
        q_cond = D * h_D + B * h_B + side_h - feed_h - q_reb
        # Diagnostic: the independently computed stage-0 (condenser) balance,
        #   V_1 h_V1 + Q_cond = (L_0 + D) h_L0.
        q_cond_stage = (L[0] + D) * hL[0] - V[1] * hV[1]
        energy_residual = abs(q_cond - q_cond_stage) / max(
            abs(q_cond), abs(q_reb), 1.0)

        distillate = Stream(
            id=f"{self.id}.distillate", components=comps,
            T=T[0], P=P_j[0], molar_flow=D, z=z_dist, H=h_D,
            phase="liquid", vapor_fraction=0.0)
        bottoms = Stream(
            id=f"{self.id}.bottoms", components=comps,
            T=T[N - 1], P=P_j[-1], molar_flow=B, z=z_bot, H=h_B,
            phase="liquid", vapor_fraction=0.0)

        self._warm = {"n_stages": N, "partial": False, "active": list(active),
                      "z": dict(z), "x": [dict(row) for row in x]}
        self.design = {
            "N": float(N - 1), "P": P_j[0], "x_D": dict(x_d), "x_B": dict(x_b),
            "T_top": T[0], "T_top_dew": T[1], "T_bottom": T[N - 1],
            "V_top": V[1], "Q_condenser": q_cond, "Q_reboiler": q_reb,
            "D": D, "B": B, "R": L[0] / D, "q": q,
            "feed_stage": feed_info[0][0] + 1,
            "partial_condenser": False, "method": "naphtali_sandholm",
            "reboiled": True,
            "feeds": [{"stage": j + 1, "F": Fi, "q": qi}
                      for j, Fi, qi in feed_info],
            "side_draws": [{"stage": stage, "phase": phase, "rate": rate}
                           for stage, phase, rate in draws],
            "distillate_flows": dict(d_flows), "bottoms_flows": dict(b_flows),
            "side_draw_flows": [dict(df) for df in draw_flows],
            "n_stages": N, "T_profile": list(T), "P_profile": list(P_j),
            "L_profile": list(L), "V_profile": list(V),
            "x_profile": [{c: row.get(c, 0.0) for c in comps} for row in x],
            "y_profile": [{c: row.get(c, 0.0) for c in comps} for row in y],
            "iterations": it, "energy_residual_rel": energy_residual,
            "prop_evals": n_prop,
        }
        out: dict[str, PortStream] = {
            "distillate": distillate, "bottoms": bottoms,
            "condenser_duty": EnergyStream(
                id=f"{self.id}.condenser_duty", duty=q_cond),
            "reboiler_duty": EnergyStream(
                id=f"{self.id}.reboiler_duty", duty=q_reb),
        }
        for i, s in enumerate(side_streams):
            out[f"side{i + 1}"] = s
        return out

    # -- integrated decanting condenser on Naphtali-Sandholm ---------------------
    def _solve_decant_ns(self, pp, comps: list[str], active: list[str],
                         N: int, P_j: list[float], F: float,
                         f: dict[str, float], z: dict[str, float],
                         fz_stage: list[dict[str, float]],
                         fh_stage: list[float],
                         feed_info: list[tuple[int, float, float]],
                         u: list[float], w: list[float],
                         draws: list[tuple[int, str, float]],
                         R: float, D: float, B: float, q: float,
                         condenser_T: float, reflux_organic: bool,
                         max_iter: int, feed_h: float,
                         ) -> dict[str, PortStream]:
        """Reboiled column with an INTEGRATED DECANTING CONDENSER on Naphtali-
        Sandholm — the fix for the heteroazeotropic entrainer column (Hameed
        §9.5.6 anhydrous-ethanol dehydration with a cyclohexane entrainer).

        An EXTERNAL decanter would feed the entrainer to the column, so the
        overhead must carry it all out → a large distillate the NS cannot
        converge. An INTERNAL decant keeps the entrainer circulating inside the
        column, so the net (aqueous) distillate stays small. Stage 0 is an
        ordinary equilibrium tray; its overhead vapour is condensed at
        ``condenser_T`` and settled into an organic layer (refluxed in full) and
        an aqueous layer (the distillate, spec ``D``). The reflux ratio
        ``R = organic/aqueous`` is an OUTPUT; both duties are recovered. See
        :meth:`_solve_ns` (``decant=True``).
        """
        # -- seed -------------------------------------------------------------
        # A Fenske seed built from the (tiny) entrainer makeup feed has ~no
        # entrainer anywhere, so the overhead does not decant and the internal
        # organic circulation cannot bootstrap (NS stalls). Instead seed an
        # ENTRAINER-RICH rectifying section so the overhead splits from the very
        # first residual and the organic reflux is present in the seed. The
        # entrainer is the smallest feed (the makeup); the bottoms product is the
        # largest feed; the remaining ("water-like") species leave overhead.
        warm = self._warm
        warm_ok = (warm is not None and warm.get("method") == "decant"
                   and warm.get("n_stages") == N and warm.get("active") == active
                   and max(abs(warm["z"].get(c, 0.0) - z.get(c, 0.0))
                           for c in active) < _WARM_Z_TOL)
        if warm_ok:
            assert warm is not None
            x = [dict(row) for row in warm["x"]]
            T = list(warm["T"])
            L, V = list(warm["L"]), list(warm["V"])
        else:
            ent = min(active, key=lambda c: f.get(c, 0.0))
            prod = max(active, key=lambda c: f.get(c, 0.0))
            mids = [c for c in active if c not in (ent, prod)]
            fmid = sum(f.get(c, 0.0) for c in mids) or 1.0
            # Concentrate the entrainer in the RECTIFYING section (above the main
            # feed) and drive it to ~0 in the stripping section, so the overhead
            # is entrainer-rich (decants from the first residual) while the
            # bottoms is the pure product. Grading over the whole column instead
            # spills entrainer into the stripping section and is a poor seed for
            # a long column (the cold solve then misses the basin).
            fmain = feed_info[max(range(len(feed_info)),
                                  key=lambda k: feed_info[k][1])][0]
            x = []
            for j in range(N):
                if j <= fmain:                 # rectifying: entrainer-rich
                    g = j / max(fmain, 1)
                    x_ent = 0.45 * (1.0 - g) + 0.05
                    x_mid = 0.20 * (1.0 - g) + 0.03
                else:                          # stripping: entrainer -> 0
                    g = (j - fmain) / max(N - 1 - fmain, 1)
                    x_ent = 0.05 * (1.0 - g) + 1e-4
                    x_mid = 0.06 * (1.0 - g) + 1e-3
                row = {ent: x_ent}
                for c in mids:
                    row[c] = x_mid * (f.get(c, 0.0) / fmid)
                row[prod] = max(1.0 - x_ent - x_mid, 1e-3)
                tot = sum(row.values())
                x.append({c: max(row.get(c, 0.0), _X_FLOOR) / tot
                          for c in active})
            L, V = self._initial_traffic(N, feed_info, R, D, F, True,
                                         [B] * N, w)
            T = None  # bubble-pointed below
        V[0] = D * (1.0 + R)
        T0, K, _hL, _hV, _nf = self._stage_bubble_points(pp, P_j, x, active, 0)
        if T is None:
            T = T0
        y = [{c: max(K[j][c] * x[j].get(c, 0.0), _X_FLOOR) for c in active}
             for j in range(N)]
        for j in range(N):
            tot = sum(y[j].values())
            y[j] = {c: v / tot for c, v in y[j].items()}
        t_bounds = (max(min(min(T), condenser_T) - 120.0, _BUBBLE_T_LO),
                    min(max(T) + 300.0, _BUBBLE_T_HI))

        (T, x, y, K, L, V, hL, hV, it, dF, n_prop, converged) = \
            self._solve_ns(
                pp, N, P_j, active, fz_stage, fh_stage, u, w, [0.0] * N, D,
                partial=False, T=T, x=x, y=y, L=L, V=V, t_bounds=t_bounds,
                max_iter=max_iter, total_condenser=False, reflux=R,
                tol=_NS_TOL_REBOILED, decant=True, condenser_T=condenser_T,
                reflux_organic=reflux_organic)
        if not converged:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: decanting-condenser "
                f"naphtali_sandholm did not converge in {it} iterations (scaled "
                f"residual {dF:.3g} vs {_NS_TOL_REBOILED}). Check the specs "
                f"(D={D:.4g} mol/s, reflux_ratio={R:.4g}, n_stages={N}, "
                f"condenser_T={condenser_T:.4g} K) or raise 'max_iter'"
            )

        # -- decant the converged overhead -> organic reflux + aqueous distillate
        V0, y0 = V[0], y[0]
        r3 = pp.flash_pt_3p(condenser_T, P_j[0], dict(y0))
        bl, bh = r3.beta_light, r3.beta_heavy
        bliq = bl + bh
        if (bl <= _DECANT_BETA_MIN or bh <= _DECANT_BETA_MIN or bliq <= 0.0
                or r3.x_light is None or r3.x_heavy is None):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: the converged overhead does not "
                f"decant into two liquids at condenser_T={condenser_T:.4g} K "
                f"(beta_light={bl:.3g}, beta_heavy={bh:.3g}) — the column has no "
                f"heteroazeotrope to split; lower condenser_T or check the feed"
            )
        if reflux_organic:
            xo, hL_o, frac_o, hL_a = r3.x_light, r3.H_light, bl / bliq, r3.H_heavy
        else:
            xo, hL_o, frac_o, hL_a = r3.x_heavy, r3.H_heavy, bh / bliq, r3.H_light
        org = V0 * frac_o
        o_flows = {c: min(org * xo.get(c, 0.0), V0 * y0.get(c, 0.0))
                   for c in active}
        org = sum(o_flows.values())
        a_flows = {c: V0 * y0.get(c, 0.0) - o_flows[c] for c in active}
        D_aq = sum(a_flows.values())

        # Bottoms by difference (the organic reflux is internal — only the
        # aqueous layer and the bottoms leave the column).
        b_flows = {c: f[c] - a_flows[c] for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: the converged aqueous distillate "
                f"carries more of a component than the feed supplies (bottoms "
                f"{neg:.3g} mol/s) — check D, reflux_ratio and condenser_T"
            )
        for c in active:
            b_flows[c] = max(b_flows[c], 0.0)
        B_eff = sum(b_flows.values())
        x_d = {c: (a_flows[c] / D_aq if D_aq > 0.0 else 0.0) for c in active}
        x_b = {c: (b_flows[c] / B_eff if B_eff > 0.0 else 0.0) for c in active}
        z_dist = {c: x_d.get(c, 0.0) for c in comps}
        z_bot = {c: x_b.get(c, 0.0) for c in comps}

        side_streams: list[Stream] = []
        side_h = 0.0
        for i, (stage, phase, rate) in enumerate(draws):
            j = stage - 1
            comp = x[j] if phase == "liquid" else y[j]
            h_draw = hL[j] if phase == "liquid" else hV[j]
            side_h += rate * h_draw
            side_streams.append(Stream(
                id=f"{self.id}.side{i + 1}", components=comps,
                T=T[j], P=P_j[j], molar_flow=rate,
                z={c: comp.get(c, 0.0) for c in comps},
                H=h_draw, phase=phase,
                vapor_fraction=0.0 if phase == "liquid" else 1.0))

        # Duties. Reboiler from the bottom-stage balance; condenser from the
        # overall balance (formation-inclusive enthalpies). The aqueous
        # distillate leaves at condenser_T with the heavy-layer enthalpy.
        h_dist, h_bot = float(hL_a), hL[N - 1]
        q_reb = B_eff * h_bot + V[N - 1] * hV[N - 1] - L[N - 2] * hL[N - 2]
        q_cond = D_aq * h_dist + B_eff * h_bot + side_h - feed_h - q_reb
        # Independent condenser balance (heat removed condensing V_0 into the two
        # cold layers, as a duty): organic returns internally, aqueous leaves.
        q_cond_stage = org * float(hL_o) + D_aq * h_dist - V0 * hV[0]
        energy_residual = abs(q_cond - q_cond_stage) / max(
            abs(q_cond), abs(q_reb), 1.0)

        distillate = Stream(
            id=f"{self.id}.distillate", components=comps,
            T=condenser_T, P=P_j[0], molar_flow=D_aq, z=z_dist, H=h_dist,
            phase="liquid", vapor_fraction=0.0)
        bottoms = Stream(
            id=f"{self.id}.bottoms", components=comps,
            T=T[N - 1], P=P_j[-1], molar_flow=B_eff, z=z_bot, H=h_bot,
            phase="liquid", vapor_fraction=0.0)

        self.design = {
            "N": float(N - 1), "P": P_j[0], "x_D": dict(x_d), "x_B": dict(x_b),
            # The integrated decant condenser cools the overhead from its dew
            # point (~the top-tray temperature T[0]) down to ``condenser_T``,
            # where the two liquids settle; expose those as the condenser hot/cold
            # ends so the economics sizer + pinch analysis can cost it (T_top =
            # condensate/distillate temperature, T_top_dew = overhead dew point).
            "T_top": condenser_T, "T_top_dew": T[0],
            "T_bottom": T[N - 1], "condenser_T": condenser_T,
            "V_top": V0, "Q_condenser": q_cond, "Q_reboiler": q_reb,
            "D": D_aq, "B": B_eff,
            "R": (org / D_aq if D_aq > 0.0 else 0.0), "q": q,
            "reflux_ratio_organic": (org / D_aq if D_aq > 0.0 else 0.0),
            "organic_reflux": org, "organic_flows": dict(o_flows),
            "x_organic": {c: (o_flows[c] / org if org > 0.0 else 0.0)
                          for c in active},
            "feed_stage": feed_info[0][0] + 1,
            "partial_condenser": False, "method": "naphtali_sandholm",
            "reboiled": True, "decant_condenser": True,
            "feeds": [{"stage": j + 1, "F": Fi, "q": qi}
                      for j, Fi, qi in feed_info],
            "side_draws": [{"stage": stage, "phase": phase, "rate": rate}
                           for stage, phase, rate in draws],
            "distillate_flows": dict(a_flows), "bottoms_flows": dict(b_flows),
            "side_draw_flows": [],
            "n_stages": N, "T_profile": list(T), "P_profile": list(P_j),
            "L_profile": list(L), "V_profile": list(V),
            "x_profile": [{c: row.get(c, 0.0) for c in comps} for row in x],
            "y_profile": [{c: row.get(c, 0.0) for c in comps} for row in y],
            "iterations": it, "energy_residual_rel": energy_residual,
            "prop_evals": n_prop,
        }
        # Warm start for recycle re-solves (a closed entrainer loop re-calls this
        # column with near-identical feeds many times): reuse the converged
        # entrainer-rich profile + traffic instead of rebuilding the seed.
        self._warm = {"n_stages": N, "method": "decant", "active": list(active),
                      "z": dict(z), "x": [dict(row) for row in x],
                      "T": list(T), "L": list(L), "V": list(V)}
        out: dict[str, PortStream] = {
            "distillate": distillate, "bottoms": bottoms,
            "condenser_duty": EnergyStream(
                id=f"{self.id}.condenser_duty", duty=q_cond),
            "reboiler_duty": EnergyStream(
                id=f"{self.id}.reboiler_duty", duty=q_reb),
        }
        for i, s in enumerate(side_streams):
            out[f"side{i + 1}"] = s
        return out

    def _store_cache(self, key: tuple, out: dict[str, PortStream]) -> None:
        """Record the exact-repeat cache entry for this solve."""
        assert self.design is not None
        self._cache_key = key
        self._cache_out = {
            name: (s.with_() if isinstance(s, Stream)
                   else EnergyStream(id=s.id, duty=s.duty))
            for name, s in out.items()
        }
        self._cache_design = {k: (list(v) if isinstance(v, list) else
                                  dict(v) if isinstance(v, dict) else v)
                              for k, v in self.design.items()}

    # -- non-reboiled mode (steam-stripped main fractionator) --------------------
    def _solve_nonreboiled(self, pp, comps: list[str], active: list[str],
                           N: int, P_j: list[float], F: float,
                           f: dict[str, float],
                           fz_stage: list[dict[str, float]],
                           fh_stage: list[float],
                           feed_info: list[tuple[int, float, float]],
                           feed_flashes: list,
                           u: list[float], w: list[float],
                           draws: list[tuple[int, str, float]],
                           duties: list[tuple[int, float]],
                           R0: float, D: float, B: float, q: float,
                           partial: bool, max_iter: int,
                           feed_h: float, method: str) -> dict[str, PortStream]:
        """MESH for the non-reboiled column (condenser, open stripping steam,
        side draws, pumparound stage duties): see the module docstring.

        Two temperature/traffic updates share the same component-balance
        (Thomas) core, selected by ``method``:

        * ``"bubble_point"`` (the robust non-reboiled method) — stage
          temperatures from K-value bubble points of the stage liquids
          (:func:`_bubble_T_kvalues`); the vapor traffic from *envelope*
          energy balances (each V_j from its own around-stage-j-down-to-the-
          bottom balance, so there is no stage-to-stage error accumulation),
          evaluated on a sensible basis by subtracting each stream's
          composition-weighted formation-enthalpy offset (``hf_mix``) — which
          keeps the envelope pivot hLs_{j-1} - hVs_j ~ -(latent heat) < 0 even
          across the steam stage, where the formation-inclusive offsets of
          water (-242 kJ/mol) and the pseudo-components (0) would otherwise
          make it sign-indefinite. The reflux follows from the condenser total
          balance, so the distillate-rate spec is the column's single
          specification. The temperature move and traffic update are both
          under-relaxed, with the damping reduced whenever the *combined*
          temperature+traffic change stalls — the open-steam bottom stage can
          otherwise sit in a period-2 limit cycle (the energy balance flipping
          the bottom traffic between two states each iteration).
        * ``"sum_rates"`` — Burningham-Otto: liquid traffic from the
          component-flow sums, temperatures from a simultaneous Newton solve
          of the stage energy balances (:func:`_newton_T_draws`), condenser
          pinned to saturation. NOT recommended: on towers where stage
          temperatures are equilibrium-dominated (every distillation-like
          section) this update limit-cycles — the documented Seader 3e ch.
          10.4 division of labor; ``bubble_point`` is the robust choice here
          and is what the crude-tower example/test use.

        In both, the condenser duty absorbs the stage-1 energy balance and
        closes the overall balance exactly; mass balances close by
        difference, machine-exact."""
        u0 = 0.0 if partial else D      # condenser liquid product draw
        v0 = D if partial else 0.0      # condenser vapor product
        floor = _V_FLOOR_FRAC * F

        # Stage duties (pumparound heat) enter the stage energy balances
        # exactly like a feed enthalpy flow. They are ramped in over the
        # first _SR_DUTY_RAMP iterations: dumping tens of MW on a stage whose
        # traffic is still the initializer's guess sends the profile wild.
        Q_stage = [0.0] * N
        for stage, duty in duties:
            Q_stage[stage - 1] += duty

        # Feeds / draws entering at-or-below each stage, for the
        # bottom-section total balances V_j = L_{j-1} + F>=j - (u+w)>=j - B.
        F_ge = [0.0] * (N + 1)
        UW_ge = [0.0] * (N + 1)
        for j in range(N - 1, -1, -1):
            F_ge[j] = F_ge[j + 1] + sum(fz_stage[j].values())
            UW_ge[j] = UW_ge[j + 1] + u[j] + w[j]

        # Per-component constant enthalpy offset (J/mol) used to recondition the
        # stage energy balances onto a sensible-only basis. Subtracting a
        # per-component constant from every molar enthalpy is exactly admissible
        # (it cancels against the conserved envelope component balances at
        # convergence) and it keeps the envelope pivot hLs_{j-1} - hVs_j ~
        # -(latent heat) < 0: on the raw formation-inclusive basis the offsets
        # span water (-242 kJ/mol) vs the crude pseudo-components (0), and that
        # disparity makes the bare pivot hL_{j-1} - hV_j sign-indefinite, which
        # wrecks the envelope form on the book's crude tower. The package's
        # formation_enthalpies() are precisely those baked-in offsets; if a
        # backend does not expose them (the method is not on the
        # PropertyPackage protocol), fall back to each pure component's vapor
        # enthalpy at a fixed reference state -- an equally admissible
        # per-component constant.
        _form_fn = getattr(pp, "formation_enthalpies", None)
        hf_c: dict[str, float]
        if callable(_form_fn):
            _hf_all = _form_fn()
            hf_c = {c: float(_hf_all.get(c, 0.0)) for c in active}
        else:
            hf_c = {c: pp.enthalpy_vapor(400.0, 1.0e4, {c: 1.0}) for c in active}

        def hf_mix(row: dict[str, float]) -> float:
            """Composition-weighted constant enthalpy offset (J/mol for a
            mole-fraction row; J/s for a molar-flow row) -- the per-component
            constant subtracted to put the energy balances on a sensible
            basis."""
            return sum(v * hf_c[c] for c, v in row.items())

        # Per-stage feed formation-enthalpy FLOW (J/s), the constant paralleling
        # the absolute feed-enthalpy flow fh_eff; their difference is the feed's
        # sensible enthalpy flow used by the envelope balances.
        fh_form = [hf_mix(fz) for fz in fz_stage]
        # ... and the equivalent top-down cumulative A_j (L_j = V_{j+1} + A_j;
        # A_{N-1} = B), shared with the reboiled bubble-point path.
        A = [0.0] * N
        run = -D
        for j in range(N):
            run += sum(fz_stage[j].values()) - u[j] - w[j]
            A[j] = run

        T, x, y, K, L, V = self._initial_sr_profiles(
            pp, N, P_j, active, f, D, R0, feed_info, feed_flashes,
            F_ge, UW_ge, B, v0, floor)
        t_bounds = (max(min(T) - 120.0, _BUBBLE_T_LO),
                    min(max(T) + 300.0, _BUBBLE_T_HI))

        if method == "inside_out":
            (T, x, y, K, L, V, hL, hV, it, dT, dF, n_prop,
             converged) = self._inside_out_loops(
                pp, N, P_j, active, fz_stage, fh_stage, u, w, A, B, v0, u0,
                floor, Q_stage, duties, hf_c, fh_form, F_ge, UW_ge,
                T, x, y, L, V, t_bounds, max_iter)
            if not converged:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: inside_out iteration did not "
                    f"converge in {it} outer passes (max |dT|={dT:.3g} K vs "
                    f"{_IO_OUTER_TOL_T} K, max traffic change={dF:.3g} vs "
                    f"{_IO_OUTER_TOL_F}). Check the specs (D={D:.4g} mol/s, "
                    f"n_stages={N}, draws={sum(u) + sum(w):.4g} mol/s) or raise "
                    f"'max_iter'"
                )
            return self._finish_nonreboiled(
                pp, comps, active, N, P_j, F, f, feed_info, u, w, draws,
                duties, Q_stage, D, B, q, partial, u0, v0, feed_h,
                method, T, x, y, K, L, V, hL, hV, it, dT, dF, 1.0, n_prop)

        if method == "naphtali_sandholm":
            if not partial:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: method='naphtali_sandholm' "
                    f"currently supports a partial condenser only (the crude "
                    f"main-fractionator form); set partial_condenser=True or "
                    f"use 'bubble_point'"
                )
            # One NS attempt for a given draw schedule (u_, w_) and optional
            # warm profile: inside-out homotopy warm start -> NS Newton ->
            # FD-Jacobian endgame fallback. Returns the _solve_ns 12-tuple.
            #
            # Homotopy warm start: the inside-out method robustly drives the
            # profile to the right PHYSICS (a hot resid bottom, steam rising)
            # even where it stalls short of a tight residual, whereas the raw
            # Fenske/CMO seed leaves a cold steam-flooded bottom that traps the
            # NS Newton in a region where the cubic-EOS heavy-vapour enthalpy is
            # degenerate. Seeding NS from the inside-out profile lets the full
            # Newton finish to machine precision.
            def _ns_attempt(u_, w_, warm, maxit=max_iter, use_fd=True):
                # Draw-dependent cumulatives for THIS schedule (the bottoms
                # B_ = F - D - draws absorbs the scaled-away side product).
                UW_ge_ = [0.0] * (N + 1)
                for jj in range(N - 1, -1, -1):
                    UW_ge_[jj] = UW_ge_[jj + 1] + u_[jj] + w_[jj]
                A_ = [0.0] * N
                run_ = -D
                for jj in range(N):
                    run_ += sum(fz_stage[jj].values()) - u_[jj] - w_[jj]
                    A_[jj] = run_
                B_ = F - D - (sum(u_) + sum(w_))
                if warm is None:
                    T0, x0, y0, _K0, L0, V0 = self._initial_sr_profiles(
                        pp, N, P_j, active, f, D, R0, feed_info, feed_flashes,
                        F_ge, UW_ge_, B_, v0, floor)
                else:
                    T0, x0, y0, L0, V0 = warm
                (Tw, xw, yw, _Kw, Lw, Vw, *_rest) = self._inside_out_loops(
                    pp, N, P_j, active, fz_stage, fh_stage, u_, w_, A_, B_, v0,
                    u0, floor, Q_stage, duties, hf_c, fh_form, F_ge, UW_ge_,
                    T0, x0, y0, L0, V0, t_bounds, min(maxit, _NS_WARM_ITERS))
                res = self._solve_ns(
                    pp, N, P_j, active, fz_stage, fh_stage, u_, w_, Q_stage, D,
                    partial, Tw, xw, yw, Lw, Vw, t_bounds, maxit)
                if not res[-1] and use_fd and hasattr(pp, "stage_derivs"):
                    # Endgame-robustness fallback. The analytic Jacobian
                    # (pp.stage_derivs) and the residual (pp.k_values /
                    # enthalpies) are SEPARATE thermo code paths; the cubic-EOS
                    # fugacity routines can differ enough across platform thermo
                    # builds (deps are unpinned, so e.g. Windows/py3.12 resolves
                    # a different thermo than Linux) that the analytic Jacobian
                    # becomes slightly inconsistent with its residual and the
                    # final quadratic plunge stalls short (observed: one CI job
                    # flat at scaled residual ~4e-3 while three reach 1e-9). The
                    # FD Jacobian is differenced from the SAME residual, so it is
                    # consistent to machine precision on every platform.
                    res2 = self._solve_ns(
                        pp, N, P_j, active, fz_stage, fh_stage, u_, w_, Q_stage,
                        D, partial, Tw, xw, yw, Lw, Vw, t_bounds, maxit,
                        force_fd=True)
                    res = res2[:10] + (res[10] + res2[10], res2[11])
                return res

            # Fast path: a direct analytic Newton from the inside-out warm start.
            # The resid main fractionator converges here in ~11 iterations; the
            # continuation and the FD fallback below never run for it. The
            # expensive FD-Jacobian endgame fallback is deferred PAST the cheap
            # draw continuation, so a draw-heavy basin stall (which the FD
            # fallback cannot fix anyway -- it is an endgame digit-recovery tool,
            # not a basin tool) is rescued by continuation without first grinding
            # out ~N(2C+1) finite-difference Jacobian builds.
            (T, x, y, K, L, V, hL, hV, it, dF, n_prop,
             converged) = _ns_attempt(u, w, None, use_fd=False)
            if not converged and (sum(u) + sum(w)) > 0.0:
                # Broader-basin fallback: warm-start DRAW CONTINUATION. A
                # draw-heavy tower (large side draws relative to its stage
                # count) can leave the damped Newton stuck at a non-root
                # stationary point of ||R|| (the LM escape oscillates without
                # closing). Homotopy on the draw rates fixes it: start from a
                # light-draw column (NS solves it easily -- it is close to a
                # plain distillate/bottoms split), then ramp the draws up to
                # their target in steps, warm-starting each NS solve from the
                # previous converged profile so every step stays in-basin. This
                # is the same warm-start-continuation trick that made the
                # reboiled+refluxed amine regenerator converge (_solve_reboiled_ns).
                # Adaptive step: grow the draw fraction toward 1.0, halving the
                # step whenever a target fails to converge and re-trying from the
                # last good profile (the s=0.85->1.0 jump is the hard one on a
                # heavily-drawn tower; a finer approach walks it in). Each solve
                # warm-starts from the previous, so the successful steps finish
                # in a handful of Newton iterations.
                # A warm-started step that is GOING to converge does so in a
                # handful of Newton iterations (~7-11); a capped budget + no FD
                # fallback keeps each intermediate step cheap (a stalled step
                # halves the draw step rather than grinding out the full budget).
                # Only the final, full-draw target gets the full iteration budget
                # and the FD endgame fallback.
                warm = None
                np_cont = 0
                s_done = 0.0
                step = _NS_CONT_STEP0
                r = None
                while True:
                    s = min(1.0, s_done + step)
                    final = s >= 1.0
                    us = [s * ui for ui in u]
                    ws = [s * wi for wi in w]
                    r = _ns_attempt(
                        us, ws, warm,
                        maxit=max_iter if final else min(max_iter, _NS_CONT_IT),
                        use_fd=final)
                    np_cont += r[10]
                    if getattr(self, "_ns_debug", False):
                        print(f"  [draw-cont] s={s:.3f} step={step:.3f} "
                              f"converged={r[-1]} resid={r[9]:.3e} it={r[8]}")
                    if r[-1]:
                        s_done = s
                        warm = (r[0], r[1], r[2], r[4], r[5])  # T, x, y, L, V
                        if final:
                            converged = True
                            break
                        step = min(step * _NS_CONT_GROW, _NS_CONT_STEP0)
                    else:
                        step *= 0.5
                        if step < _NS_CONT_STEP_MIN:
                            break
                if converged and r is not None:
                    (T, x, y, K, L, V, hL, hV, it, dF, _np_last,
                     converged) = r
                    n_prop += np_cont
            if not converged and hasattr(pp, "stage_derivs"):
                # Endgame FD-Jacobian fallback (last resort; the ONLY fallback
                # for a no-draw column, where continuation cannot run). The
                # analytic Jacobian (pp.stage_derivs) and the residual
                # (pp.k_values / enthalpies) are SEPARATE thermo code paths; the
                # cubic-EOS fugacity routines can differ enough across platform
                # thermo builds that near the root the analytic Jacobian becomes
                # slightly inconsistent with its residual and the final quadratic
                # plunge stalls short (observed: one CI job flat at scaled
                # residual ~4e-3 while three reach 1e-9). The FD Jacobian is
                # differenced from the SAME residual, so it is consistent to
                # machine precision on every platform.
                (T2, x2, y2, K2, L2, V2, hL2, hV2, it2, dF2, np2,
                 conv2) = _ns_attempt(u, w, None, use_fd=True)
                n_prop += np2
                if conv2:
                    (T, x, y, K, L, V, hL, hV, it, dF, converged) = (
                        T2, x2, y2, K2, L2, V2, hL2, hV2, it2, dF2, conv2)
            if not converged:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: naphtali_sandholm Newton did "
                    f"not converge in {it} iterations (scaled residual "
                    f"{dF:.3g} vs {_NS_TOL}). Check the specs (D={D:.4g} mol/s, "
                    f"n_stages={N}, draws={sum(u) + sum(w):.4g} mol/s) or raise "
                    f"'max_iter'"
                )
            return self._finish_nonreboiled(
                pp, comps, active, N, P_j, F, f, feed_info, u, w, draws,
                duties, Q_stage, D, B, q, partial, u0, v0, feed_h,
                method, T, x, y, K, L, V, hL, hV, it, 0.0, dF, 1.0, n_prop)

        omega = 1.0
        move_prev = math.inf
        worsened = 0
        converged = False
        n_prop = 0
        it = 0
        dT = dF = resid_rel = math.inf
        hL = [pp.enthalpy_liquid(T[j], P_j[j], x[j]) for j in range(N)]
        hV = [pp.enthalpy_vapor(T[j], P_j[j], y[j]) for j in range(N)]
        for it in range(1, max_iter + 1):
            ramp = min(1.0, it / _SR_DUTY_RAMP) if duties else 1.0
            fh_eff = [fh + ramp * qd for fh, qd in zip(fh_stage, Q_stage)]
            if method == "sum_rates":
                K = [pp.k_values(T[j], P_j[j], x[j], y[j]) for j in range(N)]
                n_prop += N

            # -- component balances (Thomas), shared by both updates -----------
            cols: dict[str, list[float]] = {}
            for c in active:
                a = [0.0] + [L[j - 1] for j in range(1, N)]
                b = [-(L[0] + u0 + v0 * K[0][c])] + \
                    [-(L[j] + u[j] + (V[j] + w[j]) * K[j][c])
                     for j in range(1, N)]
                cc = [V[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
                d = [-fz_stage[j].get(c, 0.0) for j in range(N)]
                cols[c] = _thomas(a, b, cc, d)
            sum_raw = [sum(max(cols[c][j], 0.0) for c in active)
                       for j in range(N)]
            if min(sum_raw) <= 0.0:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: the component balances "
                    f"left a stage with no liquid (iteration {it}); the "
                    f"column is infeasible as specified"
                )
            L_prev, V_prev = list(L), list(V)

            if method == "sum_rates":
                # -- sum-rates traffic (bottoms pinned at B by the spec);
                # ratio-clamped against transients --------------------------
                L_new = [L[j] * min(max(sum_raw[j], 1.0 / _SR_L_RATIO),
                                    _SR_L_RATIO) for j in range(N)]
                L = [max(lo + omega * (ln - lo), floor)
                     for lo, ln in zip(L, L_new)]
                L[N - 1] = B
                V = [v0] + [max(L[j - 1] + F_ge[j] - UW_ge[j] - B, floor)
                            for j in range(1, N)]
                for j in range(N):
                    row = {c: max(cols[c][j], _X_FLOOR) for c in active}
                    tot = sum(row.values())
                    x[j] = {c: v / tot for c, v in row.items()}
                    yr = {c: max(K[j][c] * x[j][c], _X_FLOOR) for c in active}
                    tot = sum(yr.values())
                    y[j] = {c: v / tot for c, v in yr.items()}

                # condenser pinned to saturation; the rest from the
                # simultaneous Newton energy solve
                T0_new, K0, y0 = _bubble_T_kvalues(pp, P_j[0], x[0], T[0])
                y[0] = y0
                n_prop += 6
                T_new, hL, hV, resid_rel, n_h = _newton_T_draws(
                    pp, N, P_j, x, y, L, V, fh_eff, u, w,
                    [T0_new] + T[1:], t_bounds)
                n_prop += n_h
            else:
                # -- bubble-point update: normalize the stage liquids, take
                # each stage's K-value bubble point, then drive the vapor
                # traffic from the energy balances run bottom-up from the
                # steam stage (V_{N+1} = 0; no reboiler) -----------------------
                for j in range(N):
                    row = {c: max(cols[c][j], _X_FLOOR) for c in active}
                    tot = sum(row.values())
                    x[j] = {c: v / tot for c, v in row.items()}
                T_new = list(T)
                for j in range(N):
                    tj, K[j], y[j] = _bubble_T_kvalues(pp, P_j[j], x[j], T[j])
                    n_prop += 6
                    T_new[j] = T[j] + max(min(tj - T[j], _SR_T_OUTER_MAX),
                                          -_SR_T_OUTER_MAX)
                hL = [pp.enthalpy_liquid(T_new[j], P_j[j], x[j])
                      for j in range(N)]
                hV = [pp.enthalpy_vapor(T_new[j], P_j[j], y[j])
                      for j in range(N)]
                n_prop += 2 * N

                # Envelope energy balances (around stage j down to the
                # bottom; Seader 3e sec. 10.3 section form) on the SENSIBLE
                # basis: each V_j comes from its own envelope, so there is no
                # stage-to-stage error accumulation, and the formation-offset
                # shift keeps the pivot hLs_{j-1} - hVs_j ~ -(latent heat) < 0
                # even across the steam stage. ``hf_mix`` is the composition-
                # weighted formation enthalpy (a per-component constant, so it
                # cancels against the conserved envelope component balances at
                # convergence); ``fh_form`` is the matching per-stage feed
                # formation-enthalpy flow, so fh_eff - fh_form is the feed's
                # sensible enthalpy flow (the ramped pumparound duty in fh_eff
                # is pure heat and rightly survives the subtraction).
                hLs = [hL[j] - hf_mix(x[j]) for j in range(N)]
                hVs = [hV[j] - hf_mix(y[j]) for j in range(N)]
                fhs = [fe - ff for fe, ff in zip(fh_eff, fh_form)]
                # suffix sums: draws' and feeds' sensible enthalpy flows at
                # or below stage j
                S_uw = [0.0] * (N + 1)
                S_f = [0.0] * (N + 1)
                for j in range(N - 1, 0, -1):
                    S_uw[j] = S_uw[j + 1] + u[j] * hLs[j] + w[j] * hVs[j]
                    S_f[j] = S_f[j + 1] + fhs[j]
                # The envelope pivot hLs_{j-1} - hVs_j is ~ -(latent heat) < 0
                # at the solution, but a transient iterate (especially across
                # the open-steam bottom stage, where the rising vapor is steam-
                # rich and its sensible enthalpy can momentarily approach the
                # falling liquid's) can drive it to zero or slightly positive.
                # Clamp it to a small negative floor -- keeping the physical
                # sign while bounding V_new -- rather than aborting on a
                # recoverable transient; a genuinely degenerate *converged*
                # state is caught downstream (dried-up section / unmet energy
                # residual), never silently accepted. The clamp scale is a
                # small fraction of a typical molar latent heat.
                _pivot_floor = -1.0e3       # J/mol
                V_new = list(V)
                for j in range(1, N):
                    denom = min(hLs[j - 1] - hVs[j], _pivot_floor)
                    num = (S_uw[j] + B * hLs[N - 1] - S_f[j]
                           - A[j - 1] * hLs[j - 1])
                    V_new[j] = max(num / denom, floor)
                V = [v0] + [max(vo + omega * (vn - vo), floor)
                            for vo, vn in zip(V[1:], V_new[1:])]
                L = [max(V[j + 1] + A[j], floor)
                     for j in range(N - 1)] + [B]
                resid_rel = 0.0     # the recurrence enforces the balances
                n_prop += 0

            dF = max(
                max(abs(a_ - b_) / max(a_, b_) for a_, b_ in zip(L, L_prev)),
                max(abs(a_ - b_) / max(a_, b_, floor * 10)
                    for a_, b_ in zip(V[1:], V_prev[1:])),
            )
            dT = max(abs(tn - to) for tn, to in zip(T_new, T))
            # Per-outer-iteration temperature move: damped by omega and then
            # clamped. Damping the *temperature* (not only the traffic) is what
            # lets omega break a period-2 limit cycle on an open-steam bottom
            # stage, where the energy balance flips the bottom traffic between
            # two states each iteration and the raw temperature move pins at the
            # clamp (so a dT-only oscillation test never fires).
            T = [to + max(min(omega * (tn - to), _SR_T_OUTER_MAX),
                          -_SR_T_OUTER_MAX)
                 for tn, to in zip(T_new, T)]

            # Reduce damping when EITHER the temperature or the traffic change
            # stalls/worsens: a flat-but-large dT (clamped) with a persistently
            # large dF is exactly the limit-cycle signature, so progress is
            # measured on the combined move and a lack of improvement shrinks
            # omega geometrically.
            move = max(dT / max(_SR_T_OUTER_MAX, 1.0), min(dF, 1.0))
            if move >= move_prev - 1e-3:
                worsened += 1
                if worsened >= 2:
                    omega = max(0.1, 0.5 * omega)
                    worsened = 0
            else:
                worsened = 0
            move_prev = move

            if getattr(self, "_sr_debug", False):
                dl = [abs(a_ - b_) / max(a_, b_) for a_, b_ in zip(L, L_prev)]
                jmax = dl.index(max(dl))
                print(f"  it={it:3d} dT={dT:.3e} dF={dF:.3e} "
                      f"resid={resid_rel:.1e} omega={omega} "
                      f"jL={jmax} L[jL]={L[jmax]:.2f} L0={L[0]:.2f} "
                      f"T0={T[0]:.2f} Tbot={T[-1]:.2f}")

            if (ramp >= 1.0 and dT <= _TOL_T and dF <= _SR_TOL_FLOW
                    and resid_rel <= _SR_ENERGY_TOL):
                converged = True
                break

        if not converged:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: {method} iteration did not "
                f"converge in {max_iter} iterations (max |dT|={dT:.3g} K vs "
                f"{_TOL_T} K, max traffic change={dF:.3g} vs {_SR_TOL_FLOW}, "
                f"energy residual={resid_rel:.3g} vs {_SR_ENERGY_TOL}, "
                f"damping={omega}). Check the specs (D={D:.4g} mol/s, "
                f"n_stages={N}, draws={sum(u) + sum(w):.4g} mol/s) or raise "
                f"'max_iter'"
            )
        return self._finish_nonreboiled(
            pp, comps, active, N, P_j, F, f, feed_info, u, w, draws, duties,
            Q_stage, D, B, q, partial, u0, v0, feed_h, method,
            T, x, y, K, L, V, hL, hV, it, dT, dF, omega, n_prop)

    # -- shared non-reboiled product / duty / design assembly --------------------
    def _finish_nonreboiled(
        self, pp, comps: list[str], active: list[str], N: int,
        P_j: list[float], F: float, f: dict[str, float],
        feed_info: list[tuple[int, float, float]],
        u: list[float], w: list[float],
        draws: list[tuple[int, str, float]],
        duties: list[tuple[int, float]], Q_stage: list[float],
        D: float, B: float, q: float, partial: bool,
        u0: float, v0: float, feed_h: float, method: str,
        T: list[float], x: list[dict[str, float]], y: list[dict[str, float]],
        K: list[dict[str, float]], L: list[float], V: list[float],
        hL: list[float], hV: list[float],
        it: int, dT: float, dF: float, omega: float, n_prop: int,
    ) -> dict[str, PortStream]:
        """Assemble products, duties and the design record from a *converged*
        non-reboiled stage state — shared by every non-reboiled method
        (bubble_point / inside_out / sum_rates). The mass balance closes by
        difference (machine-exact), the condenser duty closes the overall
        energy balance exactly (there is no reboiler), and the independently
        computed condenser-stage balance is reported as ``energy_residual_rel``.
        """
        floor = _V_FLOOR_FRAC * F
        if any(lj <= floor for lj in L) or any(vj <= floor for vj in V[1:]):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: a column section dried up at "
                f"the converged point (min L={min(L):.3g}, "
                f"min V={min(V[1:]):.3g} mol/s) — the specified D={D:.4g} "
                f"mol/s and side draws are infeasible for this feed"
            )

        # -- products (mass balance closed exactly by difference) ----------------
        x_d = dict(y[0]) if partial else dict(x[0])
        d_flows = {c: D * x_d[c] for c in active}
        draw_flows: list[dict[str, float]] = []
        for stage, phase, rate in draws:
            comp = x[stage - 1] if phase == "liquid" else y[stage - 1]
            draw_flows.append({c: rate * comp[c] for c in active})
        b_flows = {c: f[c] - d_flows[c]
                   - sum(df[c] for df in draw_flows) for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: converged distillate/side draws "
                f"would carry more of a component than the feed supplies "
                f"(bottoms flow {neg:.3g} mol/s) — tighten tolerances or "
                f"check the rate specs"
            )
        if neg < 0.0:
            # Same trace-negative handling as the bubble-point path: keep
            # every per-component balance machine-exact without touching D.
            donor = max(b_flows, key=lambda c: b_flows[c])
            for c, v in b_flows.items():
                if v < 0.0:
                    d_flows[c] += v
                    d_flows[donor] -= v
                    b_flows[c] = 0.0
            b_flows[donor] = (f[donor] - d_flows[donor]
                              - sum(df[donor] for df in draw_flows))
            x_d = {c: v / D for c, v in d_flows.items()}
        x_b = {c: v / B for c, v in b_flows.items()}

        z_dist = {c: x_d.get(c, 0.0) for c in comps}
        z_bot = {c: x_b.get(c, 0.0) for c in comps}
        h_D = hV[0] if partial else hL[0]
        # The bottoms leave as liquid at the bottom-stage temperature — a
        # steam-stripped liquid is below its dry-EOS bubble point, so no
        # saturation flash is imposed here.
        h_B = pp.enthalpy_liquid(T[N - 1], P_j[N - 1], x_b)

        distillate = Stream(
            id=f"{self.id}.distillate", components=comps,
            T=T[0], P=P_j[0], molar_flow=D, z=z_dist, H=h_D,
            phase="vapor" if partial else "liquid",
            vapor_fraction=1.0 if partial else 0.0,
        )
        bottoms = Stream(
            id=f"{self.id}.bottoms", components=comps,
            T=T[N - 1], P=P_j[N - 1], molar_flow=B, z=z_bot, H=h_B,
            phase="liquid", vapor_fraction=0.0,
        )
        side_streams: list[Stream] = []
        side_h = 0.0
        for i, (stage, phase, rate) in enumerate(draws):
            j = stage - 1
            comp = x[j] if phase == "liquid" else y[j]
            h_draw = hL[j] if phase == "liquid" else hV[j]
            side_h += rate * h_draw
            side_streams.append(Stream(
                id=f"{self.id}.side{i + 1}", components=comps,
                T=T[j], P=P_j[j], molar_flow=rate,
                z={c: comp.get(c, 0.0) for c in comps},
                H=h_draw, phase=phase,
                vapor_fraction=0.0 if phase == "liquid" else 1.0,
            ))

        # The condenser duty closes the overall energy balance exactly
        # (there is no reboiler); the independently computed condenser-stage
        # balance is the diagnostic residual.
        q_cond = D * h_D + B * h_B + side_h - feed_h - sum(Q_stage)
        q_cond_stage = (L[0] + u0) * hL[0] + v0 * hV[0] - V[1] * hV[1]
        e_scale = max(abs(q_cond), 1.0)
        energy_residual = abs(q_cond - q_cond_stage) / e_scale

        self._warm = {"method": "nonreboiled", "n_stages": N,
                      "active": list(active),
                      "z": {c: f[c] / F for c in active},
                      "T": list(T), "x": [dict(r) for r in x],
                      "y": [dict(r) for r in y], "K": [dict(r) for r in K],
                      "L": list(L), "V": list(V)}

        self.design = {
            # FUG-compatible sizing keys ("N" counts the trays — every stage
            # but the condenser; there is no reboiler). T_top_dew for the
            # condenser LMTD is the temperature of the vapor entering it
            # (saturated at stage-2 conditions).
            "N": float(N - 1), "P": P_j[0], "x_D": dict(x_d), "x_B": dict(x_b),
            "T_top": T[0], "T_top_dew": T[1], "T_bottom": T[N - 1],
            "V_top": V[1], "Q_condenser": q_cond, "Q_reboiler": 0.0,
            "D": D, "B": B, "R": L[0] / D, "q": q,
            "feed_stage": feed_info[0][0] + 1,
            "partial_condenser": partial,
            "method": method, "reboiled": False,
            "feeds": [{"stage": j + 1, "F": Fi, "q": qi}
                      for j, Fi, qi in feed_info],
            "side_draws": [{"stage": stage, "phase": phase, "rate": rate}
                           for stage, phase, rate in draws],
            "stage_duties": [{"stage": stage, "duty": duty}
                             for stage, duty in duties],
            "distillate_flows": dict(d_flows), "bottoms_flows": dict(b_flows),
            "side_draw_flows": [dict(df) for df in draw_flows],
            "n_stages": N,
            "T_profile": list(T),
            "P_profile": list(P_j),
            "L_profile": list(L),
            "V_profile": list(V),
            "x_profile": [{c: row.get(c, 0.0) for c in comps} for row in x],
            "y_profile": [{c: row.get(c, 0.0) for c in comps} for row in y],
            "iterations": it, "max_dT": dT, "max_dV_rel": dF,
            "damping": omega, "flash_calls": n_prop,
            "energy_residual_rel": energy_residual,
        }
        out: dict[str, PortStream] = {
            "distillate": distillate,
            "bottoms": bottoms,
            "condenser_duty": EnergyStream(id=f"{self.id}.condenser_duty",
                                           duty=q_cond),
            "reboiler_duty": EnergyStream(id=f"{self.id}.reboiler_duty",
                                          duty=0.0),
        }
        for i, s in enumerate(side_streams):
            out[f"side{i + 1}"] = s
        return out

    # -- inside-out (Boston-Sullivan) core ---------------------------------------
    def _inside_out_loops(
        self, pp, N: int, P_j: list[float], active: list[str],
        fz_stage: list[dict[str, float]], fh_stage: list[float],
        u: list[float], w: list[float], A: list[float], B: float,
        v0: float, u0: float, floor: float, Q_stage: list[float],
        duties: list[tuple[int, float]], hf_c: dict[str, float],
        fh_form: list[float], F_ge: list[float], UW_ge: list[float],
        T: list[float], x: list[dict[str, float]],
        y: list[dict[str, float]], L: list[float], V: list[float],
        t_bounds: tuple[float, float], max_iter: int,
    ):
        """Inside-out MESH with a damped-Newton inner loop (Boston & Sullivan
        1974; Russell 1983; Seader, Henley & Roper 3e secs. 10.4-10.5) for the
        non-reboiled steam-stripped main fractionator. A Newton inner loop is
        used (not successive substitution), because the strong temperature-flow
        coupling of a multi-draw / pumparound crude tower drives every
        substitution-based update into a coupled limit cycle.

        * **Outer loop** — at the current (T, x, y) profile, evaluate the
          rigorous property package once per stage and fit per-stage *inside*
          models held constant through the inner loop:
            - a base K-value ``ln Kb_j = a_j - b_j / T_j`` (Clausius-Clapeyron
              form; ``b_j`` from a finite-difference slope of the
              mole-fraction-weighted log-mean K), and relative volatilities
              ``alpha_i,j = K_i,j / Kb_j`` (frozen inside);
            - a sensible-basis enthalpy model linear in T,
              ``hLs_j(T) = hLs_j + cpL_j (T - T0_j)`` and likewise for the
              vapor, with each molar enthalpy's composition-weighted
              formation-enthalpy offset (``hf_c``) subtracted so the open-steam
              bottom stays well-conditioned (the per-mole latent pivot inverts
              sign when a heavy liquid descends past a much lighter vapour).
        * **Inner loop** — a damped Newton on the ``2N-1`` tear variables (the
          stage temperatures ``T_0..T_{N-1}`` and vapour flows ``V_1..V_{N-1}``;
          the liquid traffic follows from the total balance
          ``L_j = V_{j+1} + A_j``). The residuals are the ``N`` per-stage
          summation defects (``sum_i K_i,j x_i,j - 1``, with x from the
          component balances) and the ``N-1`` interior/bottom energy balances
          (scaled by a *fixed* reference enthalpy flow, so the Newton cannot
          drive the traffic to a spurious infinite-recycle root). The step is a
          truncated-SVD (regularized) Newton step inside a trust region, with a
          sum-of-squares backtracking line search — a Newton step inverts the
          temperature-flow coupling that successive substitution cannot. The
          Jacobian is cheap (finite differences over the frozen models, no
          rigorous calls).

        The single specification is the distillate rate (the reflux ratio is an
        output); the condenser duty absorbs the stage-1 balance downstream.

        **Scope / status.** This converges narrow-to-medium wide-boiling
        steam-stripped towers and on the light crude tower reproduces the
        bubble-point method's solution to machine precision (a cross-method
        validation; see ``tests/test_m15_crude_column``). It does NOT yet
        converge a *full resid-bearing* crude tower: at an extreme wide-boiling
        side-draw stage the lumped summation residual loses its sensitivity to
        the stage temperature (the direct K-sensitivity is cancelled by the
        liquid-composition coupling), so the reduced (T, V) Jacobian is
        rank-deficient there and the stage cannot be driven to equilibrium. The
        robust completion is the full Naphtali-Sandholm formulation with
        per-component flow variables (equilibrium enforced per component), which
        does not have this degeneracy — tracked as the next step.

        Returns ``(T, x, y, K, L, V, hL, hV, n_outer, dT, dF, n_prop,
        converged)`` with ``hL``/``hV`` the rigorous (formation-inclusive)
        molar enthalpies at the converged profile."""
        t_lo, t_hi = t_bounds

        # -- component material balances (Thomas per component) -----------------
        # Returns the normalized stage liquids and the per-stage summation
        # defect ``s_j = sum_i x_i,j`` (the unknowns solved here are the liquid
        # mole fractions; with consistent flows they sum to 1, so the defect is
        # the sum-rates correction factor for the liquid traffic).
        def component_balances(Lf, Vf, K):
            cols: dict[str, list[float]] = {}
            for c in active:
                a_ = [0.0] + [Lf[j - 1] for j in range(1, N)]
                b_ = [-(Lf[0] + u0 + v0 * K[0][c])] + \
                     [-(Lf[j] + u[j] + (Vf[j] + w[j]) * K[j][c])
                      for j in range(1, N)]
                cc = [Vf[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
                d_ = [-fz_stage[j].get(c, 0.0) for j in range(N)]
                cols[c] = _thomas(a_, b_, cc, d_)
            xs: list[dict[str, float]] = []
            defect: list[float] = []
            ok = True
            for j in range(N):
                row = {c: max(cols[c][j], _X_FLOOR) for c in active}
                tot = sum(row.values())
                if tot <= 0.0:
                    ok = False
                defect.append(tot)
                xs.append({c: vv / tot for c, vv in row.items()})
            return xs, defect, ok

        def hf_mix(row):
            return sum(vv * hf_c.get(c, 0.0) for c, vv in row.items())

        def kb_of(a_j, b_j, T_j):
            return math.exp(min(max(a_j - b_j / T_j, -50.0), 50.0))

        n_prop = 0
        converged = False
        dT = dF = math.inf
        Kb = [1.0] * N
        it = 0
        for it in range(1, max_iter + 1):
            ramp = min(1.0, it / _IO_DUTY_RAMP) if duties else 1.0
            fh_eff = [fh + ramp * qd for fh, qd in zip(fh_stage, Q_stage)]
            fhs = [fe - ff for fe, ff in zip(fh_eff, fh_form)]

            # -- OUTER: rigorous evaluation + inside-model fit at (T, x) --------
            alpha: list[dict[str, float]] = []
            a_par = [0.0] * N
            b_par = [0.0] * N
            T0 = list(T)
            hLs0 = [0.0] * N
            hVs0 = [0.0] * N
            cpL = [0.0] * N
            cpV = [0.0] * N
            for j in range(N):
                # settle the equilibrium vapor and the rigorous K at (T_j, x_j)
                yj = dict(y[j])
                Kj = y[j]
                for _ in range(3):
                    Kj = pp.k_values(T[j], P_j[j], x[j], yj)
                    s = sum(Kj[c] * x[j][c] for c in active)
                    yj = {c: max(Kj[c] * x[j][c], _X_FLOOR) / s for c in active}
                y[j] = yj
                Kj2 = pp.k_values(T[j] + _IO_DH_DT, P_j[j], x[j], yj)
                n_prop += 4
                lnkb1 = sum(x[j][c] * math.log(Kj[c]) for c in active)
                lnkb2 = sum(x[j][c] * math.log(Kj2[c]) for c in active)
                inv_diff = 1.0 / T[j] - 1.0 / (T[j] + _IO_DH_DT)
                bj = (lnkb2 - lnkb1) / inv_diff if inv_diff != 0.0 else 0.0
                bj = (min(max(bj, _IO_B_LO), _IO_B_HI) if bj > 0.0
                      else _IO_B_DEFAULT)
                kb = math.exp(lnkb1)
                a_par[j] = lnkb1 + bj / T[j]
                b_par[j] = bj
                Kb[j] = kb
                alpha.append({c: Kj[c] / kb for c in active})
                # sensible-basis enthalpy model (formation offset subtracted)
                hl0 = pp.enthalpy_liquid(T[j], P_j[j], x[j]) - hf_mix(x[j])
                hv0 = pp.enthalpy_vapor(T[j], P_j[j], yj) - hf_mix(yj)
                hl1 = (pp.enthalpy_liquid(T[j] + _IO_DH_DT, P_j[j], x[j])
                       - hf_mix(x[j]))
                hv1 = (pp.enthalpy_vapor(T[j] + _IO_DH_DT, P_j[j], yj)
                       - hf_mix(yj))
                n_prop += 4
                hLs0[j] = hl0
                hVs0[j] = hv0
                cpL[j] = (hl1 - hl0) / _IO_DH_DT
                cpV[j] = (hv1 - hv0) / _IO_DH_DT

            # -- INNER: damped Newton on the (T, V) tear variables -------------
            # Frozen alpha / base-K / enthalpy models. The 2N-1 inner unknowns
            # are the stage temperatures T_0..T_{N-1} and the vapour flows
            # V_1..V_{N-1} (V_0 = the condenser product is fixed; the liquid
            # traffic follows from the total balance L_j = V_{j+1} + A_j). The
            # 2N-1 residuals are the per-stage summation defects
            # (sum_i K_i,j x_i,j - 1, with x from the component balances) and
            # the interior/bottom-stage energy balances (sensible basis). A
            # Newton step inverts the full temperature-flow coupling that
            # successive substitution cannot, which is what tames the
            # wide-boiling, multi-draw crude tower. The Jacobian is cheap
            # (finite differences over the frozen models, no rigorous calls).
            nv = 2 * N - 1
            # FIXED reference enthalpy flow for the energy residuals — a typical
            # molar enthalpy times the total feed. Scaling each energy balance
            # by a fixed reference (not a traffic-dependent one) is essential:
            # a per-stage |L*h| scale would make the energy residual vanish at
            # infinite traffic (huge cancelling terms over a huge scale), and
            # the Newton would happily run the traffic away to that spurious
            # minimum. With a fixed scale, high traffic is properly penalised.
            e_ref = max(max(abs(hLs0[j]) for j in range(N)),
                        max(abs(hVs0[j]) for j in range(N)),
                        1.0) * max(F_ge[0], 1.0)

            def io_resid(vrs):
                Tv = [min(max(float(vrs[j]), t_lo), t_hi) for j in range(N)]
                Vv = [v0] + [max(float(vrs[N + j - 1]), floor)
                             for j in range(1, N)]
                Lv = [max(Vv[j + 1] + A[j], floor)
                      for j in range(N - 1)] + [B]
                Kbv = [kb_of(a_par[j], b_par[j], Tv[j]) for j in range(N)]
                Km = [{c: alpha[j][c] * Kbv[j] for c in active}
                      for j in range(N)]
                xv, _defect, ok = component_balances(Lv, Vv, Km)
                if not ok:
                    return None, None, None, None, False
                hLs = [hLs0[j] + cpL[j] * (Tv[j] - T0[j]) for j in range(N)]
                hVs = [hVs0[j] + cpV[j] * (Tv[j] - T0[j]) for j in range(N)]
                r = np.empty(nv)
                for j in range(N):                          # summation defects
                    r[j] = sum(Km[j][c] * xv[j][c] for c in active) - 1.0
                for k, j in enumerate(range(1, N)):         # energy balances
                    ein = (Lv[j - 1] * hLs[j - 1]
                           + (Vv[j + 1] * hVs[j + 1] if j < N - 1 else 0.0)
                           + fhs[j])
                    eout = (Lv[j] + u[j]) * hLs[j] + (Vv[j] + w[j]) * hVs[j]
                    r[N + k] = (ein - eout) / e_ref
                return r, Tv, Vv, xv, True

            # Trust region around the outer-fit point: the inside models are
            # only valid near (T, V) where they were fitted, so the inner Newton
            # is confined to a box (+-dT in temperature, a bounded ratio in
            # vapour flow). This globalizes the Newton — it cannot run away to a
            # spurious high-traffic root — while the outer loop walks the box in.
            lo = np.empty(nv)
            hi = np.empty(nv)
            for j in range(N):
                lo[j] = max(t_lo, T[j] - _IO_TR_DT)
                hi[j] = min(t_hi, T[j] + _IO_TR_DT)
            for j in range(1, N):
                lo[N + j - 1] = max(floor, V[j] / _IO_TR_VR)
                hi[N + j - 1] = max(V[j] * _IO_TR_VR, floor * 100.0)
            vrs = np.array([T[j] for j in range(N)]
                           + [V[j] for j in range(1, N)], dtype=float)
            vrs = np.minimum(np.maximum(vrs, lo), hi)
            r, Tv, Vv, x_in, ok = io_resid(vrs)
            inner_ok = ok
            inner_used = 0
            rin = float(np.max(np.abs(r))) if ok else math.inf
            for _nit in range(_IO_NEWTON_MAX):
                if not inner_ok:
                    break
                inner_used = _nit + 1
                rin = float(np.max(np.abs(r)))
                if rin < _IO_NEWTON_TOL:
                    break
                rss = float(r @ r)                  # ||R||^2 merit
                jac = np.empty((nv, nv))
                for k in range(nv):
                    dv = (_IO_DT_PERT if k < N
                          else max(_IO_DV_FRAC * abs(vrs[k]), _IO_DV_MIN))
                    vp = vrs.copy()
                    vp[k] += dv
                    rp, _, _, _, okp = io_resid(vp)
                    jac[:, k] = (rp - r) / dv if okp else 0.0
                # Truncated-SVD (regularized) Newton step: lstsq drops the
                # singular directions of the finite-difference Jacobian (an
                # insensitive stage makes J rank-deficient and a plain solve()
                # then throws an unbounded step into that null space, which the
                # line search rejects and the iteration stalls). The min-norm
                # least-squares step stays bounded and descends. Backtracking
                # on the sum-of-squares merit globalizes it.
                delta = np.linalg.lstsq(jac, -r, rcond=_IO_RCOND)[0]
                step = 1.0
                accepted = False
                for _ls in range(_IO_LS_MAX):
                    vt = np.minimum(np.maximum(vrs + step * delta, lo), hi)
                    rt, Tt, Vt, xt, okt = io_resid(vt)
                    if okt and float(rt @ rt) < (1.0 - 1e-4 * step) * rss:
                        vrs, r, Tv, Vv, x_in = vt, rt, Tt, Vt, xt
                        rin = float(np.max(np.abs(r)))
                        accepted = True
                        break
                    step *= 0.5
                if not accepted:
                    break
            if not inner_ok:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: inside_out inner balances "
                    f"left a stage with no liquid (outer pass {it}); the "
                    f"column is infeasible as specified"
                )
            T_in, V_in = list(Tv), list(Vv)
            L_in = [max(V_in[j + 1] + A[j], floor)
                    for j in range(N - 1)] + [B]
            Kb = [kb_of(a_par[j], b_par[j], T_in[j]) for j in range(N)]
            y_in = []
            for j in range(N):
                yr = {c: max(alpha[j][c] * Kb[j] * x_in[j][c], _X_FLOOR)
                      for c in active}
                tot = sum(yr.values())
                y_in.append({c: vv / tot for c, vv in yr.items()})

            # -- accept the inner (Newton) solution; outer SS over the rigorous
            # model refit. The Newton converges the frozen models exactly, so
            # the inner solution is accepted in full; the outer move is the
            # change between successive rigorous fits.
            dT = max(abs(T_in[j] - T[j]) for j in range(N))
            dF = max(abs(L_in[j] - L[j]) / max(L_in[j], L[j], floor * 10.0)
                     for j in range(N))
            rin = float(np.max(np.abs(r)))
            T, V, L, x, y = T_in, V_in, L_in, x_in, y_in

            if getattr(self, "_io_debug", False):
                km = int(np.argmax(np.abs(r)))
                rkind = (f"sum@{km}" if km < N else f"egy@{km - N + 1}")
                print(f"  io outer={it:3d} dT={dT:.2e} dF={dF:.2e} "
                      f"newton={inner_used}/{rin:.1e}({rkind}) "
                      f"T0={T[0]:.0f} Tbot={T[-1]:.0f} R={L[0]/max(V[1],1e-9):.2f} "
                      f"V1={V[1]:.0f} Vbot={V[-1]:.0f}")

            if (ramp >= 1.0 and rin < _IO_NEWTON_ACCEPT
                    and dT < _IO_OUTER_TOL_T and dF < _IO_OUTER_TOL_F):
                converged = True
                break

        # rigorous (formation-inclusive) enthalpies + K at the converged state
        hL = [pp.enthalpy_liquid(T[j], P_j[j], x[j]) for j in range(N)]
        hV = [pp.enthalpy_vapor(T[j], P_j[j], y[j]) for j in range(N)]
        n_prop += 2 * N
        K_final = [{c: alpha[j][c] * Kb[j] for c in active} for j in range(N)]
        return (T, x, y, K_final, L, V, hL, hV, it, dT, dF, n_prop, converged)

    # -- Naphtali-Sandholm simultaneous correction -------------------------------
    def _solve_ns(self, pp, N: int, P_j: list[float], active: list[str],
                  fz_stage: list[dict[str, float]], fh_stage: list[float],
                  u: list[float], w: list[float], Q_stage: list[float],
                  D: float, partial: bool,
                  T: list[float], x: list[dict[str, float]],
                  y: list[dict[str, float]], L: list[float], V: list[float],
                  t_bounds: tuple[float, float], max_iter: int,
                  force_fd: bool = False, total_condenser: bool = False,
                  reflux: float = 1.0, tol: float = _NS_TOL,
                  decant: bool = False, condenser_T: float = 0.0,
                  reflux_organic: bool = True):
        """Naphtali & Sandholm (1971) simultaneous-correction method for the
        non-reboiled, partial-condenser steam-stripped main fractionator —
        Seader, Henley & Roper 3e sec. 10.4. ALL the MESH equations are solved
        at once by a damped Newton, with per-component flow variables, so the
        method has no division-of-labour assumption and (unlike the reduced
        inside-out Newton) no draw-stage degeneracy: it converges the full
        resid-bearing crude atmospheric tower.

        Variables (``N(2C+1)`` of them): the component liquid flows ``l_i,j``,
        the component vapour flows ``v_i,j`` and the temperature ``T_j`` of each
        stage. The stage totals and mole fractions are ``L_j = sum_i l_i,j``,
        ``V_j = sum_i v_i,j``, ``x_i,j = l_i,j / L_j``, ``y_i,j = v_i,j / V_j``.

        Residuals (``N(2C+1)``):
          * **M** (component material, per stage & component) —
            ``l_i,j(1+U_j/L_j) + v_i,j(1+W_j/V_j) - l_i,j-1 - v_i,j+1 - f_i,j``;
          * **E** (equilibrium, per stage & component) —
            ``v_i,j - K_i,j (V_j/L_j) l_i,j`` (equilibrium enforced *per
            component* — this is what removes the lumped-summation degeneracy);
          * **H** (energy, per interior/bottom stage) — the stage enthalpy
            balance, scaled by a fixed reference enthalpy flow; the condenser
            energy balance is replaced by the **distillate-rate specification**
            ``sum_i v_i,0 = D`` (its duty is recovered downstream).

        The Jacobian is built by finite differences with per-stage property
        caching (perturbing a stage's variables re-evaluates only that stage's
        rigorous K and enthalpies, since the residuals couple only to the
        nearest neighbours). A bounded, line-searched Newton step keeps every
        component flow positive and the temperatures in range. Returns
        ``(T, x, y, K, L, V, hL, hV, n_iter, resid, n_prop, converged)``.
        """
        t_lo, t_hi = t_bounds
        C = len(active)
        span = 2 * C + 1
        # A total-condenser/reboiler column adds ONE global unknown -- the
        # reboiler duty -- and ONE equation (the stage-0 bubble point), so the
        # uniform per-stage layout is preserved and only the system is extended
        # by one. The extra variable lives at index N*span; the extra residual
        # at the same index.
        # An INTEGRATED DECANTING CONDENSER (decant=True) is also a reboiled
        # column with one extra global unknown (the reboiler duty) and one extra
        # equation, but the extra equation is the AQUEOUS-DISTILLATE-RATE spec
        # (sum_i a_i = D) rather than the total-condenser bubble point, and stage
        # 0 stays an ORDINARY equilibrium tray whose overhead vapour is decanted
        # at ``condenser_T`` into an organic layer (refluxed in full, the
        # liquid-from-above on stage 0) and an aqueous layer (the net
        # distillate). Both modes share the +1 layout below.
        glob = total_condenser or decant
        nv = N * span + (1 if glob else 0)
        qi_reb = N * span                                # index of the reb-duty var
        P = np.array(P_j, dtype=float)
        Fcomp = np.zeros((N, C))
        for j in range(N):
            for c, val in fz_stage[j].items():
                if c in active:
                    Fcomp[j, active.index(c)] = val
        Fh = np.array(fh_stage, dtype=float)
        U = np.array(u, dtype=float)
        W = np.array(w, dtype=float)
        Q = np.array(Q_stage, dtype=float)
        Ftot = max(float(Fcomp.sum()), 1.0)

        # -- initial NS variables from the seed profile ------------------------
        lmat = np.full((N, C), _NS_FLOOR)
        vmat = np.full((N, C), _NS_FLOOR)
        for j in range(N):
            for i, c in enumerate(active):
                lmat[j, i] = max(L[j] * x[j].get(c, 0.0), _NS_FLOOR)
                vmat[j, i] = max(V[j] * y[j].get(c, 0.0), _NS_FLOOR)
        Tv = np.array(T, dtype=float)
        n_prop = 0

        def stage_eval(j, lj, vj, Tj):
            nonlocal n_prop
            Lj = float(lj.sum())
            Vj = float(vj.sum())
            xj = {active[i]: lj[i] / Lj for i in range(C)}
            yj = {active[i]: vj[i] / Vj for i in range(C)}
            Kd = pp.k_values(Tj, P[j], xj, yj)
            Karr = np.array([Kd[active[i]] for i in range(C)])
            hLj = pp.enthalpy_liquid(Tj, P[j], xj)
            hVj = pp.enthalpy_vapor(Tj, P[j], yj)
            n_prop += 3
            return Lj, Vj, Karr, hLj, hVj

        # property cache over all stages
        Lc = np.empty(N)
        Vc = np.empty(N)
        Kc = np.empty((N, C))
        hLc = np.empty(N)
        hVc = np.empty(N)

        def refresh_all():
            for j in range(N):
                Lc[j], Vc[j], Kc[j], hLc[j], hVc[j] = stage_eval(
                    j, lmat[j], vmat[j], Tv[j])

        refresh_all()
        e_ref = max(float(np.abs(hLc).max()), float(np.abs(hVc).max()),
                    1.0) * Ftot
        # Reboiler-duty unknown (total condenser / decanting condenser): seed
        # from the bottom-stage energy balance of the initial profile.
        qreb = (Lc[N - 1] * hLc[N - 1] + Vc[N - 1] * hVc[N - 1]
                - Lc[N - 2] * hLc[N - 2]) if glob else 0.0

        P0 = float(P[0])

        def decant_overhead(vm0):
            """Decant the stage-0 overhead vapour ``vm0`` (a per-component flow
            array, total ``V_0``) at ``condenser_T``: a total condenser whose
            condensate settles into an organic and an aqueous liquid layer. The
            organic layer is refluxed in full (the cold liquid-from-above on
            stage 0); the aqueous layer is the net distillate. Returns
            ``(O, o_arr, A, hL_o, hL_a)`` -- the organic total flow, its
            per-component flows, the aqueous total flow, and the two layer molar
            enthalpies at ``condenser_T``. Component mass is conserved exactly
            (``o_arr + a = vm0``) by construction, so the column mass balance
            closes regardless of the flash's residual-vapour bookkeeping."""
            nonlocal n_prop
            V0 = float(vm0.sum())
            y0 = {active[i]: max(float(vm0[i]), _NS_FLOOR) / V0 for i in range(C)}
            r3 = pp.flash_pt_3p(condenser_T, P0, y0)
            n_prop += 1
            bl, bh = r3.beta_light, r3.beta_heavy
            bliq = bl + bh
            if (bl <= _DECANT_BETA_MIN or bh <= _DECANT_BETA_MIN
                    or bliq <= 0.0 or r3.x_light is None or r3.x_heavy is None):
                # No clean two-liquid split yet (early Newton iterate, or a
                # momentarily single-phase overhead): reflux a fixed fraction so
                # the internal traffic stays sane until the decant forms.
                frac_o = _DECANT_FALLBACK_R / (_DECANT_FALLBACK_R + 1.0)
                o_arr = frac_o * vm0
                org = float(o_arr.sum())
                hbulk = pp.enthalpy_liquid(condenser_T, P0, y0)
                return org, o_arr, V0 - org, hbulk, hbulk
            if reflux_organic:                 # organic = light (less dense)
                xo, hL_o, frac_o = r3.x_light, r3.H_light, bl / bliq
                hL_a = r3.H_heavy
            else:                              # reflux the heavy (aqueous) layer
                xo, hL_o, frac_o = r3.x_heavy, r3.H_heavy, bh / bliq
                hL_a = r3.H_light
            org = V0 * frac_o
            o_arr = np.minimum(
                np.array([org * xo.get(active[i], 0.0) for i in range(C)]), vm0)
            org = float(o_arr.sum())
            return org, o_arr, V0 - org, float(hL_o), float(hL_a)

        def residuals(lm, vm, Lq, Vq, Kq, hLq, hVq, qrb=0.0, decant_pre=None):
            R = np.empty(nv)
            # Decanting condenser: the stage-0 overhead is decanted ONCE per
            # residual evaluation (needed by the stage-0 M rows, the stage-0
            # energy row and the appended distillate-rate row). The decant
            # depends ONLY on the stage-0 vapour ``vm[0]``, so the Jacobian build
            # passes the base decant (``decant_pre``) for every perturbation that
            # does not touch a stage-0 vapour variable — sparing ~N*(2C+1) full
            # 3-phase flashes per Jacobian (the dominant cost at book scale).
            if decant:
                if decant_pre is not None:
                    Odec, o_arr, _Adec, hLo_dec, _hLa = decant_pre
                else:
                    Odec, o_arr, _Adec, hLo_dec, _hLa = decant_overhead(vm[0])
            # -- M (component material) + E (equilibrium) rows, vectorized over
            #    the (N, C) flow arrays (residuals couple only to neighbours).
            Mr = (lm * (1.0 + (U / Lq)[:, None])
                  + vm * (1.0 + (W / Vq)[:, None]) - Fcomp)
            Mr[1:, :] -= lm[:-1, :]                       # in_l = l_{i,j-1}
            Mr[:-1, :] -= vm[1:, :]                       # in_v = v_{i,j+1}
            if decant:
                # stage-0 liquid-from-above is the COLD organic reflux (internal)
                Mr[0, :] -= o_arr
            Mr /= Ftot
            Er = (vm - (Kq * (Vq / Lq)[:, None]) * lm) / Ftot
            # -- H (energy) rows, vectorized; stage 0 overwritten per condenser
            #    type. The reboiler (bottom) stage carries the global reb-duty
            #    unknown ``qrb`` (total-condenser AND decanting-condenser are
            #    reboiled).
            in_lH = np.empty(N)
            in_lH[0] = 0.0
            in_vH = np.empty(N)
            in_vH[N - 1] = 0.0
            if N > 1:
                in_lH[1:] = (Lq * hLq)[:-1]
                in_vH[:-1] = (Vq * hVq)[1:]
            qcol = Q.copy()
            if glob:
                qcol[N - 1] += qrb
            Hr = ((Lq + U) * hLq + (Vq + W) * hVq - in_lH - in_vH - Fh - qcol) / e_ref
            if total_condenser:
                # vapour pinned to zero (E rows), reflux pinned L_0 = R*D; T_0
                # set by the appended bubble-point row.
                Er[0, :] = vm[0, :] / Ftot
                Hr[0] = (Lq[0] - reflux * D) / Ftot
            elif decant:
                # stage 0 is an ordinary tray fed the cold organic reflux.
                Hr[0] = ((Lq[0] + U[0]) * hLq[0] + (Vq[0] + W[0]) * hVq[0]
                         - Odec * hLo_dec - (Vq[1] * hVq[1] if N > 1 else 0.0)
                         - Fh[0] - Q[0]) / e_ref
            else:
                Hr[0] = (Vq[0] - D) / D
            R2 = R[:N * span].reshape(N, span)
            R2[:, :C] = Mr
            R2[:, C:2 * C] = Er
            R2[:, 2 * C] = Hr
            if total_condenser:
                # Appended global equation: the condenser liquid is at its
                # bubble point, sum_i K_i,0 x_i,0 = 1 -> determines T_0.
                R[qi_reb] = float((Kq[0] * lm[0]).sum() / Lq[0] - 1.0)
            elif decant:
                # Appended global equation: the AQUEOUS distillate rate
                # sum_i a_i = V_0 - sum_i o_i = D.
                R[qi_reb] = (Vq[0] - Odec - D) / D
            return R

        def build_jac():
            jac = np.empty((nv, nv))
            # The decant flash depends only on stage-0 vapour; precompute it once
            # and reuse for every perturbation that leaves ``vmat[0]`` unchanged.
            base_dec = decant_overhead(vmat[0]) if decant else None
            for j in range(N):
                for kind in range(span):
                    col = j * span + kind
                    lm, vm, tj = lmat, vmat, Tv[j]
                    touches_v0 = (j == 0 and C <= kind < 2 * C)
                    if kind < C:                       # perturb l_{i,j}
                        dk = max(_NS_DL_FRAC * lmat[j, kind], _NS_DL_MIN)
                        lm = lmat.copy()
                        lm[j, kind] += dk
                    elif kind < 2 * C:                 # perturb v_{i,j}
                        i = kind - C
                        dk = max(_NS_DL_FRAC * vmat[j, i], _NS_DL_MIN)
                        vm = vmat.copy()
                        vm[j, i] += dk
                    else:                              # perturb T_j
                        dk = _NS_DT
                        tj = Tv[j] + dk
                    Lj, Vj, Kj, hLj, hVj = stage_eval(j, lm[j], vm[j], tj)
                    Lc2, Vc2 = Lc.copy(), Vc.copy()
                    Kc2, hLc2, hVc2 = Kc.copy(), hLc.copy(), hVc.copy()
                    Lc2[j], Vc2[j], Kc2[j] = Lj, Vj, Kj
                    hLc2[j], hVc2[j] = hLj, hVj
                    dec_pre = None if touches_v0 else base_dec
                    R1 = residuals(lm, vm, Lc2, Vc2, Kc2, hLc2, hVc2, qreb,
                                   decant_pre=dec_pre)
                    jac[:, col] = (R1 - R0) / dk
            if glob:
                # Column for the reboiler-duty unknown. The variable is scaled by
                # e_ref (the duty is ~e_ref in magnitude while the energy residual
                # is divided by e_ref, so the raw column would be ~1/e_ref^2 and
                # the endgame would stall); a unit step in the SCALED variable is
                # e_ref watts, giving an O(1) column. Perturbing the duty changes
                # only the reboiler energy row (no property re-evaluation).
                dk = _NS_DT
                R1 = residuals(lmat, vmat, Lc, Vc, Kc, hLc, hVc, qreb + dk * e_ref,
                               decant_pre=base_dec)
                jac[:, qi_reb] = (R1 - R0) / dk
            return jac

        # active-component indices in the property package's full ordering
        pp_comps = list(getattr(pp, "components", active))
        act_idx = np.array([pp_comps.index(c) for c in active])
        # The analytic Jacobian hardcodes the partial-condenser distillate-rate
        # spec and has no reboiler-duty column; the total-condenser/reboiler and
        # decanting-condenser forms use the finite-difference Jacobian (it
        # differences whatever the residual set, including the appended global
        # equation and the per-residual VLLE decant flash).
        has_analytic = (hasattr(pp, "stage_derivs") and not force_fd
                        and not glob)

        def build_jac_analytic():
            # Per-stage analytic K-value / enthalpy derivatives (Naphtali-
            # Sandholm). One rigorous evaluation per stage (vs N*(2C+1) for the
            # finite-difference Jacobian) and exact — so the Newton direction
            # and the LM gradient are noise-free, which is what makes the
            # wide-boiling resid tower converge reproducibly across platforms.
            SD = []
            for j in range(N):
                xj = {active[i]: lmat[j, i] / Lc[j] for i in range(C)}
                yj = {active[i]: vmat[j, i] / Vc[j] for i in range(C)}
                SD.append(pp.stage_derivs(Tv[j], P[j], xj, yj))
            J = np.zeros((nv, nv))
            ix = np.ix_(act_idx, act_idx)
            for j in range(N):
                d = SD[j]
                K = d["K"][act_idx]
                dlnK_dT = d["dlnK_dT"][act_idx]
                dphL = d["dlnphiL_dns"][ix]
                dphV = d["dlnphiV_dns"][ix]
                HbarL = d["HbarL"][act_idx]
                HbarV = d["HbarV"][act_idx]
                CpL, CpV, hL, hV = d["CpL"], d["CpV"], d["hL"], d["hV"]
                Lj, Vj, Uj, Wj = Lc[j], Vc[j], U[j], W[j]
                Sj = Vj / Lj
                lj, vj = lmat[j], vmat[j]
                b = j * span
                # -- M rows (material balance), scaled by 1/Ftot --------------
                Mdl = (np.diag(np.full(C, 1.0 + Uj / Lj))
                       - (Uj / Lj**2) * lj[:, None])
                Mdv = (np.diag(np.full(C, 1.0 + Wj / Vj))
                       - (Wj / Vj**2) * vj[:, None])
                J[b:b + C, b:b + C] = Mdl / Ftot
                J[b:b + C, b + C:b + 2 * C] = Mdv / Ftot
                if j > 0:
                    pb = (j - 1) * span
                    for i in range(C):
                        J[b + i, pb + i] += -1.0 / Ftot
                if j < N - 1:
                    nb = (j + 1) * span
                    for i in range(C):
                        J[b + i, nb + C + i] += -1.0 / Ftot
                # -- E rows (equilibrium), scaled by 1/Ftot -------------------
                eb = b + C
                # dE_i/dl_k = -K_i S [(l_i/L)(dphL[i,k]-1) + delta_ik]
                Edl = -(K * Sj)[:, None] * (
                    (lj / Lj)[:, None] * (dphL - 1.0) + np.eye(C))
                # dE_i/dv_k = delta_ik + K_i (l_i/L)(dphV[i,k]-1)
                Edv = np.eye(C) + (K * lj / Lj)[:, None] * (dphV - 1.0)
                EdT = -K * dlnK_dT * Sj * lj
                J[eb:eb + C, b:b + C] = Edl / Ftot
                J[eb:eb + C, b + C:b + 2 * C] = Edv / Ftot
                J[eb:eb + C, b + 2 * C] = EdT / Ftot
                # -- H row (energy balance, j>=1), scaled by 1/e_ref ----------
                if j >= 1:
                    hr = b + 2 * C
                    J[hr, b:b + C] = (hL + (Lj + Uj) * (HbarL - hL) / Lj) / e_ref
                    J[hr, b + C:b + 2 * C] = (
                        hV + (Vj + Wj) * (HbarV - hV) / Vj) / e_ref
                    J[hr, b + 2 * C] = (
                        (Lj + Uj) * CpL + (Vj + Wj) * CpV) / e_ref
                    pd = SD[j - 1]
                    pb = (j - 1) * span
                    J[hr, pb:pb + C] = -pd["HbarL"][act_idx] / e_ref
                    J[hr, pb + 2 * C] = -Lc[j - 1] * pd["CpL"] / e_ref
                    if j < N - 1:
                        nd = SD[j + 1]
                        nb = (j + 1) * span
                        J[hr, nb + C:nb + 2 * C] = -nd["HbarV"][act_idx] / e_ref
                        J[hr, nb + 2 * C] = -Vc[j + 1] * nd["CpV"] / e_ref
                else:
                    # spec row (stage 0): sum_k v_{k,0} = D
                    J[b + 2 * C, b + C:b + 2 * C] = 1.0 / D
            return J

        def cap(delta):
            # Fractional-step cap (Naphtali-Sandholm): the global step is scaled
            # so no SIGNIFICANT component flow that merely shrinks drops below
            # _NS_TRUST of its value in one step. Two kinds of flows are
            # deliberately excluded from setting that scale, so neither throttles
            # the whole step:
            #   * trace flows (cur <= thr) -- already protected by the floor;
            #   * flows the Newton step drives clean THROUGH zero (cur + dd <= 0)
            #     -- a component collapsing out of a stage (e.g. a heavy pseudo-
            #     component carried too high up the column by the warm start).
            # The second exclusion is the fix for the wide-boiling resid tower.
            # Without it, ONE such collapsing flow (its Newton step is hugely
            # negative relative to its current value) drives the global scale to
            # ~1e-3, freezing EVERY variable into a ~16-iteration crawl while
            # that flow bleeds down a factor of two per step. The crawl is a
            # marginal, near-stationary phase: the build it was tuned on breaks
            # out of it by iter ~34 and finishes, but it is fragile to small
            # numerical differences and stalls there indefinitely on others
            # (the CI failure: a flat scaled residual ~0.04-0.07). The floor in
            # ``evals`` already keeps these flows positive, so letting them
            # collapse in a single step is both correct and what removes the
            # crawl -- the solve then converges in ~11 clean Newton steps.
            smax = 1.0
            for j in range(N):
                thr_l = _NS_STEP_REL * Lc[j]
                thr_v = _NS_STEP_REL * Vc[j]
                for kind in range(2 * C):
                    cur = lmat[j, kind] if kind < C else vmat[j, kind - C]
                    thr = thr_l if kind < C else thr_v
                    dd = delta[j * span + kind]
                    if dd < 0.0 and cur > thr and cur + dd > 0.0:
                        smax = min(smax, -_NS_TRUST * cur / dd)
            return min(1.0, smax)

        def evals(delta, frac):
            d = frac * delta
            lm = lmat.copy()
            vm = vmat.copy()
            tt = Tv.copy()
            Lt, Vt = Lc.copy(), Vc.copy()
            Kt, hLt, hVt = Kc.copy(), hLc.copy(), hVc.copy()
            for j in range(N):
                b = j * span
                for i in range(C):
                    lm[j, i] = max(lmat[j, i] + d[b + i], _NS_FLOOR)
                    vm[j, i] = max(vmat[j, i] + d[b + C + i], _NS_FLOOR)
                tt[j] = min(max(Tv[j] + d[b + 2 * C], t_lo), t_hi)
                Lt[j], Vt[j], Kt[j], hLt[j], hVt[j] = stage_eval(
                    j, lm[j], vm[j], tt[j])
            # The reboiler-duty step is in the e_ref-scaled variable (see
            # build_jac): a unit Newton step is e_ref watts.
            qrb = qreb + d[qi_reb] * e_ref if glob else 0.0
            Rt = residuals(lm, vm, Lt, Vt, Kt, hLt, hVt, qrb)
            return float(Rt @ Rt), (lm, vm, tt, Lt, Vt, Kt, hLt, hVt, Rt, qrb)

        def _ns_try(delta, rss):
            ss, st = evals(delta, cap(delta))
            return (ss, st) if ss < rss else (None, None)

        def _ns_line_search(delta, rss):
            frac = cap(delta)
            for _ls in range(_NS_LS_MAX):
                ss, st = evals(delta, frac)
                if ss < (1.0 - 1e-4 * frac) * rss:
                    return ss, st
                frac *= 0.5
            return None, None

        R0 = residuals(lmat, vmat, Lc, Vc, Kc, hLc, hVc, qreb)
        converged = False
        it = 0
        rnorm = float(np.max(np.abs(R0)))
        lam = _NS_LM_INIT      # Levenberg-Marquardt damping
        for it in range(1, max_iter + 1):
            rnorm = float(np.max(np.abs(R0)))
            if rnorm < tol:
                converged = True
                break
            rss = float(R0 @ R0)
            jac = (build_jac_analytic() if has_analytic else build_jac())
            if getattr(self, "_ns_jac_check", False) and it <= 1 and not glob:
                jfd = build_jac()
                ja = build_jac_analytic()
                scale = np.maximum(np.abs(jfd), np.abs(ja)) + 1e-12
                rel = np.abs(ja - jfd) / scale
                km = int(np.argmax(rel))
                print(f"  [jac check] max rel diff {rel.max():.2e} at "
                      f"({km // nv},{km % nv}); analytic={ja.flat[km]:.4e} "
                      f"fd={jfd.flat[km]:.4e}")
            # Newton (Gauss-Newton) direction, regularized by lstsq so an
            # insensitive variable cannot throw an unbounded component into the
            # step. A line search along this direction reaches the quadratic
            # endgame (machine precision).
            try:
                delta_n = np.linalg.solve(jac, -R0)
                # Iterative refinement: the wide-boiling MESH Jacobian is
                # ill-conditioned (cond ~1e8-1e12 from the energy rows), so a
                # single LU solve loses several digits of the Newton step. On a
                # favourable BLAS/CPU build the loss is small enough to still
                # plunge to 1e-9; on an unfavourable one it caps the endgame at
                # a marginal residual (~3e-3) and the test flakes across runner
                # CPUs. One or two refinement passes recover the lost digits
                # (delta += J^-1 (-R0 - J delta)), restoring the quadratic
                # endgame reproducibly. Cheap: the iterate is near the root.
                for _ref in range(_NS_REFINE):
                    lin_res = -R0 - jac @ delta_n
                    if (float(np.max(np.abs(lin_res)))
                            <= 1e-13 * float(np.max(np.abs(R0)) + 1e-30)):
                        break
                    delta_n = delta_n + np.linalg.solve(jac, lin_res)
            except np.linalg.LinAlgError:
                delta_n = np.linalg.lstsq(jac, -R0, rcond=_IO_RCOND)[0]
            ss, state = _ns_line_search(delta_n, rss)
            used = "N"
            if ss is None:
                # the Newton direction stalled (not a descent direction at this
                # iterate — a platform-sensitive ill-conditioning). Escape with
                # Levenberg-Marquardt steps of growing damping: each is a
                # shorter, well-scaled blend toward gradient descent that still
                # reduces ||R||^2. The fixed energy scaling + flow bounds keep
                # LM from wandering to a spurious high-traffic state.
                jtj = jac.T @ jac
                jtr = jac.T @ R0
                jdiag = np.maximum(np.diag(jtj), 1e-30)
                lam = max(lam, _NS_LM_INIT)
                for _ in range(_NS_LS_MAX):
                    try:
                        delta_l = np.linalg.solve(
                            jtj + lam * np.diag(jdiag), -jtr)
                    except np.linalg.LinAlgError:
                        delta_l = -jtr / jdiag
                    ss, state = _ns_try(delta_l, rss)
                    if ss is not None:
                        used = "LM"
                        break
                    lam = min(lam * _NS_LM_UP, _NS_LM_MAX)
            if state is None:
                accepted = False
            else:
                accepted = True
                lmat, vmat, Tv, R0 = state[0], state[1], state[2], state[8]
                Lc[:], Vc[:] = state[3], state[4]
                Kc[:], hLc[:], hVc[:] = state[5], state[6], state[7]
                qreb = state[9]
                lam = max(lam / _NS_LM_DOWN, _NS_LM_MIN)
            if getattr(self, "_ns_debug", False):
                km = int(np.argmax(np.abs(R0)))
                jj, kk = km // span, km % span
                rk = "M" if kk < C else "E" if kk < 2 * C else "G"
                print(f"  ns it={it:3d} |R|={rnorm:.2e}({rk}@{jj}) {used} "
                      f"lam={lam:.0e} T0={Tv[0]:.0f} Tbot={Tv[-1]:.0f} "
                      f"R={Lc[0] / max(Vc[1], 1e-9):.2f} V1={Vc[1]:.0f} "
                      f"Vbot={Vc[-1]:.0f}")
            if not accepted:
                break

        # -- assemble the converged profile ------------------------------------
        Lout = [float(Lc[j]) for j in range(N)]
        Vout = [float(Vc[j]) for j in range(N)]
        xout = [{active[i]: float(lmat[j, i] / Lc[j]) for i in range(C)}
                for j in range(N)]
        yout = [{active[i]: float(vmat[j, i] / Vc[j]) for i in range(C)}
                for j in range(N)]
        Kout = [{active[i]: float(Kc[j, i]) for i in range(C)}
                for j in range(N)]
        Tout = [float(Tv[j]) for j in range(N)]
        hLout = [float(hLc[j]) for j in range(N)]
        hVout = [float(hVc[j]) for j in range(N)]
        return (Tout, xout, yout, Kout, Lout, Vout, hLout, hVout,
                it, rnorm, n_prop, converged)

    def _initial_sr_profiles(self, pp, N: int, P_j: list[float],
                             active: list[str], f: dict[str, float],
                             D: float, R0: float,
                             feed_info: list[tuple[int, float, float]],
                             feed_flashes: list,
                             F_ge: list[float], UW_ge: list[float],
                             B: float, v0: float, floor: float):
        """Starting profiles for the non-reboiled mode: the last converged
        profiles when the layout/feed are unchanged; else a Fenske-style
        composition ramp (K-values harvested from the already-computed feed
        flashes — no extra flash), per-stage K-value bubble temperatures
        marched outward from the feed stage, and constant-molal-overflow
        traffic seeded by ``R0``. Returns ``(T, x, y, K, L, V)``."""
        F = sum(f.values())
        z = {c: f[c] / F for c in active}
        w_ = self._warm
        if (w_ is not None and w_.get("method") == "nonreboiled"
                and w_["n_stages"] == N and w_["active"] == active
                and max(abs(w_["z"].get(c, 0.0) - z[c]) for c in active)
                < _WARM_Z_TOL):
            return (list(w_["T"]), [dict(r) for r in w_["x"]],
                    [dict(r) for r in w_["y"]], [dict(r) for r in w_["K"]],
                    list(w_["L"]), list(w_["V"]))

        # K at feed conditions, harvested from the feed PH flashes. A
        # vapor-only feed (stripping steam) contributes no K; its components
        # default to "very light" (4x the largest harvested K) — it only
        # seeds the initializer.
        K_feed: dict[str, float] = {}
        for res in feed_flashes:
            if res.x and res.y:
                for c in active:
                    xa, ya = res.x.get(c, 0.0), res.y.get(c, 0.0)
                    if xa > 0.0 and ya > 0.0 and c not in K_feed:
                        K_feed[c] = ya / xa
        if K_feed:
            kmax = max(K_feed.values())
            for c in active:
                K_feed.setdefault(c, 4.0 * kmax)
        try:
            x = fenske_profile(pp, P_j[0], f, D, N,
                               K_feed=K_feed or None)
        except _NoVLEError as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: no vapor-liquid equilibrium for "
                f"the feed at P={P_j[0]:.4g} Pa; check the pressure"
            ) from exc

        # Per-stage bubble temperatures of the initial liquid profile,
        # marching outward from the (first two-phase) feed stage so every
        # secant starts from its neighbor's solution.
        j_seed, T_seed = feed_info[0][0], float(feed_flashes[0].T)
        for (j0, _, _), res in zip(feed_info, feed_flashes):
            if res.x and res.y:
                j_seed, T_seed = j0, float(res.T)
                break
        T = [0.0] * N
        y: list[dict[str, float]] = [{} for _ in range(N)]
        K: list[dict[str, float]] = [{} for _ in range(N)]
        g = T_seed
        for j in range(j_seed, N):
            T[j], K[j], y[j] = _bubble_T_kvalues(pp, P_j[j], x[j], g)
            g = T[j]
        g = T[j_seed]
        for j in range(j_seed - 1, -1, -1):
            T[j], K[j], y[j] = _bubble_T_kvalues(pp, P_j[j], x[j], g)
            g = T[j]

        # Constant-molal-overflow-style traffic: D(1+R0) above the main feed,
        # the feeds' flashed vapor below it; L from the total balances.
        f0 = feed_info[0][0]
        V = [v0] * N
        for j in range(1, N):
            if j <= f0:
                V[j] = D * (1.0 + R0)
            else:
                vap = sum((1.0 - qk) * Fk for jk, Fk, qk in feed_info
                          if jk >= j)
                V[j] = max(vap, 0.02 * F)
        L = [floor] * N
        for j in range(1, N):
            L[j - 1] = max(V[j] - F_ge[j] + UW_ge[j] + B, floor)
        L[N - 1] = B
        return T, x, y, K, L, V

    # -- initialization ------------------------------------------------------------
    def _initial_x(self, pp, P: float, f: dict[str, float], D: float, N: int,
                   z: dict[str, float]) -> list[dict[str, float]]:
        """Initial liquid-composition profile.

        Warm path: reuse the last converged profile when the column layout is
        unchanged and the feed composition is close (the equation-oriented
        solver and recycle sweeps re-solve with near-identical feeds).
        Cold path: :func:`fenske_profile`.
        """
        active = list(f)
        w = self._warm
        if (w is not None and w.get("method", "bubble_point") == "bubble_point"
                and w["n_stages"] == N and w["active"] == active
                and max(abs(w["z"].get(c, 0.0) - z.get(c, 0.0)) for c in active)
                < _WARM_Z_TOL):
            return [dict(row) for row in w["x"]]
        try:
            return fenske_profile(pp, P, f, D, N)
        except _NoVLEError as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: no vapor-liquid equilibrium for the "
                f"feed at P={P:.4g} Pa (bubble flash returned a single phase); "
                f"check the pressure"
            ) from exc

    @staticmethod
    def _initial_traffic(N: int, feed_info: list[tuple[int, float, float]],
                         R: float, D: float, F: float, partial: bool,
                         A: list[float], w: list[float],
                         ) -> tuple[list[float], list[float]]:
        """Constant-molal-overflow traffic from R, D and the feed qualities.
        V[0] is the condenser vapor product (0 for a total condenser); V[1] is
        fixed by the condenser balance at D(1+R) throughout the iteration.
        Below each feed the vapor drops by its vapor fraction (1-q)F_k, and
        above each vapor side draw by its rate."""
        v_rect = D * (1.0 + R)
        V = [D if partial else 0.0]
        for j in range(1, N):
            vj = v_rect
            for f0_k, F_k, q_k in feed_info:
                if f0_k < j:
                    vj -= (1.0 - q_k) * F_k
            vj -= sum(w[m] for m in range(1, j))
            V.append(max(vj, 0.2 * F))
        L = [V[j + 1] + A[j] for j in range(N - 1)] + [A[N - 1]]
        return [max(lj, _V_FLOOR_FRAC * F) for lj in L], V

    # -- MESH building blocks ---------------------------------------------------
    @staticmethod
    def _component_balances(N: int, active: list[str],
                            fz_stage: list[dict[str, float]], D: float,
                            K: list[dict[str, float]], L: list[float],
                            V: list[float], partial: bool,
                            u: list[float], w: list[float],
                            ) -> list[dict[str, float]]:
        """One tridiagonal solve per component (Thomas algorithm) for the stage
        liquid compositions, then floor and renormalize each stage (the
        summation equations). Coefficients per Seader 3e eqs. 10-8..10-12:
        feeds enter the right-hand side on their stages (``fz_stage``), liquid
        and vapor side draws (``u``/``w``) the diagonal. The distillate is the
        liquid draw off stage 1 of a total condenser (U_1 = D); a partial
        condenser's vapor product is V_1 = D."""
        u1 = 0.0 if partial else D
        cols: dict[str, list[float]] = {}
        for c in active:
            a = [0.0] + [L[j - 1] for j in range(1, N)]
            b = [-(L[0] + u1 + V[0] * K[0][c])] + \
                [-(L[j] + u[j] + (V[j] + w[j]) * K[j][c]) for j in range(1, N)]
            cc = [V[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
            d = [-fz_stage[j].get(c, 0.0) for j in range(N)]
            cols[c] = _thomas(a, b, cc, d)
        x: list[dict[str, float]] = []
        for j in range(N):
            row = {c: max(cols[c][j], _X_FLOOR) for c in active}
            tot = sum(row.values())
            x.append({c: v / tot for c, v in row.items()})
        return x

    def _read_reactions(self, n_stages: int, reboiled: bool, method: str
                        ) -> tuple[list[tuple[KineticReaction, frozenset[int]]], float]:
        """Parse the optional ``reactions`` / ``tray_holdup`` params for reactive
        distillation. Returns a list of ``(KineticReaction, reactive stages as a
        0-based frozenset)`` and the per-stage liquid holdup volume (m^3).

        Each ``reactions`` entry is a JSON-friendly kinetic-reaction dict (see
        :class:`caldyr.unitops.reaction.KineticReaction` — ``stoich``/``key``/
        ``k0``/``Ea``/``orders``) plus a ``"stages": [first, last]`` 1-based,
        inclusive range naming the trays the reaction runs on. ``k0`` is on the
        SI molar basis (concentrations in mol/m^3), ``Ea`` in J/mol. The rate
        ``r = k0·exp(-Ea/RT)·Π C_i^order_i`` [mol/(m^3·s)] is turned into a stage
        extent by the liquid holdup: ``ξ_j = r_j · tray_holdup``."""
        raw = self.params.get("reactions")
        if not raw:
            return [], 0.0
        if method != "bubble_point" or not reboiled:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: reactive distillation ('reactions') "
                f"is implemented for the default reboiled bubble-point column "
                f"(method='bubble_point', reboiled=True); got method={method!r}, "
                f"reboiled={reboiled}"
            )
        try:
            holdup = float(self.params.get("tray_holdup", 1.0))
        except (TypeError, ValueError) as exc:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: 'tray_holdup' (liquid holdup volume "
                f"per reactive stage, m^3) must be numeric"
            ) from exc
        if holdup <= 0.0:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: tray_holdup={holdup} m^3 must be > 0")
        out: list[tuple[KineticReaction, frozenset[int]]] = []
        for k, entry in enumerate(raw):
            try:
                rxn = KineticReaction.from_param(entry)
            except (ValueError, KeyError, TypeError) as exc:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: reactions[{k}] is not a valid "
                    f"kinetic reaction ({exc})"
                ) from exc
            stages_spec = entry.get("stages")
            if (not isinstance(stages_spec, (list, tuple)) or len(stages_spec) != 2):
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: reactions[{k}] needs "
                    f"'stages': [first, last] (1-based, inclusive); got "
                    f"{stages_spec!r}")
            first, last = int(stages_spec[0]), int(stages_spec[1])
            if not 1 <= first <= last <= n_stages:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: reactions[{k}] stages "
                    f"[{first}, {last}] out of range 1..{n_stages}")
            out.append((rxn, frozenset(range(first - 1, last))))
        return out, holdup

    def _stage_generation(self, pp, reactions, holdup: float,
                          x: list[dict[str, float]], T: list[float],
                          P_j: list[float], L: list[float], active: list[str]
                          ) -> list[dict[str, float]]:
        """Per-stage component generation (mol/s) from the kinetic reactions on
        their reactive stages, evaluated at the current liquid profile. Negative
        for consumed reactants, positive for products; sums to the stage net mole
        change (zero for a mole-conserving reaction). Concentrations are the
        liquid molar densities C_i = x_i / v_liquid(T, P, x).

        The extent is capped so a stage cannot consume more than a fraction
        ``_RXN_MAX_CONV`` of any reactant in its liquid traffic — a numerical
        rail that keeps the per-stage flash well-conditioned while the lagged
        reaction source settles (at the converged steady state, which for a
        reversible reaction is equilibrium-bounded, the cap does not bind)."""
        N = len(x)
        g: list[dict[str, float]] = [{} for _ in range(N)]
        for rxn, stages in reactions:
            for j in stages:
                v_liq = pp.volume_liquid(T[j], P_j[j], x[j])     # m^3/mol
                if v_liq <= 0.0:
                    continue
                conc = {c: x[j].get(c, 0.0) / v_liq for c in active}
                xi = rxn.rate(conc, T[j]) * holdup               # mol/s extent
                # Cap by reactant (xi>0) or product (xi<0) availability in L_j.
                for c, nu in rxn.stoich.items():
                    consumed = -nu if xi > 0 else nu            # >0 if depleting c
                    if consumed > 0:
                        avail = _RXN_MAX_CONV * L[j] * x[j].get(c, 0.0)
                        xi = min(xi, avail / consumed) if xi > 0 else \
                            max(xi, -avail / consumed)
                for c, nu in rxn.stoich.items():
                    g[j][c] = g[j].get(c, 0.0) + nu * xi
        return g

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

    def _energy_balances(self, N: int, fh_stage: list[float], A: list[float],
                         u: list[float], w: list[float],
                         hL: list[float], hV: list[float],
                         V: list[float], v_floor: float) -> list[float]:
        """Forward recurrence for the vapor traffic from the stage energy
        balances (Seader 3e eqs. 10-28..10-30), with L eliminated via the
        section total balance L_j = V_{j+1} + A_j (A = cumulative feeds minus
        draws minus distillate). ``fh_stage`` carries the feed enthalpy flow
        F_k*h_Fk landing on each stage; liquid draws u_j leave at hL_j and
        vapor draws w_j at hV_j. V_2 stays fixed at D(1+R) by the condenser
        balance."""
        V_new = list(V)
        for j in range(1, N - 1):           # 0-based stage j (= stage j+1, 1-based)
            denom = hV[j + 1] - hL[j]
            if denom <= 0.0:
                raise RigorousColumnError(
                    f"RigorousColumn {self.id!r}: degenerate energy balance on "
                    f"stage {j + 2} (saturated vapor enthalpy below the liquid's, "
                    f"dH={denom:.4g} J/mol) — the property package returned an "
                    f"unphysical state"
                )
            num = ((A[j] + u[j]) * hL[j] - A[j - 1] * hL[j - 1]
                   - fh_stage[j] - V_new[j] * (hL[j - 1] - hV[j])
                   + w[j] * hV[j])
            V_new[j + 1] = max(num / denom, v_floor)
        return V_new
