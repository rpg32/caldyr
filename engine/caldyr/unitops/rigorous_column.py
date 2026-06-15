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
      * ``reboiled`` — default True. ``False`` removes the reboiler: stage
        ``n_stages`` becomes an ordinary (adiabatic) stage that may receive a
        feed — open stripping steam, the crude-tower configuration. The
        column then has a single specification (the distillate rate);
        ``reflux_ratio`` only seeds the initial traffic and the converged
        reflux ratio is reported in ``design['R']``. Requires
        ``method="sum_rates"``.
      * ``method`` — the MESH temperature/traffic update. Three supported
        combinations (validated by the test suite):

        - ``"bubble_point"`` + ``reboiled=True`` (the **default**): the
          Wang-Henke bubble-point method — a saturated-liquid flash per stage
          sets the temperature and the vapor traffic comes from a forward
          energy recurrence. The narrow/medium-boiling distillation workhorse.
        - ``"bubble_point"`` + ``reboiled=False``: the same per-stage K-value
          bubble temperatures, but the vapor traffic comes from *envelope*
          energy balances (each V_j from its own around-stage-j-to-bottom
          balance) run on a formation-offset-conditioned sensible basis. This
          is the **robust method for a steam-stripped main fractionator**
          (open stripping steam, side draws, pumparound stage duties); it
          converges the crude atmospheric tower (see ``examples/20_crude_tower``
          and ``tests/test_m15_crude_column``).
        - ``"sum_rates"`` + ``reboiled=False``: a Burningham-Otto variant
          (liquid traffic from the component-flow sums; temperatures from a
          simultaneous Newton solve of the stage energy balances). **Provided
          but not recommended**: on equilibrium-dominated (distillation-like)
          sections it limit-cycles — the documented Seader 3e ch. 10.4
          division of labor. Prefer ``"bubble_point"`` for non-reboiled towers.

        ``reboiled=True`` with ``"sum_rates"`` is rejected.
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
        if method not in ("bubble_point", "sum_rates"):
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: method={method!r} must be "
                f"'bubble_point' or 'sum_rates'"
            )
        reboiled = bool(self.params.get("reboiled", True))
        if reboiled and method == "sum_rates":
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: method='sum_rates' is "
                f"implemented for the non-reboiled (steam-stripped main "
                f"fractionator) form only — set reboiled=False, or use the "
                f"default bubble_point method for a reboiled column"
            )
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
        return (n_stages, feed_stages, draws, R, D, P, dP, partial, max_iter,
                method, reboiled, duties)

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
         method, reboiled, duties) = self._read_params(F, P_in)

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

        # Cumulative (feeds - draws - distillate) above-and-including stage j:
        # the total balance reads L_j = V_{j+1} + A_j (A_{N-1} = B).
        A = [0.0] * N
        run = -D
        for j in range(N):
            run += sum(fz_stage[j].values()) - u[j] - w[j]
            A[j] = run

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
        for it in range(1, max_iter + 1):
            x = self._component_balances(N, active, fz_stage, D, K, L, V,
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

            V_new = self._energy_balances(N, fh_stage, A, u, w, hL, hV, V,
                                          v_floor)
            dV = max(abs(vn - vo) for vn, vo in zip(V_new[1:], V[1:])) / max(V_new[1:])
            V = [V[0]] + [max(vo + omega * (vn - vo), v_floor)
                          for vo, vn in zip(V[1:], V_new[1:])]
            L = [max(V[j + 1] + A[j], v_floor) for j in range(N - 1)] + [B]

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
        b_flows = {c: f[c] - d_flows[c]
                   - sum(df[c] for df in draw_flows) for c in active}
        neg = min(b_flows.values())
        if neg < -1e-7 * F:
            raise RigorousColumnError(
                f"RigorousColumn {self.id!r}: converged distillate/side draws "
                f"would carry more of a component than the feed supplies "
                f"(bottoms flow {neg:.3g} mol/s) — tighten tolerances or check "
                f"the rate specs"
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
            b_flows[donor] = (f[donor] - d_flows[donor]
                              - sum(df[donor] for df in draw_flows))
            x_d = {c: v / D for c, v in d_flows.items()}
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
        q_reb = D * res_d.H + B * res_b.H + side_h - feed_h - q_cond
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
            "D": D, "B": B, "R": R, "q": q, "feed_stage": feed_stages[0],
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
        t_bounds = (max(min(T) - 60.0, _BUBBLE_T_LO),
                    min(max(T) + 200.0, _BUBBLE_T_HI))

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
