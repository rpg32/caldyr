"""Absorber / stripper and reboiled absorber: sum-rates (Burningham-Otto) MESH.

Columns *without* a condenser — gas absorbers, strippers and reboiled
absorbers — are wide-boiling: the stage temperatures are governed by the
energy balances (absorption heat, solvent sensible heat), not by bubble
points, so the bubble-point (Wang-Henke) method of
:mod:`caldyr.unitops.rigorous_column` fails on them. The standard workhorse is
the **sum-rates method** of Burningham & Otto (1967), per Seader, Henley &
Roper, *Separation Process Principles* 3e, ch. 10.4:

1. **Initialize** temperatures linearly between the two feed temperatures and
   the traffic at the feed rates (absorber) or a boilup estimate (reboiled
   absorber).
2. **Component balances**: with the current ``K_ij`` (phi-phi, evaluated at
   the stage T — *not* at saturation), ``L_j``, ``V_j``, each component's
   stage balances form a tridiagonal system solved by the Thomas algorithm.
3. **Sum-rates step**: the new total liquid flows are the *sums of the
   component liquid flows* (``L_j = sum_i l_ij``) and the vapor flows follow
   from the total mass balance around the column bottom
   (``V_j = L_{j-1} + sum_{k>=j} F_k - L_N``). This replaces the energy-based
   traffic update of the bubble-point method.
4. **Energy balances**: the N stage enthalpy balances are solved
   *simultaneously* for the new temperature profile by Newton's method with a
   tridiagonal Jacobian (the original Burningham-Otto step), using per-phase
   enthalpies at the stage temperature (``pp.enthalpy_liquid/vapor``).
5. Repeat 2-4 until the temperatures and traffic stop moving. A
   non-convergent or degenerate iteration raises :class:`AbsorberError` with
   diagnostics — never a silent wrong answer.

**Stage convention**: stages are numbered from the top, ``1 .. n_stages``.
The :class:`Absorber` has no condenser or reboiler — the liquid (solvent)
feed enters stage 1, the gas feed enters stage ``n_stages``, the vapor
product leaves stage 1 and the liquid product leaves stage ``n_stages``. A
**stripper is the same unit**: feed the volatile-rich liquid on top and the
stripping gas at the bottom (gas desorption is the mirror image of
absorption — Hameed, *Chemical Process Simulations using Aspen Hysys*, 2025,
sec. 9.2); no separate unit type is needed.

The :class:`ReboiledAbsorber` (stripping tower; Hameed 2025 sec. 9.3.5) is
the same MESH system with the bottom stage replaced by a reboiler: a single
(liquid) feed near the top, the stripping vapor generated internally by the
reboiler heat input. **Method choice**: a reboiled absorber boils a
condensable liquid feed, so its stages are narrow-boiling — exactly the
regime where the sum-rates method oscillates (its energy-balance temperature
step fights the composition-determined boiling temperatures; Seader 3e
ch. 10.4 discusses this division of labor) and the bubble-point method
shines. The ReboiledAbsorber therefore runs a Wang-Henke *bubble-point*
inner loop (stage temperatures from bubble points, vapor traffic from the
energy-balance recurrence) **driven by the overhead vapor rate**: with
``V_top`` fixed, the N-1 upper-stage energy balances determine the full
traffic and the reboiler duty follows from the overall energy balance,
exactly. The heat-input specification is one of ``vapor_rate`` (direct),
``boilup_ratio`` or ``reboiler_duty`` (the latter two close an outer Brent
iteration over ``V_top``; both are strictly monotone in it).

Mass balances close to machine precision (products by difference);
the overall energy balance is closed exactly by resolving the bottom liquid
enthalpy from the column energy balance (its temperature then comes from a
PH flash, which agrees with the converged stage-N temperature to within the
solver tolerance — the residual is reported in
``design['energy_residual_rel']``). The full converged stage profiles
(per-stage T, P, L, V, x, y) are stored on ``unit.design``.

Validation (see tests/test_m12_absorber.py and the module reports there):
the Kremser closed form — Hameed 2025 eq. (9.1) / Seader 3e ch. 5 — for
dilute near-isothermal absorption and stripping, structural monotonicity
(more stages / more solvent -> more absorption), and the worked reboiled
absorber of Hameed 2025 sec. 9.3.5 (n-pentane/n-heptane at 110 kPa, PR).
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.optimize import brentq

from ..core import EnergyStream, Port, Stream, UnitOp
from ..core.unitop import PortStream
from .base import register
from .rigorous_column import _NoVLEError, _thomas, fenske_profile

_MAX_ITER = 120          # default sum-rates iteration cap ('max_iter' overrides)
_TOL_T = 1e-4            # K, stage-temperature convergence
_TOL_FLOW = 1e-7         # relative liquid/vapor traffic convergence
_X_FLOOR = 1e-15         # floor for stage mole fractions
_FLOW_FLOOR_FRAC = 1e-9  # traffic floor, as a fraction of the total feed
_NEWTON_MAX = 60         # inner Newton iterations on the energy balances
_NEWTON_TOL = 1e-10      # relative energy-balance residual target (inner)
_ENERGY_TOL = 1e-9       # relative energy residual required at convergence
_NEWTON_STEP_MAX = 12.0  # K, per-stage Newton step clamp
_DH_DT = 0.05            # K, finite-difference step for dh/dT
_WARM_Z_TOL = 0.1        # max |dz| of the feed for warm-starting

# -- Naphtali-Sandholm (simultaneous-correction) constants, mirroring the
# proven settings in :mod:`caldyr.unitops.rigorous_column`. NS is the robust
# fallback for columns the alternating sum-rates iteration cannot converge —
# notably an amine *regenerator* (reactive desorption + condensing stripping
# steam), where the sum-rates traffic update limit-cycles.
_NS_TOL = 1e-9           # scaled residual infinity-norm convergence
_NS_FLOOR = 1e-12        # mol/s, floor on the component flows
_NS_DL_FRAC = 1e-6       # relative finite-difference step on the flows
_NS_DL_MIN = 1e-10       # mol/s, floor on the flow finite-difference step
_NS_DT = 0.01            # K, finite-difference step on the temperatures
_NS_LS_MAX = 24          # backtracking line-search halvings
_NS_TRUST = 0.5          # max fractional reduction a flow may take in one step
_NS_STEP_REL = 1e-4      # only flows above this fraction of the stage total cap
_NS_LM_INIT = 1e-3       # Levenberg-Marquardt damping: initial value
_NS_LM_MIN = 1e-10       # ... lower bound (-> Gauss-Newton, quadratic endgame)
_NS_LM_MAX = 1e10        # ... upper bound (-> short scaled gradient descent)
_NS_LM_UP = 4.0          # ... growth on a rejected trial
_NS_LM_DOWN = 3.0        # ... shrink on an accepted step
_NS_REFINE = 2           # iterative-refinement passes on the Newton linear solve
_NS_TSPEC_SCALE = 50.0   # K, scaling of a pinned-temperature (condenser) residual


class AbsorberError(ValueError):
    """Specification or convergence error in an Absorber / ReboiledAbsorber
    (bad stage counts, infeasible heat-input spec, sum-rates failure, ...)."""


def _normalized(row: dict[str, float]) -> dict[str, float]:
    tot = sum(row.values())
    return {c: v / tot for c, v in row.items()}


def _solve_sum_rates(
    unit_id: str,
    pp: Any,
    active: list[str],
    N: int,
    P_j: list[float],
    feeds: list[tuple[int, dict[str, float], float]],
    Q: list[float],
    T: list[float],
    x: list[dict[str, float]],
    y: list[dict[str, float]],
    L: list[float],
    V: list[float],
    max_iter: int,
    T_bounds: tuple[float, float],
    murphree: dict[str, float] | None = None,
    bottom_vapor_flows: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Run the sum-rates (Burningham-Otto) MESH iteration.

    ``feeds`` is a list of ``(stage_0based, component_flows mol/s, F*h W)``;
    ``Q`` is the per-stage heat input (W, positive heats the stage). ``T``,
    ``x``, ``y``, ``L``, ``V`` are starting profiles (mutated copies are
    returned in the result dict, never the inputs); ``T_bounds`` brackets the
    physically admissible stage temperatures for the Newton energy solve.

    ``murphree`` is an optional per-component vapor (Murphree) stage efficiency
    ``{component: E}`` (default 1.0 = equilibrium). A finite efficiency lets a
    stage approach equilibrium only partially — essential for the *kinetically*
    limited acid gases in amine sweetening (e.g. CO2 with a tertiary amine: book
    §15.3 uses E_CO2 ~ 0.15, E_H2S ~ 0.8). The Murphree relation
    ``y_j = E*K_j*x_j + (1-E)*y_(j+1)`` is solved *exactly* (not lagged) by
    writing the component balance of an E<1 species in its **vapor** component
    flows ``v_j`` — that keeps the (1-E)*y_(j+1) coupling tridiagonal (see
    :func:`_murphree_vapor_column`). The entering vapor below the bottom stage is
    the gas feed (``bottom_vapor_flows``). E=1 species keep the standard
    liquid-flow tridiagonal (a vapor formulation is ill-conditioned for a
    barely-volatile solvent like the amine). The two formulations are identical
    at E=1, so a mix of efficiencies is consistent.
    """
    T, L, V = list(T), list(L), list(V)
    x = [dict(row) for row in x]
    y = [dict(row) for row in y]

    f_stage: list[dict[str, float]] = [{} for _ in range(N)]
    fh_stage = [0.0] * N
    for j, flows, fh in feeds:
        for c, v in flows.items():
            f_stage[j][c] = f_stage[j].get(c, 0.0) + v
        fh_stage[j] += fh
    F_total = sum(sum(flows.values()) for _, flows, _ in feeds)
    # Total feed entering at or below stage j (for the total mass balance).
    F_below = [0.0] * (N + 1)
    for j in range(N - 1, -1, -1):
        F_below[j] = F_below[j + 1] + sum(f_stage[j].values())
    floor = _FLOW_FLOOR_FRAC * F_total

    omega = 1.0
    dT_prev = math.inf
    worsened = 0
    converged = False
    n_prop = 0
    it = 0
    dT = dF = resid_rel = math.inf
    eff = {c: float((murphree or {}).get(c, 1.0)) for c in active}
    gas_f = {c: max((bottom_vapor_flows or {}).get(c, 0.0), 0.0) for c in active}
    gas_tot = sum(gas_f.values())
    y_gas = {c: (gas_f[c] / gas_tot if gas_tot > 0.0 else 0.0) for c in active}

    for it in range(1, max_iter + 1):
        K = [pp.k_values(T[j], P_j[j], x[j], y[j]) for j in range(N)]
        n_prop += N

        # -- component balances (Thomas) + sum-rates traffic update ----------
        cols: dict[str, list[float]] = {}
        for c in active:
            if eff[c] < 1.0:
                # Kinetically limited (Murphree E<1): exact vapor-flow tridiagonal.
                cols[c] = _murphree_vapor_column(
                    c, N, K, L, V, eff[c], f_stage, gas_f, y_gas)
            else:
                # Equilibrium species: standard liquid-mole-fraction tridiagonal.
                a = [0.0] + [L[j - 1] for j in range(1, N)]
                b = [-(L[j] + V[j] * K[j][c]) for j in range(N)]
                cc = [V[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
                d = [-f_stage[j].get(c, 0.0) for j in range(N)]
                cols[c] = _thomas(a, b, cc, d)
        sum_raw = [sum(max(cols[c][j], 0.0) for c in active) for j in range(N)]
        if min(sum_raw) <= 0.0:
            raise AbsorberError(
                f"{unit_id}: sum-rates produced a stage with no liquid "
                f"(iteration {it}); the column is infeasible as specified"
            )
        L_new = [max(L[j] * sum_raw[j], floor) for j in range(N)]
        L_prev, V_prev = list(L), list(V)
        L = [lo + omega * (ln - lo) for lo, ln in zip(L, L_new)]
        V = [max((L[j - 1] if j > 0 else 0.0) + F_below[j] - L[N - 1], floor)
             for j in range(N)]
        dF = max(
            max(abs(a - b) / max(a, b) for a, b in zip(L, L_prev)),
            max(abs(a - b) / max(a, b, floor * 10) for a, b in zip(V, V_prev)),
        )

        for j in range(N):
            row = {c: max(cols[c][j], _X_FLOOR) for c in active}
            x[j] = _normalized(row)
        # Vapor compositions: equilibrium y=Kx for E=1 species, the Murphree
        # relation y = E*K*x + (1-E)*y_below for E<1 (bottom-up so y_below is
        # this iterate's value; the gas feed is the bottom stage's y_below).
        for j in range(N - 1, -1, -1):
            below = y[j + 1] if j < N - 1 else y_gas
            y_raw = {
                c: (K[j][c] * x[j][c] if eff[c] >= 1.0
                    else eff[c] * K[j][c] * x[j][c]
                    + (1.0 - eff[c]) * below.get(c, 0.0))
                for c in active
            }
            y[j] = _normalized(y_raw)

        # -- energy balances: simultaneous Newton for the new T profile ------
        # Best-effort mid-iteration (the traffic may still be inconsistent);
        # tight closure is *required* at the converged point below.
        T_new, hL, hV, resid_rel, n_h = _newton_temperatures(
            pp, N, P_j, x, y, L, V, fh_stage, Q, T, T_bounds)
        n_prop += n_h
        dT = max(abs(tn - to) for tn, to in zip(T_new, T))
        T = T_new

        # Damp the liquid update if the temperature profile oscillates.
        if dT > dT_prev:
            worsened += 1
            if worsened >= 2:
                omega = max(0.25, 0.5 * omega)
                worsened = 0
        else:
            worsened = 0
        dT_prev = dT

        if dT <= _TOL_T and dF <= _TOL_FLOW and resid_rel <= _ENERGY_TOL:
            converged = True
            break

    if not converged:
        raise AbsorberError(
            f"{unit_id}: sum-rates (Burningham-Otto) iteration did not "
            f"converge in {max_iter} iterations (max |dT|={dT:.3g} K vs "
            f"{_TOL_T} K, max traffic change={dF:.3g} vs {_TOL_FLOW}, "
            f"energy residual={resid_rel:.3g} vs {_ENERGY_TOL}, "
            f"damping={omega}). Check the specifications and feed states, or "
            f"raise 'max_iter'"
        )
    return {
        "T": T, "L": L, "V": V, "x": x, "y": y, "hL": hL, "hV": hV,
        "iterations": it, "max_dT": dT, "max_dF_rel": dF, "damping": omega,
        "prop_calls": n_prop, "energy_residual_rel": resid_rel,
        "F_total": F_total,
    }


def _solve_ns_absorber(
    unit_id: str, pp: Any, active: list[str], N: int, P_j: list[float],
    feeds: list[tuple[int, dict[str, float], float]], Q: list[float],
    T: list[float], x: list[dict[str, float]], y: list[dict[str, float]],
    L: list[float], V: list[float], max_iter: int,
    t_bounds: tuple[float, float], cond_T: float | None = None,
    reb_T: float | None = None,
) -> dict[str, Any]:
    """Naphtali-Sandholm (1971) simultaneous-correction MESH for a gas absorber,
    stripper, or amine regenerator (Seader, Henley & Roper 3e §10.4), optionally
    with a **partial condenser** as the top stage and/or a **reboiler** as the
    bottom stage. ALL the MESH equations are solved at once by a damped Newton
    over per-component flow variables, so the method has none of the alternating
    sum-rates traffic instability that limit-cycles on reactive desorption.

    This is a simplification of the crude-tower NS in
    :mod:`caldyr.unitops.rigorous_column`: no side draws and no
    distillate-rate spec — every stage carries an ordinary energy balance (the
    top stage's vapour and the bottom stage's liquid are simply the products).

    When ``cond_T`` is given, stage 0 is a **partial (reflux) condenser**: its
    energy balance is replaced by a pinned-temperature spec ``T_0 = cond_T`` and
    its duty becomes a free output recovered from that balance after
    convergence. The overhead vapour leaves at ``cond_T`` (so condensables drop
    out as internal reflux to stage 1 — drying the product). When ``reb_T`` is
    given, the bottom stage is a **reboiler**: its energy balance is likewise
    replaced by ``T_{N-1} = reb_T`` and its (positive) duty boils up the
    stripping vapour internally — no open stripping steam, so the regenerator's
    water loop closes (only a tiny makeup for the dry overhead). The liquid feed
    enters stage 1 when there is a condenser (else stage 0); a gas feed, if any,
    enters the bottom stage.

    Variables (``N(2C+1)``): component liquid flows ``l_i,j``, component vapour
    flows ``v_i,j`` and stage temperature ``T_j``. Residuals: component material
    balance ``l_i,j + v_i,j - l_i,j-1 - v_i,j+1 - f_i,j``; per-component
    equilibrium ``v_i,j - K_i,j (V_j/L_j) l_i,j``; the stage energy balance
    (scaled). The Jacobian is finite-difference with per-stage property caching;
    a fractional-step cap + 2-norm line search keep the flows positive and the
    temperatures bounded, with a Levenberg-Marquardt escape if the Newton
    direction stalls. Returns the same profile dict as :func:`_solve_sum_rates`.
    """
    t_lo, t_hi = t_bounds
    C = len(active)
    span = 2 * C + 1
    nv = N * span
    P = np.array(P_j, dtype=float)
    Fcomp = np.zeros((N, C))
    Fh = np.zeros(N)
    for j, flows, fh in feeds:
        for c, val in flows.items():
            if c in active:
                Fcomp[j, active.index(c)] += val
        Fh[j] += fh
    Qa = np.array(Q, dtype=float)
    Ftot = max(float(Fcomp.sum()), 1.0)

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

    Lc = np.empty(N)
    Vc = np.empty(N)
    Kc = np.empty((N, C))
    hLc = np.empty(N)
    hVc = np.empty(N)

    for j in range(N):
        Lc[j], Vc[j], Kc[j], hLc[j], hVc[j] = stage_eval(j, lmat[j], vmat[j], Tv[j])
    e_ref = max(float(np.abs(hLc).max()), float(np.abs(hVc).max()), 1.0) * Ftot

    def residuals(lm, vm, Lq, Vq, Kq, hLq, hVq, Tq):
        R = np.empty(nv)
        for j in range(N):
            base = j * span
            in_lh = Lq[j - 1] * hLq[j - 1] if j > 0 else 0.0
            in_vh = Vq[j + 1] * hVq[j + 1] if j < N - 1 else 0.0
            for i in range(C):
                in_l = lm[j - 1, i] if j > 0 else 0.0
                in_v = vm[j + 1, i] if j < N - 1 else 0.0
                R[base + i] = (lm[j, i] + vm[j, i] - in_l - in_v
                               - Fcomp[j, i]) / Ftot
                R[base + C + i] = (vm[j, i]
                                   - Kq[j, i] * (Vq[j] / Lq[j]) * lm[j, i]) / Ftot
            if cond_T is not None and j == 0:
                # Partial condenser: pin the overhead temperature; its duty is a
                # free output recovered from the energy balance after solving.
                R[base + 2 * C] = (Tq[0] - cond_T) / _NS_TSPEC_SCALE
            elif reb_T is not None and j == N - 1:
                # Reboiler: pin the bottoms temperature; its duty (the internal
                # boilup) is the free output, likewise recovered after solving.
                R[base + 2 * C] = (Tq[N - 1] - reb_T) / _NS_TSPEC_SCALE
            else:
                out = Lq[j] * hLq[j] + Vq[j] * hVq[j]
                R[base + 2 * C] = (out - in_lh - in_vh - Fh[j] - Qa[j]) / e_ref
        return R

    def build_jac(R0):
        jac = np.empty((nv, nv))
        for j in range(N):
            for kind in range(span):
                col = j * span + kind
                lm, vm, tj, tvec = lmat, vmat, Tv[j], Tv
                if kind < C:
                    dk = max(_NS_DL_FRAC * lmat[j, kind], _NS_DL_MIN)
                    lm = lmat.copy()
                    lm[j, kind] += dk
                elif kind < 2 * C:
                    i = kind - C
                    dk = max(_NS_DL_FRAC * vmat[j, i], _NS_DL_MIN)
                    vm = vmat.copy()
                    vm[j, i] += dk
                else:
                    dk = _NS_DT
                    tj = Tv[j] + dk
                    tvec = Tv.copy()
                    tvec[j] = tj
                Lj, Vj, Kj, hLj, hVj = stage_eval(j, lm[j], vm[j], tj)
                Lc2, Vc2 = Lc.copy(), Vc.copy()
                Kc2, hLc2, hVc2 = Kc.copy(), hLc.copy(), hVc.copy()
                Lc2[j], Vc2[j], Kc2[j] = Lj, Vj, Kj
                hLc2[j], hVc2[j] = hLj, hVj
                jac[:, col] = (residuals(lm, vm, Lc2, Vc2, Kc2, hLc2, hVc2, tvec)
                               - R0) / dk
        return jac

    def cap(delta):
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
            Lt[j], Vt[j], Kt[j], hLt[j], hVt[j] = stage_eval(j, lm[j], vm[j], tt[j])
        Rt = residuals(lm, vm, Lt, Vt, Kt, hLt, hVt, tt)
        return float(Rt @ Rt), (lm, vm, tt, Lt, Vt, Kt, hLt, hVt, Rt)

    def line_search(delta, rss):
        frac = cap(delta)
        for _ in range(_NS_LS_MAX):
            ss, st = evals(delta, frac)
            if ss < (1.0 - 1e-4 * frac) * rss:
                return ss, st
            frac *= 0.5
        return None, None

    R0 = residuals(lmat, vmat, Lc, Vc, Kc, hLc, hVc, Tv)
    converged = False
    it = 0
    rnorm = float(np.max(np.abs(R0)))
    lam = _NS_LM_INIT
    for it in range(1, max_iter + 1):
        rnorm = float(np.max(np.abs(R0)))
        if rnorm < _NS_TOL:
            converged = True
            break
        rss = float(R0 @ R0)
        jac = build_jac(R0)
        try:
            delta_n = np.linalg.solve(jac, -R0)
            for _ref in range(_NS_REFINE):
                lin = -R0 - jac @ delta_n
                if (float(np.max(np.abs(lin)))
                        <= 1e-13 * float(np.max(np.abs(R0)) + 1e-30)):
                    break
                delta_n = delta_n + np.linalg.solve(jac, lin)
        except np.linalg.LinAlgError:
            delta_n = np.linalg.lstsq(jac, -R0, rcond=1e-8)[0]
        ss, state = line_search(delta_n, rss)
        if ss is None:                  # Newton stalled: Levenberg-Marquardt
            jtj = jac.T @ jac
            jtr = jac.T @ R0
            jdiag = np.maximum(np.diag(jtj), 1e-30)
            lam = max(lam, _NS_LM_INIT)
            for _ in range(_NS_LS_MAX):
                try:
                    delta_l = np.linalg.solve(jtj + lam * np.diag(jdiag), -jtr)
                except np.linalg.LinAlgError:
                    delta_l = -jtr / jdiag
                ss2, st2 = evals(delta_l, cap(delta_l))
                if ss2 < rss:
                    ss, state = ss2, st2
                    break
                lam = min(lam * _NS_LM_UP, _NS_LM_MAX)
        if state is None:
            break
        lmat, vmat, Tv, R0 = state[0], state[1], state[2], state[8]
        Lc[:], Vc[:] = state[3], state[4]
        Kc[:], hLc[:], hVc[:] = state[5], state[6], state[7]
        lam = max(lam / _NS_LM_DOWN, _NS_LM_MIN)

    if not converged:
        raise AbsorberError(
            f"{unit_id}: Naphtali-Sandholm simultaneous-correction did not "
            f"converge in {max_iter} iterations (scaled residual {rnorm:.3g} vs "
            f"{_NS_TOL}); check the specifications, or raise 'max_iter'"
        )
    x_out = [_normalized({active[i]: max(float(lmat[j, i]), _X_FLOOR)
                          for i in range(C)}) for j in range(N)]
    y_out = [_normalized({active[i]: max(float(vmat[j, i]), _X_FLOOR)
                          for i in range(C)}) for j in range(N)]
    # Partial-condenser duty (heat ADDED at stage 0 — negative, heat removed):
    # the stage-0 energy balance with no feed and no liquid from above.
    cond_duty = (float(Lc[0] * hLc[0] + Vc[0] * hVc[0] - Vc[1] * hVc[1])
                 if cond_T is not None else None)
    # Reboiler duty (heat ADDED at the bottom stage — positive): its energy
    # balance with no vapour from below.
    reb_duty = (float(Lc[N - 1] * hLc[N - 1] + Vc[N - 1] * hVc[N - 1]
                      - Lc[N - 2] * hLc[N - 2] - Fh[N - 1])
                if reb_T is not None else None)
    return {
        "T": [float(t) for t in Tv], "L": [float(v) for v in Lc],
        "V": [float(v) for v in Vc], "x": x_out, "y": y_out,
        "hL": [float(h) for h in hLc], "hV": [float(h) for h in hVc],
        "iterations": it, "max_dT": 0.0, "max_dF_rel": 0.0, "damping": 1.0,
        "prop_calls": n_prop, "energy_residual_rel": rnorm, "F_total": Ftot,
        "condenser_duty": cond_duty, "reboiler_duty": reb_duty,
    }


def _murphree_vapor_column(
    c: str, N: int, K: list[dict[str, float]], L: list[float], V: list[float],
    e: float, f_stage: list[dict[str, float]], gas_f: dict[str, float],
    y_gas: dict[str, float],
) -> list[float]:
    """Component balance of one kinetically-limited species (Murphree E<1),
    solved EXACTLY in its vapor component flows ``v_j`` then returned as the
    liquid mole fractions the sum-rates loop consumes (``l_j / L_j``).

    Writing the tray relation ``v_j = E*S_j*l_j + (1-E)*(V_j/V_(j+1))*v_(j+1)``
    (stripping factor ``S_j = V_j*K_j/L_j``) and substituting the liquid flows
    ``l_j = [v_j - V_j*(1-E)*y_below]/(E*S_j)`` into the stage component balance
    keeps the (1-E) coupling between adjacent vapors *tridiagonal* in ``v`` — so
    the Murphree solution is obtained in one Thomas sweep, with no lag (unlike an
    effective-K successive substitution, which converges only linearly). The
    vapor entering below the bottom stage is the gas feed (``gas_f`` flows /
    ``y_gas`` composition). At E=1 this reduces to the standard liquid-flow
    balance (gas feed as the bottom entering vapor), so it is consistent with the
    equilibrium path.
    """
    S = [V[j] * K[j][c] / L[j] for j in range(N)]   # stripping factor (K>0)
    a = [0.0] * N
    b = [0.0] * N
    cc = [0.0] * N
    d = [0.0] * N
    for j in range(N):
        # Liquid-side feed of this species at stage j (the gas feed enters as the
        # bottom stage's *entering vapor*, not as a stage source).
        f_liq = f_stage[j].get(c, 0.0) - (gas_f[c] if j == N - 1 else 0.0)
        b[j] = -1.0 / (e * S[j]) - 1.0
        d[j] = -f_liq
        if j >= 1:
            a[j] = 1.0 / (e * S[j - 1])
            b[j] += -(1.0 - e) * (V[j - 1] / V[j]) / (e * S[j - 1])
        if j <= N - 2:
            cc[j] = 1.0 + (1.0 - e) * (V[j] / V[j + 1]) / (e * S[j])
        if j == N - 1:
            d[j] += -gas_f[c] - V[N - 1] * (1.0 - e) * y_gas[c] / (e * S[N - 1])
    v = _thomas(a, b, cc, d)                          # vapor component flows v_j
    cols = [0.0] * N
    for j in range(N):
        y_below = (v[j + 1] / V[j + 1]) if j < N - 1 else y_gas[c]
        l_j = (v[j] - V[j] * (1.0 - e) * y_below) / (e * S[j])
        cols[j] = l_j / L[j]                          # liquid mole fraction
    return cols


def _newton_temperatures(
    pp: Any, N: int, P_j: list[float],
    x: list[dict[str, float]], y: list[dict[str, float]],
    L: list[float], V: list[float],
    fh_stage: list[float], Q: list[float],
    T0: list[float], T_bounds: tuple[float, float],
) -> tuple[list[float], list[float], list[float], float, int]:
    """Solve the N stage energy balances for the temperature profile by
    Newton's method with a tridiagonal Jacobian (Burningham-Otto; Seader 3e
    eqs. 10-65..10-67), steps clamped and temperatures bounded to the
    physically admissible window. Best-effort: returns the last iterate with
    its residual rather than raising — the *outer* sum-rates loop requires a
    tight residual before declaring convergence. Returns
    ``(T, hL, hV, residual_rel, n_prop_calls)`` with the phase enthalpies at
    the returned temperatures."""
    T = list(T0)
    t_lo, t_hi = T_bounds
    n_prop = 0
    best: tuple[float, list[float], list[float], list[float]] | None = None
    for _ in range(_NEWTON_MAX):
        hL = [pp.enthalpy_liquid(T[j], P_j[j], x[j]) for j in range(N)]
        hV = [pp.enthalpy_vapor(T[j], P_j[j], y[j]) for j in range(N)]
        cpL = [(pp.enthalpy_liquid(T[j] + _DH_DT, P_j[j], x[j]) - hL[j]) / _DH_DT
               for j in range(N)]
        cpV = [(pp.enthalpy_vapor(T[j] + _DH_DT, P_j[j], y[j]) - hV[j]) / _DH_DT
               for j in range(N)]
        n_prop += 4 * N

        H = [0.0] * N
        scale = [1.0] * N
        for j in range(N):
            inflow = ((L[j - 1] * hL[j - 1] if j > 0 else 0.0)
                      + (V[j + 1] * hV[j + 1] if j < N - 1 else 0.0)
                      + fh_stage[j] + Q[j])
            outflow = L[j] * hL[j] + V[j] * hV[j]
            H[j] = inflow - outflow
            scale[j] = max(abs(L[j] * hL[j]), abs(V[j] * hV[j]),
                           abs(fh_stage[j]), 1.0)
        resid_rel = max(abs(H[j]) / scale[j] for j in range(N))
        if best is None or resid_rel < best[0]:
            best = (resid_rel, list(T), hL, hV)
        if resid_rel <= _NEWTON_TOL:
            return T, hL, hV, resid_rel, n_prop

        a = [0.0] + [L[j - 1] * cpL[j - 1] for j in range(1, N)]
        b = [-(L[j] * cpL[j] + V[j] * cpV[j]) for j in range(N)]
        c = [V[j + 1] * cpV[j + 1] for j in range(N - 1)] + [0.0]
        # b_j < 0 always (Cp > 0, flows > 0), so Thomas is well-posed.
        dT = _thomas(a, b, c, [-h for h in H])
        T = [min(max(tj + max(min(d, _NEWTON_STEP_MAX), -_NEWTON_STEP_MAX),
                     t_lo), t_hi)
             for tj, d in zip(T, dT)]
    assert best is not None
    resid_rel, T, hL, hV = best
    return T, hL, hV, resid_rel, n_prop


def _active_components(comps: list[str], *zs: dict[str, float]) -> list[str]:
    active = [c for c in comps if any(z.get(c, 0.0) > 0.0 for z in zs)]
    if len(active) < 2:
        raise AbsorberError(
            f"absorption needs at least two components with feed flow; got "
            f"{active}"
        )
    return active


def _pad(row: dict[str, float], comps: list[str]) -> dict[str, float]:
    return {c: row.get(c, 0.0) for c in comps}


@register("Absorber")
class Absorber(UnitOp):
    """Gas absorber / stripper: MESH without condenser or reboiler, solved by
    the sum-rates (Burningham-Otto) method. See the module docstring.

    Ports: ``gas_in`` (bottom-stage vapor feed), ``liquid_in`` (top-stage
    solvent feed), ``vapor_out`` (top), ``liquid_out`` (bottom). To use the
    unit as a **stripper**, feed the volatile-rich liquid on ``liquid_in``
    and the stripping gas on ``gas_in`` — same physics, mirrored direction.

    Params (JSON-friendly scalars; ``.flow`` round-trips):
      * ``n_stages`` — theoretical stages (no condenser/reboiler). Required,
        >= 1.
      * ``P`` — top-stage pressure, Pa (default: gas feed pressure).
      * ``dP_stage`` — linear pressure rise per stage going down, Pa
        (default 0).
      * ``method`` — MESH solver: ``"sum_rates"`` (default, Burningham-Otto) or
        ``"naphtali_sandholm"``. Sum-rates is fast and converges most absorbers
        and strippers; the simultaneous-correction Naphtali-Sandholm method is
        the robust fallback for cases where the alternating sum-rates iteration
        limit-cycles — notably an **amine regenerator** (reactive desorption +
        condensing stripping steam). NS does not support ``murphree`` yet.
      * ``murphree`` — optional Murphree vapor stage efficiency: a scalar
        (every component) or a ``{component: E}`` dict (others default to 1.0),
        each in (0, 1]. Use it for kinetically limited mass transfer — e.g.
        amine sweetening, where CO2 reacts slowly with a tertiary amine
        (book §15.3: E_CO2 ~ 0.15, E_H2S ~ 0.8) so an equilibrium stage would
        over-predict CO2 removal. Default: equilibrium stages. (sum_rates only.)
      * ``condenser_T`` — optional **partial (reflux) condenser** temperature, K
        (Naphtali-Sandholm only, ``n_stages >= 2``). Stage 0 becomes a partial
        condenser pinned to this temperature; condensables drop out as internal
        reflux to stage 1, so the overhead product leaves dried — the natural
        refinement for an amine regenerator's wet acid-gas overhead. The reported
        ``design['condenser_duty']`` (negative) closes the overall energy balance.
        The liquid feed then enters the top tray (stage 1).
      * ``reboiler_T`` — optional **reboiler** temperature, K (Naphtali-Sandholm
        only). The bottom stage boils up the stripping vapour internally (pinned
        to this temperature, ``design['reboiler_duty']`` positive), so ``gas_in``
        is then optional — an internally-boiled stripper with no open steam.
        Robust on well-conditioned systems; **combining it with ``condenser_T``,
        or using it on the reactive amine regenerator, is numerically stiff** in
        the FD-Jacobian endgame (a robust reboiled-NS amine regenerator is
        tracked follow-up work), so the validated amine path uses open stripping
        steam plus the reflux condenser.
      * ``max_iter`` — iteration cap (default 120).
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
            Port("gas_in", "inlet"),
            Port("liquid_in", "inlet"),
            Port("vapor_out", "outlet"),
            Port("liquid_out", "outlet"),
        ]

    def _read_params(self, P_in: float) -> tuple[int, float, float, int]:
        try:
            n_stages = int(self.params["n_stages"])
        except (KeyError, TypeError, ValueError) as exc:
            raise AbsorberError(
                f"Absorber {self.id!r}: integer 'n_stages' is required "
                f"(got {self.params.get('n_stages')!r})"
            ) from exc
        if n_stages < 1:
            raise AbsorberError(
                f"Absorber {self.id!r}: n_stages={n_stages} must be >= 1"
            )
        if n_stages > 1000:
            raise AbsorberError(
                f"Absorber {self.id!r}: n_stages={n_stages} exceeds the "
                f"1000-stage limit"
            )
        P = float(self.params.get("P") or P_in)
        if P <= 0.0:
            raise AbsorberError(f"Absorber {self.id!r}: P={P} Pa must be > 0")
        dP = float(self.params.get("dP_stage", 0.0))
        if dP < 0.0:
            raise AbsorberError(
                f"Absorber {self.id!r}: dP_stage={dP} Pa must be >= 0 "
                f"(pressure increases going down)"
            )
        max_iter = int(self.params.get("max_iter", _MAX_ITER))
        return n_stages, P, dP, max_iter

    def _murphree_map(self, active: list[str]) -> dict[str, float] | None:
        """Parse the optional ``murphree`` param into a per-component efficiency
        map. Accepts a scalar (applied to every component) or a
        ``{component: E}`` dict (unspecified components default to 1.0). Returns
        ``None`` for the equilibrium case (no param / all 1.0)."""
        spec = self.params.get("murphree")
        if spec is None:
            return None
        if isinstance(spec, dict):
            eff = {c: float(spec.get(c, 1.0)) for c in active}
        else:
            eff = {c: float(spec) for c in active}
        for c, e in eff.items():
            if not 0.0 < e <= 1.0:
                raise AbsorberError(
                    f"Absorber {self.id!r}: Murphree efficiency for {c!r} is "
                    f"{e} — it must lie in (0, 1]"
                )
        return eff if any(abs(e - 1.0) > 1e-12 for e in eff.values()) else None

    def _condenser_T(self, method: str, N: int) -> float | None:
        """Parse the optional ``condenser_T`` param (the reflux/overhead
        temperature, K, of a partial condenser on the top stage). Returns
        ``None`` when absent. The partial condenser dries the overhead product
        (condensables drop out as internal reflux) — the natural refinement for
        an amine regenerator's wet acid-gas product. Requires the
        Naphtali-Sandholm method and ``n_stages >= 2`` (stage 0 is the
        condenser, the top tray is stage 1)."""
        spec = self.params.get("condenser_T")
        if spec is None:
            return None
        cond_T = float(spec)
        if method != "naphtali_sandholm":
            raise AbsorberError(
                f"Absorber {self.id!r}: 'condenser_T' (a partial condenser) is "
                f"only supported with method='naphtali_sandholm'"
            )
        if N < 2:
            raise AbsorberError(
                f"Absorber {self.id!r}: a partial condenser needs n_stages >= 2 "
                f"(stage 0 is the condenser, the top tray is stage 1); got {N}"
            )
        if cond_T <= 0.0:
            raise AbsorberError(
                f"Absorber {self.id!r}: condenser_T={cond_T} K must be > 0"
            )
        return cond_T

    def _reboiler_T(self, method: str, N: int, has_condenser: bool) -> float | None:
        """Parse the optional ``reboiler_T`` param (the bottoms temperature, K,
        of a reboiler on the bottom stage). Returns ``None`` when absent. The
        reboiler boils up the stripping vapour internally — so an amine
        regenerator needs no open stripping steam and its water loop closes
        (only a tiny makeup for the dry overhead). Requires the
        Naphtali-Sandholm method; with a condenser too the column needs
        ``n_stages >= 3`` (condenser, >=1 tray, reboiler)."""
        spec = self.params.get("reboiler_T")
        if spec is None:
            return None
        reb_T = float(spec)
        if method != "naphtali_sandholm":
            raise AbsorberError(
                f"Absorber {self.id!r}: 'reboiler_T' (a reboiler) is only "
                f"supported with method='naphtali_sandholm'"
            )
        n_min = 3 if has_condenser else 2
        if N < n_min:
            raise AbsorberError(
                f"Absorber {self.id!r}: a reboiler needs n_stages >= {n_min} "
                f"(the bottom stage is the reboiler); got {N}"
            )
        if reb_T <= 0.0:
            raise AbsorberError(
                f"Absorber {self.id!r}: reboiler_T={reb_T} K must be > 0"
            )
        return reb_T

    def _reboiler_duty(self, method: str, N: int, reb_T: float | None,
                       active: list[str]) -> float | None:
        """Parse the optional ``reboiler_duty`` param (W, heat added at the bottom
        stage) — the alternative reboiler spec to ``reboiler_T``: instead of
        pinning the bottoms temperature it sets the reboiler heat input, which
        keeps the bottom-stage energy balance intact (better conditioned). This
        is the **robust amine-regenerator** spec: paired with ``condenser_T`` it
        gives a dried acid-gas product AND a closed water loop (no open steam), and
        the solve warm-starts from an open-steam proxy so the otherwise-stiff
        condenser+reboiler combination converges. Requires Naphtali-Sandholm and
        water in the system (for the steam proxy)."""
        spec = self.params.get("reboiler_duty")
        if spec is None:
            return None
        reb_duty = float(spec)
        if method != "naphtali_sandholm":
            raise AbsorberError(
                f"Absorber {self.id!r}: 'reboiler_duty' is only supported with "
                f"method='naphtali_sandholm'"
            )
        if reb_T is not None:
            raise AbsorberError(
                f"Absorber {self.id!r}: specify only one of 'reboiler_T' "
                f"(pin the bottoms temperature) or 'reboiler_duty' (set the heat)"
            )
        if N < 2:
            raise AbsorberError(
                f"Absorber {self.id!r}: a reboiler needs n_stages >= 2; got {N}"
            )
        if reb_duty <= 0.0:
            raise AbsorberError(
                f"Absorber {self.id!r}: reboiler_duty={reb_duty} W must be > 0 "
                f"(it boils up the stripping vapour)"
            )
        if "water" not in active:
            raise AbsorberError(
                f"Absorber {self.id!r}: 'reboiler_duty' needs water in the system "
                f"(the warm-start proxy boils steam); got components {active}"
            )
        return reb_duty

    def _solve_reboiled_ns(self, uid, pp, active, N, P_j, feeds, Q, reb_duty,
                           profiles, max_iter, t_bounds, cond_T, P, T_l):
        """Solve an internally-boiled (reboiler-duty) NS column via a two-stage
        continuation: (1) a robust **open-steam proxy** — the same column with a
        synthetic steam feed on the bottom stage (sized from the duty) and no
        reboiler — converges from a cold start; (2) the **real** column (the
        reboiler duty on the bottom stage, no steam) re-solves from the proxy's
        in-basin profile in a handful of Newton steps."""
        # Proxy steam ~ the boilup the duty produces (duty / water latent heat),
        # bounded by the liquid feed so the proxy column is well-posed.
        f_liq = feeds[0][1]
        F_liq = sum(f_liq.values()) or 1.0
        z_l = {c: f_liq[c] / F_liq for c in active}
        steam = min(max(reb_duty / 45.0e3, 0.1 * F_liq), F_liq)
        T_steam = min(t_bounds[1] - 5.0, max(T_l, 393.0))
        h_steam = pp.enthalpy(T_steam, P, {"water": 1.0})
        # Proxy seed: cool condenser top -> hot steam bottom, with steam traffic.
        Tp, xp, yp, Lp, Vp = self._initial_profiles(
            N, active, {"water": 1.0}, z_l, T_steam, T_l, steam, F_liq, cond_T, None)
        proxy_feeds = list(feeds) + [(N - 1, {"water": steam}, steam * h_steam)]
        proxy = _solve_ns_absorber(
            uid, pp, active, N, P_j, proxy_feeds, [0.0] * N,
            Tp, xp, yp, Lp, Vp, max_iter, t_bounds, cond_T, None)
        # Re-solve the real internally-boiled column from the proxy profile.
        return _solve_ns_absorber(
            uid, pp, active, N, P_j, feeds, Q,
            proxy["T"], proxy["x"], proxy["y"], proxy["L"], proxy["V"],
            max_iter, t_bounds, cond_T, None)

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        gas = inlets.get("gas_in")
        liq = inlets.get("liquid_in")
        if liq is None or not liq.molar_flow:
            raise AbsorberError(
                f"Absorber {self.id!r}: missing or empty inlet on 'liquid_in'"
            )
        method = str(self.params.get("method", "sum_rates"))
        if method not in ("sum_rates", "naphtali_sandholm"):
            raise AbsorberError(
                f"Absorber {self.id!r}: unknown method {method!r}; expected "
                f"'sum_rates' or 'naphtali_sandholm'"
            )
        # A reboiler supplies the stripping vapour internally, so 'gas_in' is
        # then optional (a regenerator with no open steam).
        reboiled = (self.params.get("reboiler_T") is not None
                    or self.params.get("reboiler_duty") is not None)
        gas_present = gas is not None and bool(gas.molar_flow)
        if not gas_present and not reboiled:
            raise AbsorberError(
                f"Absorber {self.id!r}: missing or empty inlet on 'gas_in' (or set "
                f"'reboiler_T'/'reboiler_duty' for an internally-boiled stripper)"
            )

        T_l, P_l, F_l = liq.require_state()
        z_l = liq.normalized_z()
        comps = list(liq.components)
        gas_key: tuple
        if gas is not None and gas_present:
            T_g, P_g, F_g = gas.require_state()
            z_g = gas.normalized_z()
            comps = list(gas.components)
            h_g = gas.H if gas.H is not None else pp.enthalpy(T_g, P_g, z_g)
            gas_key = (T_g, P_g, F_g, tuple(sorted(z_g.items())), gas.H)
        else:
            T_g, P_g, F_g, z_g, h_g = T_l, P_l, 0.0, {}, 0.0
            gas_key = (None,)
        n_stages, P, dP, max_iter = self._read_params(P_g)

        key = (repr(sorted(self.params.items())), gas_key,
               T_l, P_l, F_l, tuple(sorted(z_l.items())), liq.H)
        if key == self._cache_key and self._cache_out is not None:
            assert self._cache_design is not None
            self.design = _copy_design(self._cache_design)
            return _copy_out(self._cache_out)

        active = _active_components(comps, z_g, z_l)
        N = n_stages
        P_j = [P + j * dP for j in range(N)]
        h_l = liq.H if liq.H is not None else pp.enthalpy(T_l, P_l, z_l)
        f_gas = {c: F_g * z_g.get(c, 0.0) for c in active}
        f_liq = {c: F_l * z_l.get(c, 0.0) for c in active}
        f_tot = {c: f_gas[c] + f_liq[c] for c in active}

        murphree = self._murphree_map(active)
        cond_T = self._condenser_T(method, N)
        reb_T = self._reboiler_T(method, N, cond_T is not None)
        reb_duty = self._reboiler_duty(method, N, reb_T, active)
        # With a partial condenser, stage 0 is the condenser and the liquid feed
        # enters the top tray (stage 1); otherwise it enters stage 0.
        liq_stage = 1 if cond_T is not None else 0
        feeds = [(liq_stage, f_liq, F_l * h_l)]
        if gas_present:
            feeds.append((N - 1, f_gas, F_g * h_g))
        Q = [0.0] * N
        if reb_duty is not None:
            Q[N - 1] = reb_duty                # reboiler heat on the bottom stage
        T0, x0, y0, L0, V0 = self._initial_profiles(
            N, active, z_g, z_l, T_g, T_l, F_g, F_l, cond_T, reb_T)
        t_lo = min(T_g, T_l) - 60.0
        t_hi = max(T_g, T_l) + 90.0
        if cond_T is not None:
            t_lo = min(t_lo, cond_T - 30.0)
        if reb_T is not None:
            t_hi = max(t_hi, reb_T + 30.0)
        if reb_duty is not None:
            t_hi = max(t_hi, T_l + 60.0)
        t_bounds = (t_lo, t_hi)
        uid = f"Absorber {self.id!r}"
        if method == "naphtali_sandholm":
            if murphree is not None:
                raise AbsorberError(
                    f"Absorber {self.id!r}: the 'naphtali_sandholm' method does "
                    f"not support a Murphree efficiency yet — use 'sum_rates'"
                )
            if reb_duty is not None:
                # A reboiler duty + condenser is FD-Jacobian-stiff from a cold
                # start, so warm-start from a robust open-steam proxy (a synthetic
                # bottom steam feed sized from the duty) then re-solve the real
                # internally-boiled column from that in-basin profile.
                res = self._solve_reboiled_ns(
                    uid, pp, active, N, P_j, feeds, Q, reb_duty,
                    (T0, x0, y0, L0, V0), max_iter, t_bounds, cond_T, P, T_l)
            else:
                res = _solve_ns_absorber(
                    uid, pp, active, N, P_j, feeds, Q, T0, x0, y0, L0, V0,
                    max_iter, t_bounds, cond_T, reb_T)
        else:
            res = _solve_sum_rates(
                uid, pp, active, N, P_j, feeds, Q, T0, x0, y0, L0, V0, max_iter,
                T_bounds=t_bounds, murphree=murphree, bottom_vapor_flows=f_gas)

        q_cond = res.get("condenser_duty") if cond_T is not None else None
        q_reb = res.get("reboiler_duty") if reb_T is not None else reb_duty
        q_added = (q_cond or 0.0) + (q_reb or 0.0)
        out, design = _build_products(
            self.id, pp, comps, active, f_tot, res, P_j,
            feeds_enthalpy=F_g * h_g + F_l * h_l, q_added=q_added,
            vapor_port="vapor_out", liquid_port="liquid_out")
        if q_cond is not None:
            design["condenser_T"] = cond_T
            design["condenser_duty"] = q_cond
        if q_reb is not None:
            if reb_T is not None:
                design["reboiler_T"] = reb_T
            design["reboiler_duty"] = q_reb
        design["N"] = float(N)            # ideal stages (incl. condenser if any)
        design["absorbed"] = {
            c: (1.0 - design["vapor_flows"][c] / f_gas[c]) if f_gas[c] > 0.0
            else 0.0
            for c in active
        }
        self.design = design

        self._warm = {"n_stages": N, "active": list(active),
                      "z": dict(_normalized(f_tot)),
                      "T": list(res["T"]), "x": [dict(r) for r in res["x"]],
                      "y": [dict(r) for r in res["y"]],
                      "L": list(res["L"]), "V": list(res["V"])}
        self._cache_key = key
        self._cache_out = _copy_out(out)
        self._cache_design = _copy_design(design)
        return out

    def _initial_profiles(self, N, active, z_g, z_l, T_g, T_l, F_g, F_l,
                          cond_T=None, reb_T=None):
        w = self._warm
        z_tot = _normalized({c: F_g * z_g.get(c, 0.0) + F_l * z_l.get(c, 0.0)
                             for c in active})
        if (w is not None and w["n_stages"] == N and w["active"] == active
                and max(abs(w["z"].get(c, 0.0) - z_tot[c]) for c in active)
                < _WARM_Z_TOL):
            return w["T"], w["x"], w["y"], w["L"], w["V"]
        # End temperatures anchor the linear seed (cool top with a condenser,
        # hot bottom with a reboiler) so the Newton starts on the right slope.
        t_top = cond_T if cond_T is not None else T_l
        t_bot = reb_T if reb_T is not None else T_g
        T0 = [t_top + (t_bot - t_top) * (j / max(N - 1, 1)) for j in range(N)]
        x0 = [_normalized({c: max(z_l.get(c, 0.0), _X_FLOOR) for c in active})
              for _ in range(N)]
        # With no gas feed (internally boiled), seed the vapour from the liquid
        # feed (the boilup rises off the descending liquid).
        z_v = z_g if F_g > 0.0 else z_l
        y0 = [_normalized({c: max(z_v.get(c, 0.0), _X_FLOOR) for c in active})
              for _ in range(N)]
        boilup = F_l if reb_T is not None else F_g
        L0, V0 = [F_l] * N, [max(boilup, F_g)] * N
        if cond_T is not None:
            # Seed the partial-condenser stage cool: condensables reflux, so the
            # overhead vapour leaving is a fraction of the rising traffic.
            T0[0] = cond_T
            V0[0] = 0.5 * max(boilup, F_g)
        if reb_T is not None:
            T0[N - 1] = reb_T
        return T0, x0, y0, L0, V0


@register("ReboiledAbsorber")
class ReboiledAbsorber(UnitOp):
    """Reboiled absorber (stripping tower): an absorber whose bottom stage is
    a reboiler — the stripping vapor is generated internally by the reboiler
    heat input instead of being fed (Hameed 2025 sec. 9.3.5). Solved with a
    Wang-Henke bubble-point inner loop driven by the overhead vapor rate;
    the reboiler duty closes the overall energy balance exactly (see the
    module docstring for why bubble-point rather than sum-rates here).

    Ports: ``feed`` (liquid feed near the top), ``vapor_out`` (top),
    ``bottoms`` (reboiler liquid), ``reboiler_duty`` (energy).

    Params:
      * ``n_stages`` — theoretical stages **including the reboiler** as stage
        n_stages. Required, >= 2. (A HYSYS reboiled absorber with "N stages"
        plus its reboiler corresponds to ``n_stages = N + 1`` here.)
      * ``feed_stage`` — feed stage, 1..n_stages-1 (default 1: top).
      * exactly one of ``vapor_rate`` (overhead vapor, mol/s),
        ``boilup_ratio`` (V_reboiler / bottoms) or ``reboiler_duty`` (W).
        The latter two are met by an outer Brent iteration on the overhead
        vapor rate (both are strictly monotone in it).
      * ``P`` — top-stage pressure, Pa (default: feed pressure).
      * ``P_bottom`` — reboiler pressure, Pa (linear stage profile), or
        ``dP_stage`` — Pa per stage going down (not both; default uniform).
      * ``max_iter`` — inner-iteration cap (default 120).
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
            Port("feed", "inlet"),
            Port("vapor_out", "outlet"),
            Port("bottoms", "outlet"),
            Port("reboiler_duty", "outlet", "energy"),
        ]

    def _read_params(self, P_in: float):
        try:
            n_stages = int(self.params["n_stages"])
        except (KeyError, TypeError, ValueError) as exc:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: integer 'n_stages' is "
                f"required (got {self.params.get('n_stages')!r})"
            ) from exc
        if n_stages < 2:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: n_stages={n_stages} must be "
                f">= 2 — the count includes the reboiler as stage n_stages"
            )
        if n_stages > 1000:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: n_stages={n_stages} exceeds "
                f"the 1000-stage limit"
            )
        feed_stage = int(self.params.get("feed_stage", 1))
        if not 1 <= feed_stage <= n_stages - 1:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: feed_stage={feed_stage} is "
                f"out of range — it must lie between stage 1 (top) and stage "
                f"{n_stages - 1} (above the reboiler) for n_stages={n_stages}"
            )

        specs = {k: self.params.get(k)
                 for k in ("reboiler_duty", "boilup_ratio", "vapor_rate")}
        given = {k: v for k, v in specs.items() if v is not None}
        if len(given) != 1:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: specify exactly one of "
                f"'reboiler_duty' (W), 'boilup_ratio' or 'vapor_rate' "
                f"(mol/s); got {given or specs}"
            )
        spec_name, spec_value = next(iter(given.items()))
        spec_value = float(spec_value)
        if spec_value <= 0.0:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: {spec_name}={spec_value} "
                f"must be > 0"
            )

        P = float(self.params.get("P") or P_in)
        if P <= 0.0:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: P={P} Pa must be > 0")
        if (self.params.get("P_bottom") is not None
                and self.params.get("dP_stage") is not None):
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: give 'P_bottom' or "
                f"'dP_stage', not both"
            )
        if self.params.get("P_bottom") is not None:
            P_bot = float(self.params["P_bottom"])
            if P_bot < P:
                raise AbsorberError(
                    f"ReboiledAbsorber {self.id!r}: P_bottom={P_bot} Pa must "
                    f"be >= the top pressure P={P} Pa"
                )
            dP = (P_bot - P) / max(n_stages - 1, 1)
        else:
            dP = float(self.params.get("dP_stage", 0.0))
            if dP < 0.0:
                raise AbsorberError(
                    f"ReboiledAbsorber {self.id!r}: dP_stage={dP} Pa must "
                    f"be >= 0"
                )
        max_iter = int(self.params.get("max_iter", _MAX_ITER))
        return n_stages, feed_stage, spec_name, spec_value, P, dP, max_iter

    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        feed = inlets.get("feed")
        if feed is None or not feed.molar_flow:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: missing or empty inlet on "
                f"'feed'"
            )
        T_in, P_in, F = feed.require_state()
        z = feed.normalized_z()
        comps = list(feed.components)
        n_stages, feed_stage, spec_name, spec_value, P, dP, max_iter = \
            self._read_params(P_in)

        key = (repr(sorted(self.params.items())),
               T_in, P_in, F, tuple(sorted(z.items())), feed.H)
        if key == self._cache_key and self._cache_out is not None:
            assert self._cache_design is not None
            self.design = _copy_design(self._cache_design)
            return _copy_out(self._cache_out)

        active = _active_components(comps, z)
        N = n_stages
        f0 = feed_stage - 1
        P_j = [P + j * dP for j in range(N)]
        h_F = feed.H if feed.H is not None else pp.enthalpy(T_in, P_in, z)
        f = {c: F * z[c] for c in active}
        z_act = {c: z[c] for c in active}

        if spec_name == "vapor_rate" and not 0.0 < spec_value < F:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: vapor_rate={spec_value:.6g} "
                f"mol/s must lie strictly between 0 and the feed rate "
                f"F={F:.6g} mol/s"
            )

        def run(v_top: float) -> dict[str, Any]:
            x0 = self._initial_x(pp, N, active, z_act, P_j[f0], f, v_top)
            res = _solve_reboiled_wh(
                f"ReboiledAbsorber {self.id!r}", pp, active, N, P_j, f0,
                f, F, h_F, v_top, x0, max_iter)
            # Remember the profile: the next spec-iteration (or recycle
            # sweep) warm-starts from it.
            self._warm = {"n_stages": N, "active": list(active),
                          "z": dict(z_act),
                          "x": [dict(r) for r in res["x"]]}
            return res

        if spec_name == "vapor_rate":
            res = run(spec_value)
        else:
            res = self._meet_spec(run, spec_name, spec_value, F, h_F)

        # The reboiler duty closes the overall energy balance exactly
        # (enthalpies are absolute / formation-inclusive):
        #   F h_F + Q_reb = V_top h_V,top + B h_B
        q_reb = (res["V"][0] * res["hV"][0]
                 + res["L"][N - 1] * res["hL"][N - 1] - F * h_F)
        out, design = _build_products(
            self.id, pp, comps, active, f, res, P_j,
            feeds_enthalpy=F * h_F, q_added=q_reb,
            vapor_port="vapor_out", liquid_port="bottoms")
        # Sizing keys: "N" is the ideal-stage count excluding the reboiler
        # (it has no tray); the reboiler is costed from Q_reboiler/T_bottom.
        design["N"] = float(N - 1)
        design["Q_reboiler"] = q_reb
        design["boilup_ratio"] = res["V"][N - 1] / res["L"][N - 1]
        design["feed_stage"] = feed_stage
        self.design = design

        out["reboiler_duty"] = EnergyStream(
            id=f"{self.id}.reboiler_duty", duty=q_reb)
        self._cache_key = key
        self._cache_out = _copy_out(out)
        self._cache_design = _copy_design(design)
        return out

    def _meet_spec(self, run, spec_name: str, target: float, F: float,
                   h_F: float) -> dict[str, Any]:
        """Outer Brent iteration on the overhead vapor rate to meet a
        boilup-ratio or reboiler-duty spec (both are strictly increasing in
        V_top: more overhead vapor means more boilup and more reboiler heat)."""
        cache: dict[float, dict[str, Any]] = {}

        def g(v_top: float) -> float:
            res = cache.get(v_top)
            if res is None:
                res = run(v_top)
                cache[v_top] = res
            if spec_name == "boilup_ratio":
                return res["V"][-1] / res["L"][-1] - target
            q = (res["V"][0] * res["hV"][0]
                 + res["L"][-1] * res["hL"][-1] - F * h_F)
            return q - target

        lo, hi = 0.02 * F, 0.98 * F
        g_lo, g_hi = g(lo), g(hi)
        for _ in range(8):
            if g_lo > 0.0 and lo > 1e-4 * F:
                hi, g_hi = lo, g_lo
                lo = max(lo / 4.0, 1e-4 * F)
                g_lo = g(lo)
            elif g_hi < 0.0 and hi < (1.0 - 1e-4) * F:
                lo, g_lo = hi, g_hi
                hi = min(hi * 1.02, (1.0 - 1e-4) * F)
                g_hi = g(hi)
            else:
                break
        if g_lo > 0.0 or g_hi < 0.0:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: no overhead vapor rate in "
                f"({lo:.4g}, {hi:.4g}) mol/s meets {spec_name}={target:.6g} "
                f"(spec residuals {g_lo:.3g}..{g_hi:.3g}); the spec is "
                f"infeasible for this feed"
            )
        v_top = float(brentq(g, lo, hi, rtol=1e-10))
        return cache.get(v_top) or run(v_top)

    def _initial_x(self, pp, N, active, z, P_feed, f, v_top):
        """Stage liquid-composition starting profile: the last converged one
        when the layout/feed are unchanged, else a Fenske-style interpolation
        (:func:`caldyr.unitops.rigorous_column.fenske_profile`)."""
        w = self._warm
        if (w is not None and w["n_stages"] == N and w["active"] == active
                and max(abs(w["z"].get(c, 0.0) - z[c]) for c in active)
                < _WARM_Z_TOL):
            return [dict(r) for r in w["x"]]
        try:
            return fenske_profile(pp, P_feed, f, v_top, N)
        except _NoVLEError as exc:
            raise AbsorberError(
                f"ReboiledAbsorber {self.id!r}: no vapor-liquid equilibrium "
                f"for the feed at P={P_feed:.4g} Pa (bubble flash returned a "
                f"single phase); check the pressure"
            ) from exc


# -- reboiled-absorber inner loop (Wang-Henke, V_top-driven) -------------------
def _solve_reboiled_wh(
    unit_id: str, pp: Any, active: list[str], N: int, P_j: list[float],
    f0: int, f: dict[str, float], F: float, h_F: float, V_top: float,
    x0: list[dict[str, float]], max_iter: int,
) -> dict[str, Any]:
    """Bubble-point (Wang-Henke) MESH for a reboiled absorber with the
    overhead vapor rate fixed at ``V_top``: stage temperatures from bubble
    points, vapor traffic from the energy-balance forward recurrence (the
    N-1 upper-stage balances; the reboiler balance is the one closed by the
    duty, which the caller computes from the overall balance). Returns the
    same profile dict as :func:`_solve_sum_rates`."""
    floor = _FLOW_FLOOR_FRAC * F
    z = {c: f[c] / F for c in active}
    cum = [F if j >= f0 else 0.0 for j in range(N)]   # feed at or above stage j
    A = [cum[j] - V_top for j in range(N)]            # L_j = V_{j+1} + A_j
    B = F - V_top

    x = [dict(row) for row in x0]
    V = [V_top] * N
    L = [max(V[j + 1] + A[j], floor) for j in range(N - 1)] + [max(B, floor)]
    T, K, hL, hV, n_prop = _stage_bubble_points(unit_id, pp, P_j, x, active, 0)

    omega = 1.0
    dT_prev = math.inf
    worsened = 0
    converged = False
    it = 0
    dT = dV = math.inf
    for it in range(1, max_iter + 1):
        cols: dict[str, list[float]] = {}
        for c in active:
            a = [0.0] + [L[j - 1] for j in range(1, N)]
            b = [-(L[j] + V[j] * K[j][c]) for j in range(N)]
            cc = [V[j + 1] * K[j + 1][c] for j in range(N - 1)] + [0.0]
            d = [-F * z[c] if j == f0 else 0.0 for j in range(N)]
            cols[c] = _thomas(a, b, cc, d)
        for j in range(N):
            row = {c: max(cols[c][j], _X_FLOOR) for c in active}
            x[j] = _normalized(row)

        T_new, K, hL, hV, n_prop = _stage_bubble_points(
            unit_id, pp, P_j, x, active, n_prop)
        dT = max(abs(tn - to) for tn, to in zip(T_new, T))
        T = T_new

        if dT > dT_prev:
            worsened += 1
            if worsened >= 2:
                omega = max(0.25, 0.5 * omega)
                worsened = 0
        else:
            worsened = 0
        dT_prev = dT

        # Energy-balance recurrence for the vapor traffic, V_top fixed.
        V_new = list(V)
        for j in range(N - 1):
            denom = hV[j + 1] - hL[j]
            if denom <= 0.0:
                raise AbsorberError(
                    f"{unit_id}: degenerate energy balance on stage {j + 2} "
                    f"(saturated vapor enthalpy below the liquid's, "
                    f"dH={denom:.4g} J/mol) — the property package returned "
                    f"an unphysical state"
                )
            feed_j = F * h_F if j == f0 else 0.0
            if j == 0:
                num = A[0] * hL[0] + V_new[0] * hV[0] - feed_j
            else:
                num = (A[j] * hL[j] - A[j - 1] * hL[j - 1] - feed_j
                       + V_new[j] * (hV[j] - hL[j - 1]))
            V_new[j + 1] = max(num / denom, floor)
        dV = (max(abs(vn - vo) for vn, vo in zip(V_new[1:], V[1:]))
              / max(V_new[1:]))
        V = [V_top] + [max(vo + omega * (vn - vo), floor)
                       for vo, vn in zip(V[1:], V_new[1:])]
        L = [max(V[j + 1] + A[j], floor) for j in range(N - 1)] + \
            [max(B, floor)]

        if dT <= _TOL_T and dV <= _TOL_FLOW:
            converged = True
            break

    if not converged:
        raise AbsorberError(
            f"{unit_id}: bubble-point iteration did not converge in "
            f"{max_iter} iterations (max |dT|={dT:.3g} K vs {_TOL_T} K, "
            f"max |dV|/V={dV:.3g} vs {_TOL_FLOW}, damping={omega}) for "
            f"V_top={V_top:.4g} mol/s — check the spec, or raise 'max_iter'"
        )
    if any(lj <= floor for lj in L) or any(vj <= floor for vj in V):
        raise AbsorberError(
            f"{unit_id}: a column section dried up at the converged point "
            f"(min L={min(L):.3g}, min V={min(V):.3g} mol/s) — "
            f"V_top={V_top:.4g} mol/s is infeasible for this feed"
        )

    y = [_normalized({c: K[j][c] * x[j][c] for c in active}) for j in range(N)]
    # Diagnostic: the reboiler-stage balance vs the overall-balance duty.
    q_overall = V[0] * hV[0] + L[N - 1] * hL[N - 1] - F * h_F
    q_stage = (V[N - 1] * hV[N - 1] + L[N - 1] * hL[N - 1]
               - L[N - 2] * hL[N - 2])
    e_scale = max(abs(q_overall), abs(q_stage), 1.0)
    return {
        "T": T, "L": L, "V": V, "x": x, "y": y, "hL": hL, "hV": hV,
        "iterations": it, "max_dT": dT, "max_dF_rel": dV, "damping": omega,
        "prop_calls": n_prop,
        "energy_residual_rel": abs(q_overall - q_stage) / e_scale,
        "F_total": F,
    }


def _stage_bubble_points(unit_id: str, pp: Any, P_j: list[float],
                         x: list[dict[str, float]], active: list[str],
                         n_prop: int):
    """One saturated-liquid flash per stage: new temperatures, K-values and
    saturated phase enthalpies (the bubble-point analog of the sum-rates
    k_values + Newton-temperatures pair)."""
    T: list[float] = []
    K: list[dict[str, float]] = []
    hL: list[float] = []
    hV: list[float] = []
    for j, row in enumerate(x):
        res = pp.bubble_point(P_j[j], row)
        n_prop += 1
        if (res.y is None or res.H_liquid is None or res.H_vapor is None
                or not math.isfinite(res.T)):
            raise AbsorberError(
                f"{unit_id}: bubble-point flash failed on stage {j + 1} "
                f"(T={res.T!r}, composition {row}) at P={P_j[j]:.4g} Pa"
            )
        T.append(res.T)
        K.append({c: res.y[c] / row[c] for c in active})
        hL.append(res.H_liquid)
        hV.append(res.H_vapor)
    return T, K, hL, hV, n_prop


# -- shared product / design assembly -----------------------------------------
def _build_products(
    uid: str, pp: Any, comps: list[str], active: list[str],
    f_tot: dict[str, float], res: dict[str, Any], P_j: list[float],
    feeds_enthalpy: float, q_added: float,
    vapor_port: str, liquid_port: str,
) -> tuple[dict[str, PortStream], dict[str, Any]]:
    """Assemble the two product streams and the ``design`` dict from a
    converged sum-rates result. Component mass balances close to machine
    precision (the liquid product is the feed minus the vapor product); the
    overall energy balance closes exactly (the liquid product enthalpy is
    resolved from it, its temperature from a PH flash)."""
    N = len(P_j)
    T, L, V = res["T"], res["L"], res["V"]
    x, y, hV = res["x"], res["y"], res["hV"]
    F_total = res["F_total"]

    v_flows = {c: V[0] * y[0][c] for c in active}
    l_flows = {c: f_tot[c] - v_flows[c] for c in active}
    neg = min(l_flows.values())
    if neg < -1e-7 * F_total:
        raise AbsorberError(
            f"{uid!r}: converged overhead vapor would carry more of a "
            f"component than the feeds supply (liquid flow {neg:.3g} mol/s) "
            f"— tighten tolerances or check the specs"
        )
    if neg < 0.0:
        l_flows = {c: max(v, 0.0) for c, v in l_flows.items()}
        scale = (F_total - V[0]) / sum(l_flows.values())
        l_flows = {c: v * scale for c, v in l_flows.items()}
    L_out = sum(l_flows.values())
    x_out = {c: v / L_out for c, v in l_flows.items()}

    vapor = Stream(
        id=f"{uid}.{vapor_port}", components=comps,
        T=T[0], P=P_j[0], molar_flow=V[0], z=_pad(y[0], comps),
        H=hV[0], phase="vapor", vapor_fraction=1.0,
    )
    # Close the overall energy balance exactly: the liquid product carries
    # whatever enthalpy the balance requires; its temperature follows from a
    # PH flash (agrees with the converged stage-N temperature to within the
    # energy-balance tolerance).
    h_liq = (feeds_enthalpy + q_added - V[0] * hV[0]) / L_out
    res_l = pp.flash_ph(P_j[-1], h_liq, x_out)
    liquid = Stream(
        id=f"{uid}.{liquid_port}", components=comps,
        T=res_l.T, P=P_j[-1], molar_flow=L_out, z=_pad(x_out, comps),
        H=h_liq, phase=res_l.phase, vapor_fraction=res_l.vapor_fraction,
    )

    design: dict[str, Any] = {
        "P": P_j[0],
        "n_stages": N,
        "T_top": T[0], "T_bottom": res_l.T,
        "V_top": V[0], "L_bottom": L_out,
        "vapor_flows": dict(v_flows), "liquid_flows": dict(l_flows),
        "x_bottom": dict(x_out), "y_top": dict(y[0]),
        "T_profile": list(T), "P_profile": list(P_j),
        "L_profile": list(L), "V_profile": list(V),
        "x_profile": [_pad(row, comps) for row in x],
        "y_profile": [_pad(row, comps) for row in y],
        "iterations": res["iterations"], "max_dT": res["max_dT"],
        "max_dF_rel": res["max_dF_rel"], "damping": res["damping"],
        "prop_calls": res["prop_calls"],
        "energy_residual_rel": res["energy_residual_rel"],
    }
    out: dict[str, PortStream] = {vapor_port: vapor, liquid_port: liquid}
    return out, design


def _copy_out(out: dict[str, PortStream]) -> dict[str, PortStream]:
    return {
        name: (s.with_() if isinstance(s, Stream)
               else EnergyStream(id=s.id, duty=s.duty))
        for name, s in out.items()
    }


def _copy_design(design: dict[str, Any]) -> dict[str, Any]:
    return {k: (list(v) if isinstance(v, list) else
                dict(v) if isinstance(v, dict) else v)
            for k, v in design.items()}
