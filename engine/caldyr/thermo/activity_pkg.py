"""Activity-coefficient (gamma-phi) property package wrapping `thermo`.

The liquid is modelled with an excess-Gibbs activity-coefficient model and the
vapor as an ideal gas — the standard low-pressure gamma-phi approach. Unlike a
cubic EOS, this captures strongly non-ideal polar liquids and the azeotropes
they form (e.g. ethanol + water), and it stays robust on near-pure-water trays
where a cubic EOS loses its all-liquid root.

Three excess-Gibbs models are available, selected by the ``property_package``
string:

  * ``thermo:NRTL`` — NRTL with binary interaction parameters from the ChemSep
    dataset bundled with `thermo` (via its IPDB). Pairs without data fall back
    to an ideal-solution liquid, warned about rather than applied silently.
    Excellent VLE for fitted pairs; the bundled params are VLE-fit, so a single
    NRTL set carries no liquid-liquid miscibility gap.
  * ``thermo:UNIFAC`` — **Modified UNIFAC (Dortmund)**, a predictive
    group-contribution model fit *simultaneously* to VLE and LLE. Per-component
    UNIFAC groups are assigned automatically from the DDBST fragmentation
    database shipped with `thermo` (no per-system fitting). This is the model
    that unblocks heteroazeotropic (3-phase) distillation: it gives both the
    heteroazeotrope (VLE) and the decanter's liquid-liquid split (LLE) from one
    consistent set of parameters.
  * ``thermo:UNIFAC-LLE`` — the original Magnussen (1981) UNIFAC-LLE parameter
    table (version 0), LLE-fit. Stronger near a miscibility gap, weaker VLE than
    the Dortmund model; offered as a secondary cross-check.

Validated against the ethanol/water minimum-boiling azeotrope (~89.4 mol%
ethanol, 78.2 C at 1 atm — DePriester / Gmehling DECHEMA), which a cubic EOS
cannot represent, and (UNIFAC) the ethanol/water/cyclohexane ternary decant of
heterogeneous-azeotrope ethanol dehydration.
"""
from __future__ import annotations

import math
import warnings
from functools import lru_cache

from .base import ThreePhaseResult
from ._flasher import FlasherPackage, formation_props

# ChemSep NRTL stores tau as b_ij / T and the non-randomness as alpha_ij.
_NRTL_TABLE = "ChemSep NRTL"

# Supported excess-Gibbs liquid models (uppercase keys, as `make_package` sends).
_NRTL = "NRTL"
_UNIFAC_DMD = "UNIFAC"        # Modified UNIFAC (Dortmund), thermo version=1
_UNIFAC_LLE = "UNIFAC-LLE"    # original Magnussen LLE table, thermo version=0


@lru_cache(maxsize=1)
def _orig_to_lle_subgroup() -> dict[int, int]:
    """Remap original-UNIFAC subgroup IDs (what the DDBST 'UNIFAC' fragmentation
    returns) onto the UNIFAC-LLE table's subgroup IDs. The two numberings agree
    for the hydrocarbon backbone but diverge from subgroup 15 on (the LLE table
    inserts the P1/P2 main groups, shifting H2O from 16 to 17, etc.). We match by
    subgroup *name*, which is stable across the two tables."""
    from thermo.unifac import LLEUFSG, UFSG

    name_to_lle: dict[str, int] = {}
    for sid, entry in LLEUFSG.items():
        name_to_lle.setdefault(entry.group, sid)
    return {sid: name_to_lle[entry.group]
            for sid, entry in UFSG.items() if entry.group in name_to_lle}


def _unifac_chemgroups(cass: tuple[str, ...], ddbst_model: str) -> list[dict[int, int]]:
    """Per-component UNIFAC subgroup-count dicts from the DDBST fragmentation
    database (``ddbst_model`` is 'MODIFIED_UNIFAC' or 'UNIFAC'). Raises naming
    any component the database cannot fragment — there is no honest gamma without
    its groups."""
    from thermo.unifac import UNIFAC_group_assignment_DDBST

    groups: list[dict[int, int]] = []
    missing: list[str] = []
    for cas in cass:
        g = UNIFAC_group_assignment_DDBST(cas, ddbst_model)
        if not g:
            missing.append(cas)
        groups.append(g)
    if missing:
        raise ValueError(
            f"no DDBST {ddbst_model} group assignment for component(s) {missing}; "
            f"UNIFAC cannot model them. Use 'thermo:NRTL' (if ChemSep has the "
            f"binary parameters) or a cubic EOS instead."
        )
    return groups


def _build_unifac_ge(cass: tuple[str, ...], model: str):
    """Build a `thermo` UNIFAC excess-Gibbs model for an ordered component tuple,
    with groups assigned automatically from the DDBST database. ``model`` is
    ``_UNIFAC_DMD`` (Dortmund, version 1) or ``_UNIFAC_LLE`` (Magnussen, version
    0)."""
    from thermo.unifac import UNIFAC, LLEUFIP, LLEUFSG

    n = len(cass)
    xs = [1.0 / n] * n
    if model == _UNIFAC_DMD:
        chemgroups = _unifac_chemgroups(cass, "MODIFIED_UNIFAC")
        return UNIFAC.from_subgroups(T=298.15, xs=xs, chemgroups=chemgroups,
                                     version=1)
    # UNIFAC-LLE: DDBST returns original-UNIFAC numbering; remap to the LLE table.
    remap = _orig_to_lle_subgroup()
    orig = _unifac_chemgroups(cass, "UNIFAC")
    lle_groups: list[dict[int, int]] = []
    for cas, g in zip(cass, orig):
        out: dict[int, int] = {}
        for sid, count in g.items():
            if sid not in remap:
                raise ValueError(
                    f"component {cas} uses UNIFAC subgroup {sid}, which has no "
                    f"UNIFAC-LLE equivalent; use 'thermo:UNIFAC' (Dortmund) for "
                    f"this system."
                )
            out[remap[sid]] = out.get(remap[sid], 0) + count
        lle_groups.append(out)
    return UNIFAC.from_subgroups(T=298.15, xs=xs, chemgroups=lle_groups,
                                 version=0, subgroups=LLEUFSG,
                                 interaction_data=LLEUFIP)


# Liquid-liquid (isoactivity) flash tuning.
_LLE_MAX_IT = 300        # successive-substitution iterations per seed
_LLE_TOL = 1e-11         # max |d ln K| convergence on the distribution ratios
_LLE_TRIVIAL = 0.02      # min |x_I - x_II|_inf for a split to count as two phases
_LLE_PSI_EPS = 1e-4      # keep the phase split strictly interior (0, 1)
_LLE_G_TOL = 1e-9        # the split must lower the molar Gibbs energy by this


def _has_offdiagonal(matrix) -> bool:
    return any(abs(matrix[i][j]) > 0.0
               for i in range(len(matrix)) for j in range(len(matrix)) if i != j)


@lru_cache(maxsize=64)
def _build_flasher(components: tuple[str, ...], model: str):
    """Build (and cache) a gamma-phi thermo flasher for an ordered component
    tuple. Pure-component systems use `FlashPureVLS` (the activity model is
    irrelevant for one component: gamma = 1)."""
    from thermo import (
        ChemicalConstantsPackage,
        FlashPureVLS,
        FlashVL,
        GibbsExcessLiquid,
        IdealGas,
        NRTL,
    )
    from thermo.interaction_parameters import IPDB

    if model not in (_NRTL, _UNIFAC_DMD, _UNIFAC_LLE):
        raise ValueError(
            f"unsupported activity model {model!r}; expected one of "
            f"{(_NRTL, _UNIFAC_DMD, _UNIFAC_LLE)}"
        )

    constants, props = ChemicalConstantsPackage.from_IDs(list(components))
    gas = IdealGas(HeatCapacityGases=props.HeatCapacityGases)
    hf, gf = formation_props(components, constants.Hfgs, constants.Gfgs)
    liquid_kwargs = dict(
        VaporPressures=props.VaporPressures,
        HeatCapacityGases=props.HeatCapacityGases,
        VolumeLiquids=props.VolumeLiquids,
        equilibrium_basis="Psat",
        caloric_basis="Psat",
    )

    if len(components) == 1:
        liquid = GibbsExcessLiquid(**liquid_kwargs)
        return FlashPureVLS(constants, props, gas=gas, liquids=[liquid], solids=[]), hf, gf

    cas = constants.CASs
    if model == _NRTL:
        tau_bs = IPDB.get_ip_asymmetric_matrix(_NRTL_TABLE, cas, "bij")
        alpha_cs = IPDB.get_ip_asymmetric_matrix(_NRTL_TABLE, cas, "alphaij")
        if not _has_offdiagonal(tau_bs):
            warnings.warn(
                f"no ChemSep NRTL parameters for {list(components)}; the liquid "
                f"falls back to an ideal solution (no activity correction). Results "
                f"for polar mixtures may be inaccurate.",
                stacklevel=2,
            )
        n = len(components)
        ge_model = NRTL(T=298.15, xs=[1.0 / n] * n, tau_bs=tau_bs, alpha_cs=alpha_cs)
    else:
        ge_model = _build_unifac_ge(tuple(cas), model)

    liquid = GibbsExcessLiquid(GibbsExcessModel=ge_model, **liquid_kwargs)
    return FlashVL(constants, props, liquid=liquid, gas=gas), hf, gf


class ActivityPackage(FlasherPackage):
    """Gamma-phi property package: an excess-Gibbs liquid (NRTL or UNIFAC) plus
    an ideal-gas vapor, over a fixed component list. Selected by
    ``property_package`` strings ``"thermo:NRTL"``, ``"thermo:UNIFAC"`` (Modified
    UNIFAC / Dortmund) or ``"thermo:UNIFAC-LLE"`` (Magnussen LLE table)."""

    SUPPORTED = (_NRTL, _UNIFAC_DMD, _UNIFAC_LLE)

    def __init__(self, components: list[str], model: str = _NRTL) -> None:
        model = model.upper()
        if model not in self.SUPPORTED:
            raise ValueError(
                f"unsupported activity model {model!r}; expected one of {self.SUPPORTED}"
            )
        self.model = model
        flasher, hf, gf = _build_flasher(tuple(components), model)
        self._init(components, flasher, hf, gf)

    @classmethod
    def from_spec(cls, spec: str, components: list[str]) -> "ActivityPackage":
        """Build from a flowsheet `property_package` string like ``"thermo:NRTL"``
        or ``"thermo:UNIFAC"``."""
        backend, _, model = spec.partition(":")
        if backend != "thermo":
            raise ValueError(
                f"ActivityPackage cannot build backend {backend!r} (got {spec!r})"
            )
        return cls(components, model or _NRTL)

    # -- three-phase (VLLE) flashes via an NRTL liquid-liquid split ----------
    # The cubic-EOS FlasherPackage builds three-phase results from thermo's
    # FlashVLN (two trial EOS liquids). That route does NOT work for a
    # Gibbs-excess (activity) liquid: thermo's FlashVLN stability search leaves
    # the two GE liquids on the same root, so an NRTL miscibility gap is never
    # found. Instead we do the VLLE in two transparent steps: the existing
    # gamma-phi VL flash sets the vapor and the (combined) liquid, then an
    # isoactivity liquid-liquid flash (gamma_i^I x_i^I = gamma_i^II x_i^II,
    # Rachford-Rice on the distribution ratios K_i = gamma_i^I/gamma_i^II)
    # splits that liquid into two phases when an NRTL miscibility gap exists.
    # This is the standard low-pressure VLLE method (Seader 3e sec. 4.5 /
    # Prausnitz) and it correctly captures heterogeneous azeotropes (the
    # decant of a water/organic overhead) that the cubic EOS cannot.
    def _ge_model(self):
        liq = getattr(self._flasher, "liquid", None)
        if liq is None:
            liqs = getattr(self._flasher, "liquids", None)
            liq = liqs[0] if liqs else None
        return getattr(liq, "GibbsExcessModel", None)

    def _gammas(self, T: float, xs: list[float]) -> list[float]:
        ge = self._ge_model()
        if ge is None:
            return [1.0] * len(xs)
        return [float(g) for g in ge.to_T_xs(float(T), list(xs)).gammas()]

    def _g_mix(self, T: float, xs: list[float]) -> float:
        """Molar Gibbs energy of mixing over RT: sum_i x_i (ln x_i + ln gamma_i).
        Lower means more stable; the LLE split is accepted only if the
        two-phase value beats the single-liquid one."""
        gam = self._gammas(T, xs)
        return sum(xi * (math.log(xi) + math.log(gi))
                   for xi, gi in zip(xs, gam) if xi > 0.0)

    @staticmethod
    def _rr_psi(z: list[float], K: list[float]) -> float | None:
        """Solve the Rachford-Rice phase fraction for a liquid-liquid split
        (psi = fraction in phase II, distribution K_i = x_i^II / x_i^I). Returns
        the interior root in (0, 1) or None when the phases do not separate."""
        def f(psi: float) -> float:
            return sum(zi * (Ki - 1.0) / (1.0 + psi * (Ki - 1.0))
                       for zi, Ki in zip(z, K))
        f0, f1 = f(0.0), f(1.0)
        if f0 * f1 > 0.0:                      # no sign change -> single phase
            return None
        lo, hi = 0.0, 1.0
        for _ in range(100):                   # bisection (monotone in psi)
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if f0 * fm <= 0.0:
                hi = mid
            else:
                lo, f0 = mid, fm
        return 0.5 * (lo + hi)

    def _lle_one(self, T: float, z: list[float], K0: list[float]):
        """One isoactivity liquid-liquid flash from an initial distribution
        ``K0`` by successive substitution. Returns ``(psi, xI, xII)`` or None."""
        n = len(z)
        K = list(K0)
        psi = 0.5
        for _ in range(_LLE_MAX_IT):
            root = self._rr_psi(z, K)
            if root is None:
                return None
            psi = min(max(root, _LLE_PSI_EPS), 1.0 - _LLE_PSI_EPS)
            xI = [z[i] / (1.0 + psi * (K[i] - 1.0)) for i in range(n)]
            xII = [K[i] * xI[i] for i in range(n)]
            sI, sII = sum(xI), sum(xII)
            if sI <= 0.0 or sII <= 0.0:
                return None
            xI = [v / sI for v in xI]
            xII = [v / sII for v in xII]
            gI, gII = self._gammas(T, xI), self._gammas(T, xII)
            Knew = [gI[i] / gII[i] for i in range(n)]
            if max(abs(math.log(Knew[i] / K[i])) for i in range(n)) < _LLE_TOL:
                K = Knew
                break
            K = Knew
        root = self._rr_psi(z, K)
        if root is None:
            return None
        psi = root
        if not (_LLE_PSI_EPS < psi < 1.0 - _LLE_PSI_EPS):
            return None
        xI = [z[i] / (1.0 + psi * (K[i] - 1.0)) for i in range(n)]
        xII = [K[i] * xI[i] for i in range(n)]
        sI, sII = sum(xI), sum(xII)
        xI = [v / sI for v in xI]
        xII = [v / sII for v in xII]
        if max(abs(a - b) for a, b in zip(xI, xII)) < _LLE_TRIVIAL:
            return None                        # trivial (identical phases)
        return psi, xI, xII

    def _lle_split(self, T: float, z: list[float]):
        """Split a liquid of composition ``z`` into two liquid phases at ``T``
        if an NRTL miscibility gap exists. Tries seeds enriched in each
        component (to break the trivial symmetry) and keeps the converged split
        that most lowers the Gibbs energy. Returns ``(psi, xI, xII)`` (psi =
        fraction in phase II) or None for a single stable liquid."""
        n = len(z)
        if n < 2 or self._ge_model() is None:
            return None
        g_single = self._g_mix(T, z)
        best = None
        best_g = g_single - _LLE_G_TOL
        for k in range(n):                     # seed phase II near-pure in comp k
            xII0 = [(0.02 / (n - 1)) for _ in range(n)]
            xII0[k] = 0.98
            gI0 = self._gammas(T, z)
            gII0 = self._gammas(T, xII0)
            K0 = [gI0[i] / gII0[i] for i in range(n)]
            out = self._lle_one(T, z, K0)
            if out is None:
                continue
            psi, xI, xII = out
            g_two = (1.0 - psi) * self._g_mix(T, xI) + psi * self._g_mix(T, xII)
            if g_two < best_g:
                best_g, best = g_two, (psi, xI, xII)
        return best

    def flash_pt_3p(self, T: float, P: float,
                    z: dict[str, float]) -> ThreePhaseResult:
        base = self.flash_pt(T, P, z)
        return self._assemble_3p(T, P, z, base)

    def flash_ph_3p(self, P: float, H: float,
                    z: dict[str, float]) -> ThreePhaseResult:
        base = self.flash_ph(P, H, z)
        return self._assemble_3p(base.T, P, z, base)

    def _assemble_3p(self, T: float, P: float, z: dict[str, float],
                     base) -> ThreePhaseResult:
        """Combine the gamma-phi VL flash (``base``) with an isoactivity
        liquid-liquid split of its liquid into a :class:`ThreePhaseResult`
        (vapor + light + heavy liquid, liquids ordered by mass density)."""
        T = float(T)
        vf = float(base.vapor_fraction)
        y = base.y
        h_vapor = base.H_vapor
        liq_frac = 1.0 - vf
        x_liq = base.x if base.x is not None else self._comp(self._zs(z))
        split = self._lle_split(T, self._zs(x_liq)) if liq_frac > 1e-12 else None

        def liq_entry(xs: list[float]):
            comp = self._comp(xs)
            h = self.enthalpy_liquid(T, P, comp)
            v = self.volume_liquid(T, P, comp)
            mw = sum(xi * mi for xi, mi in zip(xs, self._mw_list()))
            rho = mw / v if v > 0.0 else None
            return comp, h, rho

        if split is None:
            # Single liquid: report it as the light phase (heavy absent), the
            # same graceful degradation as the cubic-EOS three-phase result.
            if liq_frac > 1e-12:
                xs = self._zs(x_liq)
                xl, hl, rhol = liq_entry(xs)
                bl, bh, xh, hh, rhoh = liq_frac, 0.0, None, None, None
            else:
                xl = hl = rhol = None
                bl = bh = 0.0
                xh = hh = rhoh = None
        else:
            psi, xI, xII = split
            (cI, hI, rI), (cII, hII, rII) = liq_entry(xI), liq_entry(xII)
            betaI, betaII = liq_frac * (1.0 - psi), liq_frac * psi
            # order by mass density: light = less dense
            if (rI or 0.0) <= (rII or 0.0):
                xl, hl, rhol, bl = cI, hI, rI, betaI
                xh, hh, rhoh, bh = cII, hII, rII, betaII
            else:
                xl, hl, rhol, bl = cII, hII, rII, betaII
                xh, hh, rhoh, bh = cI, hI, rI, betaI

        h_bulk = (vf * (h_vapor or 0.0) + bl * (hl or 0.0) + bh * (hh or 0.0))
        return ThreePhaseResult(
            T=T, P=float(P), H=h_bulk,
            beta_vapor=vf, beta_light=bl, beta_heavy=bh,
            y=dict(y) if y is not None else None,
            x_light=xl, x_heavy=xh,
            H_vapor=h_vapor, H_light=hl, H_heavy=hh,
            rho_light=rhol, rho_heavy=rhoh,
        )

    def _mw_list(self) -> list[float]:
        """Per-component molar mass (kg/mol) for liquid-density ordering."""
        if self._mws is None:
            mws = getattr(self._flasher.constants, "MWs", None)
            self._mws = ([float(m) / 1000.0 for m in mws] if mws is not None
                         else [0.0] * len(self.components))
        return self._mws
