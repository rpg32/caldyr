"""Kent-Eisenberg acid-gas solubility model for aqueous alkanolamine solutions.

Gas sweetening (scrubbing CO2/H2S out of a gas with an amine like DEA/MDEA)
works because the acid gas *chemically reacts* with the amine in the liquid to
form ionic species — so only the small free-molecular fraction of the acid gas
has a vapour pressure. A cubic-EOS or activity package cannot represent that
reactive equilibrium; an absorber/regenerator therefore needs a dedicated
acid-gas property method. This is the open analogue of HYSYS's "Acid Gas –
Chemical Solvents" package (Hameed 2025, *Chemical Process Simulations using
Aspen HYSYS*, §15.3).

**Model — modified Kent-Eisenberg** (Kent & Eisenberg 1976), with the constants
and apparent-constant factors of **Haji-Sulaiman, Aroua & Benamor**, *Chem. Eng.
Res. Des.* 76:961 (1998). The liquid-phase ionic equilibria are solved with
literature constants for the carbonate / sulfide / water reactions and Henry's
law; the amine protonation and (for primary/secondary amines) the carbamate
formation carry an apparent-constant factor ``F`` that lumps the solution
non-idealities as a function of CO2 partial pressure and amine concentration.

Reaction set (molarity basis, mol/L):
  * amine protonation     ``AmH+  <->  Am + H+``                  K_prot
  * carbamate formation   ``Am + HCO3-  <->  AmCOO- + H2O``       K_carb (sec. amine)
  * CO2 first ionization  ``CO2 + H2O  <->  H+ + HCO3-``          K_co2
  * bicarbonate           ``HCO3-  <->  H+ + CO3^2-``             K_hco3
  * water                 ``H2O  <->  H+ + OH-``                  K_w
  * H2S first ionization  ``H2S  <->  H+ + HS-``                  K_h2s
  * Henry (physical)      ``P_i = H_i * [i]``         for CO2, H2S
The bisulfide ``HS- <-> H+ + S^2-`` is omitted (pK>12 -> negligible at amine pH).
Everything is eliminated to one charge-balance equation in ``[H+]``; the loading
``alpha = (mol acid gas)/(mol amine)`` follows from each gas's mass balance.

Validation (``tests/test_m16_amine_kent_eisenberg.py``):
  * **CO2 in DEA** — vs Haji-Sulaiman et al. (1998) 2 M / 4 M data, ~7% AAD
    (matches the paper's own 9.2%); the carbamate caps the fast loading near
    0.5 mol/mol and bicarbonate carries it higher, as observed.
  * **CO2 in MDEA** — vs the same source (tertiary amine, no carbamate).
  * **H2S in MDEA** — a prediction from the amine protonation + the literature
    H2S ionization/Henry constants (verified vs pure-water H2S solubility);
    physically validated, quantitative H2S-MDEA data fit tracked as follow-up.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

ATM_PER_KPA = 1.0 / 101.325

# -- temperature-dependent constants: value = exp(a/T + b ln T + c T + d), T[K].
# Acid-dissociation constants are molarity-basis (mol/L); Henry's constants are
# atm.L/mol (so partial pressures convert to atm internally). K_co2/K_hco3/K_w
# and H_CO2 are the literature values used by Haji-Sulaiman 1998 (Edwards et al.
# 1978; CO2 Henry Little et al. 1990); K_h2s is the Edwards 1978 H2S first
# ionization; H_H2S is a van't Hoff fit to the NIST/Sander (2015) H2S water
# solubility. Each verified here against its known 25 C magnitude (the guard
# that caught an OCR-corrupted constant during sourcing).
_K_CO2 = (-12092.1, -36.7816, 0.0, 235.482)    # CO2 + H2O <-> H+ + HCO3-
_K_HCO3 = (-12431.7, -35.4819, 0.0, 220.067)   # HCO3- <-> H+ + CO3^2-
_K_W = (-13445.9, -22.4773, 0.0, 140.932)      # H2O <-> H+ + OH-
_K_H2S = (-12995.4, -33.5471, 0.0, 218.599)    # H2S <-> H+ + HS-
_H_CO2 = (-6789.04, -11.4519, -0.010454, 94.4914)   # Henry CO2, atm.L/mol
_H_H2S = (-2100.0, 0.0, 0.0, 9.3329)           # Henry H2S, atm.L/mol (Sander)


def _arr(p: tuple[float, float, float, float], T: float) -> float:
    a, b, c, d = p
    return math.exp(a / T + b * math.log(T) + c * T + d)


@dataclass(frozen=True)
class AmineSystem:
    """Acid-gas equilibrium parameters for one aqueous amine (Haji-Sulaiman 1998
    Tables 1 & 3).

    ``prot`` is the (a,b,c,d) of the protonation constant ``K_prot`` (acid
    dissociation of the protonated amine). ``prot_F`` = (g, k) gives the
    apparent factor ``F = exp(g ln P_CO2[kPa] + k ln[amine_total])`` so
    ``K_prot' = K_prot·F``. ``carb`` is the carbamate *formation* constant
    (``[AmCOO-] = K_carb·[HCO3-][Am]``; Aroua et al. 1997) for primary/secondary
    amines, with ``carb_F`` = (g, j) and ``F = exp(g ln P_CO2 + j/[amine])``;
    both ``None`` for tertiary amines (no carbamate).
    """
    name: str
    prot: tuple[float, float, float, float]
    prot_F: tuple[float, float]
    carb: tuple[float, float, float, float] | None = None
    carb_F: tuple[float, float] | None = None


# Protonation: DEA from Perrin 1965 (pKa 8.88), MDEA from Littel 1990 (pKa 8.6).
# Carbamate: DEA formation constant from Aroua, Ben Amor & Haji-Sulaiman 1997.
# F-factor exponents from Haji-Sulaiman 1998 Table 3.
DEA = AmineSystem(
    name="DEA",
    prot=(-3071.15, 6.776904, 0.0, -48.7594), prot_F=(-0.4559, 0.2584),
    carb=(-17067.2, -66.8007, 0.0, 439.709), carb_F=(-0.002386, 2.88))
# MDEA protonation F-factor re-fitted JOINTLY to the full CO2 (Haji-Sulaiman
# 1998, 39 pts) and H2S (Jou et al. 1982, 12 pts) sets so both are quantitative —
# H2S ~17% AAD (was a ~37% prediction from the CO2-only fit), CO2 ~24% AAD
# (unchanged band). MDEA is the H2S-selective amine, so a data-fit H2S matters.
MDEA = AmineSystem(
    name="MDEA",
    prot=(-8483.95, -13.8328, 0.0, 87.39717), prot_F=(-0.1068, -0.2445))
# MEA (primary amine, strong carbamate). Protonation: van't Hoff base
# (ΔH_diss ~ 50 kJ/mol so the 1/T coeff is −ΔH/R ≈ −6038) with the constant term
# + pressure exponent FITTED to Lawson & Garst (1976) Table IV (H2S in 15.2 wt%
# MEA, 313–373 K) to 7.4% AAD: effective pKa(298) ≈ 9.95 vs the physical 9.50 —
# the offset absorbs the lumped H2S non-idealities, the Kent-Eisenberg
# fitted-constant approach. Carbamate: a physical formation slope (≈ Aroua 1997's
# 1781·ln10 ≈ 4100 K) with the magnitude FITTED to Lawson Table V (CO2 in 15.2
# wt% MEA — sparse, high-T); the carbamate caps the loading near 0.5 then
# bicarbonate carries it past (see tests). carb_F = 1 (single-concentration data).
MEA = AmineSystem(
    name="MEA",
    prot=(-6038.0, 0.0, 0.0, -2.648), prot_F=(-0.110, 0.0),
    carb=(4100.0, 0.0, 0.0, -8.895), carb_F=(0.0, 0.0))

_SYSTEMS = {"DEA": DEA, "MDEA": MDEA, "MEA": MEA}


def amine_system(name: str) -> AmineSystem:
    """Look up a supported amine by name (case-insensitive)."""
    try:
        return _SYSTEMS[name.upper()]
    except KeyError:
        raise ValueError(
            f"unsupported amine {name!r}; Kent-Eisenberg acid-gas solubility is "
            f"available for {sorted(_SYSTEMS)}"
        ) from None


def _clamp_exp(arg: float) -> float:
    """``exp`` of an apparent-factor exponent, clamped to avoid the overflow that
    extreme *extrapolation* can produce (e.g. a near-zero amine molarity on a
    transient column iterate makes the ``jc/M`` carbamate term explode). The
    clamp is far outside the fitted range, so valid results are untouched."""
    return math.exp(min(max(arg, -700.0), 700.0))


def _apparent(sys: AmineSystem, T: float, P_prot_kpa: float, P_co2_kpa: float,
              amine_M: float) -> tuple[float, float | None]:
    """Apparent protonation and carbamate constants K' = K·F at (T, P, M).

    The protonation apparent factor is an ionic-strength correction, so it is
    driven by the **total acid-gas partial pressure** ``P_prot_kpa`` = P_CO2 +
    P_H2S (Haji-Sulaiman parameterized it with P_CO2 as the proxy because their
    data were CO2-only; using the total generalizes it to H2S and mixed systems
    and reduces to the original for CO2-only). The carbamate apparent factor
    stays CO2-specific (the carbamate forms from bicarbonate). Earlier the
    protonation factor used P_CO2 even for H2S, where the 1e-6 floor blew K_prot
    up by ~700x and made H2S loadings ~30x too low."""
    pp = max(P_prot_kpa, 1e-6)         # ionic-strength proxy: total acid gas
    pc = max(P_co2_kpa, 1e-6)          # carbamate factor is CO2-specific
    m = max(amine_M, 1e-6)             # floor: F is undefined at zero amine
    g, k = sys.prot_F
    kp = _arr(sys.prot, T) * _clamp_exp(g * math.log(pp) + k * math.log(m))
    kc = None
    if sys.carb is not None and sys.carb_F is not None:
        gc, jc = sys.carb_F
        kc = _arr(sys.carb, T) * _clamp_exp(gc * math.log(pc) + jc / m)
    return kp, kc


def speciate(T: float, P_co2: float, P_h2s: float, amine_M: float,
             amine: str | AmineSystem = "MDEA") -> dict[str, float]:
    """Solve the liquid-phase ionic equilibrium at (T[K], P_CO2[kPa], P_H2S[kPa],
    amine molarity). Returns species concentrations (mol/L) plus ``H+``."""
    sys = amine if isinstance(amine, AmineSystem) else amine_system(amine)
    if amine_M <= 0:
        raise ValueError("amine molarity must be positive")
    Kco2, Khco3, Kw, Kh2s = (_arr(p, T) for p in (_K_CO2, _K_HCO3, _K_W, _K_H2S))
    co2 = max(P_co2, 0.0) * ATM_PER_KPA / _arr(_H_CO2, T)
    h2s = max(P_h2s, 0.0) * ATM_PER_KPA / _arr(_H_H2S, T)
    Kp, Kc = _apparent(sys, T, P_co2 + P_h2s, P_co2, amine_M)

    def carb_ratio(h: float) -> float:        # [AmCOO-]/[Am] = K_carb'·[HCO3-]
        return Kc * Kco2 * co2 / h if Kc is not None else 0.0

    def charge(h: float) -> float:
        hco3 = Kco2 * co2 / h
        co3 = Kco2 * Khco3 * co2 / h / h
        hs = Kh2s * h2s / h
        n = amine_M / (1.0 + h / Kp + carb_ratio(h))
        return h + n * h / Kp - Kw / h - hco3 - 2.0 * co3 - hs - n * carb_ratio(h)

    h = brentq(charge, 1e-14, 1e-1, xtol=1e-25, rtol=1e-14, maxiter=400)
    n = amine_M / (1.0 + h / Kp + carb_ratio(h))
    return {
        "H+": h, "OH-": Kw / h,
        "CO2": co2, "HCO3-": Kco2 * co2 / h, "CO3--": Kco2 * Khco3 * co2 / h / h,
        "AmCOO-": n * carb_ratio(h),
        "H2S": h2s, "HS-": Kh2s * h2s / h,
        "amine": n, "amineH+": n * h / Kp,
    }


def acid_gas_loadings(T: float, P_co2: float, P_h2s: float, amine_M: float,
                      amine: str | AmineSystem = "MDEA") -> tuple[float, float]:
    """``(alpha_CO2, alpha_H2S)`` loadings (mol acid gas / mol amine) over CO2
    and H2S partial pressures (kPa)."""
    sp = speciate(T, P_co2, P_h2s, amine_M, amine)
    a_co2 = (sp["CO2"] + sp["HCO3-"] + sp["CO3--"] + sp["AmCOO-"]) / amine_M
    a_h2s = (sp["H2S"] + sp["HS-"]) / amine_M
    return a_co2, a_h2s


def co2_loading(T: float, P_co2: float, amine_M: float,
                amine: str | AmineSystem = "MDEA") -> float:
    """CO2 loading ``alpha`` (mol CO2 / mol amine) over a CO2 partial pressure
    (``P_co2`` in **kPa**), no H2S present."""
    return acid_gas_loadings(T, P_co2, 0.0, amine_M, amine)[0]


def h2s_loading(T: float, P_h2s: float, amine_M: float,
                amine: str | AmineSystem = "MDEA") -> float:
    """H2S loading ``alpha`` (mol H2S / mol amine) over an H2S partial pressure
    (``P_h2s`` in **kPa**), no CO2 present (the protonation F-factor then uses
    its CO2-pressure floor — H2S-only loading is a prediction, see module docs)."""
    return acid_gas_loadings(T, 0.0, P_h2s, amine_M, amine)[1]


def _invert(load_fn, T: float, alpha: float, amine_M: float,
            amine: str | AmineSystem) -> float:
    if alpha <= 0.0:
        return 0.0

    def resid(logp: float) -> float:
        return load_fn(T, math.exp(logp), amine_M, amine) - alpha

    lo, hi = math.log(1e-4), math.log(5.0e4)
    if resid(hi) < 0.0:
        return math.exp(hi)
    if resid(lo) > 0.0:
        return math.exp(lo)
    return math.exp(brentq(resid, lo, hi, xtol=1e-8, rtol=1e-10, maxiter=200))


def co2_partial_pressure(T: float, alpha: float, amine_M: float,
                         amine: str | AmineSystem = "MDEA") -> float:
    """Equilibrium CO2 partial pressure (**kPa**) over an amine solution at CO2
    loading ``alpha`` — the vapour-side driving force an absorber stage needs."""
    return _invert(co2_loading, T, alpha, amine_M, amine)


def h2s_partial_pressure(T: float, alpha: float, amine_M: float,
                         amine: str | AmineSystem = "MDEA") -> float:
    """Equilibrium H2S partial pressure (**kPa**) over an amine solution at H2S
    loading ``alpha``."""
    return _invert(h2s_loading, T, alpha, amine_M, amine)


def acid_gas_partial_pressures(
    T: float, alpha_co2: float, alpha_h2s: float, amine_M: float,
    amine: str | AmineSystem = "MDEA",
) -> tuple[float, float]:
    """Joint inverse of :func:`acid_gas_loadings`: the ``(P_CO2, P_H2S)`` partial
    pressures (**kPa**) whose equilibrium loadings equal ``(alpha_co2,
    alpha_h2s)``.

    When both acid gases are present they compete for the amine through the
    shared proton/charge balance, so the equilibrium CO2 pressure depends on the
    H2S loading and vice-versa — inverting each gas independently (with
    :func:`co2_partial_pressure` / :func:`h2s_partial_pressure`) ignores that
    coupling. This solves the 2x2 system simultaneously (a damped Newton on the
    log-pressures, seeded from the independent inversions, which is already
    close). Degenerate cases fall back to the relevant 1-D inversion. This is the
    vapour-side driving force a CO2+H2S sweetening absorber/regenerator stage
    needs."""
    import numpy as np

    if alpha_co2 <= 0.0 and alpha_h2s <= 0.0:
        return 0.0, 0.0
    if alpha_h2s <= 0.0:
        return co2_partial_pressure(T, alpha_co2, amine_M, amine), 0.0
    if alpha_co2 <= 0.0:
        return 0.0, h2s_partial_pressure(T, alpha_h2s, amine_M, amine)

    p = np.array([
        math.log(max(co2_partial_pressure(T, alpha_co2, amine_M, amine), 1e-6)),
        math.log(max(h2s_partial_pressure(T, alpha_h2s, amine_M, amine), 1e-6)),
    ])
    target = np.array([math.log(alpha_co2), math.log(alpha_h2s)])

    def f(pv: np.ndarray) -> np.ndarray:
        c, h = acid_gas_loadings(T, math.exp(pv[0]), math.exp(pv[1]), amine_M, amine)
        return np.array([math.log(max(c, 1e-300)), math.log(max(h, 1e-300))]) - target

    eps = 1e-6
    for _ in range(40):
        r = f(p)
        if float(np.max(np.abs(r))) < 1e-9:
            break
        jac = np.zeros((2, 2))
        for k in range(2):
            pk = p.copy()
            pk[k] += eps
            jac[:, k] = (f(pk) - r) / eps
        try:
            dp = np.linalg.solve(jac, -r)
        except np.linalg.LinAlgError:
            break
        p = p + dp
    return math.exp(float(p[0])), math.exp(float(p[1]))
