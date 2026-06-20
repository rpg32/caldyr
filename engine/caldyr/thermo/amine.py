"""Kent-Eisenberg acid-gas solubility model for aqueous alkanolamine solutions.

Gas sweetening (scrubbing CO2/H2S out of a gas with an amine like MEA/DEA/MDEA)
works because the acid gas *chemically reacts* with the amine in the liquid to
form ionic species — so only the small free-molecular fraction of the acid gas
has a vapour pressure. A cubic-EOS or activity package cannot represent that
reactive equilibrium; an absorber/regenerator therefore needs a dedicated
acid-gas property method. This is the open analogue of HYSYS's "Acid Gas –
Chemical Solvents" package (Hameed 2025, *Chemical Process Simulations using
Aspen HYSYS*, §15.3).

**Model — modified Kent-Eisenberg** (Kent & Eisenberg, *Hydrocarbon Processing*
55(2):87, 1976; Haji-Sulaiman, Aroua & Benamor, *Chem. Eng. Res. Des.* 76:961,
1998). The liquid-phase ionic equilibria are solved with *literature* constants
for the carbonate / sulfide / water reactions and Henry's law, and a *fitted*
apparent constant for the amine protonation that lumps the solution
non-idealities (the defining idea of the method) with a single amine-
concentration correction ``[amine]**conc_exp``.

Scope today: **CO2 and H2S in aqueous MDEA** (methyldiethanolamine — a *tertiary*
amine, so no carbamate forms and the reaction set is the minimal one). The
CO2-MDEA branch is validated against the Haji-Sulaiman et al. (1998) 2 M / 4 M
solubility data; the H2S branch is a *prediction* — the amine protonation
constant is taken straight from the CO2 fit and combined with the literature
H2S ionization and Henry constants (verified against pure-water H2S solubility)
— so it is physically grounded but not yet fitted/validated against amine-
specific H2S data (see ``tests/test_m16_amine_kent_eisenberg.py``). The structure
(an ``AmineSystem`` carrying the amine-specific constants; a general ionic-
equilibrium solver over both acid gases) is built to extend to carbamate-forming
primary/secondary amines (MEA, DEA).

Reaction set (molarity basis, mol/L):
  * amine protonation     ``AmH+  <->  Am + H+``                 (K1, fitted)
  * CO2 first ionization  ``CO2 + H2O  <->  H+ + HCO3-``         (K2)
  * bicarbonate           ``HCO3-  <->  H+ + CO3^2-``            (K3)
  * water                 ``H2O  <->  H+ + OH-``                 (K4 = Kw)
  * H2S first ionization  ``H2S  <->  H+ + HS-``                 (K5)
  * Henry (physical)      ``P_i = H_i * [i]``        for CO2, H2S
The bisulfide ``HS- <-> H+ + S^2-`` is omitted: its pK (>12) makes ``S^2-``
negligible at amine-solution pH (8-11). Unknowns are eliminated to a single
charge-balance equation in ``[H+]``, solved by a bracketed root find; each acid
gas's loading ``alpha = (mol gas)/(mol amine)`` follows from its mass balance.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

ATM_PER_KPA = 1.0 / 101.325

# -- temperature-dependent constants: value = exp(a/T + b ln T + c T + d), T[K].
# The acid-dissociation constants are molarity-basis (mol/L); the Henry's
# constants are in atm.L/mol (so partial pressures here are in atm internally).
# K2,K3,K4 and H_CO2 are the canonical literature values (Edwards, Maurer,
# Newman & Prausnitz, *AIChE J.* 24:966, 1978; CO2 Henry from Little, Bos &
# Knoop 1990); K5 is the Edwards 1978 H2S first-ionization constant; H_H2S is a
# van't Hoff fit to the NIST/Sander (2015) H2S water solubility. Every one is
# verified here against its known 25 C magnitude (K2~4.4e-7, K3~4.7e-11,
# Kw~1e-14, H_CO2~29, pKa1(H2S)~7.0, [H2S]~0.1 mol/L at 1 atm) — the guard that
# caught an OCR-corrupted constant during sourcing.
_K2 = (-12092.1, -36.7816, 0.0, 235.482)    # CO2 + H2O <-> H+ + HCO3-
_K3 = (-12431.7, -35.4819, 0.0, 220.067)    # HCO3- <-> H+ + CO3^2-
_K4 = (-13445.9, -22.4773, 0.0, 140.932)    # H2O <-> H+ + OH-  (Kw)
_K5 = (-12995.4, -33.5471, 0.0, 218.599)    # H2S <-> H+ + HS-
_H_CO2 = (-6789.04, -11.4519, -0.010454, 94.4914)   # Henry CO2, atm.L/mol
_H_H2S = (-2100.0, 0.0, 0.0, 9.3329)        # Henry H2S, atm.L/mol (Sander 2015)


def _arr(p: tuple[float, float, float, float], T: float) -> float:
    a, b, c, d = p
    return math.exp(a / T + b * math.log(T) + c * T + d)


@dataclass(frozen=True)
class AmineSystem:
    """Acid-gas equilibrium parameters for one aqueous amine.

    ``k1_a``/``k1_b`` give the amine protonation constant ``K1 = exp(k1_a +
    k1_b/T)`` (acid dissociation of the protonated amine, mol/L). The *apparent*
    protonation constant used in the charge balance is ``K1' = K1 *
    [amine]**conc_exp`` — a single amine-concentration (ionic-strength)
    correction (Haji-Sulaiman 1998), shared by both acid gases since it is a
    property of the amine, not of CO2 or H2S. ``carbamate`` flags
    primary/secondary amines (MEA, DEA) that also form a carbamate ion — not
    modelled yet, so only tertiary amines are exposed.
    """
    name: str
    k1_a: float
    k1_b: float
    conc_exp: float
    carbamate: bool = False


# K1 fitted to the Haji-Sulaiman et al. (1998) 2 M + 4 M MDEA CO2 solubility
# data (303-323 K, 0.1-100 kPa) with the carbonate/water/Henry constants held at
# their literature values. The fit gives K1(298 K)=1.7e-9 (pKa(MDEAH+)~8.8,
# consistent with the ~8.5 literature value) and reproduces the data to ~7 %
# (2 M) / ~19 % (4 M) AAD — see tests/test_m16_amine_kent_eisenberg.py.
MDEA = AmineSystem(name="MDEA", k1_a=-10.9358, k1_b=-2758.80,
                   conc_exp=0.8338, carbamate=False)

_SYSTEMS = {"MDEA": MDEA}


def amine_system(name: str) -> AmineSystem:
    """Look up a supported amine by name (case-insensitive)."""
    try:
        return _SYSTEMS[name.upper()]
    except KeyError:
        raise ValueError(
            f"unsupported amine {name!r}; Kent-Eisenberg acid-gas solubility is "
            f"available for {sorted(_SYSTEMS)} — primary/secondary amines (MEA, "
            f"DEA) need the carbamate reaction, not yet modelled"
        ) from None


def speciate(T: float, P_co2: float, P_h2s: float, amine_M: float,
             amine: str | AmineSystem = "MDEA") -> dict[str, float]:
    """Solve the liquid-phase ionic equilibrium at (T, P_CO2, P_H2S, amine M).

    ``T`` in K, partial pressures in **kPa**, ``amine_M`` the amine molarity
    (mol/L). Returns the molar concentrations (mol/L) of every species plus
    ``H+``. The apparent protonation constant depends on the free-amine
    concentration, so a short fixed-point loop wraps the ``[H+]`` root find.
    """
    sys = amine if isinstance(amine, AmineSystem) else amine_system(amine)
    if amine_M <= 0:
        raise ValueError("amine molarity must be positive")
    K2, K3, K4, K5 = (_arr(p, T) for p in (_K2, _K3, _K4, _K5))
    co2 = max(P_co2, 0.0) * ATM_PER_KPA / _arr(_H_CO2, T)   # free molecular CO2
    h2s = max(P_h2s, 0.0) * ATM_PER_KPA / _arr(_H_H2S, T)   # free molecular H2S
    free_amine = amine_M
    h = 1e-9
    for _ in range(60):
        K1p = math.exp(sys.k1_a + sys.k1_b / T) * (free_amine ** sys.conc_exp)

        def charge(hp: float) -> float:
            # [H+] + [AmH+] - [OH-] - [HCO3-] - 2[CO3^2-] - [HS-] = 0
            hco3 = K2 * co2 / hp
            co3 = K2 * K3 * co2 / hp / hp
            hs = K5 * h2s / hp
            oh = K4 / hp
            amh = amine_M * hp / (K1p + hp)
            return hp + amh - oh - hco3 - 2.0 * co3 - hs

        h = brentq(charge, 1e-14, 1e-1, xtol=1e-22, rtol=1e-13, maxiter=200)
        new_free = amine_M * K1p / (K1p + h)
        if abs(new_free - free_amine) <= 1e-10 * amine_M:
            free_amine = new_free
            break
        free_amine = new_free
    return {
        "H+": h, "OH-": K4 / h,
        "CO2": co2, "HCO3-": K2 * co2 / h, "CO3--": K2 * K3 * co2 / h / h,
        "H2S": h2s, "HS-": K5 * h2s / h,
        "amine": free_amine, "amineH+": amine_M - free_amine,
    }


def acid_gas_loadings(T: float, P_co2: float, P_h2s: float, amine_M: float,
                      amine: str | AmineSystem = "MDEA") -> tuple[float, float]:
    """``(alpha_CO2, alpha_H2S)`` loadings (mol acid gas / mol amine) of an
    aqueous amine in equilibrium with CO2 and H2S partial pressures (kPa)."""
    sp = speciate(T, P_co2, P_h2s, amine_M, amine)
    a_co2 = (sp["CO2"] + sp["HCO3-"] + sp["CO3--"]) / amine_M
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
    (``P_h2s`` in **kPa**), no CO2 present."""
    return acid_gas_loadings(T, 0.0, P_h2s, amine_M, amine)[1]


def _invert(load_fn, T: float, alpha: float, amine_M: float,
            amine: str | AmineSystem) -> float:
    """Bracketed inverse of a (monotone) loading-vs-partial-pressure relation,
    returning the partial pressure (kPa) that gives ``alpha``."""
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
