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
1998). The liquid-phase ionic equilibria are solved with *literature*
equilibrium constants for the carbonate/water reactions and a *fitted* apparent
constant for the amine protonation that lumps the solution non-idealities (the
defining idea of the method: force the fit through the one reaction whose
constant is hardest to know a-priori).

Scope today: **CO2 in aqueous MDEA** (methyldiethanolamine — a *tertiary* amine,
so no carbamate forms and the reaction set is the minimal one). Validated
against the Haji-Sulaiman et al. (1998) 2 M and 4 M MDEA solubility data in
``tests/test_m16_amine_kent_eisenberg.py`` (reproduces it to better than the
published model's accuracy). The structure (an ``AmineSystem`` carrying the
amine-specific constants, a general ionic-equilibrium solver) is built to extend
to H2S and to carbamate-forming primary/secondary amines (MEA, DEA).

Reaction set (concentration / molarity basis, mol/L):
  * amine protonation     ``AmH+  <->  Am + H+``                    (K1, fitted)
  * CO2 first ionization  ``CO2 + H2O  <->  H+ + HCO3-``            (K2)
  * bicarbonate           ``HCO3-  <->  H+ + CO3^2-``               (K3)
  * water                 ``H2O  <->  H+ + OH-``                    (K4 = Kw)
  * Henry (physical)      ``P_CO2 = H_CO2 * [CO2]``
Unknowns are eliminated to a single charge-balance equation in ``[H+]``, solved
by a bracketed root find; the CO2 loading ``alpha = (mol CO2)/(mol amine)``
follows from the carbon balance.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

ATM_PER_KPA = 1.0 / 101.325

# -- temperature-dependent constants: value = exp(a/T + b ln T + c T + d), T[K].
# K2,K3,K4 are molarity-basis (mol/L) acid-dissociation constants; H_CO2 is the
# CO2 Henry's constant in atm.L/mol. These four are the canonical literature
# values (Edwards, Maurer, Newman & Prausnitz, *AIChE J.* 24:966, 1978; CO2
# Henry from Little, Bos & Knoop, 1990) — verified here against their known
# 25 C magnitudes (K2~4.4e-7, K3~4.7e-11, Kw~1e-14, H~29 atm.L/mol).
_K2 = (-12092.1, -36.7816, 0.0, 235.482)    # CO2 + H2O <-> H+ + HCO3-
_K3 = (-12431.7, -35.4819, 0.0, 220.067)    # HCO3- <-> H+ + CO3^2-
_K4 = (-13445.9, -22.4773, 0.0, 140.932)    # H2O <-> H+ + OH-  (Kw)
_H_CO2 = (-6789.04, -11.4519, -0.010454, 94.4914)   # Henry, atm.L/mol


def _arr(p: tuple[float, float, float, float], T: float) -> float:
    a, b, c, d = p
    return math.exp(a / T + b * math.log(T) + c * T + d)


@dataclass(frozen=True)
class AmineSystem:
    """Acid-gas equilibrium parameters for one aqueous amine.

    ``k1_a``/``k1_b`` give the amine protonation constant ``K1 = exp(k1_a +
    k1_b/T)`` (acid dissociation of the protonated amine, mol/L). ``f_g``/``f_k``
    are the Haji-Sulaiman apparent-constant factor exponents: the *apparent*
    protonation constant used in the balance is ``K1' = K1 * P_CO2**f_g *
    [amine]**f_k`` (``P_CO2`` in atm, ``[amine]`` the free-amine molarity), which
    lets the lumped constant carry the partial-pressure and loading dependence.
    ``carbamate`` flags primary/secondary amines (MEA, DEA) that also form a
    carbamate ion — not modelled yet, so only tertiary amines are exposed.
    """
    name: str
    k1_a: float
    k1_b: float
    f_g: float
    f_k: float
    carbamate: bool = False


# K1 fitted to the Haji-Sulaiman et al. (1998) 2 M + 4 M MDEA CO2 solubility
# data (303-323 K, 0.1-100 kPa) with the carbonate/water/Henry constants above
# held at their literature values; f_g, f_k are the published MDEA values. The
# fit gives K1(298 K)=1.6e-9 (pKa(MDEAH+)~8.8, physically consistent with the
# ~8.5 literature value) and reproduces the data to ~6% (2 M) / ~18% (4 M) AAD —
# see tests/test_m16_amine_kent_eisenberg.py.
MDEA = AmineSystem(name="MDEA", k1_a=-9.4553, k1_b=-3221.39,
                   f_g=-0.03628, f_k=0.6262, carbamate=False)

_SYSTEMS = {"MDEA": MDEA}


def amine_system(name: str) -> AmineSystem:
    """Look up a supported amine by name (case-insensitive)."""
    try:
        return _SYSTEMS[name.upper()]
    except KeyError:
        raise ValueError(
            f"unsupported amine {name!r}; Kent-Eisenberg acid-gas solubility is "
            f"available for {sorted(_SYSTEMS)} (CO2 only) — primary/secondary "
            f"amines (MEA, DEA) need the carbamate reaction, not yet modelled"
        ) from None


def _speciate(T: float, P_co2_atm: float, amine_M: float,
              sys: AmineSystem) -> dict[str, float]:
    """Solve the liquid-phase ionic equilibrium at (T, P_CO2, amine molarity).

    Returns the molar concentrations (mol/L) of every species plus ``H+``. The
    apparent protonation constant depends on the free-amine concentration, so a
    short fixed-point loop wraps the ``[H+]`` root find.
    """
    if amine_M <= 0:
        raise ValueError("amine molarity must be positive")
    K2, K3, K4, H = (_arr(p, T) for p in (_K2, _K3, _K4, _H_CO2))
    co2 = max(P_co2_atm, 0.0) / H        # free molecular CO2, Henry's law
    free_amine = amine_M                  # initial guess for the F-factor
    h = 1e-9
    for _ in range(50):
        K1 = math.exp(sys.k1_a + sys.k1_b / T)
        f1 = (max(P_co2_atm, 1e-12) ** sys.f_g) * (max(free_amine, 1e-12) ** sys.f_k)
        K1p = K1 * f1

        def charge(hp: float) -> float:
            # [H+] + [AmH+] - [OH-] - [HCO3-] - 2[CO3^2-] = 0
            hco3 = K2 * co2 / hp
            co3 = K2 * K3 * co2 / hp / hp
            oh = K4 / hp
            amh = amine_M * hp / (K1p + hp)
            return hp + amh - oh - hco3 - 2.0 * co3

        h = brentq(charge, 1e-14, 1e-1, xtol=1e-22, rtol=1e-13, maxiter=200)
        new_free = amine_M * K1p / (K1p + h)
        if abs(new_free - free_amine) <= 1e-10 * amine_M:
            free_amine = new_free
            break
        free_amine = new_free
    hco3 = K2 * co2 / h
    co3 = K2 * K3 * co2 / h / h
    return {
        "H+": h, "OH-": K4 / h,
        "CO2": co2, "HCO3-": hco3, "CO3--": co3,
        "amine": free_amine, "amineH+": amine_M - free_amine,
    }


def co2_loading(T: float, P_co2: float, amine_M: float,
                amine: str | AmineSystem = "MDEA") -> float:
    """CO2 loading ``alpha`` (mol CO2 / mol amine) of an aqueous amine solution
    in equilibrium with a CO2 partial pressure.

    ``T`` in K, ``P_co2`` in **kPa**, ``amine_M`` the amine molarity (mol/L).
    """
    sys = amine if isinstance(amine, AmineSystem) else amine_system(amine)
    sp = _speciate(T, P_co2 * ATM_PER_KPA, amine_M, sys)
    return (sp["CO2"] + sp["HCO3-"] + sp["CO3--"]) / amine_M


def co2_partial_pressure(T: float, alpha: float, amine_M: float,
                         amine: str | AmineSystem = "MDEA") -> float:
    """Inverse of :func:`co2_loading`: the equilibrium CO2 partial pressure
    (**kPa**) over an aqueous amine solution carrying loading ``alpha`` at ``T``.

    This is the building block the absorber/regenerator stage needs — the
    vapour-side CO2 driving force at a given liquid loading. Solved by a
    bracketed root find on the (monotone) loading-vs-pressure relation.
    """
    sys = amine if isinstance(amine, AmineSystem) else amine_system(amine)
    if alpha <= 0.0:
        return 0.0

    def resid(logp: float) -> float:
        return co2_loading(T, math.exp(logp), amine_M, sys) - alpha

    # loading is strictly increasing in P_CO2; bracket on log P (kPa)
    lo, hi = math.log(1e-4), math.log(5.0e4)
    if resid(hi) < 0.0:           # loading beyond the bracket's reach
        return math.exp(hi)
    if resid(lo) > 0.0:
        return math.exp(lo)
    return math.exp(brentq(resid, lo, hi, xtol=1e-8, rtol=1e-10, maxiter=200))
