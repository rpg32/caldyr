"""M16 test: Kent-Eisenberg acid-gas solubility (CO2 and H2S in aqueous MDEA).

Reference & validation oracle
-----------------------------
Gas sweetening (Hameed 2025, *Chemical Process Simulations using Aspen HYSYS*,
§15.3) needs a reactive acid-gas property method — HYSYS uses its "Acid Gas –
Chemical Solvents" package, and the book reports no numeric slate. The
validatable physics is the underlying VLE: the equilibrium acid-gas loading of
an aqueous amine vs the acid-gas partial pressure.

* **CO2-MDEA** is validated against the published experimental data of
  Haji-Sulaiman, Aroua & Benamor (*Chem. Eng. Res. Des.* 76:961, 1998) for 2 M
  and 4 M aqueous MDEA at 303/313/323 K, as tabulated by the modified-Kent-
  Eisenberg study UTP FYP CHE 12506 (Table 5). Their model reproduces this to
  12.2 % AAD (2 M) / 29.7 % (4 M); ``caldyr.thermo.amine`` does comparably-or-
  better.
* **H2S-MDEA** is a *prediction* (the amine protonation constant comes from the
  CO2 fit, combined with the literature H2S ionization + Henry constants), so it
  is checked for physical correctness — the constants against their known
  magnitudes, and the loading isotherm for shape/monotonicity — rather than a
  quantitative data fit. (A quantitative check against Jou-Mather-Otto 1982
  H2S-MDEA data is a tracked follow-up; that data sits behind paywalls.)

The carbonate / sulfide / water / Henry constants are the literature values
(Edwards et al. 1978; Little et al. 1990; Sander 2015), verified here against
their known 25 C magnitudes; the amine-protonation constant is fitted to the
CO2 data above.
"""
import math

import pytest

from caldyr.thermo.amine import (
    MDEA,
    _arr,
    _H_CO2,
    _H_H2S,
    _K2,
    _K3,
    _K4,
    _K5,
    acid_gas_loadings,
    co2_loading,
    co2_partial_pressure,
    h2s_loading,
    h2s_partial_pressure,
)

# Haji-Sulaiman et al. (1998) CO2-MDEA solubility: (molarity, T[K], P_CO2[kPa],
# alpha_exp [mol CO2 / mol MDEA]).
_HS_DATA = [
    (2, 303, 1.064, 0.114), (2, 303, 3.130, 0.244), (2, 303, 4.802, 0.333),
    (2, 303, 10.535, 0.483), (2, 303, 29.756, 0.673), (2, 303, 48.370, 0.793),
    (2, 303, 95.830, 0.880),
    (2, 313, 1.064, 0.103), (2, 313, 3.069, 0.197), (2, 313, 5.176, 0.267),
    (2, 313, 10.029, 0.374), (2, 313, 30.349, 0.603), (2, 313, 47.520, 0.688),
    (2, 313, 93.956, 0.805),
    (2, 323, 0.997, 0.079), (2, 323, 2.938, 0.148), (2, 323, 4.761, 0.194),
    (2, 323, 9.725, 0.298), (2, 323, 28.435, 0.471), (2, 323, 44.136, 0.590),
    (2, 323, 91.514, 0.726),
    (4, 303, 0.099, 0.027), (4, 303, 0.984, 0.061), (4, 303, 4.918, 0.149),
    (4, 303, 9.853, 0.284), (4, 303, 29.509, 0.516), (4, 303, 49.100, 0.633),
    (4, 303, 98.200, 0.761),
    (4, 313, 0.095, 0.015), (4, 313, 0.954, 0.052), (4, 313, 4.762, 0.086),
    (4, 313, 9.523, 0.190), (4, 313, 28.521, 0.384), (4, 313, 47.535, 0.513),
    (4, 313, 95.234, 0.654),
    (4, 323, 0.090, 0.010), (4, 323, 0.901, 0.037), (4, 323, 4.514, 0.084),
    (4, 323, 9.028, 0.151), (4, 323, 27.084, 0.251), (4, 323, 45.139, 0.363),
    (4, 323, 90.279, 0.516),
]


def _aad(molarity: int) -> float:
    errs = [abs(co2_loading(T, P, m) - a) / a
            for m, T, P, a in _HS_DATA if m == molarity]
    return 100.0 * sum(errs) / len(errs)


# -- literature constants: physical-magnitude guards (caught an OCR error) ----
def test_carbonate_water_henry_constants_have_physical_magnitudes():
    T = 298.15
    assert _arr(_K2, T) == pytest.approx(4.4e-7, rel=0.1)   # CO2 ionization
    assert _arr(_K3, T) == pytest.approx(4.7e-11, rel=0.1)  # bicarbonate
    assert _arr(_K4, T) == pytest.approx(1.0e-14, rel=0.1)  # Kw
    assert _arr(_H_CO2, T) == pytest.approx(29.0, rel=0.15)  # Henry, atm.L/mol


def test_h2s_constants_have_physical_magnitudes():
    """H2S first ionization must give pKa1 ~ 7.0, and Henry must reproduce the
    pure-water H2S solubility (~0.1 mol/L at 1 atm, 25 C; less soluble hotter)."""
    assert -math.log10(_arr(_K5, 298.15)) == pytest.approx(7.0, abs=0.1)
    assert 1.0 / _arr(_H_H2S, 298.15) == pytest.approx(0.10, rel=0.1)  # mol/L@1atm
    assert _arr(_H_H2S, 373.15) > _arr(_H_H2S, 298.15)   # less soluble when hot


def test_amine_protonation_constant_is_physical():
    """The fitted MDEA protonation constant must land near the known
    pKa(MDEAH+) ~ 8.5 — i.e. the fit is chemistry, not an arbitrary curve."""
    pKa = -math.log10(math.exp(MDEA.k1_a + MDEA.k1_b / 298.15))
    assert 8.0 < pKa < 9.3


# -- CO2-MDEA: validated against Haji-Sulaiman (1998) experimental data -------
def test_reproduces_haji_sulaiman_2M_mdea():
    assert _aad(2) < 10.0     # published modified-KE: 12.2 %


def test_reproduces_haji_sulaiman_4M_mdea():
    assert _aad(4) < 22.0     # published modified-KE: 29.7 %


# -- H2S-MDEA: physical-behaviour checks (prediction, not a data fit) ---------
def test_h2s_loading_isotherm_is_sensible():
    """H2S loading must rise monotonically with P_H2S and reach a physically
    reasonable capacity (~0.5-1 mol/mol at ~100 kPa for 2 M MDEA)."""
    loads = [h2s_loading(313.15, P, 2.0) for P in (1.0, 5.0, 20.0, 100.0)]
    assert all(b > a for a, b in zip(loads, loads[1:]))
    assert 0.5 < loads[-1] < 1.2


# -- physics shared by absorption / regeneration -----------------------------
def test_loading_falls_with_temperature():
    """At fixed partial pressure both acid gases desorb as T rises — the basis
    of thermal regeneration in the stripper."""
    for fn in (co2_loading, h2s_loading):
        loads = [fn(T, 20.0, 2.0) for T in (303.15, 323.15, 343.15)]
        assert all(b < a for a, b in zip(loads, loads[1:]))


def test_partial_pressure_inverts_loading():
    for T in (303.15, 323.15):
        for P in (2.0, 25.0, 90.0):
            assert co2_partial_pressure(T, co2_loading(T, P, 3.0), 3.0) \
                == pytest.approx(P, rel=1e-4)
            assert h2s_partial_pressure(T, h2s_loading(T, P, 3.0), 3.0) \
                == pytest.approx(P, rel=1e-4)


def test_combined_co2_h2s_speciation_is_consistent():
    """With both gases present the single charge balance still closes and each
    loading lies between its single-gas value and zero (competition for the
    amine), so the ternary solve is self-consistent."""
    a_co2, a_h2s = acid_gas_loadings(313.15, 20.0, 20.0, 3.0)
    assert a_co2 > 0.0 and a_h2s > 0.0
    assert a_co2 < co2_loading(313.15, 20.0, 3.0)   # CO2 partly displaced
    assert a_h2s < h2s_loading(313.15, 20.0, 3.0)   # H2S partly displaced
