"""M16 test: the Kent-Eisenberg acid-gas solubility model (CO2 in aqueous MDEA).

Reference & validation oracle
-----------------------------
Gas sweetening (Hameed 2025, *Chemical Process Simulations using Aspen HYSYS*,
§15.3) needs a reactive acid-gas property method — HYSYS uses its "Acid Gas –
Chemical Solvents" package, and the book reports no numeric slate. The
validatable physics is the underlying VLE: the equilibrium CO2 loading of an
aqueous amine vs the CO2 partial pressure. We validate against the published
experimental data of **Haji-Sulaiman, Aroua & Benamor (*Chem. Eng. Res. Des.*
76:961, 1998)** for 2 M and 4 M aqueous MDEA at 303/313/323 K, as tabulated by
the modified-Kent-Eisenberg study (UTP FYP CHE 12506, Table 5). Their modified
Kent-Eisenberg model reproduces this data to 12.2 % AAD (2 M) and 29.7 % (4 M);
``caldyr.thermo.amine`` reproduces it to comparable-or-better accuracy.

The carbonate/water/Henry constants are the literature values (Edwards et al.
1978; Little et al. 1990) — checked here against their known 25 C magnitudes;
the amine-protonation constant is fitted to the data above.
"""
import math

import pytest

from caldyr.thermo.amine import (
    MDEA,
    co2_loading,
    co2_partial_pressure,
)
from caldyr.thermo.amine import _arr, _H_CO2, _K2, _K3, _K4

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


def test_literature_constants_have_physical_magnitudes():
    """The carbonate/water/Henry constants must reproduce their textbook 25 C
    values — the guard that caught an OCR-corrupted constant during sourcing."""
    T = 298.15
    assert _arr(_K2, T) == pytest.approx(4.4e-7, rel=0.1)   # CO2 ionization
    assert _arr(_K3, T) == pytest.approx(4.7e-11, rel=0.1)  # bicarbonate
    assert _arr(_K4, T) == pytest.approx(1.0e-14, rel=0.1)  # Kw
    assert _arr(_H_CO2, T) == pytest.approx(29.0, rel=0.15)  # Henry, atm.L/mol


def test_amine_protonation_constant_is_physical():
    """The fitted MDEA protonation constant must land near the known
    pKa(MDEAH+) ~ 8.5 — i.e. the fit is chemistry, not an arbitrary curve."""
    K1_298 = math.exp(MDEA.k1_a + MDEA.k1_b / 298.15)
    pKa = -math.log10(K1_298)
    assert 8.0 < pKa < 9.3


def test_reproduces_haji_sulaiman_2M_mdea():
    """2 M MDEA: reproduce the experimental CO2 loading to <= the published
    modified-Kent-Eisenberg accuracy (12.2 % AAD)."""
    assert _aad(2) < 12.0


def test_reproduces_haji_sulaiman_4M_mdea():
    """4 M MDEA: reproduce the experimental CO2 loading to <= the published
    modified-Kent-Eisenberg accuracy (29.7 % AAD)."""
    assert _aad(4) < 25.0


def test_loading_is_monotone_in_pressure_and_temperature():
    """The physics a sweetening unit relies on: at fixed T, loading rises with
    CO2 partial pressure (absorption); at fixed P, loading falls as T rises
    (the basis of thermal regeneration in the stripper)."""
    loads_P = [co2_loading(313.15, P, 2.0) for P in (1.0, 5.0, 20.0, 100.0)]
    assert all(b > a for a, b in zip(loads_P, loads_P[1:]))
    loads_T = [co2_loading(T, 20.0, 2.0) for T in (303.15, 323.15, 343.15)]
    assert all(b < a for a, b in zip(loads_T, loads_T[1:]))


def test_partial_pressure_inverts_loading():
    """co2_partial_pressure is the inverse of co2_loading (the absorber stage
    needs the vapour driving force at a given liquid loading)."""
    for T in (303.15, 323.15):
        for P in (2.0, 25.0, 90.0):
            alpha = co2_loading(T, P, 3.0)
            P_back = co2_partial_pressure(T, alpha, 3.0)
            assert P_back == pytest.approx(P, rel=1e-4)
