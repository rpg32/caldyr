"""M16 test: Kent-Eisenberg acid-gas solubility (CO2/H2S in aqueous DEA & MDEA).

Reference & validation oracle
-----------------------------
Gas sweetening (Hameed 2025, *Aspen HYSYS*, §15.3) needs a reactive acid-gas
property method — HYSYS uses "Acid Gas – Chemical Solvents" and the book gives
no numeric slate, so the oracle is the underlying VLE: equilibrium acid-gas
loading vs partial pressure.

* **CO2 in DEA & MDEA** — validated against the experimental data of
  Haji-Sulaiman, Aroua & Benamor (*Chem. Eng. Res. Des.* 76:961, 1998), Table 2,
  for 2 M / 4 M solutions at 303/313/323 K. ``caldyr.thermo.amine`` reproduces
  DEA to ~7% AAD (the paper's own model: 9.2%) and MDEA comparably. The DEA
  carbamate (formation constant from Aroua et al. 1997) caps the fast loading
  near 0.5 mol/mol, with bicarbonate carrying it higher — the signature of a
  secondary amine.
* **H2S in MDEA** — a prediction (amine protonation from the CO2 fit + literature
  H2S ionization/Henry constants, verified vs pure-water solubility); checked for
  physical correctness. (Quantitative H2S data fit vs Jou 1982 / Lawson 1976 in
  caldyr/docs/ is a tracked follow-up.)
"""
import math

import pytest

from caldyr.thermo.amine import (
    DEA,
    MDEA,
    _arr,
    _H_CO2,
    _H_H2S,
    _K_CO2,
    _K_H2S,
    _K_HCO3,
    _K_W,
    acid_gas_loadings,
    co2_loading,
    co2_partial_pressure,
    h2s_loading,
    h2s_partial_pressure,
)

# Haji-Sulaiman et al. (1998) Table 2: (amine, molarity, T[K], P_CO2[kPa], alpha).
_CO2_DATA = [
    # --- DEA ---
    ("DEA", 2, 303, 0.098, 0.183), ("DEA", 2, 303, 0.492, 0.325),
    ("DEA", 2, 303, 1.119, 0.388), ("DEA", 2, 303, 5.355, 0.521),
    ("DEA", 2, 303, 10.726, 0.591), ("DEA", 2, 303, 32.527, 0.699),
    ("DEA", 2, 303, 54.213, 0.730),
    ("DEA", 2, 313, 0.095, 0.172), ("DEA", 2, 313, 0.474, 0.278),
    ("DEA", 2, 313, 1.039, 0.320), ("DEA", 2, 313, 5.265, 0.459),
    ("DEA", 2, 313, 10.665, 0.538), ("DEA", 2, 313, 32.147, 0.597),
    ("DEA", 2, 313, 53.829, 0.662),
    ("DEA", 2, 323, 0.090, 0.133), ("DEA", 2, 323, 0.449, 0.152),
    ("DEA", 2, 323, 1.040, 0.272), ("DEA", 2, 323, 5.110, 0.398),
    ("DEA", 2, 323, 10.035, 0.473), ("DEA", 2, 323, 30.358, 0.546),
    ("DEA", 2, 323, 50.763, 0.611), ("DEA", 2, 323, 98.170, 0.688),
    ("DEA", 4, 303, 0.986, 0.309), ("DEA", 4, 303, 4.893, 0.471),
    ("DEA", 4, 303, 9.863, 0.524), ("DEA", 4, 303, 29.358, 0.588),
    ("DEA", 4, 303, 48.931, 0.633), ("DEA", 4, 303, 98.628, 0.671),
    ("DEA", 4, 313, 0.951, 0.281), ("DEA", 4, 313, 5.259, 0.441),
    ("DEA", 4, 313, 10.413, 0.499), ("DEA", 4, 313, 30.987, 0.561),
    ("DEA", 4, 313, 52.568, 0.599), ("DEA", 4, 313, 102.119, 0.639),
    ("DEA", 4, 323, 0.903, 0.193), ("DEA", 4, 323, 4.514, 0.344),
    ("DEA", 4, 323, 9.028, 0.445), ("DEA", 4, 323, 27.032, 0.498),
    ("DEA", 4, 323, 46.305, 0.517), ("DEA", 4, 323, 98.673, 0.601),
    # --- MDEA ---
    ("MDEA", 2, 303, 1.064, 0.114), ("MDEA", 2, 303, 3.130, 0.244),
    ("MDEA", 2, 303, 4.802, 0.333), ("MDEA", 2, 303, 10.535, 0.483),
    ("MDEA", 2, 303, 29.756, 0.673), ("MDEA", 2, 303, 48.370, 0.793),
    ("MDEA", 2, 303, 95.830, 0.880),
    ("MDEA", 2, 313, 1.064, 0.103), ("MDEA", 2, 313, 3.069, 0.197),
    ("MDEA", 2, 313, 5.176, 0.267), ("MDEA", 2, 313, 10.029, 0.374),
    ("MDEA", 2, 313, 30.349, 0.603), ("MDEA", 2, 313, 47.520, 0.688),
    ("MDEA", 2, 313, 93.956, 0.805),
    ("MDEA", 2, 323, 0.997, 0.079), ("MDEA", 2, 323, 2.938, 0.148),
    ("MDEA", 2, 323, 4.761, 0.194), ("MDEA", 2, 323, 9.725, 0.298),
    ("MDEA", 2, 323, 28.435, 0.471), ("MDEA", 2, 323, 44.136, 0.590),
    ("MDEA", 2, 323, 91.514, 0.726),
    ("MDEA", 4, 303, 0.984, 0.061), ("MDEA", 4, 303, 4.918, 0.149),
    ("MDEA", 4, 303, 9.853, 0.284), ("MDEA", 4, 303, 29.509, 0.516),
    ("MDEA", 4, 303, 49.100, 0.633), ("MDEA", 4, 303, 98.200, 0.761),
    ("MDEA", 4, 313, 0.954, 0.052), ("MDEA", 4, 313, 4.762, 0.086),
    ("MDEA", 4, 313, 9.523, 0.190), ("MDEA", 4, 313, 28.521, 0.384),
    ("MDEA", 4, 313, 47.535, 0.513), ("MDEA", 4, 313, 95.234, 0.654),
    ("MDEA", 4, 323, 0.901, 0.037), ("MDEA", 4, 323, 4.514, 0.084),
    ("MDEA", 4, 323, 9.028, 0.151), ("MDEA", 4, 323, 27.084, 0.251),
    ("MDEA", 4, 323, 45.139, 0.363), ("MDEA", 4, 323, 90.279, 0.516),
]


def _aad(amine: str) -> float:
    errs = [abs(co2_loading(T, P, m, amine) - a) / a
            for nm, m, T, P, a in _CO2_DATA if nm == amine]
    return 100.0 * sum(errs) / len(errs)


# -- literature constants: physical-magnitude guards (caught an OCR error) ----
def test_carbonate_water_henry_constants_have_physical_magnitudes():
    T = 298.15
    assert _arr(_K_CO2, T) == pytest.approx(4.4e-7, rel=0.1)   # CO2 ionization
    assert _arr(_K_HCO3, T) == pytest.approx(4.7e-11, rel=0.1)  # bicarbonate
    assert _arr(_K_W, T) == pytest.approx(1.0e-14, rel=0.1)     # Kw
    assert _arr(_H_CO2, T) == pytest.approx(29.0, rel=0.15)     # Henry, atm.L/mol


def test_h2s_constants_have_physical_magnitudes():
    assert -math.log10(_arr(_K_H2S, 298.15)) == pytest.approx(7.0, abs=0.1)
    assert 1.0 / _arr(_H_H2S, 298.15) == pytest.approx(0.10, rel=0.1)  # mol/L@1atm
    assert _arr(_H_H2S, 373.15) > _arr(_H_H2S, 298.15)


def test_amine_protonation_constants_are_physical():
    """Fitted/literature protonation constants near the known pKa's: DEAH+ ~8.88,
    MDEAH+ ~8.5 (the MDEA value whose OCR-corrupted form this project caught)."""
    pKa_dea = -math.log10(math.exp(DEA.prot[0] / 298.15 + DEA.prot[1]
                                   * math.log(298.15) + DEA.prot[3]))
    pKa_mdea = -math.log10(math.exp(MDEA.prot[0] / 298.15 + MDEA.prot[1]
                                    * math.log(298.15) + MDEA.prot[3]))
    assert 8.5 < pKa_dea < 9.2
    assert 8.2 < pKa_mdea < 8.9


# -- CO2 in DEA & MDEA: validated against Haji-Sulaiman (1998) data -----------
def test_reproduces_dea_co2_solubility():
    """CO2 in aqueous DEA — secondary amine with carbamate. The paper's own
    modified-KE model reaches 9.2% AAD on this data."""
    assert _aad("DEA") < 10.0


def test_reproduces_mdea_co2_solubility():
    """CO2 in aqueous MDEA — tertiary amine, no carbamate."""
    assert _aad("MDEA") < 25.0


def test_dea_carbamate_caps_then_bicarbonate_exceeds():
    """The secondary-amine signature: the carbamate limits the fast CO2 uptake to
    ~0.5 mol/mol, but at high CO2 pressure bicarbonate carries loading past 0.5
    (Aroua et al. 1997); MDEA (no carbamate) shows no such 0.5 shoulder."""
    dea_lo = co2_loading(303.15, 1.0, 2.0, "DEA")    # low P: carbamate regime
    dea_hi = co2_loading(303.15, 100.0, 2.0, "DEA")  # high P: bicarbonate regime
    assert 0.30 < dea_lo < 0.55
    assert dea_hi > 0.6


# -- H2S in MDEA: physical-behaviour checks (prediction, not a data fit) ------
def test_h2s_loading_isotherm_is_sensible():
    loads = [h2s_loading(313.15, P, 2.0) for P in (1.0, 5.0, 20.0, 100.0)]
    assert all(b > a for a, b in zip(loads, loads[1:]))
    assert 0.5 < loads[-1] < 1.2


# -- physics shared by absorption / regeneration -----------------------------
def test_loading_falls_with_temperature():
    """Both acid gases desorb as T rises — the basis of thermal regeneration."""
    for fn, am in ((co2_loading, "DEA"), (co2_loading, "MDEA"),
                   (h2s_loading, "MDEA")):
        loads = [fn(T, 20.0, 2.0, am) for T in (303.15, 323.15, 343.15)]
        assert all(b < a for a, b in zip(loads, loads[1:]))


def test_partial_pressure_inverts_loading():
    for am in ("DEA", "MDEA"):
        for P in (2.0, 25.0, 90.0):
            assert co2_partial_pressure(313.15, co2_loading(313.15, P, 3.0, am),
                                        3.0, am) == pytest.approx(P, rel=1e-3)
    for P in (2.0, 25.0, 90.0):
        assert h2s_partial_pressure(313.15, h2s_loading(313.15, P, 3.0),
                                    3.0) == pytest.approx(P, rel=1e-3)


def test_combined_co2_h2s_competes_for_amine():
    a_co2, a_h2s = acid_gas_loadings(313.15, 20.0, 20.0, 3.0, "MDEA")
    assert a_co2 > 0.0 and a_h2s > 0.0
    assert a_co2 < co2_loading(313.15, 20.0, 3.0, "MDEA")
    assert a_h2s < h2s_loading(313.15, 20.0, 3.0, "MDEA")
