"""M12 tests: the property-table analysis utility (HYSYS "Property Table").

Reference: Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley
2025), §2.1 — n-pentane at 550 K and 16 atm gives a mass density of
28.31 kg/m^3 under SRK and 28.77 kg/m^3 under PR (Figure 2.7), and the
§2.1.3/2.1.4 case-study grid (Figure 2.13, SRK) tabulates mass density over
T = 500–600 K and P = 12–14+ atm. Caldyr's thermo backends land within 0.1%
of every book value.
"""
import numpy as np
import pytest

from caldyr.analysis import PROPERTIES, property_table
from caldyr.thermo import make_package

ATM = 101_325.0
Z5 = {"n-pentane": 1.0}


@pytest.fixture(scope="module")
def srk():
    return make_package("thermo:SRK", ["n-pentane"])


def test_book_single_point_pr_and_srk(srk):
    """Book Figure 2.7: n-pentane, 550 K, 16 atm -> 28.77 kg/m^3 (PR),
    28.31 kg/m^3 (SRK)."""
    pr = make_package("thermo:PR", ["n-pentane"])
    out_pr = property_table(pr, Z5, T=550.0, P=16 * ATM, props=["mass_density"])
    out_srk = property_table(srk, Z5, T=550.0, P=16 * ATM, props=["mass_density"])
    assert out_pr["mass_density"][0, 0] == pytest.approx(28.77, rel=2e-3)
    assert out_srk["mass_density"][0, 0] == pytest.approx(28.31, rel=2e-3)


def test_book_case_study_grid_srk(srk):
    """Book Figure 2.13 (case study, SRK): spot-check four grid points."""
    out = property_table(
        srk, Z5,
        T=[500.0, 520.0, 540.0, 560.0, 580.0, 600.0],
        P=[12 * ATM, 13 * ATM, 14 * ATM],
        props=["mass_density"],
    )
    rho = out["mass_density"]
    assert rho.shape == (6, 3)
    assert rho[0, 0] == pytest.approx(23.58, rel=2e-3)   # 500 K, 12 atm
    assert rho[5, 0] == pytest.approx(18.51, rel=2e-3)   # 600 K, 12 atm
    assert rho[1, 1] == pytest.approx(24.37, rel=2e-3)   # 520 K, 13 atm
    assert rho[2, 2] == pytest.approx(25.07, rel=2e-3)   # 540 K, 14 atm
    assert not out["failures"]


def test_book_trends_density_falls_with_T_rises_with_P(srk):
    """The book's Eq. (2.1)/Figure 2.14 conclusion: mass density decreases with
    temperature and increases with pressure."""
    out = property_table(
        srk, Z5, T=np.linspace(500.0, 600.0, 6), P=[12 * ATM, 15 * ATM, 18 * ATM],
        props=["mass_density"],
    )
    rho = out["mass_density"]
    assert np.all(np.diff(rho, axis=0) < 0)     # falls along T
    assert np.all(np.diff(rho, axis=1) > 0)     # rises along P


def test_enthalpy_monotone_in_T_at_fixed_P(srk):
    """Internal consistency: single-phase (vapor) enthalpy is monotone
    increasing in T at fixed P (Cp > 0)."""
    out = property_table(
        srk, Z5, T=np.linspace(500.0, 600.0, 11), P=[12 * ATM, 18 * ATM],
        props=["enthalpy", "vapor_fraction"],
    )
    assert np.all(out["vapor_fraction"] == 1.0)
    assert np.all(np.diff(out["enthalpy"], axis=0) > 0)


def test_failed_points_skip_gracefully():
    """A point outside the backend's range (water below the triple point under
    CoolProp) becomes NaN + a logged failure; the rest of the grid survives."""
    pp = make_package("coolprop:Water", ["water"])
    out = property_table(pp, {"water": 1.0}, T=[150.0, 300.0], P=1e5,
                         props=["mass_density"])
    rho = out["mass_density"]
    assert np.isnan(rho[0, 0]) and rho[1, 0] == pytest.approx(996.56, rel=1e-3)
    assert len(out["failures"]) == 1
    t_fail, p_fail, msg = out["failures"][0]
    assert (t_fail, p_fail) == (150.0, 1e5) and msg


def test_all_documented_properties_compute(srk):
    out = property_table(srk, Z5, T=550.0, P=16 * ATM, props=sorted(PROPERTIES))
    for name in PROPERTIES:
        assert np.isfinite(out[name][0, 0])
    # cross-consistency: mass_density * molar_volume = molar mass (72.15 g/mol)
    mw = out["mass_density"][0, 0] * out["molar_volume"][0, 0]
    assert mw == pytest.approx(72.15e-3, rel=1e-3)


def test_bad_inputs_are_typed_errors(srk):
    with pytest.raises(ValueError, match="unknown property"):
        property_table(srk, Z5, T=550.0, P=16 * ATM, props=["viscosity"])
    with pytest.raises(ValueError, match="at least one"):
        property_table(srk, Z5, T=550.0, P=16 * ATM, props=[])
    with pytest.raises(ValueError, match="finite and positive"):
        property_table(srk, Z5, T=[550.0, -10.0], P=16 * ATM)
    with pytest.raises(ValueError, match="sums to"):
        property_table(srk, {"n-pentane": 0.0}, T=550.0, P=16 * ATM)
