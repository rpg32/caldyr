"""M13 tests: API 520 / API 526 relief-valve orifice sizing
(``caldyr.analysis.relief``).

References:

* **Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley 2025),
  sec. 7.3 worked example** (Figs. 7.36-7.41): a PSV on a steam line, blocked
  outlet — relieving rate 3000 kg/h, relieving temperature 224.5 C, relieving
  pressure 4.950 barG (set 4.5 barG + 10%), total backpressure 0.450 barG,
  Kd = 0.975. HYSYS computes a required orifice of **10.17 cm^2**, selects the
  API 526 **K** orifice (11.86 cm^2) and reports **85.74% capacity used**.
  ``test_book_psv_steam`` reproduces that with Z from the IAPWS-95 steam
  tables (coolprop:Water) and the ideal-gas k of steam at the relieving
  temperature; achieved vs book (asserted with headroom):
    - required area: 10.13 cm^2, -0.4%  (asserted 1%; the residual is HYSYS's
      PR-derived Z/k for steam vs IAPWS/ideal-gas values)
    - orifice letter: K (exact), capacity used 85.5% vs 85.74%.
* **API Standard 520 Part I, 9th ed. (2014), sec. 5.6.3** — critical vapor
  flow, SI form A[mm^2] = W / (C Kd P1 Kb Kc) sqrt(TZ/M) with
  C = 0.03948 sqrt(k (2/(k+1))^((k+1)/(k-1))), W kg/h, P1 kPaa, M g/mol.
  The hand check below implements this published form verbatim; it is
  algebraically identical to the engine's pure-SI mass-flux form
  (0.03948 = 3.6/sqrt(1000 R)), so they must agree to the rounding of the
  0.03948 constant (~0.02%).
* **API 520 Part I sec. 5.8** — liquid sizing for capacity-certified valves,
  SI form A[mm^2] = 11.78 Q / (Kd Kw Kc Kv) sqrt(G / dP[kPa]) with Q in L/min;
  same equivalence argument (11.78 = 1e6/(6e4 sqrt(2)) to rounding).
* **API Standard 526, 7th ed. (2017)** — standard orifice letters D-T,
  effective areas 0.110-26.0 in^2.
"""
import math

import pytest

from caldyr.analysis import (
    API_526_ORIFICES,
    ReliefResult,
    ReliefSizingError,
    relief_liquid,
    relief_vapor,
)
from caldyr.thermo import make_package

ATM = 101_325.0
R = 8.314462618
IN2 = 6.4516e-4                 # m^2 / in^2

# Book sec. 7.3 relieving conditions.
W_BOOK = 3000.0 / 3600.0        # kg/s
T_BOOK = 224.5 + 273.15         # K
P1_BOOK = 4.950e5 + ATM         # Pa abs (4.950 barG)
P2_BOOK = 0.450e5 + ATM         # Pa abs (0.450 barG)
M_H2O = 0.01801528              # kg/mol


def api520_vapor_area_mm2(W_kgh, T, Z, M_gmol, k, P1_kPa,
                          Kd=0.975, Kb=1.0, Kc=1.0) -> float:
    """API 520 Part I (9e) sec. 5.6.3 critical-flow equation, SI field units,
    implemented verbatim as the independent hand check."""
    c = 0.03948 * math.sqrt(k * (2.0 / (k + 1.0)) ** ((k + 1.0) / (k - 1.0)))
    return W_kgh / (c * Kd * P1_kPa * Kb * Kc) * math.sqrt(T * Z / M_gmol)


# -- book sec. 7.3: steam PSV, blocked outlet ----------------------------------
@pytest.fixture(scope="module")
def steam_z_and_k():
    """Z of steam at the relieving conditions from IAPWS-95 (the engine's
    coolprop:Water package) and the ideal-gas Cp/Cv of steam at 497.65 K.
    Hand values: Z = 0.974 (superheated-steam tables: v ~ 0.376 m^3/kg at
    ~6 bar / 224.5 C), k = Cp0/(Cp0 - R) = 1.309 with Cp0 ~ 35.2 J/mol/K
    (Poling et al. 5e, App. A)."""
    pp = make_package("coolprop:Water", ["water"])
    vm = pp.volume(T_BOOK, P1_BOOK, {"water": 1.0})
    z_comp = P1_BOOK * vm / (R * T_BOOK)
    from CoolProp.CoolProp import PropsSI
    cp0 = PropsSI("Cp0molar", "T", T_BOOK, "P", P1_BOOK, "Water")
    k = cp0 / (cp0 - R)
    assert z_comp == pytest.approx(0.974, abs=0.002)
    assert k == pytest.approx(1.309, abs=0.005)
    return z_comp, k


def test_book_psv_steam(steam_z_and_k):
    z_comp, k = steam_z_and_k
    res = relief_vapor(W_BOOK, T_BOOK, M_H2O, z_comp, k, P1_BOOK,
                       backpressure=P2_BOOK, Kd=0.975)
    # Book (HYSYS): calculated orifice 10.17 cm^2 (achieved 10.13, -0.4%).
    assert res.area_m2 * 1e4 == pytest.approx(10.17, rel=1e-2)
    # Book: selected orifice K, 11.86 cm^2 (API 526: 1.838 in^2).
    assert res.orifice == "K"
    assert res.orifice_area_m2 == pytest.approx(1.838 * IN2, rel=1e-12)
    assert res.orifice_area_m2 * 1e4 == pytest.approx(11.86, rel=1e-3)
    # Book: capacity used 85.74% (achieved 85.5%).
    assert res.capacity_used == pytest.approx(0.8574, abs=0.01)
    # Blocked-outlet steam at 0.45 barG backpressure is well into critical
    # flow: P_crit = P1*(2/(k+1))^(k/(k-1)) ~ 0.544*P1 >> P2.
    assert res.critical is True
    assert res.phase == "vapor"
    assert res.details["P_critical_Pa"] == pytest.approx(
        P1_BOOK * (2.0 / (k + 1.0)) ** (k / (k - 1.0)), rel=1e-12)
    assert res.details["P_critical_Pa"] > P2_BOOK


def test_book_psv_matches_api520_field_form(steam_z_and_k):
    z_comp, k = steam_z_and_k
    res = relief_vapor(W_BOOK, T_BOOK, M_H2O, z_comp, k, P1_BOOK,
                       backpressure=P2_BOOK, Kd=0.975)
    hand = api520_vapor_area_mm2(3000.0, T_BOOK, z_comp, 18.01528, k,
                                 P1_BOOK / 1000.0)
    assert res.area_m2 * 1e6 == pytest.approx(hand, rel=5e-4)


def test_api520_hand_calc_hydrocarbon():
    """Independent hand point, API 520 9e sec. 5.6.3 worked by hand:
    a heavy hydrocarbon vapor, W = 24,270 kg/h, M = 51 g/mol, k = 1.11,
    Z = 0.90, T = 348 K, P1 = 670 kPaa, Kd = 0.975:
      C    = 0.03948*sqrt(1.11*(2/2.11)^(2.11/0.11)) = 0.024889
      A    = 24270/(0.024889*0.975*670)*sqrt(348*0.90/51) = 3699 mm^2
    -> API 526 'P' orifice (6.38 in^2 = 4116 mm^2)."""
    res = relief_vapor(24_270.0 / 3600.0, 348.0, 0.051, 0.90, 1.11, 670e3,
                       backpressure=ATM, Kd=0.975)
    assert res.area_m2 * 1e6 == pytest.approx(3699.0, rel=2e-3)
    assert res.orifice == "P"
    hand = api520_vapor_area_mm2(24_270.0, 348.0, 0.90, 51.0, 1.11, 670.0)
    assert res.area_m2 * 1e6 == pytest.approx(hand, rel=5e-4)


def test_correction_factors_scale_area_inversely():
    base = relief_vapor(W_BOOK, T_BOOK, M_H2O, 0.974, 1.309, P1_BOOK,
                        backpressure=P2_BOOK)
    with_rd = relief_vapor(W_BOOK, T_BOOK, M_H2O, 0.974, 1.309, P1_BOOK,
                           backpressure=P2_BOOK, Kc=0.9)
    assert with_rd.area_m2 == pytest.approx(base.area_m2 / 0.9, rel=1e-12)


# -- orifice-letter selection (API 526) ----------------------------------------
def test_orifice_table_is_api526():
    letters = [letter for letter, _ in API_526_ORIFICES]
    areas = [a for _, a in API_526_ORIFICES]
    assert letters == ["D", "E", "F", "G", "H", "J", "K", "L", "M", "N",
                       "P", "Q", "R", "T"]
    assert areas == sorted(areas)                  # strictly increasing
    assert areas[0] == 0.110 and areas[-1] == 26.0


def test_orifice_selection_boundaries():
    # Tiny requirement -> smallest (D); between E and F -> F; the selected
    # orifice is always >= the requirement.
    tiny = relief_vapor(0.01, 400.0, 0.044, 1.0, 1.2, 10e5, backpressure=ATM)
    assert tiny.orifice == "D"
    assert tiny.orifice_area_m2 >= tiny.area_m2

    res = ReliefResult(area_m2=0.25e-3, orifice=None, orifice_area_m2=None,
                       capacity_used=None, phase="vapor", critical=True)
    from caldyr.analysis.relief import _select_orifice
    picked = _select_orifice(0.25e-3, res)         # 2.5 cm^2: F is 1.98, G 3.24
    assert picked.orifice == "G"
    assert picked.capacity_used == pytest.approx(0.25e-3 / (0.503 * IN2))


def test_requirement_beyond_t_orifice_flags_multiple_valves():
    big = relief_vapor(300.0, 400.0, 0.018, 1.0, 1.3, 3e5, backpressure=ATM)
    assert big.area_m2 > 26.0 * IN2
    assert big.orifice is None and big.orifice_area_m2 is None
    assert any("multiple valves" in n for n in big.notes)


# -- liquid sizing (API 520 sec. 5.8) -------------------------------------------
def test_liquid_sizing_matches_api520_si_form():
    """1514 L/min (400 gpm) of a G = 0.90 liquid, set pressure 17.24 barG with
    10% accumulation (relieving 19.98 barA) against 4.46 barA backpressure
    (dP = 1551.3 kPa), Kd = 0.65. API 520 sec. 5.8 SI form by hand:
      A = 11.78 * 1514 / 0.65 * sqrt(0.90 / 1551.3) = 661 mm^2
    -> API 526 'J' orifice (1.287 in^2 = 830 mm^2)."""
    rho = 900.0                                    # G = 0.90
    q_lpm = 1514.0
    w = q_lpm / 60.0e3 * rho                       # kg/s
    p1 = 17.24e5 * 1.10 + ATM                      # set 17.24 barG + 10% accum.
    p2 = 3.45e5 + ATM
    dp_kpa = (p1 - p2) / 1000.0
    res = relief_liquid(w, rho, p1, p2)            # Kd = 0.65 default
    hand_mm2 = 11.78 * q_lpm / 0.65 * math.sqrt(0.90 / dp_kpa)
    assert res.area_m2 * 1e6 == pytest.approx(hand_mm2, rel=1e-3)
    assert res.area_m2 * 1e6 == pytest.approx(661.0, rel=5e-3)
    assert res.orifice == "J"
    assert res.phase == "liquid" and res.critical is None
    assert res.details["dP_Pa"] == pytest.approx(p1 - p2)


def test_liquid_correction_factors():
    base = relief_liquid(10.0, 800.0, 20e5, 5e5)
    visc = relief_liquid(10.0, 800.0, 20e5, 5e5, Kv=0.8, Kw=0.9)
    assert visc.area_m2 == pytest.approx(base.area_m2 / (0.8 * 0.9), rel=1e-12)


# -- typed errors ----------------------------------------------------------------
def test_subcritical_backpressure_raises():
    # k = 1.3 -> critical ratio 0.546; backpressure at 0.7*P1 is subcritical.
    with pytest.raises(ReliefSizingError, match="subcritical"):
        relief_vapor(1.0, 400.0, 0.018, 1.0, 1.3, 10e5, backpressure=7e5)
    assert issubclass(ReliefSizingError, ValueError)


def test_invalid_vapor_inputs_raise():
    with pytest.raises(ReliefSizingError, match="W"):
        relief_vapor(0.0, 400.0, 0.018, 1.0, 1.3, 10e5)
    with pytest.raises(ReliefSizingError, match="k"):
        relief_vapor(1.0, 400.0, 0.018, 1.0, 1.0, 10e5)
    with pytest.raises(ReliefSizingError, match="Kd"):
        relief_vapor(1.0, 400.0, 0.018, 1.0, 1.3, 10e5, Kd=-0.5)
    with pytest.raises(ReliefSizingError, match="backpressure"):
        relief_vapor(1.0, 400.0, 0.018, 1.0, 1.3, 10e5, backpressure=-1.0)


def test_invalid_liquid_inputs_raise():
    with pytest.raises(ReliefSizingError, match="P2"):
        relief_liquid(1.0, 800.0, 5e5, 5e5)        # no driving dP
    with pytest.raises(ReliefSizingError, match="rho"):
        relief_liquid(1.0, 0.0, 20e5, 5e5)
