"""M12 tests: the psychrometric ``humidity`` utility (CoolProp HAPropsSI,
ASHRAE RP-1485 real-gas humid air).

References:
* Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley 2025), §2.4,
  Figure 2.28 — air at 30 C, 1 atm, RH = 50%: humidity 1.349e-2 kg water/kg
  dry air, dew point 18.86 C, water partial pressure 2.150 kPa. (HYSYS PR
  values; ASHRAE-grade psychrometrics lands within ~1% / 0.5 K of them.)
  The book's wet-bulb readout (20.49 C) disagrees with the psychrometric
  chart (~22 C at 30 C db / 50% RH) — we validate Twb against the chart, not
  the book; see test_chart_point.
* ASHRAE Handbook—Fundamentals, Psychrometric Chart No. 1 (sea level); also
  Cengel & Boles 8e, Fig. A-31: at 25 C db and 50% RH, w ~ 0.0099 kg/kg,
  Twb ~ 17.9 C, h ~ 50.4 kJ/kg dry air.
"""
import pytest

from caldyr.analysis import humidity

P_ATM = 101_325.0


def test_book_worked_example_30C_rh50():
    """Book §2.4 Figure 2.28 (T = 30 C, P = 1 atm, RH = 50%)."""
    h = humidity(303.15, P_ATM, rh=0.5)
    assert h["w"] == pytest.approx(1.349e-2, rel=0.015)        # book: 1.349e-2 kg/kg
    assert h["t_dp"] == pytest.approx(18.86 + 273.15, abs=0.5)  # book: 18.86 C
    assert h["p_w"] == pytest.approx(2150.0, rel=0.015)        # book: 2.150 kPa
    assert h["rh"] == pytest.approx(0.5, abs=1e-9)


def test_chart_point_25C_rh50():
    """ASHRAE psychrometric chart (Cengel & Boles Fig. A-31): 25 C, 50% RH ->
    w ~ 0.0099 kg/kg, Twb ~ 17.9 C, h ~ 50.4 kJ/kg dry air, v ~ 0.858 m^3/kg
    dry air, Tdp ~ 13.9 C."""
    h = humidity(298.15, P_ATM, rh=0.5)
    assert h["w"] == pytest.approx(0.0099, rel=0.01)
    assert h["t_wb"] == pytest.approx(17.9 + 273.15, abs=0.3)
    assert h["h"] == pytest.approx(50.4e3, rel=0.01)
    assert h["v"] == pytest.approx(0.858, rel=0.005)
    assert h["t_dp"] == pytest.approx(13.9 + 273.15, abs=0.3)


def test_wet_bulb_at_30C_rh50_matches_chart_not_book():
    """At 30 C db / 50% RH the chart wet bulb is ~22.0 C; the book's HYSYS
    screenshot shows 20.49 C — we follow ASHRAE."""
    h = humidity(303.15, P_ATM, rh=0.5)
    assert h["t_wb"] == pytest.approx(22.0 + 273.15, abs=0.3)


@pytest.mark.parametrize("key", ["w", "t_dp", "t_wb"])
def test_each_spec_inverts_back_to_the_same_state(key):
    """Specifying the state by humidity ratio / dew point / wet bulb recovers
    the same RH (round-trip through HAPropsSI)."""
    ref = humidity(303.15, P_ATM, rh=0.5)
    again = humidity(303.15, P_ATM, **{key: ref[key]})
    assert again["rh"] == pytest.approx(0.5, abs=1e-4)
    assert again["w"] == pytest.approx(ref["w"], rel=1e-4)


def test_saturated_air_dew_point_equals_dry_bulb():
    h = humidity(303.15, P_ATM, rh=1.0)
    assert h["t_dp"] == pytest.approx(303.15, abs=0.01)
    assert h["t_wb"] == pytest.approx(303.15, abs=0.01)


def test_spec_errors_are_typed():
    with pytest.raises(ValueError, match="exactly one"):
        humidity(303.15, P_ATM)
    with pytest.raises(ValueError, match="exactly one"):
        humidity(303.15, P_ATM, rh=0.5, w=0.01)
    with pytest.raises(ValueError, match="fraction"):
        humidity(303.15, P_ATM, rh=60.0)        # percent, not a fraction
    with pytest.raises(ValueError, match=">= 0"):
        humidity(303.15, P_ATM, w=-0.01)
