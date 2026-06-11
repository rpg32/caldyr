"""M11 tests: the AirCooler unit op + Turton air-cooler costing.

Physics: an AirCooler is a cooler whose sink is ambient air — the outlet can
never get closer to the design air temperature than the minimum approach
(default 10 K above 308.15 K / 35 C; GPSA Engineering Data Book; Towler &
Sinnott 2e Ch. 12), and the only utility drawn is fan electricity
(``fan_power_frac``, default 0.02 kW per kW rejected — ACHE fans draw on the
order of 1-3% of the rejected duty, cf. GPSA Section 10).

Sizing: bare-tube-equivalent area = Q / (U * LMTD) against air heating
t_air_in -> t_air_in + 10 K, with U = 40 W/m^2/K (GPSA: ~30-60 W/m^2/K
bare-tube for gas service).

Costing reference: **Turton et al., 4e, Appendix A** — air cooler, Table A.1:
K1=4.0336, K2=0.2341, K3=0.0497 (area 10-10,000 m^2). Hand-checked point: at
A = 100 m^2, log10 Cp0 = 4.0336 + 0.2341*2 + 0.0497*4 = 4.7006 -> Cp0
~ $50,200 (CEPCI 397). Bare module B1=0.96, B2=1.21 (Table A.4) -> Fbm = 2.17
at Fp = Fm = 1.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import EquipmentSize, SizingOptions, size_flowsheet
from caldyr.thermo import make_package
from caldyr.unitops import AirCooler, AirCoolerApproachError

T_IN = 450.0     # K
T_OUT = 330.0    # K
P = 10.0e5       # Pa
N = 100.0        # mol/s


def n2_air_cooler(**param_overrides) -> Flowsheet:
    params: dict = {"T_out": T_OUT}
    params.update(param_overrides)
    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    fs.add(AirCooler("AC", params))
    fs.feed("FEED", "AC:in1", T=T_IN, P=P, molar_flow=N, z={"nitrogen": 1.0})
    fs.connect("OUT", "AC:out", None)
    fs.connect("Q", "AC:duty", None)
    return fs


def test_energy_balance_and_negative_duty():
    fs = n2_air_cooler()
    rep = fs.solve()
    assert rep.converged
    out, feed = fs.streams["OUT"], fs.streams["FEED"]
    assert out.T == pytest.approx(T_OUT)
    assert rep.duties["Q"] < 0                     # heat leaves the process
    assert out.molar_flow * out.H - feed.molar_flow * feed.H == \
        pytest.approx(rep.duties["Q"], rel=1e-9)


def test_approach_violation_raises_typed_error():
    # Default air at 308.15 K + 10 K approach -> T_out must be >= 318.15 K.
    with pytest.raises(AirCoolerApproachError, match="minimum approach"):
        n2_air_cooler(T_out=315.0).solve()
    # Custom ambient/approach move the limit.
    with pytest.raises(AirCoolerApproachError):
        n2_air_cooler(T_out=330.0, t_air_in=325.0).solve()
    with pytest.raises(AirCoolerApproachError):
        n2_air_cooler(T_out=330.0, approach=25.0).solve()
    # A typed error is still a ValueError (engine convention).
    assert issubclass(AirCoolerApproachError, ValueError)
    # ...and the same spec is fine once the ambient allows it.
    rep = n2_air_cooler(T_out=330.0, t_air_in=290.0, approach=25.0).solve()
    assert rep.converged


def test_heating_spec_and_missing_t_out_raise():
    with pytest.raises(ValueError, match="only cools"):
        n2_air_cooler(T_out=500.0).solve()
    with pytest.raises(ValueError, match="T_out"):
        n2_air_cooler(T_out=None).solve()


# -- economics ---------------------------------------------------------------
def test_sizing_area_and_fan_electricity():
    fs = n2_air_cooler()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    size = sizes[0]
    assert size.equipment_type == "air_cooler"
    assert size.attribute_name == "area_m2"

    # Hand check the area: counter-current LMTD vs air 308.15 -> 318.15 K.
    q = abs(rep.duties["Q"])
    dt1, dt2 = T_IN - 318.15, T_OUT - 308.15
    lmtd = (dt1 - dt2) / math.log(dt1 / dt2)
    assert size.attribute == pytest.approx(q / (40.0 * lmtd), rel=1e-9)
    assert math.isfinite(size.attribute) and size.attribute > 0

    # No cooling water: the utility is fan electricity at 2% of the duty.
    assert size.utility == "electricity"
    assert size.utility_duty_W == pytest.approx(0.02 * q, rel=1e-12)


def test_fan_power_frac_param_is_respected():
    fs = n2_air_cooler(fan_power_frac=0.035)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    size = size_flowsheet(fs, rep, pp)[0]
    assert size.utility_duty_W == pytest.approx(0.035 * abs(rep.duties["Q"]), rel=1e-12)


def test_air_cooler_U_option_scales_the_area():
    fs = n2_air_cooler()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    a40 = size_flowsheet(fs, rep, pp)[0].attribute
    a50 = size_flowsheet(fs, rep, pp, SizingOptions(air_cooler_U=50.0))[0].attribute
    assert a50 == pytest.approx(a40 * 40.0 / 50.0, rel=1e-12)


def test_costing_is_finite_positive_end_to_end():
    fs = n2_air_cooler()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    cost = cost_equipment(size_flowsheet(fs, rep, pp)[0])
    assert math.isfinite(cost.bare_module)
    assert cost.bare_module > cost.purchased > 0


def test_purchased_cost_matches_turton_hand_point():
    """A = 100 m^2 -> log10 Cp0 = 4.7006 -> Cp0 ~ $50,200 (CEPCI 397);
    Fbm = B1 + B2 = 2.17 at low pressure in carbon steel."""
    ac = EquipmentSize("AC1", "air_cooler", 100.0, "area_m2", pressure_barg=1.0)
    c = cost_equipment(ac, year=2001)
    assert c.purchased == pytest.approx(50_200.0, rel=0.01)
    assert c.factors["Fp"] == pytest.approx(1.0)          # below 10 barg
    assert c.factors["Fbm"] == pytest.approx(0.96 + 1.21, rel=1e-9)
    assert c.bare_module == pytest.approx(2.17 * c.purchased, rel=1e-9)
