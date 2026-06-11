"""M12 tests: the Evaporator unit op (registered "Evaporator").

Reference: Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley
2025), §5.2 "Simulation of Evaporator" — 3000 kg/h of 10 wt% sucrose solution
at 30 C / 110 kPa is concentrated to 50 wt% at a shell pressure of 70 kPa,
heated by condensing saturated steam at 152 C (latent heat only). Book
results: vapor = 2400 kg/h, steam required = 2894 kg/h, "efficiency" =
vapor/steam = 0.83 (book Eq. 5.7); Figure 5.17 tabulates the steam required
vs steam temperature (120–170 C).

Reproduction strategy: PR/SRK cannot represent sucrose VLE (thermo's flash
returns non-physical vapor fractions for water+sucrose — verified), so the
book case is reproduced on ``coolprop:Water`` as the *pure-water surrogate*:
evaporate the same 2400 kg/h of water out of a 3000 kg/h feed (mass VF = 0.8,
= molar VF for a pure fluid) at 70 kPa from 30 C. The surrogate ignores the
10% sucrose, whose heat capacity (~1.2 kJ/kg/K) is far below water's
(~4.2 kJ/kg/K), so it slightly *over*-estimates the sensible-heat load:
measured steam requirement 2956.7 kg/h, +2.2% over the book's 2894 kg/h —
asserted within 3% with the bias direction checked.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.unitops import Evaporator

MW_W = 0.018015268                      # kg/mol, water (CoolProp)
FEED_KG_H = 3000.0
FEED_MOL_S = FEED_KG_H / 3600.0 / MW_W
P_SHELL = 70e3                          # Pa
T_FEED = 303.15                         # 30 C
P_FEED = 110e3


def water_evaporator(**params) -> Flowsheet:
    fs = Flowsheet(components=[Component("water")],
                   property_package="coolprop:Water")
    fs.add(Evaporator("EVAP", {"P": P_SHELL, **params}))
    fs.feed("FEED", "EVAP:in1", T=T_FEED, P=P_FEED, molar_flow=FEED_MOL_S,
            z={"water": 1.0})
    fs.connect("VAP", "EVAP:vapor", None)
    fs.connect("LIQ", "EVAP:liquid", None)
    fs.connect("Q", "EVAP:duty", None)
    return fs


def _hfg_kJ_per_kg(t_celsius: float) -> float:
    """Latent heat of saturated steam at t (C), kJ/kg, from CoolProp."""
    from CoolProp.CoolProp import PropsSI

    t = t_celsius + 273.15
    return (PropsSI("Hmass", "T", t, "Q", 1, "Water")
            - PropsSI("Hmass", "T", t, "Q", 0, "Water")) / 1e3


# -- book §5.2 reproduction (pure-water surrogate) ------------------------------
def test_book_steam_requirement_and_efficiency():
    fs = water_evaporator(vapor_fraction=0.8)       # 2400 of 3000 kg/h
    rep = fs.solve()
    assert rep.converged

    # The flash sits at Tsat(70 kPa) = 89.93 C (Cengel & Boles A-5: 89.95 C).
    assert fs.streams["VAP"].T == pytest.approx(363.08, abs=0.05)
    assert fs.streams["VAP"].molar_flow * MW_W * 3600 == pytest.approx(2400.0, rel=1e-9)
    assert fs.streams["LIQ"].molar_flow * MW_W * 3600 == pytest.approx(600.0, rel=1e-9)

    # Steam at 152 C condenses (latent heat only): m_steam = Q / h_fg(152 C).
    q_kJ_h = rep.duties["Q"] * 3.6                  # W -> kJ/h
    steam = q_kJ_h / _hfg_kJ_per_kg(152.0)          # kg/h
    assert steam == pytest.approx(2894.0, rel=0.03)  # book: 2894 kg/h
    assert steam > 2894.0                            # surrogate bias is known +2.2%

    # Book Eq. (5.7): efficiency = vapor separated / steam required = 0.83.
    assert 2400.0 / steam == pytest.approx(0.83, abs=0.025)


def test_book_steam_temperature_sensitivity():
    """Book Figure 5.17: more steam is required as its temperature rises
    because h_fg falls. Reproduced within 3% at every tabulated point, with
    the monotone-increasing trend exact."""
    book = {120: 2785.0, 130: 2817.0, 140: 2851.0,
            150: 2887.0, 160: 2926.0, 170: 2967.0}
    rep = water_evaporator(vapor_fraction=0.8).solve()
    q_kJ_h = rep.duties["Q"] * 3.6
    steam = {t: q_kJ_h / _hfg_kJ_per_kg(t) for t in book}
    for t, m_book in book.items():
        assert steam[t] == pytest.approx(m_book, rel=0.03)
    temps = sorted(book)
    assert all(steam[a] < steam[b] for a, b in zip(temps, temps[1:]))


# -- spec handling ---------------------------------------------------------------
def test_three_specs_are_equivalent_on_a_mixture():
    """VF / duty / T specs must describe the same operating point. Run the VF
    spec (brentq on T between bubble and dew) on a butane/pentane feed, then
    feed the resulting duty and temperature back as the alternative specs."""
    def build(**params):
        fs = Flowsheet(components=[Component("n-butane"), Component("n-pentane")],
                       property_package="thermo:PR")
        fs.add(Evaporator("EVAP", {"P": 5e5, **params}))
        fs.feed("FEED", "EVAP:in1", T=300.0, P=6e5, molar_flow=10.0,
                z={"n-butane": 0.5, "n-pentane": 0.5})
        fs.connect("VAP", "EVAP:vapor", None)
        fs.connect("LIQ", "EVAP:liquid", None)
        fs.connect("Q", "EVAP:duty", None)
        rep = fs.solve()
        assert rep.converged
        return fs, rep

    fs_vf, rep_vf = build(vapor_fraction=0.4)
    vap, liq, feed = fs_vf.streams["VAP"], fs_vf.streams["LIQ"], fs_vf.streams["FEED"]
    assert vap.molar_flow == pytest.approx(4.0, rel=1e-6)

    # energy balance closes: n_in*H_in + Q = n_v*H_v + n_l*H_l
    lhs = feed.molar_flow * feed.H + rep_vf.duties["Q"]
    rhs = vap.molar_flow * vap.H + liq.molar_flow * liq.H
    assert lhs == pytest.approx(rhs, rel=1e-9)
    # the vapor is enriched in the light component
    assert vap.z["n-butane"] > 0.5 > liq.z["n-butane"]

    fs_q, _ = build(duty=rep_vf.duties["Q"])
    assert fs_q.streams["VAP"].molar_flow == pytest.approx(4.0, rel=1e-6)
    assert fs_q.streams["VAP"].T == pytest.approx(vap.T, rel=1e-9)

    fs_t, rep_t = build(T=vap.T)
    assert fs_t.streams["VAP"].molar_flow == pytest.approx(4.0, rel=1e-6)
    assert rep_t.duties["Q"] == pytest.approx(rep_vf.duties["Q"], rel=1e-6)


def test_pure_fluid_vf_spec_is_exact():
    """For pure water the two-phase region is a single temperature, so the VF
    spec resolves through the saturated enthalpies (no iteration): the duty is
    exactly sensible heat to Tsat plus VF times the latent heat."""
    from caldyr.thermo import make_package

    rep = water_evaporator(vapor_fraction=0.8).solve()
    pp = make_package("coolprop:Water", ["water"])
    sat = pp.bubble_point(P_SHELL, {"water": 1.0})
    h_in = pp.enthalpy(T_FEED, P_FEED, {"water": 1.0})
    q_hand = FEED_MOL_S * ((0.2 * sat.H_liquid + 0.8 * sat.H_vapor) - h_in)
    assert rep.duties["Q"] == pytest.approx(q_hand, rel=1e-9)


def test_spec_errors_are_typed():
    with pytest.raises(ValueError, match="exactly one"):
        water_evaporator().solve()
    with pytest.raises(ValueError, match="exactly one"):
        water_evaporator(T=363.0, duty=1e6).solve()
    with pytest.raises(ValueError, match="in \\[0, 1\\]"):
        water_evaporator(vapor_fraction=1.2).solve()


# -- economics --------------------------------------------------------------------
def test_sizing_vessel_plus_heating_utility_and_costing():
    from caldyr.economics.costing import cost_equipment
    from caldyr.economics.sizing import size_flowsheet
    from caldyr.thermo import make_package

    fs = water_evaporator(vapor_fraction=0.8)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    vessel = sizes[0]
    assert vessel.equipment_type == "vessel_vertical"
    assert vessel.attribute_name == "volume_m3"
    assert vessel.attribute > 0 and math.isfinite(vessel.attribute)
    # the heating duty draws a hot utility able to reach the boiling point
    assert vessel.utility is not None
    assert vessel.utility_duty_W == pytest.approx(rep.duties["Q"], rel=1e-12)

    cost = cost_equipment(vessel)
    assert math.isfinite(cost.bare_module)
    assert cost.bare_module > cost.purchased > 0
