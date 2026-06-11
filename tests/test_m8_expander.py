"""M8 tests: the Expander (turbine) unit op — the mirror of Compressor.

Validation reference: isentropic expansion of nitrogen treated as an ideal gas
(PR at <= 10 bar is within a fraction of a percent of ideal for N2).
For a reversible adiabatic ideal-gas path (Smith, Van Ness & Abbott,
*Introduction to Chemical Engineering Thermodynamics*, 7e, Eq. 3.30b):

    T2 / T1 = (P2 / P1)^(R / Cp)

Nitrogen Cp_ig ~ 29.1 J/mol/K over 250-500 K (SVA 7e Table C.1: Cp/R = 3.280 +
0.593e-3 T -> 29.1 at the 380 K path mean). Expanding 100 mol/s from 500 K,
10 bar to 1 bar:

    T2s = 500 * (0.1)^(8.314/29.1) = 259.2 K
    W_ideal = n * Cp * (T2s - T1) = 100 * 29.1 * (259.2 - 500) = -700.7 kW

Sign convention (documented on the class): duty positive = energy added to the
process, so an expander's work duty is NEGATIVE — the extracted shaft power is
its magnitude. With eta < 1 the enthalpy drop (and so the work) scales by eta
exactly: w = eta * (h_isentropic - h_in).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.thermo import make_package

T_IN = 500.0          # K
P_IN = 10.0e5         # Pa
P_OUT = 1.0e5         # Pa
N = 100.0             # mol/s

REF_T2S = 259.2       # K, ideal-gas isentropic outlet (SVA 7e Eq. 3.30b)
REF_W_IDEAL = -700.7e3   # W, n*Cp*(T2s - T1)


def n2_expander(**param_overrides) -> Flowsheet:
    params: dict = {"P_out": P_OUT}
    params.update(param_overrides)
    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    from caldyr.unitops import Expander
    fs.add(Expander("EXP", params))
    fs.feed("FEED", "EXP:in1", T=T_IN, P=P_IN, molar_flow=N, z={"nitrogen": 1.0})
    fs.connect("OUT", "EXP:out", None)
    fs.connect("W", "EXP:work", None)
    return fs


def test_isentropic_expansion_matches_ideal_gas_reference():
    fs = n2_expander(eta=1.0)
    rep = fs.solve()
    assert rep.converged
    out = fs.streams["OUT"]
    assert out.P == pytest.approx(P_OUT)
    assert out.T == pytest.approx(REF_T2S, rel=0.02)
    # Work is EXTRACTED: negative duty, magnitude ~ n*Cp*dT.
    assert rep.duties["W"] < 0
    assert rep.duties["W"] == pytest.approx(REF_W_IDEAL, rel=0.05)


def test_efficiency_scales_the_work_exactly_and_warms_the_outlet():
    rep1 = n2_expander(eta=1.0).solve()
    fs8 = n2_expander(eta=0.8)
    rep8 = fs8.solve()
    # w = eta*(h_isentropic - h_in) -> the duty scales by eta *exactly*.
    assert rep8.duties["W"] == pytest.approx(0.8 * rep1.duties["W"], rel=1e-9)
    # Less work extracted -> the real outlet is warmer than the isentropic one.
    assert fs8.streams["OUT"].T > REF_T2S


def test_default_eta_is_080():
    rep_default = n2_expander().solve()           # no eta param at all
    rep_explicit = n2_expander(eta=0.8).solve()
    assert rep_default.duties["W"] == pytest.approx(rep_explicit.duties["W"], rel=1e-12)


def test_energy_balance_closes_exactly():
    fs = n2_expander(eta=0.8)
    rep = fs.solve()
    feed, out = fs.streams["FEED"], fs.streams["OUT"]
    assert out.molar_flow * out.H - feed.molar_flow * feed.H == \
        pytest.approx(rep.duties["W"], rel=1e-12)


def test_pressure_rise_raises():
    fs = n2_expander(P_out=20.0e5)
    with pytest.raises(ValueError, match="must be below inlet P"):
        fs.solve()


def test_bad_eta_raises():
    with pytest.raises(ValueError, match="eta"):
        n2_expander(eta=1.5).solve()
    with pytest.raises(ValueError, match="eta"):
        n2_expander(eta=0.0).solve()


def test_economics_costs_an_axial_gas_turbine():
    """Sized from |shaft power| and costed with the Turton 4e Table A.1 axial
    gas turbine correlation + direct bare-module factor (Fig. A.19, CS)."""
    fs = n2_expander(eta=0.8)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    size = sizes[0]
    assert size.equipment_type == "turbine_axial"
    assert size.attribute_name == "power_kW"
    # ~ 0.8 * 700 kW shaft power, inside the 100-4000 kW correlation range.
    assert size.attribute == pytest.approx(abs(rep.duties["W"]) / 1e3, rel=1e-12)
    assert 100.0 < size.attribute < 4000.0
    assert size.utility is None        # produces power; draws no utility

    cost = cost_equipment(size)
    assert cost.bare_module > cost.purchased > 0
    assert cost.factors["Fbm"] == pytest.approx(3.5)
    assert not cost.warnings
