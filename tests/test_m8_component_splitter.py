"""M8 tests: the ComponentSplitter — a black-box separator that splits each
component between overhead and bottoms by a specified fraction.

There is no physics to validate against (the split *is* the spec); what must
hold exactly, by construction, is the component balance (overhead rate =
split * feed rate for every component) and the energy balance (the duty port
reports the net enthalpy change, so feed enthalpy + duty = product enthalpy).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import ComponentSplitter

P_ATM = 101325.0
FEED_Z = {"benzene": 0.4, "toluene": 0.4, "p-xylene": 0.2}
N_FEED = 10.0
SPLITS = {"benzene": 0.98, "toluene": 0.50}      # p-xylene falls to default


def btx_splitter(**param_overrides) -> Flowsheet:
    params: dict = {"splits": dict(SPLITS)}
    params.update(param_overrides)
    fs = Flowsheet(
        components=[Component("benzene"), Component("toluene"), Component("p-xylene")],
        property_package="thermo:PR")
    fs.add(ComponentSplitter("CS", params))
    fs.feed("FEED", "CS:in1", T=350.0, P=P_ATM, molar_flow=N_FEED, z=dict(FEED_Z))
    fs.connect("OV", "CS:overhead", None)
    fs.connect("BT", "CS:bottoms", None)
    fs.connect("Q", "CS:duty", None)
    return fs


def test_component_balance_is_exact():
    fs = btx_splitter()
    assert fs.solve().converged
    feed, ov, bt = fs.streams["FEED"], fs.streams["OV"], fs.streams["BT"]
    for c in feed.components:
        f = N_FEED * FEED_Z[c]
        split = SPLITS.get(c, 0.0)                     # default_split = 0.0
        assert ov.molar_flow * ov.z.get(c, 0.0) == pytest.approx(split * f, abs=1e-12)
        assert bt.molar_flow * bt.z.get(c, 0.0) == pytest.approx((1 - split) * f, abs=1e-12)
    assert ov.molar_flow + bt.molar_flow == pytest.approx(N_FEED, rel=1e-12)


def test_default_split_routes_unlisted_components():
    fs = btx_splitter(default_split=1.0)               # unlisted -> all overhead
    fs.solve()
    ov = fs.streams["OV"]
    assert ov.molar_flow * ov.z["p-xylene"] == pytest.approx(N_FEED * FEED_Z["p-xylene"],
                                                             rel=1e-12)


def test_duty_closes_the_energy_balance_exactly():
    fs = btx_splitter(T_overhead=320.0, T_bottoms=380.0)
    rep = fs.solve()
    feed, ov, bt = fs.streams["FEED"], fs.streams["OV"], fs.streams["BT"]
    h_in = feed.molar_flow * feed.H + rep.duties["Q"]
    h_out = ov.molar_flow * ov.H + bt.molar_flow * bt.H
    assert h_out == pytest.approx(h_in, rel=1e-12)


def test_outlet_t_p_specs_default_to_inlet():
    fs = btx_splitter()
    fs.solve()
    for sid in ("OV", "BT"):
        assert fs.streams[sid].T == pytest.approx(350.0)
        assert fs.streams[sid].P == pytest.approx(P_ATM)

    fs2 = btx_splitter(T_overhead=320.0, P_bottoms=2 * P_ATM)
    fs2.solve()
    assert fs2.streams["OV"].T == pytest.approx(320.0)
    assert fs2.streams["OV"].P == pytest.approx(P_ATM)        # P default: inlet
    assert fs2.streams["BT"].T == pytest.approx(350.0)        # T default: inlet
    assert fs2.streams["BT"].P == pytest.approx(2 * P_ATM)


def test_everything_to_overhead_leaves_an_empty_bottoms():
    fs = btx_splitter(splits={}, default_split=1.0)
    fs.solve()
    assert fs.streams["BT"].molar_flow == pytest.approx(0.0, abs=1e-15)
    assert fs.streams["OV"].molar_flow == pytest.approx(N_FEED, rel=1e-12)
    # The empty stream still has a sane state (inlet composition).
    assert sum(fs.streams["BT"].z.values()) == pytest.approx(1.0, rel=1e-9)


def test_bad_split_fraction_raises():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        btx_splitter(splits={"benzene": 1.2}).solve()


def test_unknown_component_in_splits_raises():
    with pytest.raises(ValueError, match="not in the flowsheet"):
        btx_splitter(splits={"hexane": 0.5}).solve()


def test_flow_round_trip_is_exact():
    fs = btx_splitter(T_overhead=320.0)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))


def test_economics_sizes_a_vertical_vessel():
    """Costed as a vertical process vessel sized by residence time, exactly
    like a flash drum (a deliberate placeholder heuristic for a black box)."""
    fs = btx_splitter()
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    assert sizes[0].equipment_type == "vessel_vertical"
    assert sizes[0].attribute_name == "volume_m3"
    assert sizes[0].attribute > 0
    assert cost_equipment(sizes[0]).bare_module > 0
