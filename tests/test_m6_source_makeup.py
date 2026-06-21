"""Tests for two flowsheet-plumbing primitives:

* :class:`~caldyr.unitops.Source` — a boundary feed expressed as a unit op, so
  its rate is a *parameter* the logical ops / optimizer can drive (a plain
  ``Flowsheet.feed`` is fixed).
* :class:`~caldyr.unitops.Makeup` — an in-loop make-up controller that holds one
  component of a stream at a target molar flow by injecting a pure top-up. The
  robust way to close a recycle with a drifting solvent/water inventory (used by
  the §15.3 amine plant in ``examples/24`` / ``test_m16_amine_package``).
"""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.economics.sizing import SIZER_REGISTRY
from caldyr.thermo import make_package
from caldyr.unitops import REGISTRY, Makeup, Source

COMPS = ["nitrogen", "water", "methane"]


def _pp():
    return make_package("thermo:PR", COMPS)


# -- Source -------------------------------------------------------------------
def test_source_emits_a_resolved_stream_from_params():
    src = Source("SRC", {"molar_flow": 50.0, "T": 313.15, "P": 1.0e5,
                         "z": {c: (1.0 if c == "water" else 0.0) for c in COMPS}})
    out = src.solve({}, _pp())["out"]
    assert out.molar_flow == pytest.approx(50.0)
    assert out.T == pytest.approx(313.15)
    assert out.z["water"] == pytest.approx(1.0)
    assert out.H is not None and out.phase is not None      # fully resolved


def test_source_normalizes_composition():
    src = Source("SRC", {"molar_flow": 10.0, "T": 300.0, "P": 1.0e5,
                         "z": {"nitrogen": 3.0, "methane": 1.0}})
    out = src.solve({}, _pp())["out"]
    assert out.z["nitrogen"] == pytest.approx(0.75)
    assert out.z["methane"] == pytest.approx(0.25)


def test_source_bad_params_raise():
    with pytest.raises(ValueError, match="required"):
        Source("S", {"T": 300.0, "P": 1e5, "z": {"water": 1.0}}).solve({}, _pp())
    with pytest.raises(ValueError, match=">= 0"):
        Source("S", {"molar_flow": -1.0, "T": 300.0, "P": 1e5,
                     "z": {"water": 1.0}}).solve({}, _pp())
    with pytest.raises(ValueError, match="positive total"):
        Source("S", {"molar_flow": 1.0, "T": 300.0, "P": 1e5, "z": {}}).solve({}, _pp())


def test_source_registered_and_costless():
    assert REGISTRY["Source"] is Source
    assert SIZER_REGISTRY["Source"](Source("S", {}), None) == []


# -- Makeup -------------------------------------------------------------------
def _feed_stream(flows, T=300.0, P=1.0e5):
    from caldyr.core import Stream
    F = sum(flows.values())
    return Stream(id="s", components=COMPS, T=T, P=P, molar_flow=F,
                  z={c: flows.get(c, 0.0) / F for c in COMPS})


def test_makeup_tops_a_component_up_to_target():
    """Inlet carries 30 mol/s water; target 80 -> add 50 mol/s pure water."""
    mk = Makeup("MK", {"component": "water", "target": 80.0})
    inlet = _feed_stream({"nitrogen": 70.0, "water": 30.0})
    out = mk.solve({"in1": inlet}, _pp())["out"]
    assert mk.design["makeup_flow"] == pytest.approx(50.0)
    assert out.molar_flow == pytest.approx(150.0)            # 100 in + 50 water
    assert out.molar_flow * out.z["water"] == pytest.approx(80.0)
    assert out.molar_flow * out.z["nitrogen"] == pytest.approx(70.0)   # untouched


def test_makeup_surplus_adds_nothing():
    """If the inlet already exceeds the target, no make-up is added (a surplus
    needs a purge, not a make-up)."""
    mk = Makeup("MK", {"component": "water", "target": 20.0})
    inlet = _feed_stream({"nitrogen": 70.0, "water": 30.0})
    out = mk.solve({"in1": inlet}, _pp())["out"]
    assert mk.design["makeup_flow"] == pytest.approx(0.0)
    assert out.molar_flow == pytest.approx(100.0)


def test_makeup_bad_params_raise():
    pp = _pp()
    inlet = _feed_stream({"nitrogen": 70.0, "water": 30.0})
    with pytest.raises(ValueError, match="component"):
        Makeup("MK", {"component": "argon", "target": 80.0}).solve({"in1": inlet}, pp)
    with pytest.raises(ValueError, match="target"):
        Makeup("MK", {"component": "water"}).solve({"in1": inlet}, pp)
    with pytest.raises(ValueError, match=">= 0"):
        Makeup("MK", {"component": "water", "target": -1.0}).solve({"in1": inlet}, pp)


def test_makeup_registered_and_costless():
    assert REGISTRY["Makeup"] is Makeup
    assert SIZER_REGISTRY["Makeup"](Makeup("M", {}), None) == []


def test_source_makeup_round_trip_in_flowsheet():
    """Both units persist through `.flow` IO with their params intact."""
    fs = Flowsheet(components=[Component(c) for c in COMPS],
                   property_package="thermo:PR")
    fs.add(Source("SRC", {"molar_flow": 5.0, "T": 300.0, "P": 1e5,
                          "z": {"water": 1.0}}))
    fs.add(Makeup("MK", {"component": "water", "target": 50.0, "T": 300.0}))
    fs.connect("S", "SRC:out", "MK:in1")
    fs.connect("OUT", "MK:out", None)
    fs2 = from_dict(to_dict(fs))
    assert fs2.units["SRC"].params["molar_flow"] == 5.0
    assert fs2.units["MK"].params["target"] == 50.0
    assert type(fs2.units["SRC"]).__name__ == "Source"
    assert type(fs2.units["MK"]).__name__ == "Makeup"
