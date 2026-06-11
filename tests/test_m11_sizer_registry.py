"""M11 hygiene: the sizer registry that replaced the isinstance/type-name
dispatch chain in economics/sizing.py. The pre-existing economics suite is the
behavioral regression net; these tests pin the registry mechanics."""
import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.sizing import (
    SIZER_REGISTRY,
    EquipmentSize,
    SizerContext,
    register_sizer,
    size_flowsheet,
)
from caldyr.thermo import make_package
from caldyr.unitops import Valve


def test_registry_covers_every_previously_dispatched_unit_type():
    expected = {
        "Mixer", "Splitter",                       # negligible (no items)
        "Heater", "HeatExchanger", "Pump", "Compressor", "Expander",
        "FlashDrum", "ComponentSplitter", "ThreePhaseSeparator",
        "ConversionReactor", "EquilibriumReactor", "GibbsReactor",
        "CSTR", "PFR", "ShortcutColumn", "RigorousColumn",
        "FiredHeater", "AirCooler",                # M11 additions
    }
    assert expected <= set(SIZER_REGISTRY)


def test_unsized_unit_type_still_raises_the_same_error():
    """An unregistered unit type raises the same typed error as before the
    refactor. (Valve itself gained a negligible-cost sizer in M13, so the
    probe is a deliberately unregistered subclass.)"""
    class UnsizedThrottle(Valve):
        pass

    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    fs.add(UnsizedThrottle("V1", {"P_out": 1.0e5}))
    fs.feed("FEED", "V1:in1", T=300.0, P=5.0e5, molar_flow=10.0, z={"nitrogen": 1.0})
    fs.connect("OUT", "V1:out", None)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    with pytest.raises(ValueError, match="no sizer for unit type 'UnsizedThrottle'"):
        size_flowsheet(fs, rep, pp)


def test_register_sizer_extends_size_flowsheet(monkeypatch):
    """Registering a sizer for a new type makes size_flowsheet handle it —
    the extension point new unit ops use."""
    def _valve_sizer(unit, ctx: SizerContext) -> list[EquipmentSize]:
        return [EquipmentSize(unit_id=unit.id, equipment_type="vessel_vertical",
                              attribute=1.0, attribute_name="volume_m3")]

    monkeypatch.setitem(SIZER_REGISTRY, "Valve", _valve_sizer)
    fs = Flowsheet(components=[Component("nitrogen")], property_package="thermo:PR")
    fs.add(Valve("V1", {"P_out": 1.0e5}))
    fs.feed("FEED", "V1:in1", T=300.0, P=5.0e5, molar_flow=10.0, z={"nitrogen": 1.0})
    fs.connect("OUT", "V1:out", None)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert [s.unit_id for s in sizes] == ["V1"]


def test_register_sizer_decorator_returns_the_function():
    @register_sizer("__M11TestOnlyUnit__")
    def _sizer(unit, ctx):
        return []
    try:
        assert SIZER_REGISTRY["__M11TestOnlyUnit__"] is _sizer
    finally:
        del SIZER_REGISTRY["__M11TestOnlyUnit__"]
