"""M17 tests: economics sizers for the Claus units (ClausReactor, SulfurCondenser).

Validates that the full Claus train of ``examples/26_claus_sulfur_recovery.py``
sizes and costs end to end — before these sizers, ``size_flowsheet`` raised
``no sizer for unit type 'ClausReactor'``. Correlation basis (see the sizers):

* ``ClausReactor`` -> vertical process vessel on the reactor space-time
  (Turton 4e Table A.1, vertical process vessel), the same basis as the
  Gibbs/equilibrium reactors; the refractory furnace / sour service is
  under-costed as a CS vessel (order-of-magnitude, documented in the sizer).
* ``SulfurCondenser`` -> shell-and-tube exchanger on ``A = Q/(U·LMTD)``
  (Turton 4e Table A.1), a gas-service overall U (Towler & Sinnott 2e Ch. 12)
  against cooling water; the waste-heat-boiler steam credit is not auto-booked.
"""
import importlib.util
import pathlib

import pytest

from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.thermo import make_package

_EXAMPLE = (pathlib.Path(__file__).resolve().parents[1]
            / "examples" / "26_claus_sulfur_recovery.py")

CLAUS_REACTORS = {"FURN", "CV1", "CV2"}      # 1 adiabatic furnace + 2 converters
CONDENSERS = {"C1", "C2", "C3"}


def _load_claus_flowsheet():
    """Build the actual example-26 flowsheet (importlib by path — its numeric
    module name isn't a legal import, and __name__ != '__main__' keeps main off)."""
    spec = importlib.util.spec_from_file_location("claus_example_26", _EXAMPLE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build()


@pytest.fixture(scope="module")
def sized():
    fs = _load_claus_flowsheet()
    rep = fs.solve()
    assert rep.converged
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)       # must not raise "no sizer ..."
    return {s.unit_id: s for s in sizes}


def test_every_claus_unit_is_sized(sized):
    """The whole train sizes — no unit raises 'no sizer for unit type'."""
    for uid in CLAUS_REACTORS | CONDENSERS:
        assert uid in sized


def test_claus_reactors_are_vertical_vessels(sized):
    for uid in CLAUS_REACTORS:
        s = sized[uid]
        assert s.equipment_type == "vessel_vertical"
        assert s.attribute_name == "volume_m3"
        # metres-scale reaction volume, not micro or absurd
        assert 0.01 < s.attribute < 1e4


def test_sulfur_condensers_are_shell_and_tube(sized):
    for uid in CONDENSERS:
        s = sized[uid]
        assert s.equipment_type == "heat_exchanger"
        assert s.attribute_name == "area_m2"
        assert s.attribute > 0.0
        assert s.utility is None              # WHB steam credit not auto-booked
    # C1 quenches the ~1500 K furnace effluent -> the largest area of the three
    assert sized["C1"].attribute > sized["C3"].attribute


def test_costing_runs_end_to_end_with_plausible_magnitudes(sized):
    for uid in CLAUS_REACTORS | CONDENSERS:
        c = cost_equipment(sized[uid])
        assert c.bare_module > c.purchased > 0.0
        # sanity band for one equipment item ($1k .. $100M, CEPCI-escalated)
        assert 1e3 < c.bare_module < 1e8


def test_regression_pins_furnace_and_first_condenser(sized):
    """Pin the sizes/costs the sizers produced when they were validated, so a
    silent multiple-x sizing regression can't hide inside the sanity bands.
    (An intentional correlation or CEPCI change should update these numbers.)"""
    furn = sized["FURN"]
    c1 = sized["C1"]
    assert furn.attribute == pytest.approx(105.12, rel=1e-3)      # m^3
    assert c1.attribute == pytest.approx(549.92, rel=1e-3)        # m^2
    assert cost_equipment(furn).bare_module == pytest.approx(576_657, rel=1e-3)
    assert cost_equipment(c1).bare_module == pytest.approx(348_141, rel=1e-3)
