"""M9 tests: the component database (caldyr.core.components_db) and the
economics molar-mass fallback.

References: component identity (name/formula/CAS/MW) comes from the `chemicals`
library (Caleb Bell, MIT) — the same database the thermo layer resolves
against, so an id that validates here is an id the property packages can build.
"""
import pytest

from caldyr.core import (
    COMMON_COMPONENTS,
    Component,
    Flowsheet,
    UnknownComponentError,
    resolve_component,
)
from caldyr.core.components_db import molar_mass as chem_molar_mass
from caldyr.economics import TEAConfig, analyze
from caldyr.economics import data as econ_data
from caldyr.thermo import make_package
from caldyr.unitops import Heater


# -- resolve_component ----------------------------------------------------------
def test_resolves_by_name_cas_and_formula_to_the_same_species():
    by_name = resolve_component("water")
    by_cas = resolve_component("7732-18-5")
    by_formula = resolve_component("H2O")
    assert by_name.cas == by_cas.cas == by_formula.cas == "7732-18-5"
    assert by_name.formula == "H2O"
    assert by_name.name == "water"
    # the id is preserved exactly as given (streams key on it)
    assert by_cas.id == "7732-18-5"


def test_unknown_identifier_raises_typed_error_naming_it():
    with pytest.raises(UnknownComponentError, match="notachemicalxyz"):
        resolve_component("notachemicalxyz")
    with pytest.raises(UnknownComponentError, match="empty"):
        resolve_component("  ")


def test_molar_mass_lookup():
    assert chem_molar_mass("methanol") == pytest.approx(0.0320419, rel=1e-4)
    with pytest.raises(UnknownComponentError, match="notachemicalxyz"):
        chem_molar_mass("notachemicalxyz")


# -- the curated catalog ----------------------------------------------------------
def test_catalog_is_well_formed():
    assert 80 <= len(COMMON_COMPONENTS) <= 140
    ids = [c["id"] for c in COMMON_COMPONENTS]
    assert len(set(ids)) == len(ids)                 # unique ids
    cas = [c["cas"] for c in COMMON_COMPONENTS]
    assert len(set(cas)) == len(cas)                 # unique species
    for entry in COMMON_COMPONENTS:
        assert set(entry) == {"id", "name", "formula", "cas"}
        assert all(entry[k] for k in entry)


def test_every_catalog_id_resolves_to_its_recorded_cas():
    """Guards the author-time hardcoding against drift in the chemicals
    database (and against ids that quietly resolve to the wrong species, as
    'R-125' does — which is why its id is 'pentafluoroethane')."""
    for entry in COMMON_COMPONENTS:
        resolved = resolve_component(entry["id"])
        assert resolved.cas == entry["cas"], (
            f"catalog id {entry['id']!r} resolved to CAS {resolved.cas}, "
            f"recorded {entry['cas']}")


# -- flowsheet component validation -------------------------------------------------
def test_make_package_rejects_unresolvable_component_ids():
    with pytest.raises(UnknownComponentError, match="notachemicalxyz"):
        make_package("thermo:PR", ["water", "notachemicalxyz"])


def test_flowsheet_solve_names_the_bad_component():
    fs = Flowsheet(components=[Component("water"), Component("notachemicalxyz")],
                   property_package="thermo:PR")
    fs.add(Heater("H", {"T_out": 350.0}))
    fs.feed("F", "H:in1", T=300.0, P=101325.0, molar_flow=1.0,
            z={"water": 1.0, "notachemicalxyz": 0.0})
    fs.connect("O", "H:out", None)
    fs.connect("Q", "H:duty", None)
    with pytest.raises(UnknownComponentError, match="notachemicalxyz"):
        fs.solve()


# -- economics molar-mass fallback -----------------------------------------------------
def test_econ_molar_mass_falls_back_to_chemicals():
    # methanol is deliberately NOT in the hardcoded economics table
    assert "methanol" not in econ_data.MOLAR_MASS
    assert econ_data.molar_mass("methanol") == pytest.approx(0.0320419, rel=1e-4)
    # local table still wins where present
    assert econ_data.molar_mass("water") == econ_data.MOLAR_MASS["water"]
    with pytest.raises(ValueError, match="notachemicalxyz"):
        econ_data.molar_mass("notachemicalxyz")


def test_costing_a_flowsheet_with_a_component_outside_the_old_table():
    """End-to-end: a methanol process must size/cost/price without methanol
    ever being added to economics.data.MOLAR_MASS (the old KeyError risk)."""
    fs = Flowsheet(components=[Component("methanol")], property_package="thermo:PR")
    fs.add(Heater("H", {"T_out": 320.0}))
    fs.feed("F", "H:in1", T=300.0, P=101325.0, molar_flow=10.0, z={"methanol": 1.0})
    fs.connect("PROD", "H:out", None)
    fs.connect("Q", "H:duty", None)
    report = fs.solve()
    res = analyze(fs, report, TEAConfig(
        product_component="methanol",
        prices_per_kg={"methanol": 0.45},
    ))
    # 10 mol/s * 0.032 kg/mol * 8000 h -> ~9.2e6 kg/yr
    assert res.annual_production_kg == pytest.approx(
        10.0 * 0.0320419 * 8000 * 3600, rel=1e-3)
    assert res.annual_revenue > 0
    assert res.capital.grassroots > 0
