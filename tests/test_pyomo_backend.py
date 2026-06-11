"""Pyomo grey-box EO backend tests.

Two tiers, matching the backend's honest availability handling:

* **Model construction** (needs only ``pip install pyomo``, no NLP solver
  binary): the Pyomo model builds, variables/constraints are counted, and the
  grey-box residuals/Jacobian evaluate. The strongest construction-level check:
  the grey-box residual is ~0 at the *sequential* backend's solution — i.e. the
  Pyomo translation layer encodes exactly the same physics.

* **Live solve** (needs a grey-box-capable NLP solver, i.e. cyipopt): the
  flash-recycle flowsheet solved with ``backend="pyomo"`` agrees with
  ``backend="sequential"`` to 1e-6 — the same cross-check style as
  test_m4_equation_oriented. Skips with the backend's own unavailability
  reason when cyipopt cannot be installed (pip-only Windows: no wheels, and
  the sdist needs Ipopt dev headers; the ASL ipopt.exe cannot consume
  grey-box models). When unavailable, the typed error path is tested instead.
"""
import math

import pytest

pytest.importorskip("pyomo", reason="pyomo not installed (pip install 'caldyr[eo]')")

import numpy as np  # noqa: E402

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.solver import PyomoEOSolver, PyomoSolverUnavailable  # noqa: E402
from caldyr.solver.pyomo_backend import solver_unavailable_reason  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import FlashDrum, Mixer, Splitter  # noqa: E402

UNAVAILABLE = solver_unavailable_reason()   # None when a live solve is possible


def flash_recycle() -> Flowsheet:
    """Same flash-with-recycle flowsheet as the M4 cross-check tests."""
    fs = Flowsheet(components=[Component("n-pentane"), Component("n-octane")],
                   property_package="thermo:PR")
    fs.add(Mixer("MIX", {"dP": 0.0}))
    fs.add(FlashDrum("FL", {"T": 360.0, "P": 101325.0}))
    fs.add(Splitter("SP", {"split": 0.6}))
    fs.feed("FEED", "MIX:in1", T=330.0, P=101325.0, molar_flow=10.0,
            z={"n-pentane": 0.5, "n-octane": 0.5})
    fs.connect("MIXOUT", "MIX:out", "FL:in1")
    fs.connect("VAP", "FL:vapor", None)
    fs.connect("LIQ", "FL:liquid", "SP:in1")
    fs.connect("RECY", "SP:out1", "MIX:in2")
    fs.connect("BOT", "SP:out2", None)
    fs.connect("Q", "FL:duty", None)
    return fs


def _build_model():
    fs = flash_recycle()
    pp = make_package(fs.property_package, fs.component_ids)
    model, system = PyomoEOSolver().build_model(fs, pp)
    return fs, model, system


# -- model construction (no NLP solver binary required) ---------------------
def test_model_builds_with_expected_counts():
    fs, model, system = _build_model()
    # 5 computed material streams (MIXOUT, VAP, LIQ, RECY, BOT) x
    # [n_pentane, n_octane, T, P] = 20 unknowns and 20 equality constraints.
    assert system.n == 20
    egb = model.residuals.get_external_model()
    assert egb.n_inputs() == 20
    assert egb.n_equality_constraints() == 20
    assert egb.n_outputs() == 0
    # The Pyomo block exposes one input variable per unknown, initialized at
    # the warm start (scaled values, so all finite).
    pyomo_x = [model.residuals.inputs[name].value for name in system.names]
    assert len(pyomo_x) == 20
    assert np.isfinite(pyomo_x).all()
    # Square feasibility problem: a (zero) objective is present, making the
    # model optimization-ready (swap the objective, free design variables).
    assert model.obj.expr() == 0.0


def test_greybox_residual_and_jacobian_evaluate():
    fs, model, system = _build_model()
    egb = model.residuals.get_external_model()
    egb.set_input_values(system.x0)
    r = egb.evaluate_equality_constraints()
    assert r.shape == (20,)
    assert np.isfinite(r).all()
    jac = egb.evaluate_jacobian_equality_constraints()
    assert jac.shape == (20, 20)
    dense = jac.toarray()
    assert np.isfinite(dense).all()
    # The FD Jacobian at the warm start must be nonsingular (Newton-solvable).
    assert np.linalg.matrix_rank(dense) == 20


def test_greybox_residual_is_zero_at_sequential_solution():
    """The translation layer encodes the same equations: packing the sequential
    backend's converged streams into the grey-box variables zeroes the residual."""
    fs_seq = flash_recycle()
    rep = fs_seq.solve(backend="sequential", tol=1e-10)
    assert rep.converged

    fs, model, system = _build_model()
    x = system.layout.pack(fs_seq.streams) / system.scales
    egb = model.residuals.get_external_model()
    egb.set_input_values(x)
    r = egb.evaluate_equality_constraints()
    assert float(np.max(np.abs(r))) < 1e-6


# -- live solve (requires cyipopt) -------------------------------------------
@pytest.mark.skipif(UNAVAILABLE is not None, reason=UNAVAILABLE or "")
def test_pyomo_matches_sequential_on_flash_recycle():
    fa = flash_recycle()
    ra = fa.solve(backend="sequential", tol=1e-9)
    fb = flash_recycle()
    rb = fb.solve(backend="pyomo", tol=1e-8)
    assert ra.converged and rb.converged
    assert rb.method == "pyomo"
    assert rb.tear_streams == []                 # no tearing in the EO formulation
    for sid in ("VAP", "BOT", "RECY"):
        sa, sb = fa.streams[sid], fb.streams[sid]
        assert math.isclose(sa.molar_flow, sb.molar_flow, rel_tol=1e-6, abs_tol=1e-9)
        assert math.isclose(sa.T, sb.T, rel_tol=1e-6, abs_tol=1e-6)
        for c in sa.components:
            assert math.isclose(sa.z.get(c, 0.0), sb.z.get(c, 0.0),
                                rel_tol=1e-6, abs_tol=1e-7)
    assert math.isclose(ra.duties["Q"], rb.duties["Q"], rel_tol=1e-5, abs_tol=1.0)
    # Optimization-ready: the solved model is a square grey-box NLP with a free
    # objective slot — replace Objective(expr=0) with e.g. the flash duty and
    # unfix FL.T to optimize subject to the full flowsheet equations in IPOPT.


@pytest.mark.skipif(UNAVAILABLE is None, reason="grey-box NLP solver IS available")
def test_solve_raises_typed_error_when_no_solver():
    fs = flash_recycle()
    with pytest.raises(PyomoSolverUnavailable) as excinfo:
        fs.solve(backend="pyomo")
    # The error must say what to install.
    assert "cyipopt" in str(excinfo.value) or "pyomo" in str(excinfo.value)
