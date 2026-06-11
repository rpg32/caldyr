"""Pyomo equation-oriented backend (grey-box over the flowsheet residual).

This backend exposes the *same* scaled residual system the scipy EO backend
solves — for every computed material stream, ``x_outlet − unit.solve(x_inlets)``
in scaled ``[n_1..n_C, T, P]`` block coordinates — to Pyomo, by wrapping the
flowsheet sweep as a PyNumero :class:`ExternalGreyBoxModel` attached to a
``ConcreteModel`` through ``ExternalGreyBoxBlock``. The physics stay identical
to the other two backends (the residual calls the very same ``unit.solve`` and
``PropertyPackage``), so all three converge to the same flowsheet state. The
Jacobian of the grey box is built by forward finite differences on the scaled
residual.

Solver reality (especially on pip-only Windows): a Pyomo model containing an
``ExternalGreyBoxBlock`` has no NL-file representation, so it can **only** be
solved through Pyomo's ``cyipopt`` interface — the ASL ``ipopt.exe`` (e.g. from
``idaes get-extensions``) cannot consume it. ``pip install cyipopt`` ships no
Windows wheels and its sdist needs the Ipopt development headers, so when
cyipopt (or the PyNumero ASL library) is missing this backend raises a typed
:class:`PyomoSolverUnavailable` with install guidance. Model *construction* and
residual/Jacobian *evaluation* are pure Python and work with pyomo alone — see
:meth:`PyomoEOSolver.build_model` and the tests that exercise it.

Because the model is a square feasibility problem (zero objective), it is
optimization-ready: swap the zero objective for a real one, free design
variables, and IPOPT optimizes subject to the full flowsheet equations.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .equation_oriented import EquationOrientedSolver, _Layout
from .sequential import SolveReport

_INSTALL_PYOMO = (
    "pyomo is not installed. Install it with `pip install pyomo` "
    "(or `pip install 'caldyr[eo]'`)."
)
_INSTALL_CYIPOPT = (
    "no grey-box-capable NLP solver is available: Pyomo models with an "
    "ExternalGreyBoxBlock can only be solved through the 'cyipopt' interface, and "
    "`import cyipopt` failed. The ASL ipopt.exe (e.g. from `idaes get-extensions`) "
    "cannot consume grey-box models - they have no NL-file representation. Install "
    "cyipopt with `conda install -c conda-forge cyipopt` (recommended); `pip install "
    "cyipopt` ships no Windows wheels and needs the Ipopt development headers "
    "(an Ipopt binary release unpacked next to its setup.py plus MSVC) to build from "
    "source. Model construction still works - only the solve step needs cyipopt."
)
_INSTALL_ASL = (
    "the PyNumero ASL library is not available, so Pyomo cannot assemble the "
    "grey-box NLP. Run `idaes get-extensions` (after `pip install idaes-pse`) or "
    "`pyomo download-extensions` to install the precompiled binaries."
)


class PyomoSolverUnavailable(RuntimeError):
    """Raised when the Pyomo backend cannot solve on this machine (pyomo missing,
    or no grey-box-capable NLP solver such as cyipopt). The message says exactly
    what to install."""


def solver_unavailable_reason() -> str | None:
    """Why a live ``backend="pyomo"`` solve cannot run here, or ``None`` if it can.

    Checks, in order: pyomo importable; cyipopt importable (the only interface
    that can solve grey-box models); the PyNumero ASL library present (importing
    ``idaes`` first, if installed, puts its downloaded binaries on the search path).
    """
    try:
        import pyomo.environ  # noqa: F401
    except ImportError:
        return _INSTALL_PYOMO
    try:
        import idaes  # noqa: F401  (registers the IDAES binary dir, incl. pynumero ASL)
    except ImportError:
        pass
    try:
        import cyipopt  # noqa: F401
    except ImportError:
        return _INSTALL_CYIPOPT
    from pyomo.contrib.pynumero.asl import AmplInterface
    if not AmplInterface.available():
        return _INSTALL_ASL
    return None


class FlowsheetResidualSystem:
    """The scaled EO residual system ``r(x) = 0`` (identical to
    :mod:`.equation_oriented`'s) with a finite-difference Jacobian.

    Pure Python/numpy — no Pyomo dependency — so it is testable and usable even
    when no NLP solver binary exists. The grey box delegates to it.
    """

    def __init__(self, flowsheet: Any, pp: Any, tol: float = 1e-8,
                 fd_step: float = 1e-7) -> None:
        self._fs = flowsheet
        self._pp = pp
        self._fd_step = fd_step
        self._eo = EquationOrientedSolver(tol=tol)
        self.layout: _Layout = self._eo._layout(flowsheet, pp)

        # Boundary feeds are fixed inputs: resolve once, then warm-start the
        # unknowns with one forward sweep (tears seeded empty) — same recipe as
        # the scipy EO backend, so both start from the same point.
        for conn in flowsheet.connections:
            if conn.from_unit is None and conn.to_unit is not None:
                self._eo._resolve_feed(flowsheet.streams[conn.stream_id], pp)
        self._eo._warm_start(flowsheet, self.layout, pp)

        self.scales: np.ndarray = self.layout.scales()
        self.n: int = self.layout.size
        labels = [f"n[{c}]" for c in self.layout.components] + ["T", "P"]
        self.names: list[str] = [f"{sid}.{lab}" for sid in self.layout.stream_ids
                                 for lab in labels]
        self.x0: np.ndarray = self.layout.pack(flowsheet.streams) / self.scales
        self.n_residual_evals = 0

    def residual(self, x_scaled: np.ndarray) -> np.ndarray:
        """``(computed − target)/scales`` — exactly the scipy EO residual."""
        lay, fs = self.layout, self._fs
        x = np.asarray(x_scaled, dtype=float) * self.scales
        for i, sid in enumerate(lay.stream_ids):
            fs.streams[sid] = lay.unpack_one(x, i)
        computed = np.empty_like(x)
        for i, sid in enumerate(lay.stream_ids):
            computed[i * lay.block:(i + 1) * lay.block] = \
                self._eo._computed_block(fs, sid, lay, self._pp)
        self.n_residual_evals += 1
        return (computed - x) / self.scales

    def jacobian(self, x_scaled: np.ndarray) -> np.ndarray:
        """Dense forward-finite-difference Jacobian of :meth:`residual`."""
        x = np.asarray(x_scaled, dtype=float)
        r0 = self.residual(x)
        jac = np.empty((self.n, self.n))
        for j in range(self.n):
            h = self._fd_step * max(1.0, abs(x[j]))
            xj = x.copy()
            xj[j] += h
            jac[:, j] = (self.residual(xj) - r0) / h
        return jac

    def commit(self, x_scaled: np.ndarray) -> dict[str, float]:
        """Write the solution into the flowsheet; final pass for duties/H."""
        x = np.asarray(x_scaled, dtype=float) * self.scales
        for i, sid in enumerate(self.layout.stream_ids):
            self._fs.streams[sid] = self.layout.unpack_one(x, i)
        return self._eo._finalize(self._fs, self._pp)


def _make_greybox(system: FlowsheetResidualSystem) -> Any:
    """Wrap a residual system as a PyNumero ExternalGreyBoxModel (lazy import so
    this module stays importable without pyomo)."""
    from pyomo.contrib.pynumero.interfaces.external_grey_box import ExternalGreyBoxModel
    from scipy.sparse import coo_matrix

    class _FlowsheetGreyBox(ExternalGreyBoxModel):  # type: ignore[misc]
        def __init__(self, sys: FlowsheetResidualSystem) -> None:
            self._sys = sys
            self._x = sys.x0.copy()

        def input_names(self) -> list[str]:
            return list(self._sys.names)

        def equality_constraint_names(self) -> list[str]:
            return [f"r[{name}]" for name in self._sys.names]

        def set_input_values(self, input_values: Any) -> None:
            self._x = np.asarray(input_values, dtype=float)

        def evaluate_equality_constraints(self) -> np.ndarray:
            return self._sys.residual(self._x)

        def evaluate_jacobian_equality_constraints(self) -> Any:
            return coo_matrix(self._sys.jacobian(self._x))

        def finalize_block_construction(self, pyomo_block: Any) -> None:
            # Initialize the Pyomo input variables at the warm start.
            for name, val in zip(self._sys.names, self._sys.x0):
                pyomo_block.inputs[name].set_value(float(val))

    return _FlowsheetGreyBox(system)


class PyomoEOSolver:
    """Equation-oriented backend on the Pyomo/PyNumero grey-box stack.

    ``solve`` matches the other backends' contract (returns a
    :class:`~caldyr.solver.sequential.SolveReport`, mutates ``flowsheet.streams``)
    and is selected with ``fs.solve(backend="pyomo")``. Raises
    :class:`PyomoSolverUnavailable` when no grey-box-capable NLP solver exists.
    """

    def __init__(self, tol: float = 1e-8, max_iter: int = 200) -> None:
        self.tol = tol
        self.max_iter = max_iter

    def build_model(self, flowsheet: Any, pp: Any) -> tuple[Any, FlowsheetResidualSystem]:
        """Build the Pyomo model: a square feasibility problem (zero objective)
        whose only constraints are the grey-box flowsheet residuals. Needs pyomo
        but no solver binary — usable for inspection and tests on any machine."""
        try:
            from pyomo.environ import ConcreteModel, Objective
            from pyomo.contrib.pynumero.interfaces.external_grey_box import (
                ExternalGreyBoxBlock,
            )
        except ImportError as exc:
            raise PyomoSolverUnavailable(_INSTALL_PYOMO) from exc

        system = FlowsheetResidualSystem(flowsheet, pp, tol=self.tol)
        model = ConcreteModel(name="caldyr_flowsheet_eo")
        model.residuals = ExternalGreyBoxBlock(external_model=_make_greybox(system))
        model.obj = Objective(expr=0.0)   # square system: find the feasible point
        return model, system

    def solve(self, flowsheet: Any, pp: Any) -> SolveReport:
        reason = solver_unavailable_reason()
        if reason is not None:
            raise PyomoSolverUnavailable(reason)

        model, system = self.build_model(flowsheet, pp)

        from pyomo.environ import SolverFactory, value
        from pyomo.opt import check_optimal_termination

        solver = SolverFactory("cyipopt")
        # The grey box provides no Hessians -> quasi-Newton in IPOPT.
        solver.config.options["hessian_approximation"] = "limited-memory"
        solver.config.options["tol"] = self.tol
        solver.config.options["max_iter"] = self.max_iter
        results = solver.solve(model, tee=False)

        x = np.array([value(model.residuals.inputs[name]) for name in system.names])
        res = float(np.max(np.abs(system.residual(x))))
        duties = system.commit(x)
        term = str(results.solver.termination_condition)
        converged = bool(check_optimal_termination(results)) and \
            res < max(self.tol * 100, 1e-6)
        return SolveReport(
            converged=converged, iterations=system.n_residual_evals, residual=res,
            tol=self.tol, method="pyomo",
            order=list(flowsheet.units), tear_streams=[],
            duties=duties,
            messages=[
                f"pyomo grey-box NLP ({system.n} vars, {system.n} equality "
                f"constraints) solved with cyipopt: {term}; "
                f"{system.n_residual_evals} residual evals",
            ],
        )
