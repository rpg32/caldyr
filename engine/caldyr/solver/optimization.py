"""Flowsheet optimization: minimize an objective over design variables.

This is the optimization half of the equation-oriented story (the thing
commercial sequential simulators do worst). Design variables are unit parameters
(a heater's target temperature, a splitter's split, ...); the objective and
constraints are callables of the *solved* flowsheet. The flowsheet is re-solved
inside the optimizer (scipy SLSQP), so any solve backend and the full property
physics are available — and warm starts make repeated solves cheap.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize


@dataclass
class DesignVar:
    unit_id: str
    param: str
    lower: float
    upper: float
    initial: float | None = None        # defaults to the midpoint

    def start(self) -> float:
        return self.initial if self.initial is not None else 0.5 * (self.lower + self.upper)


@dataclass
class OptimizeResult:
    success: bool
    objective: float
    design: dict[str, float]            # "unit.param" -> value
    n_solves: int
    message: str
    report: object = field(default=None)


def optimize(fs, objective, design_vars: list[DesignVar], *,
             constraints=(), backend: str = "sequential", solve_kwargs: dict | None = None,
             maxiter: int = 50) -> OptimizeResult:
    """Minimize ``objective(fs, report)`` over ``design_vars`` subject to
    inequality ``constraints`` (each ``g(fs, report) >= 0``). Returns the optimum
    with the flowsheet left in its optimal solved state."""
    solve_kwargs = solve_kwargs or {}
    lows = np.array([d.lower for d in design_vars])
    highs = np.array([d.upper for d in design_vars])
    span = highs - lows
    counter = {"n": 0}

    def apply_and_solve(x_unit: np.ndarray):
        values = lows + x_unit * span                # map [0,1] -> [lower, upper]
        for d, v in zip(design_vars, values):
            fs.units[d.unit_id].params[d.param] = float(v)
        report = fs.solve(backend=backend, **solve_kwargs)
        counter["n"] += 1
        return report

    x0 = np.array([(d.start() - lo) / s if s else 0.5
                   for d, lo, s in zip(design_vars, lows, span)])

    # Normalize the objective by its magnitude at the start so SLSQP is robust to
    # the caller's units (W vs kW vs $); the reported objective is unscaled.
    f0 = float(objective(fs, apply_and_solve(x0)))
    obj_scale = max(abs(f0), 1e-12)

    def f(x):
        report = apply_and_solve(x)
        return float(objective(fs, report)) / obj_scale

    cons = [{"type": "ineq", "fun": lambda x, g=g: float(g(fs, apply_and_solve(x)))}
            for g in constraints]

    sol = minimize(f, x0, method="SLSQP", bounds=[(0.0, 1.0)] * len(design_vars),
                   constraints=cons, options={"maxiter": maxiter, "ftol": 1e-7})

    report = apply_and_solve(sol.x)                  # leave fs at the optimum
    values = lows + sol.x * span
    return OptimizeResult(
        success=bool(sol.success), objective=float(objective(fs, report)),
        design={f"{d.unit_id}.{d.param}": float(v) for d, v in zip(design_vars, values)},
        n_solves=counter["n"], message=str(sol.message), report=report,
    )
