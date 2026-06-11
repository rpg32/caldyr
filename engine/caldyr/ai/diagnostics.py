"""Flowsheet explanation & convergence diagnosis — pure-data tools.

These return structured dicts (plus a readable ``summary``) so they work three
ways: as AI tools the LLM verbalizes, as direct REST fallbacks when no local
LLM is running, and in tests. No LLM calls happen here.
"""
from __future__ import annotations

from typing import Any


def describe_flowsheet(fs, report: Any = None) -> dict:
    """Narrate the topology and (if solved) the material story of a flowsheet."""
    units = [
        {"id": uid, "type": type(u).__name__, "params": dict(u.params)}
        for uid, u in fs.units.items()
    ]
    feeds: list[dict] = []
    products: list[dict] = []
    internal: list[dict] = []
    for conn in fs.connections:
        entry = {
            "stream": conn.stream_id,
            "from": f"{conn.from_unit}:{conn.from_port}" if conn.from_unit else None,
            "to": f"{conn.to_unit}:{conn.to_port}" if conn.to_unit else None,
        }
        s = fs.streams.get(conn.stream_id)
        if s is not None and s.molar_flow is not None:
            entry["state"] = {
                "T": s.T, "P": s.P, "molar_flow": s.molar_flow,
                "phase": s.phase, "z": s.normalized_z(),
            }
        if conn.from_unit is None:
            feeds.append(entry)
        elif conn.to_unit is None:
            internal_port = next(
                (p for p in fs.units[conn.from_unit].ports
                 if p.name == conn.from_port), None)
            (internal if internal_port is not None and internal_port.kind == "energy"
             else products).append(entry)
        else:
            internal.append(entry)

    solved = report is not None or any(
        s.H is not None for s in fs.streams.values())
    lines = [
        f"Flowsheet: {len(units)} units ({', '.join(u['id'] for u in units)}), "
        f"{len(feeds)} feed(s), {len(products)} product(s); "
        f"property package {fs.property_package}; "
        f"components {list(fs.component_ids)}.",
    ]
    if fs.logical:
        kinds = [op.get("type") for op in fs.logical]
        lines.append(f"Logical ops: {len(fs.logical)} ({', '.join(kinds)}).")
    if solved:
        for p in products:
            st = p.get("state")
            if st:
                top = max(st["z"], key=st["z"].get) if st["z"] else "?"
                lines.append(
                    f"Product {p['stream']}: {st['molar_flow']:.3g} mol/s at "
                    f"{st['T']:.1f} K, mostly {top} "
                    f"({st['z'].get(top, 0):.1%})."
                )
    else:
        lines.append("Not solved yet — stream states unknown.")

    return {
        "summary": " ".join(lines),
        "units": units,
        "feeds": feeds,
        "products": products,
        "streams": internal,
        "logical": list(fs.logical),
        "property_package": fs.property_package,
        "components": list(fs.component_ids),
        "solved": solved,
    }


def explain_convergence(fs, report) -> dict:
    """Diagnose a solve: what happened, why, and what to try next."""
    if report is None:
        return {
            "summary": "No solve has been run yet — call solve first.",
            "converged": None, "advice": ["Run solve."],
        }

    history = list(getattr(report, "history", []) or [])
    advice: list[str] = []
    trend = "n/a"
    if len(history) >= 3:
        drops = sum(1 for a, b in zip(history, history[1:]) if b < a)
        ratio = drops / (len(history) - 1)
        if ratio > 0.8:
            trend = "monotone decrease"
        elif ratio > 0.5:
            trend = "decreasing with oscillation"
        else:
            trend = "oscillating / stalled"

    if report.converged:
        if len(history) > 30:
            advice.append(
                "Converged but slowly — seed the tear streams via solver_hints "
                "tear_guesses (e.g. from this solved state) to start warm next time.")
        if trend == "decreasing with oscillation":
            advice.append(
                "Oscillatory approach — Wegstein is already damping it; the "
                "equation-oriented backend may converge in fewer evaluations.")
    else:
        advice.append(
            "Not converged. First check the spec set: every unit needs enough "
            "parameters (e.g. Heater needs T_out or Q; Flash needs P and "
            "optionally T).")
        if trend == "oscillating / stalled":
            advice.append(
                "Residual oscillates — try method='direct' at a smaller step, "
                "provide tear_guesses in solver_hints, or switch to the "
                "equation_oriented backend (no tearing).")
        if history and history[-1] > 1.0:
            advice.append(
                "Residual is still large — the recycle may be physically "
                "unstable (e.g. accumulating inerts without a purge).")
        advice.append("Increase max_iter only after the above; it rarely helps alone.")

    adjust_msgs = [m for m in report.messages if m.startswith(("adjust", "set"))]
    parts = [
        f"{'Converged' if report.converged else 'NOT converged'} via "
        f"{report.method} in {report.iterations} iteration(s)"
        + (f", final residual {report.residual:.3g} (tol {report.tol:g})"
           if report.residual is not None else "") + ".",
    ]
    if report.tear_streams:
        parts.append(f"Tear stream(s): {', '.join(report.tear_streams)}; "
                     f"residual trend: {trend}.")
    if adjust_msgs:
        parts.append("Logical ops: " + "; ".join(adjust_msgs) + ".")
    if advice:
        parts.append("Advice: " + " ".join(advice))

    return {
        "summary": " ".join(parts),
        "converged": report.converged,
        "method": report.method,
        "iterations": report.iterations,
        "residual": report.residual,
        "tolerance": report.tol,
        "tear_streams": list(report.tear_streams),
        "residual_history": history,
        "trend": trend,
        "messages": list(report.messages),
        "advice": advice,
    }
