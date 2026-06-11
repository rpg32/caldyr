"""Post-solve mass & energy closure report.

``balance_report(fs)`` audits a *solved* flowsheet: for every unit and for the
flowsheet boundary it totals mass in/out (kg/s) and energy in/out (W, enthalpy
flows plus energy-port duties) and reports the imbalance, worst offenders
first. A healthy flowsheet closes to ~solver tolerance (machine precision for
acyclic sweeps); a big imbalance pinpoints a modeling bug or an unwired port.

Conventions:
  * energy-port duties follow the engine sign convention (positive = energy
    added to the process), so they are booked on the *in* side of each unit and
    of the overall balance;
  * relative imbalance is ``|out - in| / max(|in|, |out|)`` (mass), and for
    energy is normalized by the total *enthalpy throughput* rather than the
    net, because formation-inclusive enthalpies legitimately sit far from zero.

Molar masses come from the `chemicals` database (CAS preferred, else the
component id) — the same identifier path the property packages use.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

_TINY = 1e-300


@lru_cache(maxsize=256)
def _molar_mass(key: str) -> float:
    """kg/mol from the `chemicals` identifier database."""
    from chemicals.identifiers import search_chemical
    found = search_chemical(key)
    if found is None:
        raise ValueError(f"balance_report: unknown chemical identifier {key!r}")
    return float(found.MW) / 1000.0


def _mw_map(fs) -> dict[str, float]:
    out: dict[str, float] = {}
    for c in fs.components:
        try:
            out[c.id] = _molar_mass(c.cas or c.id)
        except Exception as exc:
            raise ValueError(
                f"balance_report: cannot resolve a molar mass for component "
                f"{c.id!r}: {exc}"
            ) from exc
    return out


def _mass_kg_s(stream, mw: dict[str, float]) -> float:
    n = stream.molar_flow or 0.0
    if n <= 0.0:
        return 0.0
    return n * sum(frac * mw[c] for c, frac in stream.normalized_z().items())


def _enthalpy_W(stream) -> float:
    n = stream.molar_flow or 0.0
    if n <= 0.0:
        return 0.0
    if stream.H is None:
        raise ValueError(
            f"balance_report: stream {stream.id!r} has no resolved enthalpy; "
            f"solve the flowsheet first"
        )
    return n * stream.H


def _unit_duty(fs, report, unit) -> float:
    """Signed duty (W) over the unit's energy ports, wired or not (an unwired
    port falls back to the engine's default ``unit.port`` stream id)."""
    total = 0.0
    for p in unit.ports:
        if p.kind != "energy":
            continue
        conn = next((c for c in fs.connections
                     if c.from_unit == unit.id and c.from_port == p.name), None)
        sid = conn.stream_id if conn else f"{unit.id}.{p.name}"
        total += report.duties.get(sid, 0.0)
    return total


def balance_report(fs, report=None) -> dict[str, Any]:
    """Mass/energy closure for a solved flowsheet. ``report`` defaults to the
    flowsheet's last solve report (``fs.last_report``, set by ``fs.solve()``).

    Returns ``{"overall": {...}, "units": [...], "warnings": [...]}`` where
    ``units`` entries are sorted worst-offender-first by relative imbalance.
    """
    report = report if report is not None else getattr(fs, "last_report", None)
    if report is None:
        raise ValueError(
            "balance_report needs a solved flowsheet: call fs.solve() first "
            "(or pass the SolveReport explicitly)"
        )
    mw = _mw_map(fs)
    warnings: list[str] = []

    units: list[dict[str, Any]] = []
    for uid, unit in fs.units.items():
        m_in = m_out = h_in = h_out = 0.0
        wired_out = set()
        for conn in fs.connections:
            s = fs.streams.get(conn.stream_id)
            if conn.to_unit == uid and s is not None:
                m_in += _mass_kg_s(s, mw)
                h_in += _enthalpy_W(s)
            if conn.from_unit == uid:
                wired_out.add(conn.from_port)
                if s is not None:
                    m_out += _mass_kg_s(s, mw)
                    h_out += _enthalpy_W(s)
        for p in unit.ports:
            if p.kind == "material" and p.direction == "outlet" and p.name not in wired_out:
                warnings.append(
                    f"{uid}: material outlet port {p.name!r} is unwired; its flow "
                    f"is missing from the balance"
                )
        duty = _unit_duty(fs, report, unit)
        e_in = h_in + duty
        m_imb, e_imb = m_out - m_in, h_out - e_in
        units.append({
            "unit_id": uid,
            "mass_in_kg_s": m_in, "mass_out_kg_s": m_out,
            "mass_imbalance_kg_s": m_imb,
            "mass_rel": abs(m_imb) / max(m_in, m_out, _TINY),
            "duty_W": duty,
            "energy_in_W": e_in, "energy_out_W": h_out,
            "energy_imbalance_W": e_imb,
            "energy_rel": abs(e_imb) / max(abs(h_in), abs(h_out), _TINY),
        })
    units.sort(key=lambda u: max(u["mass_rel"], u["energy_rel"]), reverse=True)

    m_in = m_out = h_in = h_out = 0.0
    for conn in fs.connections:
        s = fs.streams.get(conn.stream_id)
        if s is None:
            continue
        if conn.from_unit is None and conn.to_unit is not None:      # feed
            m_in += _mass_kg_s(s, mw)
            h_in += _enthalpy_W(s)
        elif conn.to_unit is None and conn.from_unit is not None:    # product
            m_out += _mass_kg_s(s, mw)
            h_out += _enthalpy_W(s)
    duty_total = sum(report.duties.values())
    e_in = h_in + duty_total
    m_imb, e_imb = m_out - m_in, h_out - e_in
    overall = {
        "mass_in_kg_s": m_in, "mass_out_kg_s": m_out,
        "mass_imbalance_kg_s": m_imb,
        "mass_rel": abs(m_imb) / max(m_in, m_out, _TINY),
        "duty_W": duty_total,
        "energy_in_W": e_in, "energy_out_W": h_out,
        "energy_imbalance_W": e_imb,
        "energy_rel": abs(e_imb) / max(abs(h_in), abs(h_out), _TINY),
    }
    return {"overall": overall, "units": units, "warnings": warnings}
