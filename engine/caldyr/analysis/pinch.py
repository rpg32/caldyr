"""Pinch analysis / heat-integration targeting.

Implements the classic **problem-table (temperature-interval) algorithm** of
Linnhoff & Flower (AIChE J., 1978) as presented in Kemp, *Pinch Analysis and
Process Integration*, 2e (2007), Ch. 2, and Smith, *Chemical Process Design and
Integration* (2005), Ch. 16: given the process's hot streams (releasing heat)
and cold streams (absorbing heat) and a minimum approach ``dt_min``, compute the
minimum hot- and cold-utility targets, the pinch temperature, and the composite
curves — *before* designing any exchanger network.

Two entry points:

* :func:`pinch_from_streams` — operate on a plain list of stream specs
  (``{"Tin": ..., "Tout": ..., "Q": ...}`` dicts or :class:`ThermalStream`),
  so targeting is usable and testable without a flowsheet.
* :func:`pinch_analysis` — extract the thermal streams from a *solved*
  flowsheet and target those. Extraction reads the **duty-carrying units**
  (each utility heating/cooling demand defines one stream): every unit with an
  energy ``duty`` port and a nonzero duty contributes a stream spanning its
  material inlet -> outlet temperatures, with constant CP = Q/|Tin - Tout|;
  distillation columns contribute their condenser (hot) and reboiler (cold)
  duties from ``unit.design``. Near-isothermal duties (reboilers, condensers
  at ~constant T, isothermal flash/reactor duties) are widened to 1 K segments
  (``ISOTHERMAL_WIDTH_K``) so they cascade cleanly — a standard problem-table
  device (Kemp 2e, Sec. 3.5: latent-heat streams as small-dT segments).

  Process-process HeatExchanger units carry no utility duty and are treated as
  *already integrated*: the targets returned are over the residual utility
  demands ("remaining problem"), so ``heat_recovery_potential`` is the further
  recovery available on top of any existing matches.

Shaft-work ports (``work`` on pumps/compressors/expanders) are not thermal
duties and are excluded.
"""
from __future__ import annotations

from dataclasses import dataclass, field

#: Width given to (near-)isothermal duties so they enter the problem table as
#: finite-CP segments (see module docstring).
ISOTHERMAL_WIDTH_K = 1.0

#: Duties smaller than this (W) are ignored during flowsheet extraction.
MIN_DUTY_W = 1e-6


class PinchExtractionError(ValueError):
    """A duty-carrying unit could not be turned into a thermal stream (usually
    an unsolved flowsheet or a missing inlet/outlet stream)."""


@dataclass(frozen=True)
class ThermalStream:
    """One heating or cooling demand: ``Q`` watts released (hot) or absorbed
    (cold) over ``T_in -> T_out``, at constant CP = Q/|T_in - T_out|."""
    T_in: float                  # K, supply temperature
    T_out: float                 # K, target temperature
    Q: float                     # W, magnitude (always > 0)
    kind: str                    # "hot" (releases heat) | "cold" (absorbs heat)
    unit_id: str = ""

    @property
    def CP(self) -> float:
        """Heat-capacity flow rate, W/K."""
        return self.Q / abs(self.T_in - self.T_out)


@dataclass
class PinchResult:
    """Heat-integration targets plus plot-ready composite curves.

    Composite-curve points are ``(H, T)`` pairs (cumulative enthalpy in W,
    actual — not shifted — temperature in K); the cold curve is offset so it
    starts at ``H = qc_min``, the standard alignment at which the closest
    vertical approach between the curves equals ``dt_min``.
    """
    dt_min: float
    qh_min: float                          # minimum hot utility, W
    qc_min: float                          # minimum cold utility, W
    pinch_T_shifted: float | None          # shifted (interval) pinch T, K
    pinch_T_hot: float | None              # = shifted + dt_min/2
    pinch_T_cold: float | None             # = shifted - dt_min/2
    current_hot_utility: float             # actual heating in use, W
    current_cold_utility: float            # actual cooling in use, W
    heat_recovery_potential: float         # current hot utility - qh_min, W
    hot_composite: list[tuple[float, float]] = field(default_factory=list)
    cold_composite: list[tuple[float, float]] = field(default_factory=list)
    streams: list[ThermalStream] = field(default_factory=list)


# -- stream normalization ----------------------------------------------------
def _normalize(spec, index: int) -> ThermalStream | None:
    """Turn a dict spec or ThermalStream into a canonical ThermalStream.

    Dict keys: ``Tin``/``Tout`` (or ``T_in``/``T_out``) and ``Q``. The stream
    direction decides hot vs cold (Tin > Tout releases heat); for isothermal
    entries the engine's duty sign convention decides (Q > 0 = heat absorbed
    by the process = cold stream). Returns None for negligible duties.
    """
    if isinstance(spec, ThermalStream):
        t_in, t_out, q, unit_id = spec.T_in, spec.T_out, spec.Q, spec.unit_id
        kind = spec.kind
        if kind not in ("hot", "cold"):
            raise ValueError(f"stream {index}: kind must be 'hot' or 'cold', got {kind!r}")
        q = abs(q)
    else:
        try:
            t_in = float(spec["Tin"] if "Tin" in spec else spec["T_in"])
            t_out = float(spec["Tout"] if "Tout" in spec else spec["T_out"])
            q_signed = float(spec["Q"])
        except (KeyError, TypeError) as exc:
            raise ValueError(
                f"stream {index}: expected keys Tin/Tout/Q (or T_in/T_out/Q), got {spec!r}"
            ) from exc
        unit_id = str(spec.get("unit_id", f"stream{index}"))
        if t_in > t_out:
            kind = "hot"
        elif t_in < t_out:
            kind = "cold"
        else:
            kind = "cold" if q_signed > 0 else "hot"
        q = abs(q_signed)

    if q <= MIN_DUTY_W:
        return None
    # Widen near-isothermal duties to a 1 K segment so CP stays finite: a hot
    # stream condenses downwards, a cold stream boils upwards in temperature.
    if abs(t_in - t_out) < ISOTHERMAL_WIDTH_K:
        t_out = t_in - ISOTHERMAL_WIDTH_K if kind == "hot" else t_in + ISOTHERMAL_WIDTH_K
    # Canonical orientation (hot: high -> low, cold: low -> high). Only the
    # span [Tlo, Thi], CP and kind enter the problem table, so a stream whose
    # duty sign disagrees with its temperature direction (e.g. an
    # isothermal-spec reactor fed off-temperature) is simply reoriented.
    lo, hi = sorted((t_in, t_out))
    t_in, t_out = (hi, lo) if kind == "hot" else (lo, hi)
    return ThermalStream(T_in=t_in, T_out=t_out, Q=q, kind=kind, unit_id=unit_id)


# -- the problem-table algorithm ---------------------------------------------
def _shifted_range(s: ThermalStream, dt_min: float) -> tuple[float, float]:
    """(low, high) shifted-temperature span: hot streams -dt/2, cold +dt/2."""
    shift = -dt_min / 2.0 if s.kind == "hot" else +dt_min / 2.0
    lo, hi = sorted((s.T_in + shift, s.T_out + shift))
    return lo, hi


def _composite(streams: list[ThermalStream]) -> list[tuple[float, float]]:
    """Composite curve for one kind: (cumulative H from the cold end, T)."""
    if not streams:
        return []
    temps = sorted({t for s in streams for t in (s.T_in, s.T_out)})
    points = [(0.0, temps[0])]
    h = 0.0
    for lo, hi in zip(temps[:-1], temps[1:]):
        cp = sum(s.CP for s in streams
                 if min(s.T_in, s.T_out) <= lo and max(s.T_in, s.T_out) >= hi)
        h += cp * (hi - lo)
        points.append((h, hi))
    return points


def pinch_from_streams(streams, dt_min: float = 10.0) -> PinchResult:
    """Run the problem-table algorithm on a plain list of thermal streams.

    ``streams`` is a sequence of :class:`ThermalStream` or dicts with keys
    ``Tin``, ``Tout``, ``Q`` (temperatures in K, duty in W; see
    :func:`_normalize` for the hot/cold convention). Validated against the
    four-stream example of Kemp, *Pinch Analysis and Process Integration*, 2e
    (2007), Table 2.1 (QHmin = 20 kW, QCmin = 60 kW, pinch at 85 C shifted).
    """
    if dt_min < 0.0:
        raise ValueError(f"dt_min must be >= 0, got {dt_min}")
    norm = [s for i, spec in enumerate(streams)
            if (s := _normalize(spec, i)) is not None]
    hot = [s for s in norm if s.kind == "hot"]
    cold = [s for s in norm if s.kind == "cold"]
    current_hot = sum(s.Q for s in cold)     # heating the process buys today
    current_cold = sum(s.Q for s in hot)     # cooling the process buys today

    if not norm:
        return PinchResult(dt_min=dt_min, qh_min=0.0, qc_min=0.0,
                           pinch_T_shifted=None, pinch_T_hot=None, pinch_T_cold=None,
                           current_hot_utility=0.0, current_cold_utility=0.0,
                           heat_recovery_potential=0.0, streams=norm)

    # Temperature intervals on the shifted scale.
    spans = {s: _shifted_range(s, dt_min) for s in norm}
    bounds = sorted({t for span in spans.values() for t in span}, reverse=True)

    # Problem-table cascade: surplus per interval, cascaded from the top.
    surpluses: list[float] = []
    for hi, lo in zip(bounds[:-1], bounds[1:]):
        net_cp = 0.0
        for s in norm:
            s_lo, s_hi = spans[s]
            if s_lo <= lo + 1e-12 and s_hi >= hi - 1e-12:
                net_cp += s.CP if s.kind == "hot" else -s.CP
        surpluses.append(net_cp * (hi - lo))

    cascade = [0.0]
    for q in surpluses:
        cascade.append(cascade[-1] + q)
    qh_min = max(0.0, -min(cascade))
    flows = [qh_min + c for c in cascade]    # heat flow at each boundary
    qc_min = flows[-1]

    # Pinch: interior boundary where the cascaded heat flow is (numerically)
    # zero. A threshold problem (zero flow only at an end, or nowhere) has none.
    scale = max(qh_min, qc_min, max(s.Q for s in norm))
    pinch_shifted = next(
        (bounds[i] for i in range(1, len(bounds) - 1) if flows[i] <= 1e-9 * scale),
        None,
    )

    return PinchResult(
        dt_min=dt_min,
        qh_min=qh_min,
        qc_min=qc_min,
        pinch_T_shifted=pinch_shifted,
        pinch_T_hot=None if pinch_shifted is None else pinch_shifted + dt_min / 2.0,
        pinch_T_cold=None if pinch_shifted is None else pinch_shifted - dt_min / 2.0,
        current_hot_utility=current_hot,
        current_cold_utility=current_cold,
        heat_recovery_potential=current_hot - qh_min,
        hot_composite=_composite(hot),
        cold_composite=[(h + qc_min, t) for h, t in _composite(cold)],
        streams=norm,
    )


# -- flowsheet extraction ----------------------------------------------------
def _unit_streams(fs, uid):
    ins, outs = {}, {}
    for conn in fs.connections:
        if conn.to_unit == uid and conn.stream_id in fs.streams:
            ins[conn.to_port] = fs.streams[conn.stream_id]
        if conn.from_unit == uid and conn.stream_id in fs.streams:
            outs[conn.from_port] = fs.streams[conn.stream_id]
    return ins, outs


def _port_duty(fs, report, unit, port_name: str) -> float:
    """Duty (W) on one energy port, wired (its stream id) or not (unit.port)."""
    conn = next((c for c in fs.connections
                 if c.from_unit == unit.id and c.from_port == port_name), None)
    sid = conn.stream_id if conn else f"{unit.id}.{port_name}"
    return report.duties.get(sid, 0.0)


def extract_thermal_streams(fs, report) -> list[ThermalStream]:
    """Pull the per-unit thermal demands out of a solved flowsheet (see the
    module docstring for what counts as a thermal stream)."""
    streams: list[ThermalStream] = []
    for uid, unit in fs.units.items():
        design = getattr(unit, "design", None)
        port_names = {p.name for p in unit.ports if p.kind == "energy"}

        # Distillation columns: condenser + reboiler from the design results
        # (their duty ports carry the same numbers, but the design dict also
        # carries the temperatures the segments need).
        if design and "Q_condenser" in design and "Q_reboiler" in design:
            q_cond = abs(design["Q_condenser"])
            if q_cond > MIN_DUTY_W:
                streams.append(ThermalStream(
                    T_in=max(design["T_top_dew"], design["T_top"]),
                    T_out=min(design["T_top_dew"], design["T_top"]),
                    Q=q_cond, kind="hot", unit_id=f"{uid}.condenser"))
            q_reb = abs(design["Q_reboiler"])
            if q_reb > MIN_DUTY_W:
                streams.append(ThermalStream(
                    T_in=design["T_bottom"], T_out=design["T_bottom"],
                    Q=q_reb, kind="cold", unit_id=f"{uid}.reboiler"))
            continue

        # Generic duty-carrying units (Heater, FiredHeater, AirCooler, Flash,
        # reactors, ...): thermal `duty` ports only — `work` is shaft power.
        thermal_ports = sorted(port_names - {"work"})
        if not thermal_ports:
            continue
        q = sum(_port_duty(fs, report, unit, p) for p in thermal_ports)
        if abs(q) <= MIN_DUTY_W:
            continue
        ins, outs = _unit_streams(fs, uid)
        t_ins = [s.T for s in ins.values() if s.T is not None]
        t_outs = [s.T for s in outs.values() if s.T is not None]
        if not t_ins or not t_outs:
            raise PinchExtractionError(
                f"unit {uid!r} has duty {q:.3g} W but no resolved inlet/outlet "
                f"stream temperatures — solve the flowsheet before pinch analysis"
            )
        streams.append(ThermalStream(
            T_in=t_ins[0], T_out=t_outs[0], Q=abs(q),
            kind="cold" if q > 0 else "hot", unit_id=uid))

    # Re-normalize: widens isothermal segments and orients each span.
    out: list[ThermalStream] = []
    for i, s in enumerate(streams):
        norm = _normalize(s, i)
        if norm is not None:
            out.append(norm)
    return out


def pinch_analysis(fs, report, dt_min: float = 10.0) -> PinchResult:
    """Heat-integration targets for a *solved* flowsheet: extract the thermal
    streams from its duty-carrying units, then run :func:`pinch_from_streams`.

    ``heat_recovery_potential`` is how much of the current hot-utility duty
    could be displaced by process-process heat exchange at ``dt_min``.
    """
    return pinch_from_streams(extract_thermal_streams(fs, report), dt_min)
