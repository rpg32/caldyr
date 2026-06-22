// Analysis tools panel: HYSYS-style Property Table, relief-valve sizing
// (API 520/526), and heat-integration pinch targeting. Thin clients of the
// engine's /property-table, /relief and /pinch endpoints — no physics here.
import { Activity, FlaskConical, Gauge, Play } from "lucide-react";
import { useState } from "react";
import {
  CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { api } from "../api";
import { useStore } from "../store";
import type { FlowDoc, PinchResponse, PropertyTableResponse, ReliefResponse } from "../types";
import { Badge, Button, Hint, PanelTitle } from "./ui";

const num =
  "w-[80px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text";
const sel = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-text min-w-0";

const PROPS = ["mass_density", "molar_volume", "enthalpy", "entropy", "vapor_fraction"];

function linspace(a: number, b: number, n: number): number[] {
  if (n <= 1) return [a];
  return Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));
}

export function AnalysisPanel() {
  return (
    <div>
      <PropertyTableTool />
      <PinchTool />
      <ReliefTool />
    </div>
  );
}

// -- Property table ----------------------------------------------------------
function PropertyTableTool() {
  const edges = useStore((s) => s.edges);
  const solveRes = useStore((s) => s.solveRes);
  const toFlowDoc = useStore((s) => s.toFlowDoc);
  const toast = useStore((s) => s.toast);

  const streamIds = edges.map((e) => e.id);
  const [stream, setStream] = useState("");
  const [tMin, setTMin] = useState(300);
  const [tMax, setTMax] = useState(500);
  const [tN, setTN] = useState(6);
  const [pKpa, setPKpa] = useState("101.325, 500");
  const [prop, setProp] = useState("mass_density");
  const [res, setRes] = useState<PropertyTableResponse | null>(null);
  const [busy, setBusy] = useState(false);

  const run = async () => {
    setBusy(true);
    try {
      const T = linspace(tMin, tMax, Math.max(1, tN));
      const P = pKpa.split(",").map((s) => parseFloat(s.trim()) * 1000).filter(Number.isFinite);
      if (!P.length) throw new Error("enter one or more pressures (kPa)");
      const doc: FlowDoc = { ...toFlowDoc(), solved: solveRes?.streams };
      setRes(await api.propertyTable(doc, stream, T, P, PROPS));
      setProp((p) => p);
    } catch (e) {
      toast("error", `Property table: ${(e as Error).message}`);
    } finally {
      setBusy(false);
    }
  };

  const chartData =
    res && res.values[prop]
      ? res.T.map((t, i) => {
          const row: Record<string, number | null> = { T: t };
          res.P.forEach((p, j) => { row[`${(p / 1000).toFixed(0)} kPa`] = res.values[prop][i][j]; });
          return row;
        })
      : [];
  const series = res ? res.P.map((p) => `${(p / 1000).toFixed(0)} kPa`) : [];
  const colors = ["#3b82f6", "#f97316", "#10b981", "#a855f7", "#ef4444"];

  return (
    <section className="mb-4">
      <PanelTitle><FlaskConical size={11} className="mr-1 inline" /> Property table</PanelTitle>
      <div className="flex flex-wrap items-center gap-1.5">
        <select className={sel} value={stream} aria-label="Stream"
          onChange={(e) => setStream(e.target.value)}>
          <option value="">— stream —</option>
          {streamIds.map((s) => <option key={s} value={s}>{s}</option>)}
        </select>
      </div>
      <div className="mt-1.5 flex flex-wrap items-center gap-1.5 text-[11px] text-muted">
        T/K
        <input className={num} type="number" step="any" value={tMin} aria-label="T min"
          onChange={(e) => setTMin(parseFloat(e.target.value))} />
        →
        <input className={num} type="number" step="any" value={tMax} aria-label="T max"
          onChange={(e) => setTMax(parseFloat(e.target.value))} />
        ×
        <input className="w-[48px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text"
          type="number" min={1} max={40} value={tN} aria-label="T steps"
          onChange={(e) => setTN(parseInt(e.target.value) || 1)} />
      </div>
      <div className="mt-1.5 flex flex-wrap items-center gap-1.5 text-[11px] text-muted">
        P/kPa
        <input className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-1.5 py-1 text-text"
          value={pKpa} aria-label="Pressures (kPa, comma-separated)"
          onChange={(e) => setPKpa(e.target.value)} />
      </div>
      <div className="mt-2 flex items-center gap-2">
        <Button variant="primary" icon={<Play size={13} />} busy={busy}
          disabled={!stream} onClick={() => void run()}>Compute</Button>
        {res && (
          <select className={sel} value={prop} aria-label="Property"
            onChange={(e) => setProp(e.target.value)}>
            {res.props.map((p) => <option key={p} value={p}>{p}</option>)}
          </select>
        )}
      </div>
      {!stream && (
        <Hint>Pick a stream (solve first for internal streams) and a T×P grid —
          the HYSYS Property Table over the engine's property package.</Hint>
      )}
      {res && chartData.length > 1 && (
        <div className="mt-2" style={{ height: 180 }}>
          <ResponsiveContainer>
            <LineChart data={chartData} margin={{ top: 6, right: 10, left: 0, bottom: 0 }}>
              <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
              <XAxis dataKey="T" type="number" domain={["auto", "auto"]}
                tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                tickFormatter={(v: number) => v.toFixed(0)} />
              <YAxis tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                width={56} domain={["auto", "auto"]}
                tickFormatter={(v: number) => Number(v).toPrecision(3)} />
              <Tooltip contentStyle={{
                background: "var(--panel)", border: "1px solid var(--line)",
                borderRadius: 6, fontSize: 11 }}
                labelFormatter={(v) => `T = ${Number(v).toFixed(1)} K`} />
              {series.map((s, i) => (
                <Line key={s} type="monotone" dataKey={s} stroke={colors[i % colors.length]}
                  strokeWidth={1.6} dot={{ r: 2 }} connectNulls isAnimationActive={false} />
              ))}
            </LineChart>
          </ResponsiveContainer>
          <div className="flex flex-wrap gap-3 text-[10px] text-muted">
            {series.map((s, i) => (
              <span key={s}><i className="inline-block h-2 w-2 rounded-sm"
                style={{ background: colors[i % colors.length] }} /> {s}</span>
            ))}
          </div>
          {res.failures.length > 0 && (
            <div className="mt-1 text-[11px] text-warn">
              {res.failures.length} grid point(s) failed to flash (gaps).
            </div>
          )}
        </div>
      )}
    </section>
  );
}

// -- Pinch -------------------------------------------------------------------
function PinchTool() {
  const toFlowDoc = useStore((s) => s.toFlowDoc);
  const backend = useStore((s) => s.backend);
  const toast = useStore((s) => s.toast);
  const [dtMin, setDtMin] = useState(10);
  const [res, setRes] = useState<PinchResponse | null>(null);
  const [busy, setBusy] = useState(false);

  const run = async () => {
    setBusy(true);
    try {
      setRes(await api.pinch(toFlowDoc(), backend, dtMin));
    } catch (e) {
      toast("error", `Pinch: ${(e as Error).message}`);
    } finally {
      setBusy(false);
    }
  };

  const curve = (pts: [number, number][]) =>
    pts.map(([H, T]) => ({ H: H / 1e6, T }));   // MW, K

  return (
    <section className="mb-4 border-t border-line pt-3">
      <PanelTitle><Activity size={11} className="mr-1 inline" /> Heat integration (pinch)</PanelTitle>
      <div className="flex items-center gap-2 text-[11px] text-muted">
        ΔT min / K
        <input className={num} type="number" step="any" min={0.1} value={dtMin} aria-label="dt_min"
          onChange={(e) => setDtMin(parseFloat(e.target.value))} />
        <Button variant="primary" icon={<Play size={13} />} busy={busy}
          onClick={() => void run()}>Target</Button>
      </div>
      {!res && (
        <Hint>Solve targets from the heaters/coolers/exchangers in the flowsheet:
          minimum hot &amp; cold utility and the pinch temperature.</Hint>
      )}
      {res && (
        <>
          <table className="data-table mt-2">
            <tbody>
              <tr><td>min hot utility</td><td>{(res.qh_min / 1e6).toFixed(3)} MW</td></tr>
              <tr><td>min cold utility</td><td>{(res.qc_min / 1e6).toFixed(3)} MW</td></tr>
              <tr><td>recovery potential</td><td>{(res.heat_recovery_potential / 1e6).toFixed(3)} MW</td></tr>
              <tr><td>pinch T (hot / cold)</td><td>
                {res.pinch_T_hot != null ? `${res.pinch_T_hot.toFixed(1)} / ${res.pinch_T_cold?.toFixed(1)} K` : "— (threshold)"}
              </td></tr>
            </tbody>
          </table>
          {(res.hot_composite.length > 1 || res.cold_composite.length > 1) && (
            <div className="mt-2" style={{ height: 180 }}>
              <ResponsiveContainer>
                <LineChart margin={{ top: 6, right: 10, left: 0, bottom: 14 }}>
                  <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
                  <XAxis dataKey="H" type="number" domain={["auto", "auto"]}
                    tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                    label={{ value: "H / MW", fill: "var(--muted)", fontSize: 10, position: "insideBottom", dy: 12 }}
                    tickFormatter={(v: number) => v.toFixed(1)} allowDuplicatedCategory={false} />
                  <YAxis dataKey="T" type="number" domain={["auto", "auto"]} width={52}
                    tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                    label={{ value: "T / K", fill: "var(--muted)", fontSize: 10, angle: -90, position: "insideLeft" }}
                    tickFormatter={(v: number) => v.toFixed(0)} />
                  <Tooltip contentStyle={{
                    background: "var(--panel)", border: "1px solid var(--line)",
                    borderRadius: 6, fontSize: 11 }} />
                  <Line data={curve(res.hot_composite)} name="hot" dataKey="T" stroke="#ef4444"
                    strokeWidth={1.6} dot={false} isAnimationActive={false} />
                  <Line data={curve(res.cold_composite)} name="cold" dataKey="T" stroke="#3b82f6"
                    strokeWidth={1.6} dot={false} isAnimationActive={false} />
                </LineChart>
              </ResponsiveContainer>
              <div className="flex gap-3 text-[10px] text-muted">
                <span><i className="inline-block h-2 w-2 rounded-sm" style={{ background: "#ef4444" }} /> hot composite</span>
                <span><i className="inline-block h-2 w-2 rounded-sm" style={{ background: "#3b82f6" }} /> cold composite</span>
              </div>
            </div>
          )}
        </>
      )}
    </section>
  );
}

// -- Relief valve ------------------------------------------------------------
function ReliefTool() {
  const toast = useStore((s) => s.toast);
  const [phase, setPhase] = useState<"vapor" | "liquid">("vapor");
  const [W, setW] = useState(10);
  const [P1kpa, setP1kpa] = useState(1200);   // kPa abs
  const [T, setT] = useState(500);
  const [M, setM] = useState(0.029);          // kg/mol
  const [Z, setZ] = useState(0.95);
  const [k, setK] = useState(1.3);
  const [rho, setRho] = useState(800);
  const [P2kpa, setP2kpa] = useState(101.325);
  const [res, setRes] = useState<ReliefResponse | null>(null);
  const [busy, setBusy] = useState(false);

  const run = async () => {
    setBusy(true);
    try {
      const body =
        phase === "vapor"
          ? { phase, W, T, M, Z, k, P1: P1kpa * 1000 }
          : { phase, W, rho, P1: P1kpa * 1000, P2: P2kpa * 1000 };
      setRes(await api.relief(body));
    } catch (e) {
      toast("error", `Relief: ${(e as Error).message}`);
    } finally {
      setBusy(false);
    }
  };

  const F = (label: string, v: number, set: (n: number) => void, unit = "") => (
    <label className="flex items-center justify-between gap-2 text-[11px] text-muted">
      <span>{label}</span>
      <span className="flex items-center gap-1">
        <input className={num} type="number" step="any" value={v} aria-label={label}
          onChange={(e) => set(parseFloat(e.target.value))} />
        <span className="w-10">{unit}</span>
      </span>
    </label>
  );

  return (
    <section className="border-t border-line pt-3">
      <PanelTitle><Gauge size={11} className="mr-1 inline" /> Relief valve (API 520/526)</PanelTitle>
      <div className="mb-1.5 flex items-center gap-2 text-[11px] text-muted">
        phase
        <select className={sel} value={phase} aria-label="Relief phase"
          onChange={(e) => { setPhase(e.target.value as "vapor" | "liquid"); setRes(null); }}>
          <option value="vapor">vapor</option>
          <option value="liquid">liquid</option>
        </select>
      </div>
      <div className="grid grid-cols-1 gap-1">
        {F("W flow", W, setW, "kg/s")}
        {F("P1 (abs)", P1kpa, setP1kpa, "kPa")}
        {phase === "vapor" ? (
          <>
            {F("T", T, setT, "K")}
            {F("M", M, setM, "kg/mol")}
            {F("Z", Z, setZ, "")}
            {F("k=Cp/Cv", k, setK, "")}
          </>
        ) : (
          <>
            {F("ρ", rho, setRho, "kg/m³")}
            {F("P2 (abs)", P2kpa, setP2kpa, "kPa")}
          </>
        )}
      </div>
      <div className="mt-2">
        <Button variant="primary" icon={<Play size={13} />} busy={busy}
          onClick={() => void run()}>Size</Button>
      </div>
      {res && (
        <table className="data-table mt-2">
          <tbody>
            <tr><td>required area</td><td>{res.area_cm2.toFixed(2)} cm²</td></tr>
            <tr><td>API 526 orifice</td><td>
              {res.orifice ? <Badge ok>{res.orifice}</Badge> : <Badge ok={false}>none — parallel valves</Badge>}
            </td></tr>
            {res.capacity_used != null && (
              <tr><td>capacity used</td><td>{(res.capacity_used * 100).toFixed(1)} %</td></tr>
            )}
            {res.critical != null && (
              <tr><td>flow regime</td><td>{res.critical ? "critical (choked)" : "subcritical"}</td></tr>
            )}
          </tbody>
        </table>
      )}
      {res?.notes.map((n) => (
        <div key={n} className="mt-1 text-[11px] text-warn">{n}</div>
      ))}
    </section>
  );
}
