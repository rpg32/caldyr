// Case study: sweep one unit parameter over a range, solving at each point
// client-side, and chart any metric against it.
import { Download, Play, Square } from "lucide-react";
import { useRef, useState } from "react";
import {
  CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip, XAxis, YAxis,
} from "recharts";
import { api } from "../api";
import { downloadCsv } from "../lib/csv";
import { dimFor, dimForMetric } from "../lib/params";
import { defaultUnit, fmtDim, toDisplay } from "../lib/units";
import { useStore } from "../store";
import type { FlowDoc, MetricSpec, SolveResponse } from "../types";
import { MetricEditor, metricLabel, metricValid } from "./MetricEditor";
import { Button, DimField, Hint, PanelTitle } from "./ui";

const sel = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-text min-w-0";
const num = "w-[76px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text";

interface Point { x: number; y: number | null }

function evalMetric(res: SolveResponse, m: MetricSpec): number | null {
  if (m.type === "duty") return res.report.duties[m.stream] ?? null;
  const s = res.streams[m.stream];
  if (!s || s.molar_flow == null) return null;
  if (m.type === "flow") return s.molar_flow;
  const zsum = Object.values(s.z).reduce((a, b) => a + b, 0) || 1;
  return s.molar_flow * ((s.z[m.component ?? ""] ?? 0) / zsum);
}

/** Deep-copied flow doc with one unit param overridden. */
function withParam(flow: FlowDoc, unit: string, param: string, value: number): FlowDoc {
  const doc = structuredClone(flow) as FlowDoc & {
    units: { id: string; params: Record<string, unknown> }[];
  };
  const u = doc.units.find((x) => x.id === unit);
  if (u) u.params = { ...u.params, [param]: value };
  return doc;
}

export function StudyPanel() {
  const nodes = useStore((s) => s.nodes);
  const toFlowDoc = useStore((s) => s.toFlowDoc);
  const backend = useStore((s) => s.backend);
  const toast = useStore((s) => s.toast);
  const unitSet = useStore((s) => s.unitSet);

  const [unit, setUnit] = useState("");
  const [param, setParam] = useState("");
  const [from, setFrom] = useState(0);
  const [to, setTo] = useState(1);
  const [steps, setSteps] = useState(8);
  const [metric, setMetric] = useState<MetricSpec>({ type: "flow", stream: "" });
  const [points, setPoints] = useState<Point[]>([]);
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState("");
  const abortRef = useRef(false);

  const unitIds = nodes.filter((n) => n.data.kind === "unit").map((n) => n.id);
  const paramsOf = (u: string): string[] => {
    const n = nodes.find((x) => x.id === u);
    return Object.entries(n?.data.params ?? {})
      .filter(([, v]) => typeof v === "number").map(([k]) => k);
  };

  const valid = unit && param && metricValid(metric)
    && Number.isFinite(from) && Number.isFinite(to) && from !== to
    && steps >= 2 && steps <= 50;

  // The swept param is stored/solved in SI; convert only for display.
  const dim = dimFor(param);
  const xUnit = dim ? ` ${defaultUnit(dim, unitSet)}` : "";
  const fmtX = (x: number, digits = 4) =>
    (dim ? fmtDim(dim, x, unitSet, digits) : x.toPrecision(digits)) + xUnit;
  // The watched metric also has a dimension (mol/s for flows, W for duties).
  const yDim = dimForMetric(metric.type);
  const fmtY = (y: number, digits = 4) =>
    yDim ? fmtDim(yDim, y, unitSet, digits) : y.toPrecision(digits);
  const yUnit = yDim ? ` / ${defaultUnit(yDim, unitSet)}` : "";

  const run = async () => {
    setRunning(true);
    setPoints([]);
    abortRef.current = false;
    const base = toFlowDoc();
    const out: Point[] = [];
    let failures = 0;
    for (let i = 0; i < steps; i++) {
      if (abortRef.current) break;
      const x = from + ((to - from) * i) / (steps - 1);
      setProgress(`${i + 1} / ${steps}  (${param} = ${fmtX(x)})`);
      try {
        const res = await api.solve(withParam(base, unit, param, x), backend);
        out.push({ x, y: res.report.converged ? evalMetric(res, metric) : null });
        if (!res.report.converged) failures++;
      } catch {
        out.push({ x, y: null });
        failures++;
      }
      setPoints([...out]);
    }
    setProgress("");
    setRunning(false);
    if (failures) toast("info", `${failures} of ${steps} sweep points failed to solve (shown as gaps).`);
  };

  const exportCsv = () => downloadCsv(
    [[`${unit}.${param}${dim ? ` (${defaultUnit(dim, unitSet)})` : ""}`, metricLabel(metric)],
     ...points.map((p) => [dim ? toDisplay(dim, p.x, unitSet) : p.x, p.y])],
    "case_study.csv",
  );

  return (
    <div>
      <PanelTitle>Vary</PanelTitle>
      <div className="space-y-1.5 text-[12px]">
        <label className="flex items-center gap-2">
          <span className="w-[68px] shrink-0 text-muted">Unit</span>
          <select className={`${sel} flex-1`} value={unit} aria-label="Unit"
            onChange={(e) => { setUnit(e.target.value); setParam(""); }}>
            <option value="">— unit —</option>
            {unitIds.map((u) => <option key={u} value={u}>{u}</option>)}
          </select>
        </label>
        <label className="flex items-center gap-2">
          <span className="w-[68px] shrink-0 text-muted">Parameter</span>
          <select className={`${sel} flex-1`} value={param} aria-label="Parameter"
            onChange={(e) => setParam(e.target.value)}>
            <option value="">— param —</option>
            {paramsOf(unit).map((p) => <option key={p} value={p}>{p}</option>)}
          </select>
        </label>
        <div className="flex items-center gap-2">
          <span className="w-[68px] shrink-0 text-muted">Range</span>
          <DimField dim={dim} set={unitSet} className={num} value={from} aria-label="From" onChange={setFrom} />
          <span className="text-muted">→</span>
          <DimField dim={dim} set={unitSet} className={num} value={to} aria-label="To" onChange={setTo} />
          {dim && <span className="text-[11px] text-muted">{defaultUnit(dim, unitSet)}</span>}
        </div>
        <label className="flex items-center gap-2">
          <span className="w-[68px] shrink-0 text-muted">Steps</span>
          <input className="w-[64px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text"
            type="number" min={2} max={50} value={steps} aria-label="Steps"
            onChange={(e) => setSteps(parseInt(e.target.value) || 2)} />
        </label>
      </div>

      <PanelTitle>Watch</PanelTitle>
      <MetricEditor value={metric} onChange={setMetric} />

      <div className="mt-3 flex items-center gap-2">
        {!running ? (
          <Button variant="primary" icon={<Play size={13} />} disabled={!valid}
            onClick={() => void run()}>
            Run sweep
          </Button>
        ) : (
          <Button icon={<Square size={13} />} onClick={() => { abortRef.current = true; }}>
            Stop
          </Button>
        )}
        {points.length > 0 && !running && (
          <Button variant="ghost" icon={<Download size={13} />} onClick={exportCsv}>CSV</Button>
        )}
        <span className="text-[11px] text-muted">{progress}</span>
      </div>
      {!valid && (
        <Hint>Pick a unit parameter, a numeric range, and a metric to watch.
          Each step is a full engine solve.</Hint>
      )}

      {points.length > 1 && (
        <>
          <PanelTitle>{metricLabel(metric)} vs {unit}.{param}{dim ? ` / ${defaultUnit(dim, unitSet)}` : ""}</PanelTitle>
          <div style={{ height: 180 }}>
            <ResponsiveContainer>
              <LineChart data={points.map((p) => ({ x: dim ? toDisplay(dim, p.x, unitSet) : p.x, y: p.y }))}
                margin={{ top: 6, right: 10, left: 0, bottom: 0 }}>
                <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
                <XAxis dataKey="x" type="number" domain={["auto", "auto"]}
                  tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                  tickFormatter={(v: number) => v.toPrecision(3)} />
                <YAxis tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
                  width={56} domain={["auto", "auto"]}
                  tickFormatter={(v: number) => Number(v).toPrecision(3)} />
                <Tooltip
                  contentStyle={{
                    background: "var(--panel)", border: "1px solid var(--line)",
                    borderRadius: 6, fontSize: 11,
                  }}
                  labelFormatter={(v) => `${param} = ${Number(v).toPrecision(5)}${xUnit}`}
                />
                <Line type="monotone" dataKey="y" stroke="var(--accent)" strokeWidth={1.6}
                  dot={{ r: 2.5 }} connectNulls={false} isAnimationActive={false} />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </>
      )}

      {points.length > 0 && (
        <>
          <PanelTitle>Results</PanelTitle>
          <div className="max-h-56 overflow-auto rounded-md border border-line">
            <table className="data-table">
              <thead>
                <tr>
                  <th>{param || "x"}{dim ? ` / ${defaultUnit(dim, unitSet)}` : ""}</th>
                  <th>{metricLabel(metric)}{yUnit}</th>
                </tr>
              </thead>
              <tbody>
                {points.map((p, i) => (
                  <tr key={i}>
                    <td>{dim ? fmtDim(dim, p.x, unitSet, 4) : p.x.toPrecision(4)}</td>
                    <td>{p.y == null ? <span className="text-muted">—</span> : fmtY(p.y)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </div>
  );
}
