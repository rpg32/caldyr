import { Dices, SlidersHorizontal } from "lucide-react";
import { useMemo, useState } from "react";
import {
  Bar, BarChart, CartesianGrid, ReferenceLine, ResponsiveContainer, Tooltip,
  XAxis, YAxis,
} from "recharts";
import { useStore } from "../store";
import { AssumptionsEditor, type AssumptionsSeed } from "./AssumptionsEditor";
import { Button, Hint, PanelTitle, StaleNotice } from "./ui";

const fmt = (x: number | null | undefined, d = 2) =>
  x == null ? "—" : x.toLocaleString(undefined, { maximumFractionDigits: d });
const money = (x: number) => "$" + Math.round(x).toLocaleString();

const TOOLTIP_STYLE = {
  background: "var(--panel)", border: "1px solid var(--line)",
  borderRadius: 6, fontSize: 11,
} as const;

/** Interactive tornado: bars span low->high LCOP around the base value. */
function TornadoChart({ bars, base }: {
  bars: { variable: string; low_lcop: number; high_lcop: number; swing: number }[];
  base: number;
}) {
  const data = useMemo(() =>
    [...bars].sort((a, b) => b.swing - a.swing).map((b) => ({
      variable: b.variable,
      // floating bar: [start, end] relative span
      span: [Math.min(b.low_lcop, b.high_lcop), Math.max(b.low_lcop, b.high_lcop)],
      low: b.low_lcop, high: b.high_lcop, swing: b.swing,
    })), [bars]);
  return (
    <div style={{ height: 40 + data.length * 30 }}>
      <ResponsiveContainer>
        <BarChart data={data} layout="vertical" margin={{ top: 4, right: 12, left: 8, bottom: 0 }}>
          <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" horizontal={false} />
          <XAxis type="number" domain={["auto", "auto"]}
            tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
            tickFormatter={(v: number) => `$${v.toFixed(2)}`} />
          <YAxis type="category" dataKey="variable" width={104}
            tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)" />
          <Tooltip
            contentStyle={TOOLTIP_STYLE}
            formatter={(_v, _n, item) => {
              const p = item.payload as typeof data[number];
              return [`$${p.low.toFixed(3)} → $${p.high.toFixed(3)} (swing $${p.swing.toFixed(3)})`, "LCOP"];
            }}
          />
          <Bar dataKey="span" fill="var(--accent)" radius={3} isAnimationActive={false} />
          <ReferenceLine x={base} stroke="var(--warn)" strokeDasharray="4 3"
            label={{ value: "base", fill: "var(--warn)", fontSize: 10, position: "top" }} />
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}

/** LCOP distribution from Monte-Carlo samples. */
function McHistogram({ samples, p10, p50, p90 }: {
  samples: number[]; p10: number; p50: number; p90: number;
}) {
  const data = useMemo(() => {
    if (samples.length < 10) return [];
    const lo = Math.min(...samples);
    const hi = Math.max(...samples);
    const nBins = 24;
    const w = (hi - lo) / nBins || 1;
    const bins = Array.from({ length: nBins }, (_, i) => ({
      x: lo + (i + 0.5) * w, count: 0,
    }));
    for (const s of samples) {
      const i = Math.min(nBins - 1, Math.floor((s - lo) / w));
      bins[i].count++;
    }
    return bins;
  }, [samples]);
  if (!data.length) return null;
  return (
    <>
      <div className="mb-1 flex gap-3 text-[11px] text-muted">
        <span>P10 <b className="text-text">${p10.toFixed(3)}</b></span>
        <span>P50 <b className="text-text">${p50.toFixed(3)}</b></span>
        <span>P90 <b className="text-text">${p90.toFixed(3)}</b></span>
      </div>
      <div style={{ height: 120 }}>
        <ResponsiveContainer>
          <BarChart data={data} margin={{ top: 4, right: 8, left: 0, bottom: 0 }} barCategoryGap={1}>
            <XAxis dataKey="x" tick={{ fill: "var(--muted)", fontSize: 10 }}
              stroke="var(--line)" tickFormatter={(v: number) => `$${v.toFixed(2)}`} />
            <YAxis tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)" width={30} />
            <Tooltip contentStyle={TOOLTIP_STYLE}
              labelFormatter={(v) => `LCOP ≈ $${Number(v).toFixed(3)}/kg`} />
            <Bar dataKey="count" fill="var(--accent)" isAnimationActive={false} />
            <ReferenceLine x={p50} stroke="var(--warn)" strokeDasharray="4 3" />
          </BarChart>
        </ResponsiveContainer>
      </div>
    </>
  );
}

/** "Assumptions in this solve" — the prices, factors and heuristics that drove
 *  this result, editable in place (seeded with the values actually used). The
 *  Settings dialog offers the same editor seeded from defaults, available any
 *  time (even before a cost). */
function AssumptionsSection() {
  const res = useStore((s) => s.costRes);
  const toggleSettings = useStore((s) => s.toggleSettings);
  const [open, setOpen] = useState(false);
  const a = res?.assumptions;
  if (!a) return null;

  const seed: AssumptionsSeed = {
    financial: a.config, prices: a.prices_per_kg, utilities: a.utility_prices,
    sizing: a.sizing, factors: a.factors, citations: a.citations,
  };
  return (
    <div className="mt-3">
      <button className="flex w-full items-center gap-1.5 text-[11px] uppercase tracking-wider text-muted hover:text-text"
        onClick={() => setOpen(!open)}>
        <SlidersHorizontal size={12} />
        Assumptions in this solve {open ? "▾" : "▸"}
      </button>
      {open && (
        <div className="mt-1.5">
          <button className="mb-2 text-[11px] text-accent hover:underline"
            onClick={toggleSettings}>Open Settings (edit any time) →</button>
          <AssumptionsEditor seed={seed} />
        </div>
      )}
    </div>
  );
}

export function EconomicsPanel() {
  const res = useStore((s) => s.costRes);
  const stale = useStore((s) => s.resultsStale);
  const busy = useStore((s) => s.busy);
  const cost = useStore((s) => s.cost);
  if (!res) return <Hint>Press Cost for a techno-economic analysis.</Hint>;
  return (
    <div>
      <StaleNotice stale={stale} />
      <div className="mb-2.5 grid grid-cols-2 gap-2">
        {([
          ["LCOP", `$${res.profitability.lcop.toFixed(3)}/kg`],
          ["TCI", money(res.capital.tci)],
          ["OPEX/yr", money(res.opex.total)],
          ["NPV", money(res.profitability.npv)],
        ] as const).map(([k, v]) => (
          <div key={k} className="rounded-lg border border-line bg-panel2 p-2">
            <span className="block text-[11px] text-muted">{k}</span>
            <b className="text-[15px]">{v}</b>
          </div>
        ))}
      </div>
      <PanelTitle>Equipment (installed)</PanelTitle>
      <table className="data-table">
        <thead><tr><th>unit</th><th>type</th><th>size</th><th>Cbm</th></tr></thead>
        <tbody>
          {res.equipment.map((e, i) => (
            <tr key={`${e.unit_id}-${i}`}>
              <td>{e.unit_id}</td><td>{e.type}</td>
              <td>{fmt(e.attribute)} {e.attribute_name.split("_")[0]}</td>
              <td>{money(e.bare_module)}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {res.tornado && res.tornado.length > 0 && (
        <>
          <PanelTitle>LCOP sensitivity</PanelTitle>
          <TornadoChart bars={res.tornado} base={res.profitability.lcop} />
        </>
      )}
      <PanelTitle>Uncertainty (Monte-Carlo)</PanelTitle>
      {res.monte_carlo ? (
        <McHistogram
          samples={res.monte_carlo.lcop_samples}
          p10={res.monte_carlo.lcop.p10}
          p50={res.monte_carlo.lcop.p50}
          p90={res.monte_carlo.lcop.p90}
        />
      ) : (
        <Button icon={<Dices size={13} />} busy={busy === "cost"}
          onClick={() => void cost(500)}>
          Run 500 samples
        </Button>
      )}
      <AssumptionsSection />
    </div>
  );
}
