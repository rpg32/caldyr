// The HYSYS-Spreadsheet analog: named expressions over the solved flowsheet,
// re-evaluated live after every solve. Rows can reference earlier rows.
import { Plus, Trash2 } from "lucide-react";
import { useMemo } from "react";
import { evaluateCalcs } from "../lib/calc";
import { useStore } from "../store";
import { Hint, PanelTitle, StaleNotice } from "./ui";

const fmtVal = (v: number): string =>
  Math.abs(v) >= 1e6 || (Math.abs(v) < 1e-3 && v !== 0)
    ? v.toExponential(4)
    : v.toLocaleString(undefined, { maximumFractionDigits: 4 });

export function CalcPanel() {
  const calcs = useStore((s) => s.calcs);
  const setCalcs = useStore((s) => s.setCalcs);
  const solveRes = useStore((s) => s.solveRes);
  const stale = useStore((s) => s.resultsStale);

  const values = useMemo(
    () => evaluateCalcs(calcs, solveRes), [calcs, solveRes]);

  const update = (i: number, patch: Partial<{ name: string; expr: string }>) =>
    setCalcs(calcs.map((r, j) => (j === i ? { ...r, ...patch } : r)));

  return (
    <div>
      <StaleNotice stale={stale} />
      <PanelTitle>Calculations</PanelTitle>
      {calcs.length === 0 && (
        <Hint>
          Named formulas over the solved flowsheet, like the HYSYS Spreadsheet.
          Accessors: <code>T("S1")</code>, <code>P("S1")</code>, <code>n("S1")</code>{" "}
          (mol/s), <code>VF("S1")</code>, <code>z("S1","benzene")</code>,{" "}
          <code>duty("RXN_duty")</code>, plus <code>Math.*</code> and earlier
          rows by name. Example: <code>n("FEED") * z("FEED","benzene")</code>.
        </Hint>
      )}
      {calcs.map((row, i) => {
        const v = values[i];
        return (
          <div key={i} className="my-1.5 rounded-md border border-line bg-panel2/50 p-2">
            <div className="flex items-center gap-1.5">
              <input
                className="w-[110px] rounded-md border border-line bg-panel2 px-1.5 py-1 font-semibold text-text"
                value={row.name}
                placeholder="name"
                aria-label={`Calculation ${i + 1} name`}
                onChange={(e) => update(i, { name: e.target.value })}
              />
              <span className="text-muted">=</span>
              <span className={`ml-auto text-right font-semibold ${
                v?.error ? "text-bad" : "text-accent"}`}>
                {v?.error ? "—" : v?.value != null ? fmtVal(v.value) : "—"}
              </span>
              <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
                onClick={() => setCalcs(calcs.filter((_, j) => j !== i))}
                aria-label={`Remove ${row.name || "row"}`}>
                <Trash2 size={12} />
              </button>
            </div>
            <input
              className={`mt-1 w-full rounded-md border bg-panel2 px-1.5 py-1 font-mono text-[12px] text-text ${
                v?.error ? "border-bad" : "border-line"}`}
              value={row.expr}
              placeholder='e.g. duty("RXN_duty") / n("PRODUCT")'
              aria-label={`Calculation ${i + 1} expression`}
              onChange={(e) => update(i, { expr: e.target.value })}
            />
            {v?.error && (
              <div className="mt-0.5 text-[11px] text-bad">{v.error}</div>
            )}
          </div>
        );
      })}
      <button
        className="mt-1 cursor-pointer rounded-md border border-line bg-panel2 px-2 py-1 text-[12px] text-muted hover:border-accent hover:text-text"
        onClick={() => setCalcs([...calcs, { name: `calc${calcs.length + 1}`, expr: "" }])}
      >
        <Plus size={11} className="mr-1 inline" aria-hidden />formula
      </button>
    </div>
  );
}
