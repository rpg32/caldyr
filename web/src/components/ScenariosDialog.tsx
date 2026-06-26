// Cases / scenarios: save named snapshots of this flowsheet's parameters + cost
// assumptions ("base", "cheaper-N2", "high-capacity"…), apply one onto the
// canvas, or compare them side by side (re-costs each without touching the
// canvas). Persisted per-flowsheet in meta.ui.
import { Check, Plus, Trash2, X } from "lucide-react";
import { useEffect, useState } from "react";
import { useStore } from "../store";
import type { ScenarioResult } from "../types";
import { Button, Hint, PanelTitle } from "./ui";

const money = (x?: number) => (x == null ? "—" : "$" + Math.round(x).toLocaleString());
const lcop = (x?: number) => (x == null ? "—" : "$" + x.toFixed(3));

export function ScenariosDialog() {
  const open = useStore((s) => s.scenariosOpen);
  const toggle = useStore((s) => s.toggleScenarios);
  const scenarios = useStore((s) => s.scenarios);
  const saveScenario = useStore((s) => s.saveScenario);
  const applyScenario = useStore((s) => s.applyScenario);
  const removeScenario = useStore((s) => s.removeScenario);
  const compareScenarios = useStore((s) => s.compareScenarios);
  const product = useStore((s) => s.product);
  const [name, setName] = useState("");
  const [results, setResults] = useState<ScenarioResult[] | null>(null);
  const [busy, setBusy] = useState(false);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, toggle]);

  if (!open) return null;

  const save = () => { saveScenario(name); setName(""); setResults(null); };
  const compare = async () => {
    setBusy(true);
    try { setResults(await compareScenarios()); } finally { setBusy(false); }
  };

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40" onClick={toggle}>
      <div role="dialog" aria-modal="true" aria-label="Cases and scenarios"
        className="max-h-[85vh] w-[560px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-1 flex items-center">
          <b className="text-[15px]">Cases &amp; scenarios</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close cases"><X size={16} /></button>
        </div>
        <p className="mb-2 text-[11px] text-muted">
          A case snapshots this flowsheet's unit/feed parameters and cost assumptions.
          Save the current setup, switch with Apply, or Compare them side by side.
        </p>

        <PanelTitle>Save current as a case</PanelTitle>
        <div className="flex items-center gap-1.5">
          <input
            className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-2 py-1 text-text"
            placeholder="case name (e.g. cheaper-N2)" value={name}
            aria-label="Case name"
            onChange={(e) => setName(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && save()} />
          <Button icon={<Plus size={13} />} disabled={!name.trim()} onClick={save}>Save</Button>
        </div>

        {scenarios.length > 0 ? (
          <>
            <PanelTitle>Saved cases</PanelTitle>
            {scenarios.map((sc) => (
              <div key={sc.name}
                className="my-1 flex items-center gap-2 rounded-md border border-line bg-panel2/60 px-2 py-1 text-[12px]">
                <span className="min-w-0 flex-1 truncate">{sc.name}</span>
                <button className="flex items-center gap-1 text-muted hover:text-accent"
                  onClick={() => applyScenario(sc.name)} title="Load this case onto the canvas">
                  <Check size={12} /> apply
                </button>
                <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
                  onClick={() => { removeScenario(sc.name); setResults(null); }}
                  aria-label={`Delete ${sc.name}`}><Trash2 size={12} /></button>
              </div>
            ))}

            <div className="mt-3 flex items-center gap-2">
              <Button variant="primary" busy={busy} disabled={!product}
                onClick={() => void compare()}>Compare costs</Button>
              {!product && <span className="text-[11px] text-warn">pick a product (Flowsheet panel) to compare</span>}
            </div>

            {results && (
              <div className="mt-2 overflow-x-auto">
                <table className="data-table">
                  <thead><tr><th>case</th><th>LCOP</th><th>TCI</th><th>OPEX/yr</th><th>NPV</th></tr></thead>
                  <tbody>
                    {results.map((r) => (
                      <tr key={r.name}>
                        <td>{r.name}</td>
                        {r.ok ? (
                          <>
                            <td>{lcop(r.lcop)}/kg</td>
                            <td>{money(r.tci)}</td>
                            <td>{money(r.opex)}</td>
                            <td>{money(r.npv)}</td>
                          </>
                        ) : (
                          <td colSpan={4} className="text-bad" title={r.error}>failed: {r.error}</td>
                        )}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </>
        ) : (
          <Hint>No cases yet — set up the flowsheet (and cost assumptions), then save it as a case.</Hint>
        )}
      </div>
    </div>
  );
}
