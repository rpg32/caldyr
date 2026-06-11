// Optimization builder over POST /optimize: objective + design variables +
// constraints, run, inspect, and apply the optimum back onto the canvas.
import { Check, Play, Plus, Trash2 } from "lucide-react";
import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import type { ConstraintSpec, DesignVarSpec, MetricSpec, OptimizeResponse } from "../types";
import { MetricEditor, metricValid } from "./MetricEditor";
import { Button, Hint, PanelTitle } from "./ui";

const sel = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-text min-w-0";
const num = "w-[72px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text";

const EMPTY_METRIC: MetricSpec = { type: "flow", stream: "" };

function DesignVarRow({ dv, units, paramsOf, onChange, onRemove }: {
  dv: DesignVarSpec;
  units: string[];
  paramsOf: (unit: string) => string[];
  onChange: (dv: DesignVarSpec) => void;
  onRemove: () => void;
}) {
  return (
    <div className="my-1 flex flex-wrap items-center gap-1.5">
      <select className={sel} value={dv.unit_id} aria-label="Unit"
        onChange={(e) => onChange({ ...dv, unit_id: e.target.value, param: "" })}>
        <option value="">— unit —</option>
        {units.map((u) => <option key={u} value={u}>{u}</option>)}
      </select>
      <select className={sel} value={dv.param} aria-label="Parameter"
        onChange={(e) => onChange({ ...dv, param: e.target.value })}>
        <option value="">— param —</option>
        {paramsOf(dv.unit_id).map((p) => <option key={p} value={p}>{p}</option>)}
      </select>
      <input className={num} type="number" step="any" value={dv.lower} placeholder="min"
        aria-label="Lower bound"
        onChange={(e) => onChange({ ...dv, lower: parseFloat(e.target.value) })} />
      <span className="text-muted">…</span>
      <input className={num} type="number" step="any" value={dv.upper} placeholder="max"
        aria-label="Upper bound"
        onChange={(e) => onChange({ ...dv, upper: parseFloat(e.target.value) })} />
      <button className="cursor-pointer p-1 text-muted hover:text-bad" onClick={onRemove}
        aria-label="Remove design variable"><Trash2 size={12} /></button>
    </div>
  );
}

export function OptimizePanel() {
  const nodes = useStore((s) => s.nodes);
  const toFlowDoc = useStore((s) => s.toFlowDoc);
  const backend = useStore((s) => s.backend);
  const setParam = useStore((s) => s.setParam);
  const toast = useStore((s) => s.toast);

  const [sense, setSense] = useState<"min" | "max">("min");
  const [objective, setObjective] = useState<MetricSpec>(EMPTY_METRIC);
  const [designVars, setDesignVars] = useState<DesignVarSpec[]>([]);
  const [constraints, setConstraints] = useState<ConstraintSpec[]>([]);
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<OptimizeResponse | null>(null);

  const unitIds = nodes.filter((n) => n.data.kind === "unit").map((n) => n.id);
  const paramsOf = (unit: string): string[] => {
    const n = nodes.find((x) => x.id === unit);
    return Object.entries(n?.data.params ?? {})
      .filter(([, v]) => typeof v === "number").map(([k]) => k);
  };

  const valid = metricValid(objective)
    && designVars.length > 0
    && designVars.every((d) => d.unit_id && d.param && d.lower < d.upper)
    && constraints.every((c) => metricValid(c.metric) && Number.isFinite(c.value));

  const run = async () => {
    setRunning(true);
    setResult(null);
    try {
      const res = await api.optimize({
        flow: toFlowDoc(), backend,
        objective: { sense, metric: objective },
        design_vars: designVars,
        constraints,
      });
      setResult(res);
      if (!res.success) toast("error", `Optimizer finished without success: ${res.message}`);
    } catch (e) {
      toast("error", `Optimization failed: ${(e as Error).message}`);
    } finally {
      setRunning(false);
    }
  };

  const apply = () => {
    if (!result) return;
    for (const [key, value] of Object.entries(result.design)) {
      const dot = key.indexOf(".");
      if (dot > 0) setParam(key.slice(0, dot), key.slice(dot + 1), value);
    }
    toast("success", "Optimized parameters applied to the flowsheet — re-solve to refresh results.");
  };

  return (
    <div>
      <PanelTitle>Objective</PanelTitle>
      <div className="flex flex-wrap items-center gap-1.5">
        <select className={sel} value={sense} aria-label="Sense"
          onChange={(e) => setSense(e.target.value as "min" | "max")}>
          <option value="min">minimize</option>
          <option value="max">maximize</option>
        </select>
        <MetricEditor value={objective} onChange={setObjective} />
      </div>

      <PanelTitle>Design variables</PanelTitle>
      {designVars.map((dv, i) => (
        <DesignVarRow key={i} dv={dv} units={unitIds} paramsOf={paramsOf}
          onChange={(next) => setDesignVars(designVars.map((d, j) => (j === i ? next : d)))}
          onRemove={() => setDesignVars(designVars.filter((_, j) => j !== i))} />
      ))}
      <Button variant="ghost" icon={<Plus size={12} />}
        onClick={() => setDesignVars([...designVars, { unit_id: "", param: "", lower: 0, upper: 1 }])}>
        design variable
      </Button>

      <PanelTitle>Constraints</PanelTitle>
      {constraints.map((c, i) => (
        <div key={i} className="my-1 flex flex-wrap items-center gap-1.5">
          <MetricEditor value={c.metric}
            onChange={(m) => setConstraints(constraints.map((x, j) => (j === i ? { ...x, metric: m } : x)))} />
          <select className={sel} value={c.op} aria-label="Operator"
            onChange={(e) => setConstraints(constraints.map((x, j) =>
              (j === i ? { ...x, op: e.target.value as ">=" | "<=" } : x)))}>
            <option value=">=">≥</option>
            <option value="<=">≤</option>
          </select>
          <input className={num} type="number" step="any" value={c.value} aria-label="Constraint value"
            onChange={(e) => setConstraints(constraints.map((x, j) =>
              (j === i ? { ...x, value: parseFloat(e.target.value) } : x)))} />
          <button className="cursor-pointer p-1 text-muted hover:text-bad"
            onClick={() => setConstraints(constraints.filter((_, j) => j !== i))}
            aria-label="Remove constraint"><Trash2 size={12} /></button>
        </div>
      ))}
      <Button variant="ghost" icon={<Plus size={12} />}
        onClick={() => setConstraints([...constraints, { metric: EMPTY_METRIC, op: ">=", value: 0 }])}>
        constraint
      </Button>

      <div className="mt-3">
        <Button variant="primary" icon={<Play size={13} />} busy={running} disabled={!valid}
          onClick={() => void run()}>
          Optimize
        </Button>
      </div>
      {!valid && (
        <Hint>Pick an objective metric and at least one design variable (with min &lt; max).
          Duty metrics need a prior solve so the duty list is known.</Hint>
      )}

      {result && (
        <div className="mt-2">
          <PanelTitle>Result</PanelTitle>
          <table className="data-table">
            <tbody>
              <tr><td>success</td><td>{String(result.success)}</td></tr>
              {/* the API reports the minimization objective; un-negate for "max" */}
              <tr><td>objective</td>
                <td>{(sense === "max" ? -result.objective : result.objective).toPrecision(6)}</td></tr>
              <tr><td>engine solves</td><td>{result.n_solves}</td></tr>
              {Object.entries(result.design).map(([k, v]) => (
                <tr key={k}><td>{k}</td><td>{v.toPrecision(6)}</td></tr>
              ))}
            </tbody>
          </table>
          <Button icon={<Check size={13} />} onClick={apply}>
            Apply to flowsheet
          </Button>
        </div>
      )}
    </div>
  );
}
