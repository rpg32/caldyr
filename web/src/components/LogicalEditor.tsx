// Editor for flowsheet-level logical ops (engine: caldyr.solver.logical):
//   Set    {type:"set", target:[unit,param], source:[unit,param], multiplier, offset}
//   Adjust {type:"adjust", vary:[unit,param], bounds:[lo,hi], spec:{...}, value, tolerance}
import { Link2, Plus, Target, Trash2 } from "lucide-react";
import { useStore } from "../store";
import { PanelTitle } from "./ui";

const sel = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-text min-w-0";
const num = "w-[70px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-text";

type Op = Record<string, unknown>;

const METRIC_TYPES = ["T", "P", "molar_flow", "component_rate", "duty"] as const;

function UnitParamPicker({ value, onChange, label }: {
  value: [string, string];
  onChange: (v: [string, string]) => void;
  label: string;
}) {
  const nodes = useStore((s) => s.nodes);
  const unitIds = nodes.filter((n) => n.data.kind === "unit").map((n) => n.id);
  const params = Object.entries(
    nodes.find((n) => n.id === value[0])?.data.params ?? {})
    .filter(([, v]) => typeof v === "number").map(([k]) => k);
  return (
    <>
      <select className={sel} value={value[0]} aria-label={`${label} unit`}
        onChange={(e) => onChange([e.target.value, ""])}>
        <option value="">— unit —</option>
        {unitIds.map((u) => <option key={u} value={u}>{u}</option>)}
      </select>
      <select className={sel} value={value[1]} aria-label={`${label} param`}
        onChange={(e) => onChange([value[0], e.target.value])}>
        <option value="">— param —</option>
        {params.map((p) => <option key={p} value={p}>{p}</option>)}
      </select>
    </>
  );
}

function SetEditor({ op, onChange }: { op: Op; onChange: (op: Op) => void }) {
  const target = (op.target as [string, string]) ?? ["", ""];
  const source = (op.source as [string, string]) ?? ["", ""];
  return (
    <div className="flex flex-wrap items-center gap-1.5">
      <UnitParamPicker label="target" value={target}
        onChange={(v) => onChange({ ...op, target: v })} />
      <span className="text-muted">=</span>
      <input className={num} type="number" step="any" aria-label="Multiplier"
        value={Number(op.multiplier ?? 1)}
        onChange={(e) => onChange({ ...op, multiplier: parseFloat(e.target.value) })} />
      <span className="text-muted">×</span>
      <UnitParamPicker label="source" value={source}
        onChange={(v) => onChange({ ...op, source: v })} />
      <span className="text-muted">+</span>
      <input className={num} type="number" step="any" aria-label="Offset"
        value={Number(op.offset ?? 0)}
        onChange={(e) => onChange({ ...op, offset: parseFloat(e.target.value) })} />
    </div>
  );
}

function AdjustEditor({ op, onChange }: { op: Op; onChange: (op: Op) => void }) {
  const edges = useStore((s) => s.edges);
  const components = useStore((s) => s.components);
  const solveRes = useStore((s) => s.solveRes);
  const vary = (op.vary as [string, string]) ?? ["", ""];
  const bounds = (op.bounds as [number, number]) ?? [0, 1];
  const spec = (op.spec as Record<string, unknown>) ?? { type: "T", stream: "" };
  const dutyKeys = Object.keys(solveRes?.report.duties ?? {});
  const setSpec = (patch: Record<string, unknown>) =>
    onChange({ ...op, spec: { ...spec, ...patch } });
  return (
    <div className="flex flex-wrap items-center gap-1.5">
      <span className="text-[11px] text-muted">vary</span>
      <UnitParamPicker label="vary" value={vary}
        onChange={(v) => onChange({ ...op, vary: v })} />
      <input className={num} type="number" step="any" aria-label="Lower bound"
        value={bounds[0]}
        onChange={(e) => onChange({ ...op, bounds: [parseFloat(e.target.value), bounds[1]] })} />
      <span className="text-muted">…</span>
      <input className={num} type="number" step="any" aria-label="Upper bound"
        value={bounds[1]}
        onChange={(e) => onChange({ ...op, bounds: [bounds[0], parseFloat(e.target.value)] })} />
      <span className="basis-full" />
      <span className="text-[11px] text-muted">until</span>
      <select className={sel} value={String(spec.type ?? "T")} aria-label="Spec type"
        onChange={(e) => setSpec({ type: e.target.value, stream: "" })}>
        {METRIC_TYPES.map((t) => <option key={t} value={t}>{t}</option>)}
      </select>
      {spec.type === "duty" ? (
        <select className={sel} value={String(spec.stream ?? "")} aria-label="Spec duty"
          onChange={(e) => setSpec({ stream: e.target.value })}>
          <option value="">— duty —</option>
          {dutyKeys.map((k) => <option key={k} value={k}>{k}</option>)}
        </select>
      ) : (
        <select className={sel} value={String(spec.stream ?? "")} aria-label="Spec stream"
          onChange={(e) => setSpec({ stream: e.target.value })}>
          <option value="">— stream —</option>
          {edges.map((e) => <option key={e.id} value={e.id}>{e.id}</option>)}
        </select>
      )}
      {spec.type === "component_rate" && (
        <select className={sel} value={String(spec.component ?? "")} aria-label="Spec component"
          onChange={(e) => setSpec({ component: e.target.value })}>
          <option value="">— component —</option>
          {components.map((c) => <option key={c} value={c}>{c}</option>)}
        </select>
      )}
      <span className="text-muted">=</span>
      <input className={num} type="number" step="any" aria-label="Spec value"
        value={Number(op.value ?? 0)}
        onChange={(e) => onChange({ ...op, value: parseFloat(e.target.value) })} />
    </div>
  );
}

export function LogicalEditor() {
  const logical = useStore((s) => s.logical);
  const setLogical = useStore((s) => s.setLogical);

  const update = (i: number, op: Op) =>
    setLogical(logical.map((x, j) => (j === i ? op : x)));
  const remove = (i: number) => setLogical(logical.filter((_, j) => j !== i));

  return (
    <div>
      <PanelTitle>Logical ops (Set / Adjust)</PanelTitle>
      {logical.map((op, i) => (
        <div key={i} className="my-1.5 rounded-md border border-line bg-panel2/50 p-2">
          <div className="mb-1 flex items-center gap-1.5 text-[11px] uppercase tracking-wide text-muted">
            {op.type === "set" ? <Link2 size={11} aria-hidden /> : <Target size={11} aria-hidden />}
            {String(op.type)}
            <button className="ml-auto cursor-pointer p-0.5 text-muted hover:text-bad"
              onClick={() => remove(i)} aria-label={`Remove ${String(op.type)}`}>
              <Trash2 size={12} />
            </button>
          </div>
          {op.type === "set"
            ? <SetEditor op={op} onChange={(o) => update(i, o)} />
            : <AdjustEditor op={op} onChange={(o) => update(i, o)} />}
        </div>
      ))}
      <div className="flex gap-1.5">
        <button className="cursor-pointer rounded-md border border-line bg-panel2 px-2 py-1 text-[12px] text-muted hover:border-accent hover:text-text"
          onClick={() => setLogical([...logical,
            { type: "set", target: ["", ""], source: ["", ""], multiplier: 1, offset: 0 }])}>
          <Plus size={11} className="mr-1 inline" aria-hidden />Set
        </button>
        <button className="cursor-pointer rounded-md border border-line bg-panel2 px-2 py-1 text-[12px] text-muted hover:border-accent hover:text-text"
          onClick={() => setLogical([...logical,
            { type: "adjust", vary: ["", ""], bounds: [0, 1],
              spec: { type: "T", stream: "" }, value: 0, tolerance: 1e-6 }])}>
          <Plus size={11} className="mr-1 inline" aria-hidden />Adjust
        </button>
      </div>
    </div>
  );
}
