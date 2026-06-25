import { ChartSpline, Plus, X } from "lucide-react";
import { useState } from "react";
import {
  CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip as ChartTooltip,
  XAxis, YAxis,
} from "recharts";
import { api } from "../api";
import { compositionRows, fmtFrac, streamMassFlow } from "../lib/composition";
import { compositionSum, metaFor, paramApplies, validateParam } from "../lib/params";
import { defaultUnit, fmtDim, unitsForDim, type Dim, type UnitSet } from "../lib/units";
import { useStore, type Tab } from "../store";
import type { EnvelopeResponse, ParamSchema, StreamState } from "../types";
import { AnalysisPanel } from "./AnalysisPanel";
import { CalcPanel } from "./CalcPanel";
import { DesignPanel } from "./DesignPanel";
import { EconomicsPanel } from "./EconomicsPanel";
import { LogicalEditor } from "./LogicalEditor";
import { OptimizePanel } from "./OptimizePanel";
import { StreamTable } from "./StreamTable";
import { StudyPanel } from "./StudyPanel";
import { Button, Hint, NumberInput, PanelTitle, QuantityInput } from "./ui";

const TABS: { id: Tab; label: string }[] = [
  { id: "params", label: "Params" },
  { id: "streams", label: "Streams" },
  { id: "economics", label: "Econ" },
  { id: "optimize", label: "Opt" },
  { id: "study", label: "Study" },
  { id: "calc", label: "Calc" },
  { id: "tools", label: "Tools" },
];

export function Inspector() {
  const tab = useStore((s) => s.tab);
  const setTab = useStore((s) => s.setTab);
  return (
    <aside className="flex min-h-0 min-w-0 flex-col border-l border-line bg-panel">
      <div className="flex border-b border-line" role="tablist">
        {TABS.map((t) => (
          <button
            key={t.id}
            role="tab"
            aria-selected={tab === t.id}
            className={`flex-1 cursor-pointer border-0 bg-transparent p-2 transition-colors ${
              tab === t.id ? "text-accent shadow-[inset_0_-2px_0_var(--accent)]" : "text-muted hover:text-text"
            }`}
            onClick={() => setTab(t.id)}
          >
            {t.label}
          </button>
        ))}
      </div>
      <div className="min-h-0 flex-1 overflow-auto p-2.5">
        {tab === "params" && <ParamsTab />}
        {tab === "streams" && <StreamTable />}
        {tab === "economics" && <EconomicsPanel />}
        {tab === "optimize" && <OptimizePanel />}
        {tab === "study" && <StudyPanel />}
        {tab === "calc" && <CalcPanel />}
        {tab === "tools" && <AnalysisPanel />}
      </div>
    </aside>
  );
}

function ParamsTab() {
  const selected = useStore((s) => s.selected);
  const nodes = useStore((s) => s.nodes);
  if (selected?.kind === "node") {
    const node = nodes.find((n) => n.id === selected.id);
    if (node) return <NodePanel nodeId={node.id} />;
  }
  if (selected?.kind === "edge") return <StreamPanel edgeId={selected.id} />;
  return <FlowsheetPanel />;
}

/** Numeric input with unit (display-unit-aware) and bounds validation.
 *  Dimensioned params (T/P/flow/duty/UA) edit in the active unit system; values
 *  are always stored/validated in SI. */
function NumField({
  paramKey, value, onChange, unitKey,
}: { paramKey: string; value: number; onChange: (v: number) => void; unitKey?: string }) {
  const meta = metaFor(paramKey);
  const unitSet = useStore((s) => s.unitSet);
  const overrides = useStore((s) => s.unitOverrides);
  const setUnitOverride = useStore((s) => s.setUnitOverride);
  const error = validateParam(paramKey, value); // value + bounds are SI
  const inputClass = `w-[110px] rounded-md border bg-panel2 px-2 py-1 text-right text-text ${
    error ? "border-bad" : "border-line"
  }`;
  return (
    <label className="my-1.5 flex items-center justify-between gap-2">
      <span className="text-muted" title={meta.label}>{meta.label}</span>
      {meta.dim ? (
        <QuantityInput
          dim={meta.dim} set={unitSet} className={inputClass} value={value}
          unit={unitKey ? overrides[unitKey] : undefined}
          units={unitKey ? unitsForDim(meta.dim) : undefined}
          onUnitChange={unitKey ? (u) => setUnitOverride(unitKey, u) : undefined}
          title={error ?? meta.label} aria-label={meta.label} aria-invalid={!!error}
          onChange={onChange}
        />
      ) : (
        <span className="flex items-center gap-1.5">
          <NumberInput
            className={inputClass} value={value} min={meta.min} max={meta.max}
            title={error ?? `${meta.label}${meta.unit !== "–" && meta.unit ? ` in ${meta.unit}` : ""}`}
            aria-label={meta.label} aria-invalid={!!error} onChange={onChange}
          />
          <span className="w-10 text-[11px] text-muted">{meta.unit}</span>
        </span>
      )}
    </label>
  );
}

/** Checkbox for a boolean param (e.g. decant_condenser, reboiled). */
function BoolField({
  paramKey, value, onChange,
}: { paramKey: string; value: boolean; onChange: (v: boolean) => void }) {
  const meta = metaFor(paramKey);
  return (
    <label className="my-1.5 flex items-center justify-between gap-2">
      <span className="text-muted" title={meta.hint ?? meta.label}>{meta.label}</span>
      <input
        type="checkbox"
        className="h-4 w-4 accent-accent"
        checked={!!value}
        aria-label={meta.label}
        onChange={(e) => onChange(e.target.checked)}
      />
    </label>
  );
}

/** Dropdown for an enumerated param (e.g. method, reflux_layer). */
function SelectField({
  paramKey, value, onChange,
}: { paramKey: string; value: string; onChange: (v: string) => void }) {
  const meta = metaFor(paramKey);
  const opts = meta.options ?? [];
  return (
    <label className="my-1.5 flex items-center justify-between gap-2">
      <span className="text-muted" title={meta.hint ?? meta.label}>{meta.label}</span>
      <select
        className="w-[120px] rounded-md border border-line bg-panel2 px-2 py-1 text-text"
        value={String(value ?? opts[0] ?? "")}
        aria-label={meta.label}
        onChange={(e) => onChange(e.target.value)}
      >
        {opts.map((o) => <option key={o} value={o}>{o}</option>)}
      </select>
    </label>
  );
}

/** Render a single unit param with the right widget for its type. */
function ParamField({
  paramKey, value, onChange, unitKey,
}: { paramKey: string; value: unknown; onChange: (v: unknown) => void; unitKey?: string }) {
  const meta = metaFor(paramKey);
  if (meta.type === "boolean" || typeof value === "boolean")
    return <BoolField paramKey={paramKey} value={!!value} onChange={onChange} />;
  if (meta.type === "select")
    return <SelectField paramKey={paramKey} value={String(value ?? "")} onChange={onChange} />;
  if (typeof value === "number")
    return <NumField paramKey={paramKey} value={value} onChange={onChange} unitKey={unitKey} />;
  // unknown non-scalar (feeds, reactions, z, ...): read-only JSON.
  return (
    <div className="my-1.5 flex items-center justify-between gap-2 text-muted">
      <span>{meta.label}</span>
      <code className="max-w-[200px] overflow-hidden text-ellipsis whitespace-nowrap text-[11px] text-text">
        {JSON.stringify(value)}
      </code>
    </div>
  );
}

// -- guided unit-parameter editor (schema-driven, from /unit-types) -----------

/** Effective value of a schema param: explicit value if set, else its default. */
function effective(name: string, params: Record<string, unknown>, schema: ParamSchema[]): unknown {
  if (name in params) return params[name];
  return schema.find((s) => s.name === name)?.default;
}

/** Whether a schema param applies given the current params (its `requires`). */
function schemaApplies(e: ParamSchema, params: Record<string, unknown>, schema: ParamSchema[]): boolean {
  if (!e.requires) return true;
  return Object.entries(e.requires).every(([k, v]) => effective(k, params, schema) === v);
}

const ROW_INPUT = "w-[110px] rounded-md border bg-panel2 px-2 py-1 text-right text-text";
const ROW_SEL = "rounded-md border border-line bg-panel2 px-2 py-1 text-text";

/** One schema-described parameter: the right widget for its type, showing the
 *  engine default until the user sets it, with a reset-to-default control. */
function SchemaParamRow({ schema, value, set, unitOverride, onUnitChange, onSet, onReset }: {
  schema: ParamSchema;
  value: unknown;                       // explicit value, or undefined (= default)
  set: UnitSet;
  unitOverride?: string;
  onUnitChange: (u: string | null) => void;
  onSet: (v: unknown) => void;
  onReset: () => void;
}) {
  const isSet = value !== undefined;
  const eff = isSet ? value : schema.default;
  const missing = !!schema.required && !isSet && schema.default === undefined;
  const numEff = typeof eff === "number" ? eff : NaN;
  const outOfBounds = typeof eff === "number"
    && ((schema.min !== undefined && eff < schema.min)
      || (schema.max !== undefined && eff > schema.max));
  const border = missing || outOfBounds ? "border-bad" : "border-line";
  const dim = schema.dim as Dim | undefined;

  let widget;
  if (schema.complex || schema.type === "json") {
    widget = (
      <code className="truncate text-right text-[11px] text-text"
        title={isSet ? JSON.stringify(value) : undefined}>
        {isSet ? JSON.stringify(value) : "— set via .flow / Copilot"}
      </code>
    );
  } else if (schema.type === "boolean") {
    widget = <input type="checkbox" className="h-4 w-4 accent-accent" checked={!!eff}
      aria-label={schema.label} onChange={(e) => onSet(e.target.checked)} />;
  } else if (schema.type === "select") {
    const opts = schema.options ?? [];
    widget = (
      <select className={`${ROW_SEL} w-full`} value={String(eff ?? opts[0] ?? "")} aria-label={schema.label}
        onChange={(e) => onSet(e.target.value)}>
        {opts.map((o) => <option key={o} value={o}>{o}</option>)}
      </select>
    );
  } else if (schema.type === "string") {
    widget = <input type="text" className={`${ROW_SEL} w-full text-left`} value={String(eff ?? "")}
      aria-label={schema.label} onChange={(e) => onSet(e.target.value)} />;
  } else if (dim) {
    widget = (
      <QuantityInput dim={dim} set={set} className={`${ROW_INPUT} ${border}`} value={numEff}
        unit={unitOverride} units={unitsForDim(dim)} onUnitChange={onUnitChange}
        aria-label={schema.label} aria-invalid={missing || outOfBounds} onChange={onSet} />
    );
  } else {
    widget = (
      <>
        <NumberInput className={`${ROW_INPUT} ${border}`} value={numEff} min={schema.min} max={schema.max}
          aria-label={schema.label} aria-invalid={missing || outOfBounds} onChange={onSet} />
        <span className="w-[58px] shrink-0 text-[11px] text-muted">{schema.unit ?? ""}</span>
      </>
    );
  }

  return (
    <label className="my-1.5 flex items-center gap-2">
      <span className={`min-w-0 flex-1 truncate ${isSet ? "text-text" : "text-muted"}`}
        title={schema.help ?? schema.label}>
        {schema.label}
        {schema.required && <span className="text-bad" title="required">&nbsp;*</span>}
      </span>
      <span className="flex w-[174px] shrink-0 items-center justify-end gap-1.5">
        {widget}
      </span>
      {isSet && !schema.required ? (
        <button className="shrink-0 cursor-pointer p-0.5 text-muted hover:text-accent"
          title="Reset to default" aria-label={`Reset ${schema.label}`} onClick={onReset}>
          <X size={11} />
        </button>
      ) : <span className="w-[15px] shrink-0" aria-hidden />}
    </label>
  );
}

/** Schema-driven parameter form for a unit op: shows every available parameter
 *  (with its default, units and help), validates, and lets you reset overrides
 *  and remove unknown extras. Falls back to the free-text editor if a unit has
 *  no schema. */
function UnitParamsEditor({ nodeId, unitType, params }: {
  nodeId: string; unitType: string; params: Record<string, unknown>;
}) {
  const schema = useStore((s) => s.unitTypes.find((t) => t.type === unitType)?.params ?? []);
  const setParam = useStore((s) => s.setParam);
  const unsetParam = useStore((s) => s.unsetParam);
  const unitSet = useStore((s) => s.unitSet);
  const overrides = useStore((s) => s.unitOverrides);
  const setUnitOverride = useStore((s) => s.setUnitOverride);

  if (schema.length === 0) return <LegacyParams nodeId={nodeId} params={params} />;

  const schemaNames = new Set(schema.map((s) => s.name));
  const visible = schema.filter((e) => schemaApplies(e, params, schema));
  const extras = Object.keys(params).filter((k) => !schemaNames.has(k));

  return (
    <>
      {visible.map((e) => (
        <SchemaParamRow key={e.name} schema={e}
          value={e.name in params ? params[e.name] : undefined}
          set={unitSet} unitOverride={overrides[`${nodeId}:${e.name}`]}
          onUnitChange={(u) => setUnitOverride(`${nodeId}:${e.name}`, u)}
          onSet={(v) => setParam(nodeId, e.name, v)}
          onReset={() => unsetParam(nodeId, e.name)} />
      ))}
      {extras.length > 0 && (
        <>
          <PanelTitle>Extra parameters</PanelTitle>
          {extras.map((k) => (
            <div key={k}
              className="my-1 flex items-center gap-2 rounded-md border border-warn/40 bg-warn/5 px-2 py-1">
              <span className="text-[12px] text-warn" title="Not a known parameter for this unit type">{k}</span>
              <code className="ml-auto max-w-[130px] truncate text-[11px] text-text">{JSON.stringify(params[k])}</code>
              <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
                onClick={() => unsetParam(nodeId, k)} aria-label={`Remove ${k}`}>
                <X size={12} />
              </button>
            </div>
          ))}
        </>
      )}
    </>
  );
}

/** Legacy free-text editor — only for unit types without a published schema. */
function LegacyParams({ nodeId, params }: { nodeId: string; params: Record<string, unknown> }) {
  const setParam = useStore((s) => s.setParam);
  const unsetParam = useStore((s) => s.unsetParam);
  return (
    <>
      {Object.entries(params).filter(([k]) => paramApplies(k, params)).map(([k, v]) => (
        <div key={k} className="flex items-center gap-1">
          <div className="flex-1"><ParamField paramKey={k} value={v} unitKey={`${nodeId}:${k}`}
            onChange={(nv) => setParam(nodeId, k, nv)} /></div>
          <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
            onClick={() => unsetParam(nodeId, k)} aria-label={`Remove ${k}`}><X size={12} /></button>
        </div>
      ))}
      {Object.keys(params).length === 0 && (
        <Hint>No parameters set — this unit will use engine defaults.</Hint>
      )}
      <AddParam nodeId={nodeId} existing={Object.keys(params)} />
    </>
  );
}

function NodePanel({ nodeId }: { nodeId: string }) {
  const node = useStore((s) => s.nodes.find((n) => n.id === nodeId))!;
  const components = useStore((s) => s.components);
  const setParam = useStore((s) => s.setParam);
  const inletEdge = useStore((s) => s.edges.find((e) => e.target === nodeId));
  const inletState = useStore((s) =>
    inletEdge ? s.solveRes?.streams[inletEdge.id] : undefined);
  const p = node.data.params;

  return (
    <div>
      <PanelTitle>
        {node.data.label}{" "}
        <em className="normal-case">
          ({node.data.kind === "unit" ? node.data.unitType : node.data.kind})
        </em>
      </PanelTitle>

      {node.data.kind === "feed" && (
        <>
          {(["T", "P", "molar_flow"] as const).map((k) => (
            <NumField key={k} paramKey={k} value={Number(p[k] ?? 0)}
              unitKey={`${node.id}:${k}`}
              onChange={(v) => setParam(node.id, k, v)} />
          ))}
          <CompositionEditor nodeId={node.id} z={(p.z as Record<string, number>) ?? {}}
            components={components} />
        </>
      )}

      {node.data.kind === "unit" && (
        <>
          <UnitParamsEditor nodeId={node.id} unitType={node.data.unitType as string} params={p} />
          <DesignPanel unitId={node.id} />
        </>
      )}

      {node.data.kind === "product" && (
        inletState ? (
          <>
            <div className="mb-1 text-[11px] text-muted">
              feed stream: <span className="text-text">{inletEdge!.id}</span>
            </div>
            <StreamReadout state={inletState} />
          </>
        ) : (
          <Hint>
            Product sink. {inletEdge
              ? "Press Solve to see the stream feeding it."
              : "Connect a stream to this product, then Solve."}
          </Hint>
        )
      )}
    </div>
  );
}

function CompositionEditor({
  nodeId, z, components,
}: { nodeId: string; z: Record<string, number>; components: string[] }) {
  const setParam = useStore((s) => s.setParam);
  const sum = compositionSum(z);
  const ok = Math.abs(sum - 1) < 1e-6;
  if (!components.length) {
    return <Hint>Add components in the Flowsheet panel (deselect everything) first.</Hint>;
  }
  return (
    <>
      <div className="mt-2.5 flex items-center justify-between text-[11px]">
        <span className="text-muted">composition (mole frac)</span>
        <span className={ok ? "text-ok" : "text-warn"} title="Mole fractions should sum to 1">
          Σ = {sum.toFixed(4)}
        </span>
      </div>
      {components.map((c) => (
        <label key={c} className="my-1.5 flex items-center justify-between gap-2">
          <span className="text-muted">{c}</span>
          <NumberInput
            min={0} max={1}
            className={`w-[110px] rounded-md border bg-panel2 px-2 py-1 text-right text-text ${
              ok ? "border-line" : "border-warn"
            }`}
            value={z[c] ?? 0}
            aria-label={`Mole fraction of ${c}`}
            onChange={(v) => setParam(nodeId, "z", { ...z, [c]: v })}
          />
        </label>
      ))}
    </>
  );
}

function AddParam({ nodeId, existing }: { nodeId: string; existing: string[] }) {
  const setParam = useStore((s) => s.setParam);
  const [key, setKey] = useState("");
  const add = () => {
    const k = key.trim();
    if (!k || existing.includes(k)) return;
    setParam(nodeId, k, 0);
    setKey("");
  };
  return (
    <div className="mt-2 flex items-center gap-1.5">
      <input
        className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-2 py-1 text-text"
        placeholder="add parameter (e.g. dP)"
        value={key}
        aria-label="New parameter name"
        onChange={(e) => setKey(e.target.value)}
        onKeyDown={(e) => e.key === "Enter" && add()}
      />
      <button className="cursor-pointer rounded-md border border-line bg-panel2 p-1.5 text-muted hover:border-accent"
        onClick={add} aria-label="Add parameter" title="Add parameter">
        <Plus size={13} />
      </button>
    </div>
  );
}

function EnvelopeChart({ env }: { env: EnvelopeResponse }) {
  const bubble = env.points.map((p) => ({ T: p.T_bubble, P: p.P / 1000 }));
  const dew = env.points.map((p) => ({ T: p.T_dew, P: p.P / 1000 }));
  return (
    <div style={{ height: 170 }}>
      <ResponsiveContainer>
        <LineChart margin={{ top: 6, right: 10, left: 0, bottom: 14 }}>
          <CartesianGrid stroke="var(--line)" strokeDasharray="3 3" />
          <XAxis dataKey="T" type="number" domain={["auto", "auto"]}
            tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
            label={{ value: "T / K", fill: "var(--muted)", fontSize: 10, position: "insideBottom", dy: 12 }}
            tickFormatter={(v: number) => v.toFixed(0)} allowDuplicatedCategory={false} />
          <YAxis dataKey="P" type="number" domain={["auto", "auto"]} width={52}
            tick={{ fill: "var(--muted)", fontSize: 10 }} stroke="var(--line)"
            label={{ value: "P / kPa", fill: "var(--muted)", fontSize: 10, angle: -90, position: "insideLeft" }}
            tickFormatter={(v: number) => v.toFixed(0)} />
          <ChartTooltip
            contentStyle={{
              background: "var(--panel)", border: "1px solid var(--line)",
              borderRadius: 6, fontSize: 11,
            }}
            formatter={(v) => [`${Number(v).toFixed(1)} kPa`, "P"]}
            labelFormatter={(v) => `T = ${Number(v).toFixed(1)} K`}
          />
          <Line data={bubble} name="bubble" dataKey="P" stroke="#3b82f6"
            strokeWidth={1.6} dot={false} isAnimationActive={false} />
          <Line data={dew} name="dew" dataKey="P" stroke="#f97316"
            strokeWidth={1.6} dot={false} isAnimationActive={false} />
        </LineChart>
      </ResponsiveContainer>
      <div className="flex gap-3 text-[10px] text-muted">
        <span><i className="inline-block h-2 w-2 rounded-sm" style={{ background: "#3b82f6" }} /> bubble</span>
        <span><i className="inline-block h-2 w-2 rounded-sm" style={{ background: "#f97316" }} /> dew</span>
      </div>
    </div>
  );
}

/** Solved state (T/P/flow/phase/VF) + composition for one stream, in the active
 *  unit system. Mass flow / mass fractions use the /solve molar-mass map. */
function StreamReadout({ state }: { state: StreamState }) {
  const unitSet = useStore((s) => s.unitSet);
  const mw = useStore((s) => s.solveRes?.molar_mass);
  const massFlow = streamMassFlow(state.z, state.molar_flow, mw);
  return (
    <>
      <table className="data-table">
        <tbody>
          <tr><td>T</td><td>{fmtDim("temperature", state.T, unitSet, 2)} {defaultUnit("temperature", unitSet)}</td></tr>
          <tr><td>P</td><td>{fmtDim("pressure", state.P, unitSet, 3)} {defaultUnit("pressure", unitSet)}</td></tr>
          <tr><td>molar flow</td><td>{fmtDim("molar_flow", state.molar_flow, unitSet, 4)} {defaultUnit("molar_flow", unitSet)}</td></tr>
          <tr><td>mass flow</td><td>{fmtDim("mass_flow", massFlow, unitSet, 4)} {defaultUnit("mass_flow", unitSet)}</td></tr>
          <tr><td>phase</td><td>{state.phase ?? "—"}</td></tr>
          <tr><td>vapor frac</td><td>{state.vapor_fraction?.toFixed(3) ?? "—"}</td></tr>
        </tbody>
      </table>
      <CompositionTable state={state} />
    </>
  );
}

function CompositionTable({ state }: { state: StreamState }) {
  const unitSet = useStore((s) => s.unitSet);
  const mw = useStore((s) => s.solveRes?.molar_mass);
  const rows = compositionRows(state.z, state.molar_flow, mw);
  if (!rows.length) return null;
  const flowUnit = defaultUnit("molar_flow", unitSet);
  const haveMass = rows.some((r) => r.massFrac != null);
  return (
    <div className="mt-2">
      <PanelTitle>Composition</PanelTitle>
      <table className="data-table">
        <thead>
          <tr>
            <th>component</th>
            <th title="mole fraction">mole frac</th>
            {haveMass && <th title="mass fraction">mass frac</th>}
            <th title={`molar flow (${flowUnit})`}>{flowUnit}</th>
          </tr>
        </thead>
        <tbody>
          {rows.map((r) => (
            <tr key={r.comp}>
              <td>{r.comp}</td>
              <td>{fmtFrac(r.frac)}</td>
              {haveMass && <td>{r.massFrac != null ? fmtFrac(r.massFrac) : "—"}</td>}
              <td>{fmtDim("molar_flow", r.flow, unitSet, 4)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function StreamPanel({ edgeId }: { edgeId: string }) {
  const renameEdge = useStore((s) => s.renameEdge);
  const solveRes = useStore((s) => s.solveRes);
  const toFlowDoc = useStore((s) => s.toFlowDoc);
  const toast = useStore((s) => s.toast);
  const [draft, setDraft] = useState(edgeId);
  const [env, setEnv] = useState<EnvelopeResponse | null>(null);
  const [envBusy, setEnvBusy] = useState(false);
  const state = solveRes?.streams[edgeId];

  const fetchEnvelope = async () => {
    if (!solveRes) return;
    setEnvBusy(true);
    try {
      // Attach the solved cache so the engine can reconstruct internal streams.
      const doc = { ...toFlowDoc(), solved: solveRes.streams };
      setEnv(await api.envelope(doc, edgeId, 30));
    } catch (e) {
      toast("error", `Envelope failed: ${(e as Error).message}`);
    } finally {
      setEnvBusy(false);
    }
  };

  return (
    <div>
      <PanelTitle>Stream</PanelTitle>
      <label className="my-1.5 flex items-center justify-between gap-2">
        <span className="text-muted">name</span>
        <input
          className="w-[150px] rounded-md border border-line bg-panel2 px-2 py-1 text-text"
          value={draft}
          aria-label="Stream name"
          onChange={(e) => setDraft(e.target.value)}
          onBlur={() => draft !== edgeId && renameEdge(edgeId, draft)}
          onKeyDown={(e) => e.key === "Enter" && (e.target as HTMLInputElement).blur()}
        />
      </label>
      {state ? (
        <>
          <StreamReadout state={state} />
          {env && env.stream === edgeId ? (
            <EnvelopeChart env={env} />
          ) : (
            <Button icon={<ChartSpline size={13} />} busy={envBusy}
              onClick={() => void fetchEnvelope()}>
              Phase envelope
            </Button>
          )}
        </>
      ) : (
        <Hint>Solve the flowsheet to see this stream's state.</Hint>
      )}
    </div>
  );
}

function GroupsList() {
  const groups = useStore((s) => s.groups);
  const toggle = useStore((s) => s.toggleGroupCollapse);
  const ungroup = useStore((s) => s.ungroup);
  if (!groups.length) return null;
  return (
    <div>
      <PanelTitle>Groups</PanelTitle>
      {groups.map((g) => (
        <div key={g.id}
          className="my-1 flex items-center gap-2 rounded-md border border-line bg-panel2/60 px-2 py-1 text-[12px]">
          <span className="min-w-0 flex-1 truncate">{g.label}</span>
          <span className="text-[11px] text-muted">{g.members.length} units</span>
          <button className="cursor-pointer text-muted hover:text-accent"
            onClick={() => toggle(g.id)}>
            {g.collapsed ? "expand" : "collapse"}
          </button>
          <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
            onClick={() => ungroup(g.id)} aria-label={`Ungroup ${g.label}`}>
            <X size={11} />
          </button>
        </div>
      ))}
    </div>
  );
}

function FlowsheetPanel() {
  const components = useStore((s) => s.components);
  const catalog = useStore((s) => s.componentCatalog);
  const addComponent = useStore((s) => s.addComponent);
  const removeComponent = useStore((s) => s.removeComponent);
  const packages = useStore((s) => s.packages);
  const propertyPackage = useStore((s) => s.propertyPackage);
  const setPropertyPackage = useStore((s) => s.setPropertyPackage);
  const product = useStore((s) => s.product);
  const setProduct = useStore((s) => s.setProduct);
  const [draft, setDraft] = useState("");

  const add = () => {
    if (draft.trim()) {
      addComponent(draft);
      setDraft("");
    }
  };

  return (
    <div>
      <PanelTitle>Flowsheet</PanelTitle>

      <label className="my-1.5 flex items-center justify-between gap-2">
        <span className="text-muted">property package</span>
        <select
          className="w-[150px] rounded-md border border-line bg-panel2 px-2 py-1 text-text"
          value={propertyPackage}
          onChange={(e) => setPropertyPackage(e.target.value)}
        >
          {packages.map((p) => (
            <option key={p.id} value={p.id} title={p.use}>{p.id}</option>
          ))}
        </select>
      </label>

      <label className="my-1.5 flex items-center justify-between gap-2">
        <span className="text-muted" title="Component whose product streams set the LCOP basis">
          product (costing)
        </span>
        <select
          className="w-[150px] rounded-md border border-line bg-panel2 px-2 py-1 text-text"
          value={product}
          onChange={(e) => setProduct(e.target.value)}
        >
          <option value="">— pick —</option>
          {components.map((c) => <option key={c} value={c}>{c}</option>)}
        </select>
      </label>

      <PanelTitle>Components</PanelTitle>
      <div className="flex flex-wrap gap-1.5">
        {components.map((c) => (
          <span key={c}
            className="inline-flex items-center gap-1 rounded-full border border-line bg-panel2 px-2 py-0.5 text-[12px]">
            {c}
            <button className="cursor-pointer text-muted hover:text-bad"
              onClick={() => removeComponent(c)} aria-label={`Remove ${c}`} title={`Remove ${c}`}>
              <X size={11} />
            </button>
          </span>
        ))}
        {!components.length && <span className="text-[12px] text-muted">none yet</span>}
      </div>
      <div className="mt-2 flex items-center gap-1.5">
        <input
          className="min-w-0 flex-1 rounded-md border border-line bg-panel2 px-2 py-1 text-text"
          placeholder="add component (e.g. ammonia)"
          value={draft}
          aria-label="New component name"
          list="component-catalog"
          onChange={(e) => setDraft(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && add()}
        />
        <datalist id="component-catalog">
          {catalog.map((c) => (
            <option key={c.id} value={c.id}>{c.name} ({c.formula})</option>
          ))}
        </datalist>
        <button className="cursor-pointer rounded-md border border-line bg-panel2 p-1.5 text-muted hover:border-accent"
          onClick={add} aria-label="Add component" title="Add component">
          <Plus size={13} />
        </button>
      </div>

      <GroupsList />
      <LogicalEditor />

      <Hint>
        Select a unit or stream on the canvas to edit it here. Double-click a node
        to rename it; double-click a stream to pin its data on the canvas;
        hover a stream for a live readout. Shift+drag for marquee select.
        Keyboard: Ctrl+Z undo, Ctrl+C/V copy/paste, Ctrl+A select all,
        Ctrl+S save, Del delete.
      </Hint>
    </div>
  );
}
