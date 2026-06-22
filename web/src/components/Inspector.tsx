import { ChartSpline, Plus, X } from "lucide-react";
import { useState } from "react";
import {
  CartesianGrid, Line, LineChart, ResponsiveContainer, Tooltip as ChartTooltip,
  XAxis, YAxis,
} from "recharts";
import { api } from "../api";
import { compositionSum, metaFor, validateParam } from "../lib/params";
import { useStore, type Tab } from "../store";
import type { EnvelopeResponse } from "../types";
import { AnalysisPanel } from "./AnalysisPanel";
import { CalcPanel } from "./CalcPanel";
import { DesignPanel } from "./DesignPanel";
import { EconomicsPanel } from "./EconomicsPanel";
import { LogicalEditor } from "./LogicalEditor";
import { OptimizePanel } from "./OptimizePanel";
import { StreamTable } from "./StreamTable";
import { StudyPanel } from "./StudyPanel";
import { Button, Hint, PanelTitle } from "./ui";

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

/** Numeric input with unit suffix and bounds validation. */
function NumField({
  paramKey, value, onChange,
}: { paramKey: string; value: number; onChange: (v: number) => void }) {
  const meta = metaFor(paramKey);
  const error = validateParam(paramKey, value);
  return (
    <label className="my-1.5 flex items-center justify-between gap-2">
      <span className="text-muted" title={meta.label}>{meta.label}</span>
      <span className="flex items-center gap-1.5">
        <input
          type="number"
          step="any"
          className={`w-[110px] rounded-md border bg-panel2 px-2 py-1 text-right text-text ${
            error ? "border-bad" : "border-line"
          }`}
          value={Number.isNaN(value) ? "" : value}
          min={meta.min}
          max={meta.max}
          title={error ?? `${meta.label}${meta.unit !== "–" && meta.unit ? ` in ${meta.unit}` : ""}`}
          aria-label={meta.label}
          aria-invalid={!!error}
          onChange={(e) => onChange(parseFloat(e.target.value))}
        />
        <span className="w-10 text-[11px] text-muted">{meta.unit}</span>
      </span>
    </label>
  );
}

function NodePanel({ nodeId }: { nodeId: string }) {
  const node = useStore((s) => s.nodes.find((n) => n.id === nodeId))!;
  const components = useStore((s) => s.components);
  const setParam = useStore((s) => s.setParam);
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
              onChange={(v) => setParam(node.id, k, v)} />
          ))}
          <CompositionEditor nodeId={node.id} z={(p.z as Record<string, number>) ?? {}}
            components={components} />
        </>
      )}

      {node.data.kind === "unit" && (
        <>
          {Object.entries(p).map(([k, v]) =>
            typeof v === "number" ? (
              <NumField key={k} paramKey={k} value={v} onChange={(nv) => setParam(node.id, k, nv)} />
            ) : (
              <div key={k} className="my-1.5 flex items-center justify-between gap-2 text-muted">
                <span>{k}</span>
                <code className="max-w-[200px] overflow-hidden text-ellipsis whitespace-nowrap text-[11px] text-text">
                  {JSON.stringify(v)}
                </code>
              </div>
            ))}
          {Object.keys(p).length === 0 && (
            <Hint>No parameters set — this unit will use engine defaults.</Hint>
          )}
          <AddParam nodeId={node.id} existing={Object.keys(p)} />
          <DesignPanel unitId={node.id} />
        </>
      )}

      {node.data.kind === "product" && (
        <Hint>Product sink. The stream feeding it is reported by the engine after a solve.</Hint>
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
          <input
            type="number" step="any" min={0} max={1}
            className={`w-[110px] rounded-md border bg-panel2 px-2 py-1 text-right text-text ${
              ok ? "border-line" : "border-warn"
            }`}
            value={z[c] ?? 0}
            aria-label={`Mole fraction of ${c}`}
            onChange={(e) => setParam(nodeId, "z", { ...z, [c]: parseFloat(e.target.value) || 0 })}
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
          <table className="data-table">
            <tbody>
              <tr><td>T</td><td>{state.T?.toFixed(2) ?? "—"} K</td></tr>
              <tr><td>P</td><td>{state.P != null ? (state.P / 1000).toFixed(1) : "—"} kPa</td></tr>
              <tr><td>flow</td><td>{state.molar_flow?.toFixed(3) ?? "—"} mol/s</td></tr>
              <tr><td>phase</td><td>{state.phase ?? "—"}</td></tr>
              <tr><td>vapor frac</td><td>{state.vapor_fraction?.toFixed(3) ?? "—"}</td></tr>
            </tbody>
          </table>
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
