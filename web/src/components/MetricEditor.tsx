// Small shared editor for a MetricSpec (used by Optimize + Study panels).
import { useStore } from "../store";
import type { MetricSpec } from "../types";

export function MetricEditor({ value, onChange }: {
  value: MetricSpec;
  onChange: (m: MetricSpec) => void;
}) {
  const edges = useStore((s) => s.edges);
  const components = useStore((s) => s.components);
  const solveRes = useStore((s) => s.solveRes);
  const dutyKeys = Object.keys(solveRes?.report.duties ?? {});
  const sel = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-text min-w-0";

  return (
    <div className="flex flex-wrap items-center gap-1.5">
      <select className={sel} value={value.type} aria-label="Metric type"
        onChange={(e) => onChange({ ...value, type: e.target.value as MetricSpec["type"], stream: "" })}>
        <option value="flow">stream flow</option>
        <option value="component_rate">component rate</option>
        <option value="duty">duty</option>
      </select>
      {value.type === "duty" ? (
        <select className={sel} value={value.stream} aria-label="Duty"
          onChange={(e) => onChange({ ...value, stream: e.target.value })}>
          <option value="">— duty —</option>
          {dutyKeys.map((k) => <option key={k} value={k}>{k}</option>)}
        </select>
      ) : (
        <select className={sel} value={value.stream} aria-label="Stream"
          onChange={(e) => onChange({ ...value, stream: e.target.value })}>
          <option value="">— stream —</option>
          {edges.map((e) => <option key={e.id} value={e.id}>{e.id}</option>)}
        </select>
      )}
      {value.type === "component_rate" && (
        <select className={sel} value={value.component ?? ""} aria-label="Component"
          onChange={(e) => onChange({ ...value, component: e.target.value })}>
          <option value="">— component —</option>
          {components.map((c) => <option key={c} value={c}>{c}</option>)}
        </select>
      )}
    </div>
  );
}

export const metricValid = (m: MetricSpec): boolean =>
  !!m.stream && (m.type !== "component_rate" || !!m.component);

export const metricLabel = (m: MetricSpec): string =>
  m.type === "duty" ? `duty ${m.stream}`
  : m.type === "flow" ? `flow of ${m.stream}`
  : `${m.component} rate in ${m.stream}`;
