// Editable techno-economic assumptions, driven by a "seed" of effective values
// (engine defaults / catalog, or the values used in the last solve). Edits become
// per-flowsheet overrides in store.costConfig (sent to /cost); each row shows the
// default until overridden, with a reset. Used by the Settings dialog (seeded
// from defaults, always available) and the Econ tab (seeded from the last cost).
import { RotateCcw, X } from "lucide-react";
import { useStore } from "../store";
import type { CostConfigOverrides } from "../types";
import { Button, NumberInput, PanelTitle } from "./ui";

export interface AssumptionsSeed {
  financial: { discount_rate: number; operating_hours: number; project_years: number };
  prices: Record<string, number>;          // component -> $/kg
  utilities: Record<string, number>;        // utility -> $/GJ
  sizing: Record<string, number | string>;
  factors: Record<string, number>;
  citations: { topic: string; source: string }[];
}

const pretty = (k: string) => k.replace(/_/g, " ");
type GroupKey = "prices_per_kg" | "utility_prices" | "sizing" | "factors";

function EditRow({ label, unit, value, overridden, onChange, onReset }: {
  label: string; unit?: string; value: number; overridden: boolean;
  onChange: (v: number) => void; onReset: () => void;
}) {
  return (
    <label className="my-1 flex items-center gap-2 text-[12px]">
      <span className={`min-w-0 flex-1 truncate ${overridden ? "text-accent" : "text-muted"}`}
        title={label}>{label}</span>
      <NumberInput
        className="w-[92px] rounded-md border border-line bg-panel2 px-2 py-0.5 text-right text-text"
        value={Number.isFinite(value) ? Number(value.toPrecision(6)) : NaN}
        aria-label={label} onChange={onChange} />
      <span className="w-12 shrink-0 text-[11px] text-muted">{unit ?? ""}</span>
      {overridden ? (
        <button className="shrink-0 cursor-pointer p-0.5 text-muted hover:text-accent"
          title="Reset to default" aria-label={`Reset ${label}`} onClick={onReset}>
          <X size={11} />
        </button>
      ) : <span className="w-[15px] shrink-0" aria-hidden />}
    </label>
  );
}

function GroupRows({ values, group, unit }: {
  values: Record<string, number | string>; group: GroupKey; unit?: string;
}) {
  const cfg = useStore((s) => s.costConfig);
  const setCostConfig = useStore((s) => s.setCostConfig);
  const ov = (cfg[group] ?? {}) as Record<string, number>;
  const rows = Object.entries(values).filter(([, v]) => typeof v === "number");
  if (!rows.length) return null;
  return (
    <>
      {rows.map(([k, def]) => {
        const overridden = k in ov;
        return (
          <EditRow key={k} label={pretty(k)} unit={unit}
            value={overridden ? ov[k] : (def as number)} overridden={overridden}
            onChange={(v) => setCostConfig({ ...cfg, [group]: { ...ov, [k]: v } })}
            onReset={() => { const g = { ...ov }; delete g[k]; setCostConfig({ ...cfg, [group]: g }); }} />
        );
      })}
    </>
  );
}

export function AssumptionsEditor({ seed }: { seed: AssumptionsSeed }) {
  const cfg = useStore((s) => s.costConfig);
  const setCostConfig = useStore((s) => s.setCostConfig);
  const cost = useStore((s) => s.cost);
  const busy = useStore((s) => s.busy);
  const dirty = Object.keys(cfg).length > 0;

  const scalar = (key: "discount_rate" | "operating_hours" | "project_years",
                  def: number, unit: string) => (
    <EditRow label={pretty(key)} unit={unit}
      value={cfg[key] ?? def} overridden={cfg[key] != null}
      onChange={(v) => setCostConfig({ ...cfg, [key]: v })}
      onReset={() => { const c = { ...cfg }; delete c[key]; setCostConfig(c); }} />
  );

  return (
    <div>
      {dirty && (
        <div className="mb-2 flex items-center gap-2">
          <Button onClick={() => void cost()} busy={busy === "cost"}>Re-cost</Button>
          <button className="flex items-center gap-1 text-[11px] text-muted hover:text-accent"
            onClick={() => setCostConfig({} as CostConfigOverrides)}>
            <RotateCcw size={11} /> reset all
          </button>
          <span className="text-[11px] text-warn">edited — re-cost to apply</span>
        </div>
      )}
      <PanelTitle>Financial</PanelTitle>
      {scalar("discount_rate", seed.financial.discount_rate, "")}
      {scalar("operating_hours", seed.financial.operating_hours, "h/yr")}
      {scalar("project_years", seed.financial.project_years, "yr")}

      {Object.keys(seed.prices).length > 0 && <PanelTitle>Component prices</PanelTitle>}
      <GroupRows values={seed.prices} group="prices_per_kg" unit="$/kg" />

      {Object.keys(seed.utilities).length > 0 && <PanelTitle>Utility prices</PanelTitle>}
      <GroupRows values={seed.utilities} group="utility_prices" unit="$/GJ" />

      <PanelTitle>Sizing heuristics</PanelTitle>
      <GroupRows values={seed.sizing} group="sizing" />

      <PanelTitle>Cost factors (Turton COM + capital)</PanelTitle>
      <GroupRows values={seed.factors} group="factors" />

      <PanelTitle>Sources</PanelTitle>
      <ul className="ml-1 list-none text-[10.5px] leading-relaxed text-muted">
        {seed.citations.map((c) => (
          <li key={c.topic}><span className="text-text">{c.topic}:</span> {c.source}</li>
        ))}
      </ul>
    </div>
  );
}
