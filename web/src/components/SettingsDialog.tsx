// Settings dialog: techno-economic assumptions editable ANY TIME (no solve
// required), seeded from the engine defaults (/cost-defaults) + price catalog
// (/prices). Edits persist per-flowsheet (store.costConfig, in meta.ui) and feed
// the next Cost. The same editor appears in the Econ tab seeded from the last solve.
import { X } from "lucide-react";
import { useEffect } from "react";
import { useStore } from "../store";
import { AssumptionsEditor, type AssumptionsSeed } from "./AssumptionsEditor";
import { Hint } from "./ui";

export function SettingsDialog() {
  const open = useStore((s) => s.settingsOpen);
  const toggle = useStore((s) => s.toggleSettings);
  const defaults = useStore((s) => s.costDefaults);
  const priceCatalog = useStore((s) => s.priceCatalog);
  const components = useStore((s) => s.components);
  const product = useStore((s) => s.product);

  useEffect(() => {
    if (!open) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && toggle();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open, toggle]);

  if (!open) return null;

  // Seed prices with the flowsheet's components (+ product); catalog default when
  // known, else NaN so the field shows blank/editable (e.g. propane has no price).
  const compNames = [...new Set([...components, product].filter(Boolean))];
  const prices: Record<string, number> = {};
  for (const c of compNames) prices[c] = priceCatalog.prices_per_kg[c] ?? NaN;

  const seed: AssumptionsSeed | null = defaults && {
    financial: {
      discount_rate: defaults.config.discount_rate,
      operating_hours: defaults.config.operating_hours,
      project_years: defaults.config.project_years,
    },
    prices,
    utilities: priceCatalog.utility_prices,
    sizing: defaults.sizing,
    factors: defaults.factors,
    citations: defaults.citations,
  };

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40"
      onClick={toggle}>
      <div role="dialog" aria-modal="true" aria-label="Cost assumptions settings"
        className="max-h-[85vh] w-[460px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-1 flex items-center">
          <b className="text-[15px]">Cost assumptions</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={toggle} aria-label="Close settings"><X size={16} /></button>
        </div>
        <p className="mb-2 text-[11px] text-muted">
          Defaults for every flowsheet; edits are saved with this flowsheet and used
          on the next Cost. The engine stays SI internally.
        </p>
        {seed ? <AssumptionsEditor seed={seed} />
          : <Hint>Assumption defaults unavailable — is the engine running?</Hint>}
      </div>
    </div>
  );
}
