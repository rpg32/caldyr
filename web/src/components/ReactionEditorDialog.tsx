// Reusable reaction editor — one modal for every reaction-bearing unit op
// (EquilibriumReactor, ConversionReactor, CSTR, PFR). The unit's param schema
// carries an `editor: "reaction"` capability descriptor (editor_opts) that drives
// which controls appear; the editor always serializes to the engine's list form
// `reactions` (and clears the singular `reaction`/`conversion` keys on save).
// Pure model/validation logic lives in lib/reactions.ts; this is a thin shell.
import { Plus, Trash2, X } from "lucide-react";
import { useEffect, useMemo, useState } from "react";
import {
  atomBalance, emptyDraft, emptyRow, formatBalanceDelta, normalizeReactions,
  previewReaction, serializeReactions, validateReactions,
  type DraftReaction, type DraftRow,
} from "../lib/reactions";
import { useStore } from "../store";
import type { ReactionEditorOpts } from "../types";
import { Button, NumberInput, PanelTitle } from "./ui";

const SEL = "rounded-md border border-line bg-panel2 px-1.5 py-1 text-[12px] text-text";
const NUM = "w-[72px] rounded-md border border-line bg-panel2 px-1.5 py-1 text-right text-[12px] text-text";

export function ReactionEditorDialog() {
  const target = useStore((s) => s.reactionEditorTarget);
  const close = useStore((s) => s.closeReactionEditor);
  const setReactions = useStore((s) => s.setReactions);
  const components = useStore((s) => s.components);
  const catalog = useStore((s) => s.componentCatalog);
  const node = useStore((s) => (target ? s.nodes.find((n) => n.id === target.nodeId) : undefined));

  // component id -> Hill formula (for the advisory atom-balance hint).
  const formulaById = useMemo(
    () => Object.fromEntries(catalog.map((c) => [c.id, c.formula])),
    [catalog],
  );

  const opts = (target?.schema.editor_opts ?? {
    kind: "stoichiometric", multiple: false, conversion: false, key_required: false,
  }) as ReactionEditorOpts;

  const [drafts, setDrafts] = useState<DraftReaction[]>([emptyDraft()]);

  // Reset the working drafts whenever a new unit's editor is opened.
  useEffect(() => {
    if (target && node) setDrafts(normalizeReactions(node.data.params));
  }, [target?.nodeId]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (!target) return;
    const onKey = (e: KeyboardEvent) => e.key === "Escape" && close();
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [target, close]);

  const errors = useMemo(
    () => (target ? validateReactions(drafts, opts, components) : []),
    [drafts, opts, components, target],
  );

  if (!target || !node) return null;

  const patch = (i: number, fn: (d: DraftReaction) => DraftReaction) =>
    setDrafts((ds) => ds.map((d, idx) => (idx === i ? fn(d) : d)));

  const save = () => {
    if (errors.length > 0) return;
    setReactions(target.nodeId, serializeReactions(drafts, opts));
    close();
  };

  return (
    <div className="fixed inset-0 z-40 flex items-center justify-center bg-black/40" onClick={close}>
      <div role="dialog" aria-modal="true" aria-label="Edit reactions"
        className="max-h-[88vh] w-[680px] overflow-auto rounded-xl border border-line bg-panel p-4 shadow-2xl"
        onClick={(e) => e.stopPropagation()}>
        <div className="mb-1 flex items-center">
          <b className="text-[15px]">Edit reactions — {node.data.unitType}</b>
          <button className="ml-auto cursor-pointer text-muted hover:text-text"
            onClick={close} aria-label="Close reaction editor"><X size={16} /></button>
        </div>
        <p className="mb-2 text-[11px] text-muted">
          Enter positive coefficients in each column; the editor signs reactants negative.
          {opts.kind === "kinetic"
            ? " Rate = k₀·exp(−Eₐ/RT)·∏Cᵢ^orderᵢ (mol/m³·s); Eₐ in J/mol."
            : opts.conversion
              ? " Each reaction applies its fractional conversion of the key reactant, in order."
              : " Extent is solved from the equilibrium constant at the reactor temperature."}
        </p>

        {components.length === 0 && (
          <p className="mb-2 rounded-md border border-warn/40 bg-warn/5 px-2 py-1 text-[11px] text-warn">
            No components defined yet — add them in the Flowsheet panel first.
          </p>
        )}

        {drafts.map((d, i) => (
          <ReactionCard key={i} d={d} index={i} total={drafts.length} opts={opts}
            components={components} formulaById={formulaById}
            onChange={(fn) => patch(i, fn)}
            onRemove={drafts.length > 1 ? () => setDrafts((ds) => ds.filter((_, idx) => idx !== i)) : undefined} />
        ))}

        {opts.multiple && (
          <Button icon={<Plus size={13} />} className="mt-1"
            onClick={() => setDrafts((ds) => [...ds, emptyDraft()])}>Add reaction</Button>
        )}

        {errors.length > 0 && (
          <ul className="mt-3 list-disc rounded-md border border-bad/40 bg-bad/5 py-1 pl-5 pr-2 text-[11px] text-bad">
            {errors.map((e, i) => <li key={i}>{e}</li>)}
          </ul>
        )}

        <div className="mt-3 flex items-center gap-2">
          <Button variant="primary" disabled={errors.length > 0} onClick={save}>Save</Button>
          <Button onClick={close}>Cancel</Button>
        </div>
      </div>
    </div>
  );
}

function ReactionCard({ d, index, total, opts, components, formulaById, onChange, onRemove }: {
  d: DraftReaction; index: number; total: number; opts: ReactionEditorOpts;
  components: string[]; formulaById: Record<string, string>;
  onChange: (fn: (d: DraftReaction) => DraftReaction) => void;
  onRemove?: () => void;
}) {
  const reactantOpts = d.reactants.map((r) => r.comp).filter(Boolean);
  const keyNeeded = opts.key_required || opts.kind === "kinetic";
  // Advisory atom balance — only meaningful once both sides have a component.
  const bothSides = d.reactants.some((r) => r.comp) && d.products.some((r) => r.comp);
  const bal = bothSides ? atomBalance(d, formulaById) : null;

  return (
    <div className="my-2 rounded-lg border border-line bg-panel2/40 p-2.5">
      {total > 1 && (
        <div className="mb-1 flex items-center">
          <span className="text-[11px] font-semibold text-muted">Reaction {index + 1}</span>
          {onRemove && (
            <button className="ml-auto cursor-pointer p-0.5 text-muted hover:text-bad"
              onClick={onRemove} aria-label={`Remove reaction ${index + 1}`}><Trash2 size={13} /></button>
          )}
        </div>
      )}

      <div className="grid grid-cols-2 gap-3">
        <StoichColumn title="Reactants" rows={d.reactants} components={components}
          onChange={(rows) => onChange((cur) => ({ ...cur, reactants: rows }))} />
        <StoichColumn title="Products" rows={d.products} components={components}
          onChange={(rows) => onChange((cur) => ({ ...cur, products: rows }))} />
      </div>

      <div className="mt-2 rounded-md bg-panel px-2 py-1 text-center text-[12px] text-text">
        {previewReaction(d)}
      </div>

      {bal && bal.status === "unbalanced" && (
        <div className="mt-1 text-center text-[11px] text-warn"
          title="Advisory only — the solver does not require a balanced reaction. Lumped or pseudo-components may not balance.">
          ⚠ Not atom-balanced — net {formatBalanceDelta(bal.deltas)} (products − reactants)
        </div>
      )}
      {bal && bal.status === "balanced" && (
        <div className="mt-1 text-center text-[11px] text-ok">✓ Atom-balanced</div>
      )}
      {bal && bal.status === "unknown" && (
        <div className="mt-1 text-center text-[11px] text-muted">
          Atom balance unchecked — no formula for {bal.missing.join(", ")}
        </div>
      )}

      <div className="mt-2 flex flex-wrap items-center gap-x-4 gap-y-1.5">
        {/* Key reactant is always offered (recommended for equilibrium, required
            for conversion/kinetic); the asterisk marks when it's required. */}
          <label className="flex items-center gap-1.5 text-[12px] text-muted">
            Key reactant{keyNeeded && <span className="text-bad" title="required">*</span>}
            <select className={SEL} value={d.key} aria-label="Key reactant"
              onChange={(e) => onChange((cur) => ({ ...cur, key: e.target.value }))}>
              <option value="">—</option>
              {reactantOpts.map((c) => <option key={c} value={c}>{c}</option>)}
            </select>
          </label>
          {opts.conversion && (
            <label className="flex items-center gap-1.5 text-[12px] text-muted">
              Conversion
              <NumberInput className={NUM} value={d.conversion} min={0} max={1} aria-label="Conversion"
                onChange={(v) => onChange((cur) => ({ ...cur, conversion: v }))} />
            </label>
          )}
        </div>

      {opts.kind === "kinetic" && (
        <KineticBlock d={d} opts={opts} components={components} onChange={onChange} />
      )}
    </div>
  );
}

function StoichColumn({ title, rows, components, onChange }: {
  title: string; rows: DraftRow[]; components: string[];
  onChange: (rows: DraftRow[]) => void;
}) {
  const setRow = (i: number, patch: Partial<DraftRow>) =>
    onChange(rows.map((r, idx) => (idx === i ? { ...r, ...patch } : r)));
  return (
    <div>
      <PanelTitle>{title}</PanelTitle>
      {rows.map((r, i) => (
        <div key={i} className="my-1 flex items-center gap-1">
          <NumberInput className={`${NUM} w-[52px]`} value={r.coeff} min={0} aria-label={`${title} coefficient`}
            onChange={(v) => setRow(i, { coeff: v })} />
          <select className={`${SEL} min-w-0 flex-1`} value={r.comp} aria-label={`${title} component`}
            onChange={(e) => setRow(i, { comp: e.target.value })}>
            <option value="">—</option>
            {components.map((c) => <option key={c} value={c}>{c}</option>)}
          </select>
          <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
            onClick={() => onChange(rows.filter((_, idx) => idx !== i))}
            aria-label={`Remove ${title} row`} disabled={rows.length === 1}><X size={12} /></button>
        </div>
      ))}
      <button className="mt-0.5 flex items-center gap-1 text-[11px] text-muted hover:text-accent"
        onClick={() => onChange([...rows, emptyRow()])}><Plus size={11} /> add</button>
    </div>
  );
}

function KineticBlock({ d, opts, components, onChange }: {
  d: DraftReaction; opts: ReactionEditorOpts; components: string[];
  onChange: (fn: (d: DraftReaction) => DraftReaction) => void;
}) {
  return (
    <div className="mt-2 border-t border-line pt-2">
      <div className="flex flex-wrap items-center gap-x-4 gap-y-1.5">
        <label className="flex items-center gap-1.5 text-[12px] text-muted">
          k₀<span className="text-bad" title="required">*</span>
          <NumberInput className={NUM} value={d.k0} min={0} aria-label="k0"
            onChange={(v) => onChange((cur) => ({ ...cur, k0: v }))} />
        </label>
        <label className="flex items-center gap-1.5 text-[12px] text-muted">
          Eₐ (J/mol)<span className="text-bad" title="required">*</span>
          <NumberInput className={`${NUM} w-[96px]`} value={d.Ea} min={0} aria-label="Ea"
            onChange={(v) => onChange((cur) => ({ ...cur, Ea: v }))} />
        </label>
      </div>

      <OrderList label="Reaction orders" placeholder="default: key¹"
        rows={d.orders} components={components}
        onChange={(rows) => onChange((cur) => ({ ...cur, orders: rows }))} />

      {opts.reversible && (
        <>
          <label className="mt-2 flex items-center gap-1.5 text-[12px] text-muted">
            <input type="checkbox" className="h-4 w-4 accent-accent" checked={d.reversible}
              onChange={(e) => onChange((cur) => ({ ...cur, reversible: e.target.checked }))} />
            Reversible (adds a reverse rate term)
          </label>
          {d.reversible && (
            <div className="mt-1.5 rounded-md border border-line bg-panel p-2">
              <div className="flex flex-wrap items-center gap-x-4 gap-y-1.5">
                <label className="flex items-center gap-1.5 text-[12px] text-muted">
                  k₀,rev
                  <NumberInput className={NUM} value={d.k0_rev} min={0} aria-label="reverse k0"
                    onChange={(v) => onChange((cur) => ({ ...cur, k0_rev: v }))} />
                </label>
                <label className="flex items-center gap-1.5 text-[12px] text-muted">
                  Eₐ,rev (J/mol)
                  <NumberInput className={`${NUM} w-[96px]`} value={d.Ea_rev} min={0} aria-label="reverse Ea"
                    onChange={(v) => onChange((cur) => ({ ...cur, Ea_rev: v }))} />
                </label>
              </div>
              <OrderList label="Reverse orders" placeholder="default: product stoich"
                rows={d.orders_rev} components={components}
                onChange={(rows) => onChange((cur) => ({ ...cur, orders_rev: rows }))} />
            </div>
          )}
        </>
      )}
    </div>
  );
}

function OrderList({ label, placeholder, rows, components, onChange }: {
  label: string; placeholder: string; rows: DraftRow[]; components: string[];
  onChange: (rows: DraftRow[]) => void;
}) {
  const setRow = (i: number, patch: Partial<DraftRow>) =>
    onChange(rows.map((r, idx) => (idx === i ? { ...r, ...patch } : r)));
  return (
    <div className="mt-1.5">
      <div className="flex items-center gap-2">
        <span className="text-[11px] text-muted">{label}</span>
        {rows.length === 0 && <span className="text-[11px] text-muted italic">({placeholder})</span>}
      </div>
      {rows.map((r, i) => (
        <div key={i} className="my-1 flex items-center gap-1">
          <select className={`${SEL} min-w-0 flex-1`} value={r.comp} aria-label="Order component"
            onChange={(e) => setRow(i, { comp: e.target.value })}>
            <option value="">—</option>
            {components.map((c) => <option key={c} value={c}>{c}</option>)}
          </select>
          <NumberInput className={NUM} value={r.coeff} aria-label="Order exponent"
            onChange={(v) => setRow(i, { coeff: v })} />
          <button className="cursor-pointer p-0.5 text-muted hover:text-bad"
            onClick={() => onChange(rows.filter((_, idx) => idx !== i))}
            aria-label="Remove order row"><X size={12} /></button>
        </div>
      ))}
      <button className="mt-0.5 flex items-center gap-1 text-[11px] text-muted hover:text-accent"
        onClick={() => onChange([...rows, emptyRow()])}><Plus size={11} /> add order</button>
    </div>
  );
}
