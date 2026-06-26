// Pure helpers for the reusable reaction editor (ReactionEditorDialog).
//
// The editor edits a list of `DraftReaction`s (two columns of positive-magnitude
// rows: reactants | products, plus key/conversion/kinetics). These functions
// convert between that draft model and the engine's param shapes
// (`Reaction`/`KineticReaction`.from_param — see engine/caldyr/unitops/reaction.py)
// and validate a draft before save. Kept pure + framework-free so they're unit
// tested directly; the component is a thin shell over them.
import type { KineticReactionSpec, ReactionEditorOpts, ReactionSpec } from "../types";

export interface DraftRow {
  comp: string;
  coeff: number; // positive magnitude (the editor applies signs)
}

export interface DraftReaction {
  reactants: DraftRow[];
  products: DraftRow[];
  key: string; // "" = unset
  conversion: number; // used only when opts.conversion
  // kinetic block (used only when opts.kind === "kinetic")
  k0: number;
  Ea: number;
  orders: DraftRow[]; // {comp, order}; empty = engine default {key: 1}
  reversible: boolean;
  k0_rev: number;
  Ea_rev: number;
  orders_rev: DraftRow[];
}

export function emptyRow(comp = ""): DraftRow {
  return { comp, coeff: 1 };
}

export function emptyDraft(): DraftReaction {
  return {
    reactants: [emptyRow()],
    products: [emptyRow()],
    key: "",
    conversion: 0.8,
    k0: 1,
    Ea: 0,
    orders: [],
    reversible: false,
    k0_rev: 0,
    Ea_rev: 0,
    orders_rev: [],
  };
}

// Split a signed stoich map into positive-magnitude reactant/product rows.
function rowsFromStoich(stoich: Record<string, number>): { reactants: DraftRow[]; products: DraftRow[] } {
  const reactants: DraftRow[] = [];
  const products: DraftRow[] = [];
  for (const [comp, nu] of Object.entries(stoich)) {
    const n = Number(nu);
    if (n < 0) reactants.push({ comp, coeff: -n });
    else if (n > 0) products.push({ comp, coeff: n });
  }
  if (reactants.length === 0) reactants.push(emptyRow());
  if (products.length === 0) products.push(emptyRow());
  return { reactants, products };
}

function ordersToRows(orders?: Record<string, number>): DraftRow[] {
  if (!orders) return [];
  return Object.entries(orders).map(([comp, order]) => ({ comp, coeff: Number(order) }));
}

// Build one draft from an engine reaction spec (stoichiometric or kinetic).
function draftFromSpec(spec: ReactionSpec | KineticReactionSpec): DraftReaction {
  const d = emptyDraft();
  const { reactants, products } = rowsFromStoich(spec.stoich ?? {});
  d.reactants = reactants;
  d.products = products;
  if (spec.key) d.key = spec.key;
  const r = spec as ReactionSpec;
  if (typeof r.conversion === "number") d.conversion = r.conversion;
  const k = spec as KineticReactionSpec;
  if (typeof k.k0 === "number") d.k0 = k.k0;
  if (typeof k.Ea === "number") d.Ea = k.Ea;
  d.orders = ordersToRows(k.orders);
  if (typeof k.k0_rev === "number" && k.k0_rev > 0) {
    d.reversible = true;
    d.k0_rev = k.k0_rev;
    d.Ea_rev = typeof k.Ea_rev === "number" ? k.Ea_rev : 0;
    d.orders_rev = ordersToRows(k.orders_rev);
  }
  return d;
}

// Normalize a unit's stored reaction params (singular `reaction` (+ unit-level
// `conversion`) OR list `reactions`) into the editor's draft list. Returns one
// empty draft when nothing is set yet so the editor always has a card to fill.
export function normalizeReactions(params: Record<string, unknown>): DraftReaction[] {
  const list = params.reactions;
  if (Array.isArray(list) && list.length > 0) {
    return list.map((r) => draftFromSpec(r as ReactionSpec));
  }
  const single = params.reaction;
  if (single && typeof single === "object") {
    const spec = { ...(single as ReactionSpec) };
    if (typeof params.conversion === "number" && typeof spec.conversion !== "number") {
      spec.conversion = params.conversion as number;
    }
    return [draftFromSpec(spec)];
  }
  return [emptyDraft()];
}

// Collapse draft rows to a signed stoich map. Rows with no component are
// dropped; duplicate components on the same side sum.
function stoichFromDraft(d: DraftReaction): Record<string, number> {
  const stoich: Record<string, number> = {};
  for (const row of d.reactants) {
    if (!row.comp) continue;
    stoich[row.comp] = (stoich[row.comp] ?? 0) - Math.abs(row.coeff);
  }
  for (const row of d.products) {
    if (!row.comp) continue;
    stoich[row.comp] = (stoich[row.comp] ?? 0) + Math.abs(row.coeff);
  }
  return stoich;
}

function rowsToOrders(rows: DraftRow[]): Record<string, number> {
  const out: Record<string, number> = {};
  for (const row of rows) if (row.comp) out[row.comp] = row.coeff;
  return out;
}

// Serialize the draft list to the list-form `reactions` value the engine
// accepts for all four reactor types. (EquilibriumReactor uses only the first.)
export function serializeReactions(
  drafts: DraftReaction[],
  opts: ReactionEditorOpts,
): (ReactionSpec | KineticReactionSpec)[] {
  return drafts.map((d) => {
    const stoich = stoichFromDraft(d);
    if (opts.kind === "kinetic") {
      const spec: KineticReactionSpec = { stoich, key: d.key, k0: d.k0, Ea: d.Ea };
      if (d.orders.length > 0) spec.orders = rowsToOrders(d.orders);
      if (opts.reversible && d.reversible) {
        spec.k0_rev = d.k0_rev;
        spec.Ea_rev = d.Ea_rev;
        if (d.orders_rev.length > 0) spec.orders_rev = rowsToOrders(d.orders_rev);
      }
      return spec;
    }
    const spec: ReactionSpec = { stoich };
    if (d.key) spec.key = d.key;
    if (opts.conversion) spec.conversion = d.conversion;
    return spec;
  });
}

// Human-readable preview line, e.g. "nitrogen + 3 hydrogen → 2 ammonia".
export function previewReaction(d: DraftReaction): string {
  const side = (rows: DraftRow[]) =>
    rows
      .filter((r) => r.comp)
      .map((r) => (Math.abs(r.coeff) === 1 ? r.comp : `${r.coeff} ${r.comp}`))
      .join(" + ");
  const lhs = side(d.reactants) || "…";
  const rhs = side(d.products) || "…";
  return `${lhs} → ${rhs}`;
}

// Validate a draft list; returns one human-readable error string per problem
// (empty = OK to save). `components` is the flowsheet's declared species.
export function validateReactions(
  drafts: DraftReaction[],
  opts: ReactionEditorOpts,
  components: string[],
): string[] {
  const errors: string[] = [];
  if (drafts.length === 0) errors.push("Add at least one reaction.");
  const known = new Set(components);
  drafts.forEach((d, i) => {
    const tag = drafts.length > 1 ? `Reaction ${i + 1}: ` : "";
    const reactComps = d.reactants.filter((r) => r.comp);
    const prodComps = d.products.filter((r) => r.comp);
    if (reactComps.length === 0) errors.push(`${tag}needs at least one reactant.`);
    if (prodComps.length === 0) errors.push(`${tag}needs at least one product.`);

    const stoich = stoichFromDraft(d);
    const reactantIds = new Set(reactComps.map((r) => r.comp));
    for (const [comp, nu] of Object.entries(stoich)) {
      if (!known.has(comp)) errors.push(`${tag}"${comp}" is not a flowsheet component.`);
      if (nu === 0) errors.push(`${tag}"${comp}" appears on both sides and cancels.`);
    }
    for (const row of [...d.reactants, ...d.products]) {
      if (row.comp && !(row.coeff > 0)) errors.push(`${tag}coefficient for "${row.comp}" must be > 0.`);
    }

    const keyNeeded = opts.key_required || opts.kind === "kinetic";
    if (keyNeeded && !d.key) errors.push(`${tag}choose a key reactant.`);
    if (d.key && !reactantIds.has(d.key)) errors.push(`${tag}key "${d.key}" must be a reactant.`);

    if (opts.conversion) {
      if (!(d.conversion > 0 && d.conversion <= 1)) errors.push(`${tag}conversion must be in (0, 1].`);
    }

    if (opts.kind === "kinetic") {
      if (!(d.k0 > 0)) errors.push(`${tag}k0 must be > 0.`);
      if (!(d.Ea >= 0)) errors.push(`${tag}Ea must be ≥ 0.`);
      for (const row of d.orders) {
        if (row.comp && !known.has(row.comp)) errors.push(`${tag}order references unknown component "${row.comp}".`);
      }
      if (opts.reversible && d.reversible) {
        if (!(d.k0_rev > 0)) errors.push(`${tag}reverse k0 must be > 0 when reversible.`);
        if (!(d.Ea_rev >= 0)) errors.push(`${tag}reverse Ea must be ≥ 0.`);
        for (const row of d.orders_rev) {
          if (row.comp && !known.has(row.comp)) errors.push(`${tag}reverse order references unknown component "${row.comp}".`);
        }
      }
    }
  });
  return errors;
}

// ---- atom-balance hint (non-blocking) ----------------------------------
//
// Components carry a Hill-system formula (e.g. "N2", "H3N", "C6H12") from
// GET /components. We parse it to element counts and check that Σ νᵢ·(elements)
// is zero across the signed stoichiometry. This is advisory only — the engine
// is balance-agnostic and lumped/pseudo-components legitimately won't balance.

export type BalanceStatus = "balanced" | "unbalanced" | "unknown";

export interface BalanceResult {
  status: BalanceStatus;
  deltas: Record<string, number>; // element -> (products − reactants); nonzero only
  missing: string[]; // components with no parseable formula (status === "unknown")
}

// Parse a Hill formula into element counts; null if it contains anything we
// can't interpret (charges, hydrate dots, R-groups), so the caller treats the
// reaction as "balance unknown" rather than silently wrong.
export function parseFormula(formula: string | undefined): Record<string, number> | null {
  if (!formula || !/^[A-Za-z0-9()]+$/.test(formula)) return null;
  const stack: Record<string, number>[] = [{}];
  let i = 0;
  while (i < formula.length) {
    const ch = formula[i];
    if (ch === "(") {
      stack.push({});
      i++;
    } else if (ch === ")") {
      i++;
      let num = "";
      while (i < formula.length && /\d/.test(formula[i])) num += formula[i++];
      const mult = num ? parseInt(num, 10) : 1;
      const top = stack.pop();
      if (!top || stack.length === 0) return null;
      const below = stack[stack.length - 1];
      for (const [el, n] of Object.entries(top)) below[el] = (below[el] ?? 0) + n * mult;
    } else if (/[A-Z]/.test(ch)) {
      let sym = ch;
      i++;
      while (i < formula.length && /[a-z]/.test(formula[i])) sym += formula[i++];
      let num = "";
      while (i < formula.length && /\d/.test(formula[i])) num += formula[i++];
      const cnt = num ? parseInt(num, 10) : 1;
      const top = stack[stack.length - 1];
      top[sym] = (top[sym] ?? 0) + cnt;
    } else {
      return null;
    }
  }
  return stack.length === 1 ? stack[0] : null;
}

// Element imbalance of one draft reaction given a component-id → formula map.
export function atomBalance(d: DraftReaction, formulaById: Record<string, string>): BalanceResult {
  const stoich = stoichFromDraft(d);
  const comps = Object.keys(stoich).filter((c) => stoich[c] !== 0);
  const parsed: Record<string, Record<string, number>> = {};
  const missing: string[] = [];
  for (const c of comps) {
    const els = parseFormula(formulaById[c]);
    if (els) parsed[c] = els;
    else missing.push(c);
  }
  if (missing.length > 0) return { status: "unknown", deltas: {}, missing };
  const deltas: Record<string, number> = {};
  for (const c of comps) {
    for (const [el, n] of Object.entries(parsed[c])) deltas[el] = (deltas[el] ?? 0) + stoich[c] * n;
  }
  const nonzero: Record<string, number> = {};
  for (const [el, v] of Object.entries(deltas)) if (Math.abs(v) > 1e-9) nonzero[el] = v;
  return { status: Object.keys(nonzero).length > 0 ? "unbalanced" : "balanced", deltas: nonzero, missing: [] };
}

// "H: −2, N: +1" from a delta map (signed, products − reactants).
export function formatBalanceDelta(deltas: Record<string, number>): string {
  return Object.entries(deltas)
    .map(([el, v]) => `${el}: ${v > 0 ? "+" : "−"}${Math.abs(Math.round(v * 1e6) / 1e6)}`)
    .join(", ");
}

// One-line summary for the Inspector row (e.g. "1 reaction: nitrogen + 3 hydrogen → 2 ammonia").
export function summarizeReactions(params: Record<string, unknown>): string {
  const hasAny =
    (Array.isArray(params.reactions) && params.reactions.length > 0) ||
    (params.reaction && typeof params.reaction === "object");
  if (!hasAny) return "— none —";
  const drafts = normalizeReactions(params);
  const first = previewReaction(drafts[0]);
  if (drafts.length === 1) return `1 reaction: ${first}`;
  return `${drafts.length} reactions: ${first} …`;
}
