// Shared types mirroring the engine API (see api/models.py).

export interface Port {
  name: string;
  direction: "inlet" | "outlet";
  kind: "material" | "energy";
}

// Per-unit parameter schema served by /unit-types (engine source of truth).
export interface ParamSchema {
  name: string;
  label: string;
  type: "number" | "int" | "boolean" | "select" | "string" | "json";
  default?: unknown;
  required?: boolean;
  unit?: string;                       // SI unit suffix
  dim?: string;                        // temperature|pressure|molar_flow|power|UA
  min?: number;
  max?: number;
  options?: string[];
  requires?: Record<string, unknown>;  // show only when these params match
  complex?: boolean;                   // list/dict value — not a scalar widget
  editor?: "reaction";                 // dedicated structured editor for this param
  editor_opts?: ReactionEditorOpts;    // descriptor driving the reaction editor
  help?: string;
}

// Capability descriptor served per reaction-bearing param (see
// api/param_schemas.py). Drives the reusable ReactionEditorDialog so the web
// renders whatever the API declares — no per-unit logic in the client.
export interface ReactionEditorOpts {
  kind: "stoichiometric" | "kinetic";
  multiple: boolean;       // Add-reaction allowed (more than one reaction)
  conversion: boolean;     // per-reaction conversion field (ConversionReactor)
  key_required: boolean;   // a key reactant must be chosen
  reversible?: boolean;    // kinetic reverse params available (kinetic only)
}

// Stoichiometric reaction param value (Reaction.from_param in the engine).
// `stoich` maps component id -> signed coefficient (reactants negative).
export interface ReactionSpec {
  stoich: Record<string, number>;
  key?: string;
  conversion?: number;     // ConversionReactor list form only (0..1)
}

// Kinetic reaction param value (KineticReaction.from_param in the engine).
export interface KineticReactionSpec {
  stoich: Record<string, number>;
  key: string;
  k0: number;
  Ea: number;
  orders?: Record<string, number>;
  k0_rev?: number;
  Ea_rev?: number;
  orders_rev?: Record<string, number>;
}

export interface UnitType {
  type: string;
  doc: string;
  description?: string;   // first docstring paragraph — palette tooltip
  ports: Port[];
  params?: ParamSchema[];
}

export interface PropertyPackage {
  id: string;
  name: string;
  use: string;
}

// Canvas node payload. `kind` distinguishes engine units from the UI-only
// boundary nodes (feeds carry a spec; products are sinks).
export interface NodeData {
  kind: "unit" | "feed" | "product";
  label: string;
  unitType?: string;
  ports: Port[];
  params: Record<string, unknown>;
  [key: string]: unknown;
}

export interface StreamState {
  id: string;
  T: number | null;
  P: number | null;
  molar_flow: number | null;
  z: Record<string, number>;
  H: number | null;
  phase: string | null;
  vapor_fraction: number | null;
}

export interface SolveReport {
  converged: boolean;
  iterations: number;
  residual: number | null;
  method: string;
  tear_streams: string[];
  duties: Record<string, number>;
  messages: string[];
  history: number[];
}

export interface SolveResponse {
  report: SolveReport;
  streams: Record<string, StreamState>;
  // per-unit design results (column profiles, FUG numbers, fuel duty, ...)
  designs?: Record<string, Record<string, unknown>>;
  // per-component molar mass (kg/mol) for mass-flow / mass-fraction display
  molar_mass?: Record<string, number>;
}

export interface CostResponse {
  report: SolveReport;
  annual_production_kg: number;
  annual_revenue: number;
  capital: { isbl: number; osbl: number; grassroots: number; working_capital: number; tci: number };
  opex: { raw_materials: number; utilities: number; fixed: number; total: number };
  profitability: { lcop: number; npv: number; irr: number | null; payback_years: number | null };
  equipment: {
    unit_id: string; type: string; attribute: number; attribute_name: string;
    bare_module: number; utility: string | null;
  }[];
  tornado?: { variable: string; low_lcop: number; high_lcop: number; swing: number }[];
  monte_carlo?: {
    n: number;
    lcop: { p10: number; p50: number; p90: number; mean: number; std: number };
    npv: { p10: number; p50: number; p90: number; mean: number; std: number };
    lcop_samples: number[];
  };
  assumptions?: CostAssumptions;
}

// The numbers + correlations that drove a cost result (from /cost).
export interface CostAssumptions {
  config: {
    year: number; operating_hours: number; discount_rate: number;
    project_years: number; product_component: string; product_min_fraction: number;
  };
  prices_per_kg: Record<string, number>;
  utility_prices: Record<string, number>;
  sizing: Record<string, number | string>;
  factors: Record<string, number>;
  citations: { topic: string; source: string }[];
}

// Default price catalog (from /prices) for pre-filling cost assumptions.
export interface PriceCatalog {
  prices_per_kg: Record<string, number>;
  utility_prices: Record<string, number>;
  prices_source: string;
}

// Default TEA assumptions (from /cost-defaults) — seeds the Settings editor
// before any cost has run.
export interface CostDefaults {
  config: {
    year: number; operating_hours: number; discount_rate: number;
    project_years: number; product_min_fraction: number;
  };
  sizing: Record<string, number | string>;
  factors: Record<string, number>;
  citations: { topic: string; source: string }[];
}

// A named case: a snapshot of the flowsheet's unit/feed parameters + the cost
// assumptions, so users can keep "base", "cheaper-N2", "high-capacity", … and
// compare them. Persisted per-flowsheet in meta.ui.
export interface Scenario {
  name: string;
  params: Record<string, Record<string, unknown>>; // nodeId -> params snapshot
  costConfig: CostConfigOverrides;
}

export interface ScenarioResult {
  name: string;
  ok: boolean;
  lcop?: number;
  tci?: number;
  opex?: number;
  npv?: number;
  error?: string;
}

// Editable cost-config overrides sent to /cost (all optional; engine SI/defaults).
export interface CostConfigOverrides {
  prices_per_kg?: Record<string, number>;
  utility_prices?: Record<string, number>;
  sizing?: Record<string, number>;
  factors?: Record<string, number>;
  discount_rate?: number;
  operating_hours?: number;
  project_years?: number;
  product_component?: string;
}

// ---- optimization (mirrors api/models.py) ----
export interface MetricSpec {
  type: "duty" | "flow" | "component_rate";
  stream: string;
  component?: string | null;
}
export interface DesignVarSpec {
  unit_id: string;
  param: string;
  lower: number;
  upper: number;
  initial?: number | null;
}
export interface ConstraintSpec {
  metric: MetricSpec;
  op: ">=" | "<=";
  value: number;
}
export interface OptimizeRequest {
  flow: FlowDoc;
  backend?: string;
  objective: { sense: "min" | "max"; metric: MetricSpec };
  design_vars: DesignVarSpec[];
  constraints: ConstraintSpec[];
}
export interface OptimizeResponse {
  success: boolean;
  objective: number;
  design: Record<string, number>;
  n_solves: number;
  message: string;
  report: SolveReport;
  streams: Record<string, StreamState>;
}

export interface EnvelopeResponse {
  stream: string;
  z: Record<string, number>;
  points: { P: number; T_bubble: number; T_dew: number }[];
}

// ---- analysis tools: property table / relief / pinch ----
export interface PropertyTableResponse {
  T: number[];
  P: number[];
  props: string[];
  z: Record<string, number>;
  values: Record<string, (number | null)[][]>;   // [name][i_T][j_P]
  failures: { T: number; P: number; error: string }[];
}

export interface ReliefResponse {
  area_m2: number;
  area_cm2: number;
  orifice: string | null;
  orifice_area_m2: number | null;
  capacity_used: number | null;
  phase: string;
  critical: boolean | null;
  details: Record<string, number>;
  notes: string[];
}

export interface PinchResponse {
  report: SolveReport;
  dt_min: number;
  qh_min: number;
  qc_min: number;
  pinch_T_hot: number | null;
  pinch_T_cold: number | null;
  pinch_T_shifted: number | null;
  current_hot_utility: number;
  current_cold_utility: number;
  heat_recovery_potential: number;
  hot_composite: [number, number][];     // (H W, T K)
  cold_composite: [number, number][];
  streams: { id: string; T_in: number; T_out: number; Q: number; kind: string }[];
}

// ---- balance diagnostics (mirrors caldyr.solver.balance_report) ----
export interface BalanceRow {
  unit_id?: string;
  mass_in_kg_s: number;
  mass_out_kg_s: number;
  mass_rel: number;
  duty_W: number;
  energy_rel: number;
}
export interface BalanceResult {
  overall: BalanceRow;
  units: BalanceRow[];
  warnings: string[];
}

// The `.flow` document (engine io schema). Loosely typed on purpose.
export type FlowDoc = Record<string, unknown>;
