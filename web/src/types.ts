// Shared types mirroring the engine API (see api/models.py).

export interface Port {
  name: string;
  direction: "inlet" | "outlet";
  kind: "material" | "energy";
}

export interface UnitType {
  type: string;
  doc: string;
  ports: Port[];
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
