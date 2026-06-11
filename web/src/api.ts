// Thin client over the engine API. All physics lives server-side.
import type {
  BalanceResult, CostResponse, EnvelopeResponse, FlowDoc, OptimizeRequest,
  OptimizeResponse, Port, PropertyPackage, SolveResponse, UnitType,
} from "./types";

const BASE = "/api";

async function post<T>(path: string, body: unknown): Promise<T> {
  const res = await fetch(`${BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const detail = await res.json().catch(() => ({}));
    throw new Error(detail.detail ?? `${path} failed (${res.status})`);
  }
  return res.json() as Promise<T>;
}

async function get<T>(path: string): Promise<T> {
  const res = await fetch(`${BASE}${path}`);
  if (!res.ok) throw new Error(`${path} failed (${res.status})`);
  return res.json() as Promise<T>;
}

export const api = {
  unitTypes: () => get<UnitType[]>("/unit-types"),
  propertyPackages: () => get<PropertyPackage[]>("/property-packages"),
  components: () =>
    get<{ id: string; name: string; formula: string; cas: string }[]>("/components"),
  ports: (type: string, params: Record<string, unknown>) =>
    post<Port[]>("/ports", { type, params }),
  solve: (flow: FlowDoc, backend: string) =>
    post<SolveResponse>("/solve", { flow, backend }),
  cost: (flow: FlowDoc, product_component: string, monteCarlo = 0) =>
    post<CostResponse>("/cost", {
      flow, config: { product_component }, monte_carlo: monteCarlo,
    }),
  optimize: (req: OptimizeRequest) => post<OptimizeResponse>("/optimize", req),
  envelope: (flow: FlowDoc, stream: string, n = 30) =>
    post<EnvelopeResponse>("/envelope", { flow, stream, n }),
  balance: (flow: FlowDoc, backend: string) =>
    post<{ balance: BalanceResult }>("/balance", { flow, backend }),
  aiTool: (name: string, flow: FlowDoc, args: Record<string, unknown> = {}) =>
    post<Record<string, unknown>>("/ai/tool", { name, flow, args }),
};
