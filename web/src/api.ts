// Thin client over the engine API. All physics lives server-side.
import type {
  BalanceResult, CostConfigOverrides, CostDefaults, CostResponse, EnvelopeResponse,
  FlowDoc, OptimizeRequest, OptimizeResponse, PinchResponse, Port, PriceCatalog,
  PropertyPackage, PropertyTableResponse, ReliefResponse, SolveResponse, UnitType,
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
  cost: (flow: FlowDoc, product_component: string, monteCarlo = 0,
         overrides?: CostConfigOverrides) =>
    post<CostResponse>("/cost", {
      flow, config: { product_component, ...(overrides ?? {}) }, monte_carlo: monteCarlo,
    }),
  prices: () => get<PriceCatalog>("/prices"),
  costDefaults: () => get<CostDefaults>("/cost-defaults"),
  optimize: (req: OptimizeRequest) => post<OptimizeResponse>("/optimize", req),
  envelope: (flow: FlowDoc, stream: string, n = 30) =>
    post<EnvelopeResponse>("/envelope", { flow, stream, n }),
  balance: (flow: FlowDoc, backend: string) =>
    post<{ balance: BalanceResult }>("/balance", { flow, backend }),
  propertyTable: (
    flow: FlowDoc, stream: string, T: number[], P: number[], props: string[],
  ) => post<PropertyTableResponse>("/property-table", { flow, stream, T, P, props }),
  relief: (body: Record<string, unknown>) => post<ReliefResponse>("/relief", body),
  pinch: (flow: FlowDoc, backend: string, dt_min: number) =>
    post<PinchResponse>("/pinch", { flow, backend, dt_min }),
  aiTool: (name: string, flow: FlowDoc, args: Record<string, unknown> = {}) =>
    post<Record<string, unknown>>("/ai/tool", { name, flow, args }),
};
