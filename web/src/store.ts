// Single source of truth for the app: canvas graph, flowsheet config, results,
// undo/redo history, clipboard, theme, and async solve/cost actions. Components
// stay thin; everything stateful funnels through here.
import {
  addEdge, applyEdgeChanges, applyNodeChanges,
  type Connection, type Edge, type EdgeChange, type NodeChange,
  type NodePositionChange,
} from "@xyflow/react";
import { create } from "zustand";
import { api } from "./api";
import { canvasToFlow, flowToCanvas, type CaldyrNode } from "./flow";
import { diffFlows, mergePositions, type FlowDiff } from "./lib/diff";
import { applyHelperLines, type HelperLines } from "./lib/helperLines";
import { autoLayout } from "./lib/layout";
import { ID_RE } from "./lib/params";
import { ChatSocket, requestOverWs } from "./lib/ws";
import {
  autosave, clearAutosave, deleteProject, downloadFlow, listProjects,
  loadAutosave, loadPref, loadTheme, pickFlowFile, savePref, saveProject,
  saveTheme, type SavedProject, type Theme,
} from "./lib/persist";
import { UNIT_SETS, type UnitSet } from "./lib/units";
import type {
  BalanceResult, CostConfigOverrides, CostResponse, FlowDoc, PropertyPackage,
  SolveResponse, UnitType,
} from "./types";

export type Tab = "params" | "streams" | "economics" | "optimize" | "study" | "calc" | "tools";
export type Busy = "solve" | "cost" | null;
export type ColorMode = "none" | "phase" | "temperature";
export type ViewMode = "bfd" | "pfd" | "pid";

export interface Group {
  id: string;
  label: string;
  members: string[];
  collapsed: boolean;
  xy: [number, number];
}
export interface Selection { kind: "node" | "edge"; id: string }
export interface Toast { id: number; kind: "info" | "success" | "error"; msg: string }

export interface ChatMsg {
  role: "user" | "assistant" | "event" | "error";
  text: string;
}


interface Snapshot {
  nodes: CaldyrNode[];
  edges: Edge[];
  components: string[];
  propertyPackage: string;
  product: string;
}
interface HistoryEntry { snap: Snapshot; tag?: string }

interface Clipboard { nodes: CaldyrNode[]; edges: Edge[] }

const HISTORY_CAP = 100;
const FEED_PORTS = [{ name: "out", direction: "outlet", kind: "material" } as const];
const PRODUCT_PORTS = [{ name: "in", direction: "inlet", kind: "material" } as const];
const DEFAULT_FEED_PARAMS = { T: 298.15, P: 101325, molar_flow: 1.0, z: {} };

let toastId = 0;
const chatSocket = new ChatSocket();
const portRefreshAt = new Map<string, number>();

interface State {
  // catalog (fetched once)
  unitTypes: UnitType[];
  packages: PropertyPackage[];
  componentCatalog: { id: string; name: string; formula: string; cas: string }[];
  engineReady: boolean;
  // flowsheet
  nodes: CaldyrNode[];
  edges: Edge[];
  components: string[];
  propertyPackage: string;
  product: string;
  backend: string;
  // results
  solveRes: SolveResponse | null;
  costRes: CostResponse | null;
  resultsStale: boolean;
  // ui
  selected: Selection | null;
  tab: Tab;
  status: string;
  busy: Busy;
  toasts: Toast[];
  theme: Theme;
  colorMode: ColorMode;
  pinnedStreams: string[];
  helperLines: HelperLines;
  fitNonce: number;
  unitSet: UnitSet;
  // per-field display-unit overrides, keyed by `${nodeId}:${param}` (a display
  // preference; the stored value stays SI). Falls back to the unitSet default.
  unitOverrides: Record<string, string>;
  // per-flowsheet techno-economic assumption overrides (prices, factors, sizing,
  // rates) sent to /cost; rides in meta.ui. Empty = engine defaults.
  costConfig: CostConfigOverrides;
  inspectorWidth: number;
  viewMode: ViewMode;
  groups: Group[];
  calcs: { name: string; expr: string }[];
  // flowsheet-level logical ops + solver hints (round-tripped via .flow)
  logical: Record<string, unknown>[];
  solverHints: Record<string, unknown>;
  // chat / AI
  chatOpen: boolean;
  chatMessages: ChatMsg[];
  chatBusy: boolean;
  pendingFlow: FlowDoc | null;
  pendingDiff: FlowDiff | null;
  // balance diagnostics
  balance: BalanceResult | null;
  balanceBusy: boolean;
  // projects dialog
  projectsOpen: boolean;
  projects: SavedProject[];
  // history & clipboard
  past: HistoryEntry[];
  future: HistoryEntry[];
  clipboard: Clipboard | null;
  dragging: boolean;

  // actions
  init: () => Promise<void>;
  // NodeChange over the base Node type: the canvas mixes CaldyrNodes with
  // synthesized group blocks; non-store ids are routed to group handling.
  onNodesChange: (changes: NodeChange[]) => void;
  onEdgesChange: (changes: EdgeChange[]) => void;
  onConnect: (c: Connection) => void;
  setSelection: (sel: Selection | null) => void;
  setTab: (t: Tab) => void;
  setBackend: (b: string) => void;
  setPropertyPackage: (p: string) => void;
  setProduct: (p: string) => void;
  addComponent: (id: string) => void;
  removeComponent: (id: string) => void;
  addUnit: (kind: "unit" | "feed" | "product", unitType?: UnitType) => void;
  setParam: (nodeId: string, key: string, value: unknown) => void;
  unsetParam: (nodeId: string, key: string) => void;
  refreshPorts: (nodeId: string) => void;
  renameNode: (id: string, next: string) => boolean;
  renameEdge: (id: string, next: string) => boolean;
  undo: () => void;
  redo: () => void;
  copySelection: () => void;
  paste: () => void;
  duplicateSelection: () => void;
  newFlowsheet: () => void;
  loadFlowDoc: (doc: FlowDoc, name: string) => void;
  saveFile: () => void;
  openFile: () => Promise<void>;
  toFlowDoc: () => FlowDoc;
  autosaveNow: () => void;
  solve: () => Promise<void>;
  cost: (monteCarlo?: number) => Promise<void>;
  toast: (kind: Toast["kind"], msg: string) => void;
  dismissToast: (id: number) => void;
  toggleTheme: () => void;
  setColorMode: (m: ColorMode) => void;
  togglePin: (edgeId: string) => void;
  selectAll: () => void;
  runAutoLayout: () => Promise<void>;
  setUnitSet: (u: UnitSet) => void;
  setUnitOverride: (key: string, unit: string | null) => void;
  setCostConfig: (c: CostConfigOverrides) => void;
  setInspectorWidth: (w: number) => void;
  setViewMode: (m: ViewMode) => void;
  groupSelection: () => void;
  ungroup: (groupId: string) => void;
  toggleGroupCollapse: (groupId: string) => void;
  moveGroup: (groupId: string, dx: number, dy: number) => void;
  setCalcs: (rows: { name: string; expr: string }[]) => void;
  setLogical: (ops: Record<string, unknown>[]) => void;
  toggleChat: () => void;
  sendChat: (text: string) => Promise<void>;
  quickAction: (kind: "describe" | "diagnose") => Promise<void>;
  acceptPending: () => void;
  rejectPending: () => void;
  runBalance: () => Promise<void>;
  toggleProjects: () => void;
  saveCurrentProject: (name: string) => void;
  loadProject: (name: string) => void;
  removeProject: (name: string) => void;
}

const takeSnapshot = (s: State): Snapshot => ({
  nodes: s.nodes,
  edges: s.edges,
  components: s.components,
  propertyPackage: s.propertyPackage,
  product: s.product,
});

export const useStore = create<State>((set, get) => {
  /** Push the current state onto the undo stack. Entries with the same `tag`
   * as the top of the stack coalesce (so typing in a field is one undo step). */
  const commit = (tag?: string) => {
    const s = get();
    if (tag && s.past.length && s.past[s.past.length - 1].tag === tag) return;
    const past = [...s.past, { snap: takeSnapshot(s), tag }].slice(-HISTORY_CAP);
    set({ past, future: [] });
  };

  const markStale = () => {
    if (get().solveRes || get().costRes) set({ resultsStale: true });
  };

  const uniqueNodeId = (base: string): string => {
    const ids = new Set(get().nodes.map((n) => n.id));
    if (!ids.has(base)) return base;
    let i = 2;
    while (ids.has(`${base}_${i}`)) i++;
    return `${base}_${i}`;
  };

  const nextStreamId = (): string => {
    const ids = new Set(get().edges.map((e) => e.id));
    let i = 1;
    while (ids.has(`S${i}`)) i++;
    return `S${i}`;
  };

  const applySnapshot = (snap: Snapshot) =>
    set({ ...snap, selected: null, resultsStale: true });

  return {
    unitTypes: [],
    packages: [],
    componentCatalog: [],
    engineReady: false,
    nodes: [],
    edges: [],
    components: [],
    propertyPackage: "thermo:PR",
    product: "",
    backend: "sequential",
    solveRes: null,
    costRes: null,
    resultsStale: false,
    selected: null,
    tab: "params",
    status: "Loading engine...",
    busy: null,
    toasts: [],
    theme: loadTheme(),
    colorMode: "none",
    pinnedStreams: [],
    helperLines: {},
    fitNonce: 0,
    unitSet: (UNIT_SETS as string[]).includes(loadPref("units") ?? "")
      ? (loadPref("units") as UnitSet) : "SI",
    unitOverrides: {},  // per-flowsheet; restored from meta.ui on load
    costConfig: {},     // per-flowsheet; restored from meta.ui on load
    inspectorWidth: Math.min(640, Math.max(300, Number(loadPref("panelw")) || 360)),
    viewMode: "pfd",
    groups: [],
    calcs: [],
    logical: [],
    solverHints: {},
    chatOpen: false,
    chatMessages: [],
    chatBusy: false,
    pendingFlow: null,
    pendingDiff: null,
    balance: null,
    balanceBusy: false,
    projectsOpen: false,
    projects: listProjects(),
    past: [],
    future: [],
    clipboard: null,
    dragging: false,

    init: async () => {
      try {
        const [types, pkgs, catalog] = await Promise.all([
          api.unitTypes(), api.propertyPackages(),
          api.components().catch(() => []),  // older API: autocomplete just stays empty
        ]);
        set({ unitTypes: types, packages: pkgs, componentCatalog: catalog, engineReady: true });
        const saved = loadAutosave();
        if (saved) {
          get().loadFlowDoc(saved, "autosaved session");
          set({ past: [] });
        } else {
          set({ status: "Ready — add units from the palette or load an example" });
        }
      } catch (e) {
        set({ status: `API unreachable — is the engine running? (${(e as Error).message})` });
      }
    },

    onNodesChange: (changes) => {
      const s = get();
      // Synthesized group-block nodes: route their drags to moveGroup, their
      // deletion to ungroup; they never live in store.nodes.
      const groupIds = new Set(s.groups.map((g) => g.id));
      const groupChanges = changes.filter((c) =>
        "id" in c && typeof c.id === "string" && groupIds.has(c.id));
      for (const c of groupChanges) {
        if (c.type === "position" && c.position) {
          const g = s.groups.find((x) => x.id === c.id)!;
          get().moveGroup(g.id, c.position.x - g.xy[0], c.position.y - g.xy[1]);
        } else if (c.type === "remove") {
          get().ungroup(c.id);
        }
      }
      if (changes.some((c) => c.type === "remove")) {
        commit();
        markStale();
        const removed = new Set(changes.filter((c) => c.type === "remove").map((c) => c.id));
        if (s.selected?.kind === "node" && removed.has(s.selected.id)) set({ selected: null });
        // drop deleted nodes from any group membership
        if (get().groups.some((g) => g.members.some((m) => removed.has(m)))) {
          set({
            groups: get().groups
              .map((g) => ({ ...g, members: g.members.filter((m) => !removed.has(m)) }))
              .filter((g) => g.members.length > 1),
          });
        }
      }
      const dragStart = changes.some((c) => c.type === "position" && c.dragging);
      if (dragStart && !s.dragging) {
        commit();
        set({ dragging: true });
      }
      if (changes.some((c) => c.type === "position" && c.dragging === false)) {
        set({ dragging: false, helperLines: {} });
      }
      // Alignment guides: when exactly one node is being dragged, snap it to
      // neighbours' edges/centers and surface the guide lines.
      const moving = changes.filter(
        (c): c is NodePositionChange => c.type === "position" && c.dragging === true);
      if (moving.length === 1) {
        set({ helperLines: applyHelperLines(moving[0], get().nodes) });
      } else if (moving.length > 1) {
        set({ helperLines: {} });
      }
      set({ nodes: applyNodeChanges(changes as NodeChange<CaldyrNode>[], get().nodes) });
    },

    onEdgesChange: (changes) => {
      const s = get();
      if (changes.some((c) => c.type === "remove")) {
        commit();
        markStale();
        const removed = new Set(changes.filter((c) => c.type === "remove").map((c) => c.id));
        if (s.selected?.kind === "edge" && removed.has(s.selected.id)) set({ selected: null });
      }
      set({ edges: applyEdgeChanges(changes, get().edges) });
    },

    onConnect: (c) => {
      const s = get();
      const src = s.nodes.find((n) => n.id === c.source);
      const tgt = s.nodes.find((n) => n.id === c.target);
      if (!src || !tgt) return;
      const srcPort = src.data.ports.find((p) => p.name === (c.sourceHandle ?? ""));
      const tgtPort = tgt.data.ports.find((p) => p.name === (c.targetHandle ?? ""));
      if (srcPort && tgtPort && srcPort.kind !== tgtPort.kind) {
        get().toast("error", `Cannot connect ${srcPort.kind} port to ${tgtPort.kind} port.`);
        return;
      }
      const srcTaken = s.edges.some(
        (e) => e.source === c.source && (e.sourceHandle ?? null) === (c.sourceHandle ?? null));
      const tgtTaken = s.edges.some(
        (e) => e.target === c.target && (e.targetHandle ?? null) === (c.targetHandle ?? null));
      if (srcTaken || tgtTaken) {
        get().toast("error", "Port already connected — delete the existing stream first.");
        return;
      }
      commit();
      markStale();
      const id = nextStreamId();
      set({ edges: addEdge({ ...c, id, label: id }, get().edges) });
    },

    setSelection: (sel) => set({ selected: sel }),
    setTab: (t) => set({ tab: t }),
    setBackend: (b) => set({ backend: b }),
    setPropertyPackage: (p) => {
      commit();
      markStale();
      set({ propertyPackage: p });
    },
    setProduct: (p) => set({ product: p }),

    addComponent: (id) => {
      const c = id.trim().toLowerCase();
      if (!c) return;
      if (get().components.includes(c)) {
        get().toast("info", `"${c}" is already in the component list.`);
        return;
      }
      commit();
      markStale();
      set({ components: [...get().components, c] });
    },

    removeComponent: (id) => {
      commit();
      markStale();
      set({ components: get().components.filter((c) => c !== id) });
    },

    addUnit: (kind, unitType) => {
      commit();
      const s = get();
      const n = s.nodes.length;
      const position = { x: 140 + (n % 8) * 40, y: 80 + (n % 8) * 30 };
      let node: CaldyrNode;
      if (kind === "unit" && unitType) {
        const id = uniqueNodeId(`${unitType.type.slice(0, 4).toUpperCase()}${n + 1}`);
        node = {
          id, type: "caldyr", position,
          data: { kind: "unit", label: id, unitType: unitType.type, ports: unitType.ports, params: {} },
        };
      } else if (kind === "feed") {
        const id = uniqueNodeId(`FEED${n + 1}`);
        node = {
          id, type: "caldyr", position,
          data: { kind: "feed", label: id, ports: [...FEED_PORTS], params: { ...DEFAULT_FEED_PARAMS } },
        };
      } else {
        const id = uniqueNodeId(`PROD${n + 1}`);
        node = {
          id, type: "caldyr", position,
          data: { kind: "product", label: id, ports: [...PRODUCT_PORTS], params: {} },
        };
      }
      set({ nodes: [...s.nodes, node], selected: { kind: "node", id: node.id }, tab: "params" });
    },

    setParam: (nodeId, key, value) => {
      commit(`param:${nodeId}:${key}`);
      markStale();
      set({
        nodes: get().nodes.map((n) =>
          n.id === nodeId
            ? { ...n, data: { ...n.data, params: { ...n.data.params, [key]: value } } }
            : n),
      });
      get().refreshPorts(nodeId);
    },

    unsetParam: (nodeId, key) => {
      commit(`unset:${nodeId}:${key}`);
      markStale();
      set({
        nodes: get().nodes.map((n) => {
          if (n.id !== nodeId) return n;
          const params = { ...n.data.params };
          delete params[key];
          return { ...n, data: { ...n.data, params } };
        }),
      });
      get().refreshPorts(nodeId);
    },

    // Some unit types derive their ports from params (multi-feed columns,
    // Balance n_inlets, ...). Re-derive from the engine after param edits.
    refreshPorts: (nodeId) => {
      const node = get().nodes.find((n) => n.id === nodeId);
      if (!node || node.data.kind !== "unit" || !node.data.unitType) return;
      const requestedAt = Date.now();
      portRefreshAt.set(nodeId, requestedAt);
      void api.ports(node.data.unitType, node.data.params)
        .then((ports) => {
          if (portRefreshAt.get(nodeId) !== requestedAt) return; // stale reply
          const cur = get().nodes.find((n) => n.id === nodeId);
          if (!cur) return;
          if (JSON.stringify(cur.data.ports) === JSON.stringify(ports)) return;
          set({
            nodes: get().nodes.map((n) =>
              n.id === nodeId ? { ...n, data: { ...n.data, ports } } : n),
          });
        })
        .catch(() => { /* params mid-edit can be invalid; keep old ports */ });
    },

    renameNode: (id, next) => {
      const name = next.trim();
      if (name === id) return true;
      if (!ID_RE.test(name)) {
        get().toast("error", "Names may only use letters, digits, _ and -.");
        return false;
      }
      if (get().nodes.some((n) => n.id === name)) {
        get().toast("error", `A node named "${name}" already exists.`);
        return false;
      }
      commit();
      markStale();
      set({
        nodes: get().nodes.map((n) =>
          n.id === id ? { ...n, id: name, data: { ...n.data, label: name } } : n),
        edges: get().edges.map((e) => ({
          ...e,
          source: e.source === id ? name : e.source,
          target: e.target === id ? name : e.target,
        })),
        selected: get().selected?.id === id ? { kind: "node", id: name } : get().selected,
      });
      return true;
    },

    renameEdge: (id, next) => {
      const name = next.trim();
      if (name === id) return true;
      if (!ID_RE.test(name)) {
        get().toast("error", "Names may only use letters, digits, _ and -.");
        return false;
      }
      if (get().edges.some((e) => e.id === name)) {
        get().toast("error", `A stream named "${name}" already exists.`);
        return false;
      }
      commit();
      markStale();
      set({
        edges: get().edges.map((e) => (e.id === id ? { ...e, id: name, label: name } : e)),
        selected: get().selected?.id === id ? { kind: "edge", id: name } : get().selected,
      });
      return true;
    },

    undo: () => {
      const s = get();
      if (!s.past.length) return;
      const entry = s.past[s.past.length - 1];
      set({
        past: s.past.slice(0, -1),
        future: [...s.future, { snap: takeSnapshot(s) }],
      });
      applySnapshot(entry.snap);
    },

    redo: () => {
      const s = get();
      if (!s.future.length) return;
      const entry = s.future[s.future.length - 1];
      set({
        future: s.future.slice(0, -1),
        past: [...s.past, { snap: takeSnapshot(s) }],
      });
      applySnapshot(entry.snap);
    },

    copySelection: () => {
      const s = get();
      let picked = s.nodes.filter((n) => n.selected);
      if (!picked.length && s.selected?.kind === "node") {
        picked = s.nodes.filter((n) => n.id === s.selected!.id);
      }
      if (!picked.length) return;
      const ids = new Set(picked.map((n) => n.id));
      const inner = s.edges.filter((e) => ids.has(e.source) && ids.has(e.target));
      set({
        clipboard: {
          nodes: structuredClone(picked),
          edges: structuredClone(inner),
        },
      });
      set({ status: `Copied ${picked.length} node${picked.length > 1 ? "s" : ""}` });
    },

    paste: () => {
      const s = get();
      if (!s.clipboard) return;
      commit();
      const idMap = new Map<string, string>();
      const newNodes = s.clipboard.nodes.map((n) => {
        const id = uniqueNodeId(n.id);
        idMap.set(n.id, id);
        return {
          ...structuredClone(n),
          id,
          selected: true,
          position: { x: n.position.x + 48, y: n.position.y + 48 },
          data: { ...structuredClone(n.data), label: id },
        };
      });
      const existingEdges = new Set(get().edges.map((e) => e.id));
      const newEdges = s.clipboard.edges.map((e, i) => {
        let id = `S${existingEdges.size + i + 1}`;
        while (existingEdges.has(id)) id = `${id}c`;
        existingEdges.add(id);
        return {
          ...structuredClone(e),
          id,
          label: id,
          source: idMap.get(e.source)!,
          target: idMap.get(e.target)!,
          selected: false,
        };
      });
      set({
        nodes: [...get().nodes.map((n) => ({ ...n, selected: false })), ...newNodes],
        edges: [...get().edges, ...newEdges],
        status: `Pasted ${newNodes.length} node${newNodes.length > 1 ? "s" : ""}`,
      });
    },

    duplicateSelection: () => {
      get().copySelection();
      get().paste();
    },

    newFlowsheet: () => {
      commit();
      clearAutosave();
      set({
        nodes: [], edges: [], components: [], product: "",
        solveRes: null, costRes: null, resultsStale: false,
        selected: null, tab: "params", unitOverrides: {}, costConfig: {},
        status: "New flowsheet — add components in the panel, then units from the palette",
      });
    },

    loadFlowDoc: (doc, name) => {
      const s = get();
      if (!s.unitTypes.length) {
        get().toast("error", "Engine catalog not loaded yet — try again in a moment.");
        return;
      }
      commit();
      try {
        const c = flowToCanvas(doc, s.unitTypes);
        set({
          nodes: c.nodes,
          edges: c.edges,
          components: c.components,
          propertyPackage: c.propertyPackage,
          product: c.ui.product ?? get().product,
          backend: c.ui.backend ?? get().backend,
          colorMode: (c.ui.color_mode as ColorMode) ?? get().colorMode,
          pinnedStreams: c.ui.pinned_streams ?? [],
          viewMode: (c.ui.view_mode as ViewMode) ?? "pfd",
          unitOverrides: c.ui.unit_overrides ?? {},
          costConfig: (c.ui.cost_config as CostConfigOverrides) ?? {},
          calcs: c.ui.calcs ?? [],
          groups: (c.ui.groups as Group[]) ?? [],
          logical: c.extras.logical ?? [],
          solverHints: c.extras.solver_hints ?? {},
          solveRes: null, costRes: null, resultsStale: false,
          balance: null,
          selected: null,
          status: `Loaded ${name} — ${c.nodes.length} nodes`,
        });
      } catch (e) {
        get().toast("error", `Could not load flowsheet: ${(e as Error).message}`);
      }
    },

    toFlowDoc: () => {
      const s = get();
      return canvasToFlow(s.nodes, s.edges, s.components, s.propertyPackage, {
        product: s.product,
        backend: s.backend,
        color_mode: s.colorMode,
        pinned_streams: s.pinnedStreams,
        view_mode: s.viewMode,
        unit_overrides: s.unitOverrides,
        cost_config: s.costConfig as Record<string, unknown>,
        calcs: s.calcs,
        groups: s.groups,
      }, { logical: s.logical, solver_hints: s.solverHints });
    },

    saveFile: () => {
      downloadFlow(get().toFlowDoc());
      set({ status: "Saved flowsheet.flow" });
    },

    openFile: async () => {
      try {
        const doc = await pickFlowFile();
        if (doc) get().loadFlowDoc(doc, "file");
      } catch (e) {
        get().toast("error", `Could not open file: ${(e as Error).message}`);
      }
    },

    autosaveNow: () => {
      const s = get();
      if (s.nodes.length || s.components.length) autosave(s.toFlowDoc());
    },

    solve: async () => {
      if (get().busy) return;
      set({ busy: "solve", status: "Solving..." });
      try {
        // WebSocket first (live per-iteration progress); REST as fallback.
        let res: SolveResponse;
        try {
          interface SolveEvent {
            type: string;
            iteration?: number;
            residual?: number;
            detail?: string;
            report?: SolveResponse["report"];
            streams?: SolveResponse["streams"];
            designs?: SolveResponse["designs"];
            molar_mass?: SolveResponse["molar_mass"];
          }
          const final = await requestOverWs<SolveEvent>(
            "/ws/solve",
            { flow: get().toFlowDoc(), backend: get().backend },
            (e) => {
              if (e.type === "iteration") {
                set({ status: `Solving — iteration ${e.iteration}, residual ${e.residual?.toExponential(2)}` });
              }
            },
            (e) => e.type === "result" || e.type === "error",
          );
          if (final.type === "error") throw new Error(final.detail ?? "solve failed");
          res = { report: final.report!, streams: final.streams!, designs: final.designs,
                  molar_mass: final.molar_mass };
        } catch (wsErr) {
          if ((wsErr as Error).message.includes("WebSocket")) {
            res = await api.solve(get().toFlowDoc(), get().backend);  // fallback
          } else {
            throw wsErr;
          }
        }
        set({
          solveRes: res,
          resultsStale: false,
          balance: null,
          tab: "streams",
          status: res.report.converged
            ? `Solved (${res.report.method}, ${res.report.iterations} iterations)`
            : "Did not converge",
        });
        if (!res.report.converged) {
          get().toast("error", `Solver did not converge after ${res.report.iterations} iterations.`);
        }
      } catch (e) {
        set({ status: "Solve failed" });
        get().toast("error", `Solve failed: ${(e as Error).message}`);
      } finally {
        set({ busy: null });
      }
    },

    cost: async (monteCarlo = 0) => {
      if (get().busy) return;
      if (!get().product) {
        get().toast("error", "Pick a product component (Flowsheet panel) before costing.");
        set({ tab: "params", selected: null });
        return;
      }
      set({ busy: "cost", status: monteCarlo ? `Costing + ${monteCarlo} MC samples...` : "Costing..." });
      try {
        const res = await api.cost(get().toFlowDoc(), get().product, monteCarlo,
                                   get().costConfig);
        set({
          costRes: res,
          resultsStale: false,
          tab: "economics",
          status: `Costed — LCOP $${res.profitability.lcop.toFixed(3)}/kg`,
        });
      } catch (e) {
        set({ status: "Cost failed" });
        get().toast("error", `Cost failed: ${(e as Error).message}`);
      } finally {
        set({ busy: null });
      }
    },

    toast: (kind, msg) => {
      const t = { id: ++toastId, kind, msg };
      set({ toasts: [...get().toasts, t] });
    },

    dismissToast: (id) => set({ toasts: get().toasts.filter((t) => t.id !== id) }),

    toggleTheme: () => {
      const next: Theme = get().theme === "dark" ? "light" : "dark";
      saveTheme(next);
      document.documentElement.dataset.theme = next;
      set({ theme: next });
    },

    setColorMode: (m) => set({ colorMode: m }),

    setUnitSet: (u) => {
      savePref("units", u);
      set({ unitSet: u });
    },

    setUnitOverride: (key, unit) => {
      const next = { ...get().unitOverrides };
      if (unit == null) delete next[key];
      else next[key] = unit;
      set({ unitOverrides: next }); // persisted via meta.ui on autosave
    },

    setCostConfig: (c) => set({ costConfig: c }), // persisted via meta.ui on autosave

    setInspectorWidth: (w) => {
      const clamped = Math.min(640, Math.max(300, Math.round(w)));
      savePref("panelw", String(clamped));
      set({ inspectorWidth: clamped });
    },

    togglePin: (edgeId) => {
      const pins = get().pinnedStreams;
      set({
        pinnedStreams: pins.includes(edgeId)
          ? pins.filter((p) => p !== edgeId)
          : [...pins, edgeId],
      });
    },

    selectAll: () => set({
      nodes: get().nodes.map((n) => ({ ...n, selected: true })),
      edges: get().edges.map((e) => ({ ...e, selected: true })),
    }),

    setViewMode: (m) => set({ viewMode: m }),

    setCalcs: (rows) => set({ calcs: rows }),

    groupSelection: () => {
      const s = get();
      const members = s.nodes.filter((n) => n.selected).map((n) => n.id);
      if (members.length < 2) {
        get().toast("info", "Select at least two nodes (Shift+drag) to group them.");
        return;
      }
      commit();
      const picked = s.nodes.filter((n) => members.includes(n.id));
      const cx = picked.reduce((a, n) => a + n.position.x, 0) / picked.length;
      const cy = picked.reduce((a, n) => a + n.position.y, 0) / picked.length;
      const id = uniqueNodeId(`GRP${s.groups.length + 1}`);
      set({
        groups: [...s.groups, {
          id, label: id, members, collapsed: true,
          xy: [Math.round(cx), Math.round(cy)],
        }],
        nodes: get().nodes.map((n) => ({ ...n, selected: false })),
        status: `Grouped ${members.length} nodes into ${id} — double-click it to expand`,
      });
    },

    ungroup: (groupId) => {
      commit();
      set({ groups: get().groups.filter((g) => g.id !== groupId) });
    },

    toggleGroupCollapse: (groupId) => {
      set({
        groups: get().groups.map((g) =>
          g.id === groupId ? { ...g, collapsed: !g.collapsed } : g),
      });
    },

    moveGroup: (groupId, dx, dy) => {
      const s = get();
      const g = s.groups.find((x) => x.id === groupId);
      if (!g) return;
      set({
        groups: s.groups.map((x) =>
          x.id === groupId ? { ...x, xy: [x.xy[0] + dx, x.xy[1] + dy] } : x),
        // members ride along so expanding lands near the new location
        nodes: s.nodes.map((n) =>
          g.members.includes(n.id)
            ? { ...n, position: { x: n.position.x + dx, y: n.position.y + dy } }
            : n),
      });
    },

    setLogical: (ops) => {
      commit();
      markStale();
      set({ logical: ops });
    },

    toggleChat: () => set({ chatOpen: !get().chatOpen }),

    toggleProjects: () => set({ projectsOpen: !get().projectsOpen }),

    saveCurrentProject: (name) => {
      const trimmed = name.trim();
      if (!trimmed) return;
      set({ projects: saveProject(trimmed, get().toFlowDoc()) });
      get().toast("success", `Saved project "${trimmed}".`);
    },

    loadProject: (name) => {
      const p = get().projects.find((x) => x.name === name);
      if (!p) return;
      get().loadFlowDoc(p.doc, `project "${name}"`);
      set({ projectsOpen: false });
    },

    removeProject: (name) => set({ projects: deleteProject(name) }),

    sendChat: async (text) => {
      const msg = text.trim();
      if (!msg || get().chatBusy) return;
      set({
        chatBusy: true,
        chatMessages: [...get().chatMessages, { role: "user", text: msg }],
      });
      const push = (m: ChatMsg) =>
        set({ chatMessages: [...get().chatMessages, m] });
      try {
        if (!chatSocket.connected) {
          await chatSocket.connect(
            (e) => {
              if (e.type === "text" && e.text) {
                push({ role: "assistant", text: String(e.text) });
              } else if (e.type === "tool_call") {
                push({ role: "event", text: `→ ${e.name}` });
              } else if (e.type === "tool_result" && !e.ok) {
                push({ role: "event", text: `✗ ${e.name}: ${e.summary}` });
              } else if (e.type === "done") {
                const proposed = e.flow as FlowDoc | null;
                if (proposed) {
                  const diff = diffFlows(get().toFlowDoc(), proposed);
                  if (diff.count > 0) {
                    set({
                      pendingFlow: mergePositions(get().toFlowDoc(), proposed),
                      pendingDiff: diff,
                    });
                  }
                }
                set({ chatBusy: false });
              } else if (e.type === "error") {
                push({ role: "error", text: String(e.detail) });
                set({ chatBusy: false });
              }
            },
            () => set({ chatBusy: false }),
          );
        }
        chatSocket.send({ text: msg, flow: get().toFlowDoc() });
      } catch (err) {
        push({
          role: "error",
          text: `Chat unavailable (${(err as Error).message}). Is the API running `
            + `and a local LLM (Ollama) up? Quick actions still work without one.`,
        });
        set({ chatBusy: false });
      }
    },

    quickAction: async (kind) => {
      // No-LLM path: run the diagnostic tool directly and show its summary.
      const name = kind === "describe" ? "describe_flowsheet" : "explain_convergence";
      set({
        chatOpen: true,
        chatMessages: [...get().chatMessages, {
          role: "user",
          text: kind === "describe" ? "Describe this flowsheet." : "Diagnose the last solve.",
        }],
      });
      try {
        const out = await api.aiTool(name, get().toFlowDoc());
        set({
          chatMessages: [...get().chatMessages,
            { role: "assistant", text: String(out.summary ?? "(no summary)") }],
        });
      } catch (e) {
        set({
          chatMessages: [...get().chatMessages,
            { role: "error", text: (e as Error).message }],
        });
      }
    },

    acceptPending: () => {
      const doc = get().pendingFlow;
      if (!doc) return;
      get().loadFlowDoc(doc, "AI edit");
      set({
        pendingFlow: null, pendingDiff: null,
        status: "AI edit applied — re-solve to refresh results",
      });
    },

    rejectPending: () => set({ pendingFlow: null, pendingDiff: null }),

    runBalance: async () => {
      if (get().balanceBusy) return;
      set({ balanceBusy: true });
      try {
        const out = await api.balance(get().toFlowDoc(), get().backend);
        set({ balance: out.balance, tab: "streams" });
      } catch (e) {
        get().toast("error", `Balance check failed: ${(e as Error).message}`);
      } finally {
        set({ balanceBusy: false });
      }
    },

    runAutoLayout: async () => {
      const s = get();
      if (!s.nodes.length) return;
      try {
        const pos = await autoLayout(s.nodes, s.edges);
        commit();
        set({
          nodes: get().nodes.map((n) => {
            const p = pos.get(n.id);
            return p ? { ...n, position: p } : n;
          }),
          fitNonce: get().fitNonce + 1,
          status: "Auto-arranged flowsheet",
        });
      } catch (e) {
        get().toast("error", `Auto-layout failed: ${(e as Error).message}`);
      }
    },
  };
});
