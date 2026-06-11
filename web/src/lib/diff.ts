// Diff two .flow documents (current canvas vs an AI-proposed edit) into a
// human-reviewable change list, and merge canvas positions into a proposed
// doc so accepted changes don't scramble the layout.
import type { FlowDoc } from "../types";

interface FlowUnit { id: string; type: string; params: Record<string, unknown>; xy?: number[] }
interface FlowStream { id: string; from: string | null; to: string | null }

const unitsOf = (d: FlowDoc): FlowUnit[] => (d.units as FlowUnit[] | undefined) ?? [];
const streamsOf = (d: FlowDoc): FlowStream[] => (d.streams as FlowStream[] | undefined) ?? [];
const compsOf = (d: FlowDoc): string[] =>
  (((d.components as { id: string }[] | undefined) ?? []).map((c) => c.id));

export interface FlowDiff {
  addedUnits: { id: string; type: string }[];
  removedUnits: string[];
  changedParams: { unit: string; key: string; from: unknown; to: unknown }[];
  addedStreams: string[];
  removedStreams: string[];
  rewiredStreams: string[];
  addedComponents: string[];
  removedComponents: string[];
  packageChanged: { from: string; to: string } | null;
  logicalChanged: boolean;
  count: number;
}

export function diffFlows(current: FlowDoc, proposed: FlowDoc): FlowDiff {
  const curUnits = new Map(unitsOf(current).map((u) => [u.id, u]));
  const newUnits = new Map(unitsOf(proposed).map((u) => [u.id, u]));
  const curStreams = new Map(streamsOf(current).map((s) => [s.id, s]));
  const newStreams = new Map(streamsOf(proposed).map((s) => [s.id, s]));

  const addedUnits = [...newUnits.values()]
    .filter((u) => !curUnits.has(u.id))
    .map((u) => ({ id: u.id, type: u.type }));
  const removedUnits = [...curUnits.keys()].filter((id) => !newUnits.has(id));

  const changedParams: FlowDiff["changedParams"] = [];
  for (const [id, nu] of newUnits) {
    const cu = curUnits.get(id);
    if (!cu) continue;
    const keys = new Set([...Object.keys(cu.params ?? {}), ...Object.keys(nu.params ?? {})]);
    for (const k of keys) {
      const a = cu.params?.[k];
      const b = nu.params?.[k];
      if (JSON.stringify(a) !== JSON.stringify(b)) {
        changedParams.push({ unit: id, key: k, from: a, to: b });
      }
    }
  }

  const addedStreams = [...newStreams.keys()].filter((id) => !curStreams.has(id));
  const removedStreams = [...curStreams.keys()].filter((id) => !newStreams.has(id));
  const rewiredStreams = [...newStreams.values()]
    .filter((s) => {
      const c = curStreams.get(s.id);
      return c && (c.from !== s.from || c.to !== s.to);
    })
    .map((s) => s.id);

  const curComps = compsOf(current);
  const newComps = compsOf(proposed);
  const addedComponents = newComps.filter((c) => !curComps.includes(c));
  const removedComponents = curComps.filter((c) => !newComps.includes(c));

  const pkgFrom = (current.property_package as string) ?? "";
  const pkgTo = (proposed.property_package as string) ?? "";
  const packageChanged = pkgFrom !== pkgTo ? { from: pkgFrom, to: pkgTo } : null;

  const logicalChanged =
    JSON.stringify(current.logical ?? []) !== JSON.stringify(proposed.logical ?? []);

  const count =
    addedUnits.length + removedUnits.length + changedParams.length
    + addedStreams.length + removedStreams.length + rewiredStreams.length
    + addedComponents.length + removedComponents.length
    + (packageChanged ? 1 : 0) + (logicalChanged ? 1 : 0);

  return {
    addedUnits, removedUnits, changedParams, addedStreams, removedStreams,
    rewiredStreams, addedComponents, removedComponents, packageChanged,
    logicalChanged, count,
  };
}

/** Carry canvas unit positions (and UI meta) into a proposed doc; stagger any
 * brand-new units below the existing bounding box so they're visible. */
export function mergePositions(current: FlowDoc, proposed: FlowDoc): FlowDoc {
  const out = structuredClone(proposed) as FlowDoc & { units: FlowUnit[] };
  const curXy = new Map(unitsOf(current).map((u) => [u.id, u.xy]));
  const known = unitsOf(current).filter((u) => u.xy);
  const maxY = known.length ? Math.max(...known.map((u) => u.xy![1])) : 60;
  const minX = known.length ? Math.min(...known.map((u) => u.xy![0])) : 120;

  let n = 0;
  for (const u of out.units ?? []) {
    const xy = curXy.get(u.id);
    if (xy) {
      u.xy = xy;
    } else if (!u.xy) {
      u.xy = [minX + n * 190, maxY + 140];
      n++;
    }
  }
  // keep the canvas's UI meta (theme-independent: product, colors, pins, boundaries)
  const curMeta = (current.meta as Record<string, unknown> | undefined)?.ui;
  if (curMeta) {
    const meta = ((out.meta as Record<string, unknown>) ??= {});
    meta.ui = { ...(curMeta as object), ...((meta.ui as object) ?? {}) };
  }
  return out;
}
