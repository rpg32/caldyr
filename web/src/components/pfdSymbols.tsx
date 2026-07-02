// Equipment-scale PFD symbols (Turton / ISO 10628 flavour): the symbol *is* the
// node, drawn at a meaningful size with physically-placed nozzles. Strokes use
// currentColor so they recolor with theme/selection; closed vessel bodies take a
// subtle themed fill (.pfd-body) so they read as solid equipment over the grid.
//
// Two exports drive the PFD renderer:
//   symbolFor(kind, unitType)   -> { w, h, body } natural-size SVG
//   portAnchors(kind, type, ports) -> per-port side + fractional offset
// Anything without a custom entry falls back to evenly-spaced left/right nozzles
// (the classic card behaviour), so new engine unit types still render sanely.
import { Position } from "@xyflow/react";
import type { ReactNode } from "react";
import type { Port } from "../types";

const S = {
  fill: "none",
  stroke: "currentColor",
  strokeWidth: 1.8,
  strokeLinecap: "round" as const,
  strokeLinejoin: "round" as const,
};

export interface PfdSymbol {
  w: number;
  h: number;
  body: ReactNode;
}

const sym = (w: number, h: number, body: ReactNode): PfdSymbol => ({ w, h, body });
const svg = (w: number, h: number, children: ReactNode): ReactNode => (
  <svg width={w} height={h} viewBox={`0 0 ${w} ${h}`} aria-hidden {...S}>{children}</svg>
);

// ---- shared builders -------------------------------------------------------

/** A vertical tower shell with optional trays (solid lines) or packing (dashed),
 *  plus optional integral condenser/reboiler circles for rigorous columns. */
function column(opts: { trays?: number; packed?: boolean; condenser?: boolean; reboiler?: boolean }): PfdSymbol {
  const w = 64, h = 130;
  const x = 17, ww = 30, top = 10, bot = 120;
  const lines: ReactNode[] = [];
  const n = opts.trays ?? 0;
  for (let i = 1; i <= n; i++) {
    const y = top + ((bot - top) * i) / (n + 1);
    lines.push(<path key={`t${i}`} d={`M${x + 2} ${y} h${ww - 4}`}
      strokeDasharray={opts.packed ? "3 2.5" : undefined} />);
  }
  return sym(w, h, svg(w, h, <>
    <rect className="pfd-body" x={x} y={top} width={ww} height={bot - top} rx={8} />
    {lines}
    {opts.condenser && <circle cx={w - 6} cy={top + 6} r={5} />}
    {opts.reboiler && <circle cx={6} cy={bot - 6} r={5} />}
  </>));
}

// ---- per-type symbols ------------------------------------------------------

const SYMBOLS: Record<string, PfdSymbol> = {
  Mixer: sym(64, 64, svg(64, 64, <>
    <circle className="pfd-body" cx={32} cy={32} r={24} />
    <path d="M20 24 L34 32 L20 40 M34 32 h10" />
  </>)),
  Splitter: sym(64, 64, svg(64, 64, <>
    <circle className="pfd-body" cx={32} cy={32} r={24} />
    <path d="M14 32 h10 L38 24 M24 32 L38 40" />
  </>)),
  Heater: sym(64, 64, svg(64, 64, <>
    <circle className="pfd-body" cx={32} cy={32} r={26} />
    <path d="M14 32 h8 l4 -9 l8 18 l4 -9 h8" />
  </>)),
  FiredHeater: sym(90, 100, svg(90, 100, <>
    <path className="pfd-body" d="M10 92 L10 40 L45 10 L80 40 L80 92 Z" />
    <path d="M45 78 c-8 0 -11 -8 -5.5 -13 c1 3 2.8 4.2 4.6 4.6 c-1.9 -6 1 -11 6 -14
             c-1 6 2.2 8 4 11.5 c3 6 -1.5 11 -9 11 Z" />
    <path d="M10 84 h70" />
  </>)),
  AirCooler: sym(72, 60, svg(72, 60, <>
    <rect className="pfd-body" x={8} y={20} width={56} height={26} rx={4} />
    <path d="M24 33 l9 -6 M24 33 l9 6 M48 33 l-9 -6 M48 33 l-9 6" />
    <path d="M18 12 L24 20 M36 10 L36 20 M54 12 L48 20" />
  </>)),
  HeatExchanger: sym(72, 72, svg(72, 72, <>
    <circle className="pfd-body" cx={36} cy={36} r={30} />
    <path d="M6 36 h10 c8 0 8 -14 18 -14 s10 28 18 28 c5 0 6 -7 10 -7" />
  </>)),
  Flash: sym(56, 96, svg(56, 96, <>
    <rect className="pfd-body" x={12} y={8} width={32} height={80} rx={16} />
    <path d="M15 56 h26" strokeDasharray="3 2.5" />
  </>)),
  Evaporator: sym(56, 96, svg(56, 96, <>
    <rect className="pfd-body" x={12} y={8} width={32} height={80} rx={12} />
    <path d="M15 62 c4 -3 7 3 11 0 s7 3 11 0" strokeDasharray="0" />
    <path d="M20 46 v-10 M28 48 v-12 M36 46 v-10" strokeDasharray="3 2.5" />
  </>)),
  Valve: sym(56, 36, svg(56, 36, <>
    <path className="pfd-body" d="M6 6 L6 30 L28 18 Z" />
    <path className="pfd-body" d="M50 6 L50 30 L28 18 Z" />
  </>)),
  Pump: sym(64, 64, svg(64, 64, <>
    <circle className="pfd-body" cx={30} cy={38} r={24} />
    <path d="M22 28 L48 38 L22 48 Z" />
    <path d="M30 14 L58 14" />
  </>)),
  Compressor: sym(80, 64, svg(80, 64, <>
    <path className="pfd-body" d="M12 6 L12 58 L68 44 L68 20 Z" />
  </>)),
  Expander: sym(80, 64, svg(80, 64, <>
    <path className="pfd-body" d="M12 20 L12 44 L68 58 L68 6 Z" />
  </>)),
  ConversionReactor: sym(72, 100, svg(72, 100, <>
    <rect className="pfd-body" x={16} y={8} width={40} height={84} rx={12} />
    <path d="M28 36 L44 60 M44 36 L28 60" />
  </>)),
  EquilibriumReactor: sym(72, 100, svg(72, 100, <>
    <rect className="pfd-body" x={16} y={8} width={40} height={84} rx={12} />
    <path d="M26 42 h20 M26 54 h20" />
    <path d="M42 37 L47 42 L42 47 M30 49 L25 54 L30 59" />
  </>)),
  GibbsReactor: sym(72, 100, svg(72, 100, <>
    <rect className="pfd-body" x={16} y={8} width={40} height={84} rx={12} />
    <path d="M28 58 c0 -7 5 -11 8 -11 s8 4 8 11" />
    <path d="M36 47 V33 M30 38 L36 32 L42 38" />
  </>)),
  CSTR: sym(72, 100, svg(72, 100, <>
    <rect className="pfd-body" x={16} y={16} width={40} height={76} rx={12} />
    <path d="M36 6 V56 M26 56 a10 5 0 0 0 20 0 M26 56 a10 5 0 0 1 20 0" />
  </>)),
  PFR: sym(112, 48, svg(112, 48, <>
    <rect className="pfd-body" x={8} y={12} width={96} height={24} rx={12} />
    <path d="M28 12 v24 M46 12 v24 M64 12 v24 M82 12 v24" strokeDasharray="2 3" />
  </>)),
  ComponentSplitter: sym(60, 112, svg(60, 112, <>
    <rect className="pfd-body" x={14} y={10} width={32} height={92} rx={8} />
    <path d="M18 40 h24 M18 56 h24 M18 72 h24" strokeDasharray="2 2.5" />
  </>)),
  ThreePhaseSeparator: sym(120, 64, svg(120, 64, <>
    <rect className="pfd-body" x={8} y={14} width={104} height={36} rx={18} />
    <path d="M14 34 h92" strokeDasharray="3 2.5" />
    <path d="M14 42 h92" strokeDasharray="1.5 2.5" />
    <path d="M86 14 v36" />
  </>)),
  ShortcutColumn: column({ trays: 5, condenser: false, reboiler: false }),
  RigorousColumn: column({ trays: 6, condenser: true, reboiler: true }),
  Absorber: column({ trays: 4, packed: true }),
  ReboiledAbsorber: column({ trays: 3, packed: true, reboiler: true }),
  ExtractionColumn: sym(64, 130, svg(64, 130, <>
    <rect className="pfd-body" x={17} y={10} width={30} height={110} rx={8} />
    <circle cx={27} cy={34} r={2} /><circle cx={37} cy={50} r={2} />
    <circle cx={27} cy={66} r={2} /><circle cx={37} cy={82} r={2} />
    <circle cx={27} cy={98} r={2} />
  </>)),
  Cyclone: sym(60, 104, svg(60, 104, <>
    <path className="pfd-body" d="M14 8 L46 8 L46 40 L34 96 L26 96 L14 40 Z" />
    <path d="M22 30 c6 5 10 -5 16 0" strokeDasharray="2 2" />
    <path d="M22 44 c6 5 10 -5 16 0" strokeDasharray="2 2" />
  </>)),
  RotaryVacuumFilter: sym(84, 72, svg(84, 72, <>
    <circle className="pfd-body" cx={42} cy={34} r={26} />
    <circle cx={42} cy={34} r={11} strokeDasharray="3 2.5" />
    <path d="M14 62 h56" />
  </>)),
  BaghouseFilter: sym(72, 96, svg(72, 96, <>
    <rect className="pfd-body" x={12} y={10} width={48} height={76} rx={4} />
    <path d="M24 18 v54 M36 18 v54 M48 18 v54" />
    <path d="M14 80 h44" />
  </>)),
  Balance: sym(76, 56, svg(76, 56, <>
    <rect className="pfd-body" x={8} y={8} width={60} height={40} rx={6} />
    <path d="M38 14 v28 M20 42 h36" />
  </>)),
};

const FEED = sym(40, 30, svg(40, 30, <>
  <path d="M4 15 h22 M18 7 L28 15 L18 23" />
  <path d="M34 4 L34 26" />
</>));
const PRODUCT = sym(40, 30, svg(40, 30, <>
  <path d="M6 4 L6 26" />
  <path d="M12 15 h22 M24 7 L34 15 L24 23" />
</>));
const FALLBACK = sym(76, 54, svg(76, 54, <rect className="pfd-body" x={6} y={8} width={64} height={38} rx={6} />));

/** Natural-size symbol for a node. Boundary feed/product get stream-arrow stubs. */
export function symbolFor(kind: string, unitType?: string): PfdSymbol {
  if (kind === "feed") return FEED;
  if (kind === "product") return PRODUCT;
  return SYMBOLS[unitType ?? ""] ?? FALLBACK;
}

// ---- port anchors ----------------------------------------------------------

export interface Anchor { position: Position; offset: number }

const L = (offset: number): Anchor => ({ position: Position.Left, offset });
const R = (offset: number): Anchor => ({ position: Position.Right, offset });
const T = (offset: number): Anchor => ({ position: Position.Top, offset });
const B = (offset: number): Anchor => ({ position: Position.Bottom, offset });

// Per-type, per-port-name anchors. Physically meaningful: columns feed mid-left,
// overhead near top / bottoms near bottom, duties top (condenser) / bottom
// (reboiler); separators split vapor up / liquid down; rotating machinery in
// left, out right, work below. Ports not listed here fall back automatically.
const ANCHORS: Record<string, Record<string, Anchor>> = {
  Mixer: { in1: L(0.3), in2: L(0.7), out: R(0.5) },
  Splitter: { in1: L(0.5), out1: R(0.3), out2: R(0.7) },
  Heater: { in1: L(0.5), out: R(0.5), duty: B(0.5) },
  FiredHeater: { in1: L(0.5), out: R(0.5), duty: B(0.5) },
  AirCooler: { in1: L(0.5), out: R(0.5), duty: B(0.5) },
  // process stream horizontal (left/right), utility stream vertical (top/bottom)
  HeatExchanger: { hot_in: L(0.35), hot_out: R(0.35), cold_in: B(0.5), cold_out: T(0.5) },
  Flash: { in1: L(0.5), vapor: T(0.5), liquid: B(0.35), duty: B(0.75) },
  Evaporator: { in1: L(0.6), vapor: T(0.5), liquid: B(0.35), duty: B(0.75) },
  ThreePhaseSeparator: {
    in1: L(0.5), vapor: T(0.5), liquid_light: R(0.5), liquid_heavy: B(0.4), duty: B(0.75),
  },
  Pump: { in1: L(0.6), out: R(0.5), work: B(0.5) },
  Compressor: { in1: L(0.5), out: R(0.5), work: B(0.5) },
  Expander: { in1: L(0.5), out: R(0.5), work: B(0.5) },
  Valve: { in1: L(0.5), out: R(0.5) },
  PipeSegment: { in1: L(0.5), out: R(0.5) },
  ConversionReactor: { in1: L(0.35), out: R(0.5), duty: B(0.5) },
  EquilibriumReactor: { in1: L(0.35), out: R(0.5), duty: B(0.5) },
  GibbsReactor: { in1: L(0.35), out: R(0.5), duty: B(0.5) },
  CSTR: { in1: L(0.3), out: B(0.5), duty: R(0.8) },
  PFR: { in1: L(0.5), out: R(0.5), duty: B(0.5) },
  ComponentSplitter: { in1: L(0.5), overhead: R(0.12), bottoms: R(0.88), duty: B(0.5) },
  ShortcutColumn: {
    in1: L(0.5), distillate: R(0.1), bottoms: R(0.9),
    condenser_duty: T(0.5), reboiler_duty: B(0.5),
  },
  Absorber: { gas_in: L(0.9), liquid_in: L(0.1), vapor_out: T(0.5), liquid_out: B(0.5) },
  ReboiledAbsorber: { feed: L(0.4), vapor_out: T(0.5), bottoms: B(0.5), reboiler_duty: B(0.8) },
  ExtractionColumn: {
    feed_in: L(0.8), solvent_in: L(0.2), extract_out: R(0.8), raffinate_out: R(0.2), duty: B(0.5),
  },
  Cyclone: { gas_in: L(0.3), gas_out: T(0.5), solids_out: B(0.5), duty: R(0.8) },
  RotaryVacuumFilter: { slurry_in: L(0.5), filtrate_out: B(0.4), cake_out: R(0.5), duty: B(0.8) },
  BaghouseFilter: { gas_in: L(0.7), gas_out: T(0.5), solids_out: B(0.5), duty: R(0.85) },
};

/** RigorousColumn has a variable number of feeds (in1..inN) and side draws, so
 *  its anchors are computed from the actual port list rather than tabulated. */
function rigorousColumnAnchors(ports: Port[]): Record<string, Anchor> {
  const feeds = ports.filter((p) => /^in\d+$/.test(p.name));
  const sides = ports.filter((p) => /^side\d+$/.test(p.name));
  const out: Record<string, Anchor> = {
    distillate: R(0.08), bottoms: R(0.92),
    condenser_duty: T(0.5), reboiler_duty: B(0.5),
  };
  feeds.forEach((p, i) => { out[p.name] = L((i + 1) / (feeds.length + 1)); });
  sides.forEach((p, i) => { out[p.name] = R(0.3 + (0.4 * (i + 1)) / (sides.length + 1)); });
  return out;
}

const ANCHOR_FNS: Record<string, (ports: Port[]) => Record<string, Anchor>> = {
  RigorousColumn: rigorousColumnAnchors,
};

/** Evenly-spaced left(inlet)/right(material outlet)/bottom(energy) nozzles — the
 *  classic card behaviour, and the safety net for any un-tabulated port. */
function fallbackAnchors(ports: Port[]): Map<string, Anchor> {
  const inlets = ports.filter((p) => p.direction === "inlet");
  const outlets = ports.filter((p) => p.direction === "outlet" && p.kind !== "energy");
  const energy = ports.filter((p) => p.kind === "energy");
  const m = new Map<string, Anchor>();
  inlets.forEach((p, i) => m.set(p.name, L((i + 1) / (inlets.length + 1))));
  outlets.forEach((p, i) => m.set(p.name, R((i + 1) / (outlets.length + 1))));
  energy.forEach((p, i) => m.set(p.name, B((i + 1) / (energy.length + 1))));
  return m;
}

/** Resolve every port to a side + fractional offset. Custom table entries win;
 *  anything they omit keeps its evenly-spaced fallback so the node stays wired. */
export function portAnchors(kind: string, unitType: string | undefined, ports: Port[]): Map<string, Anchor> {
  const map = fallbackAnchors(ports);
  if (kind === "feed") { map.set("out", R(0.5)); return map; }
  if (kind === "product") { map.set("in", L(0.5)); return map; }
  const table = unitType
    ? (ANCHOR_FNS[unitType]?.(ports) ?? ANCHORS[unitType])
    : undefined;
  if (table) {
    for (const p of ports) {
      if (table[p.name]) map.set(p.name, table[p.name]);
    }
  }
  return map;
}

/** Types with a bespoke anchor table or function (for docs/tests). */
export const ANCHORED_TYPES = [...Object.keys(ANCHORS), ...Object.keys(ANCHOR_FNS)];
