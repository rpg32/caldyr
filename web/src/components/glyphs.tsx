// PFD-style equipment glyphs, one per unit type. Simple strokes inheriting
// currentColor so they recolor with theme/selection. 28x28 viewBox.
import type { ReactNode } from "react";

const S = {
  fill: "none",
  stroke: "currentColor",
  strokeWidth: 1.6,
  strokeLinecap: "round" as const,
  strokeLinejoin: "round" as const,
};

const svg = (children: ReactNode) => (
  <svg width="26" height="26" viewBox="0 0 28 28" aria-hidden {...S}>
    {children}
  </svg>
);

export const GLYPHS: Record<string, ReactNode> = {
  Mixer: svg(<>
    <path d="M3 7 L12 12 M3 21 L12 16" />
    <path d="M11 6 L11 22 L22 17 L22 11 Z" />
    <path d="M22 14 L26 14" />
  </>),
  Splitter: svg(<>
    <path d="M2 14 L6 14" />
    <path d="M6 11 L6 17 L17 22 L17 6 Z" />
    <path d="M17 9 L25 6 M17 19 L25 22" />
  </>),
  Heater: svg(<>
    <circle cx="14" cy="14" r="10" />
    <path d="M6 14 h4 l2 -4 l4 8 l2 -4 h4" />
  </>),
  FiredHeater: svg(<>
    <path d="M6 25 L6 12 L14 4 L22 12 L22 25 Z" />
    <path d="M14 21 c-2.5 0 -3.5 -2.6 -1.8 -4.3 c0.3 1 0.9 1.4 1.5 1.5 c-0.6 -2 0.3 -3.6 2 -4.6 c-0.3 2 0.7 2.6 1.3 3.8 c1 2 -0.5 3.6 -3 3.6 Z" />
  </>),
  AirCooler: svg(<>
    <rect x="3" y="10" width="22" height="9" rx="2" />
    <path d="M9 14.5 l3.5 -2.5 M9 14.5 l3.5 2.5 M19 14.5 l-3.5 -2.5 M19 14.5 l-3.5 2.5" />
    <path d="M7 7 L9 9.5 M14 6 L14 9.5 M21 7 L19 9.5" />
  </>),
  HeatExchanger: svg(<>
    <circle cx="14" cy="14" r="10" />
    <path d="M2 14 h5 c3 0 3 -5 7 -5 s4 10 7 10 c2 0 3 -2.5 5 -2.5" />
  </>),
  Flash: svg(<>
    <rect x="9" y="3" width="10" height="22" rx="5" />
    <path d="M11 11 h6 M11 17 h6" strokeDasharray="2.5 2" />
  </>),
  Valve: svg(<>
    <path d="M4 9 L14 14 L4 19 Z" />
    <path d="M24 9 L14 14 L24 19 Z" />
  </>),
  Pump: svg(<>
    <circle cx="13" cy="15" r="9" />
    <path d="M9 11 L19 15 L9 19 Z" />
    <path d="M13 6 L24 6" />
  </>),
  Compressor: svg(<>
    <path d="M5 6 L5 22 L23 17 L23 11 Z" />
    <path d="M2 14 L5 14 M23 14 L26 14" />
  </>),
  ConversionReactor: svg(<>
    <rect x="7" y="3" width="14" height="22" rx="6" />
    <path d="M11 10 L17 18 M17 10 L11 18" />
  </>),
  EquilibriumReactor: svg(<>
    <rect x="7" y="3" width="14" height="22" rx="6" />
    <path d="M10.5 12 h7 M10.5 16 h7" />
    <path d="M15.5 9.5 L17.5 12 L15.5 14.5 M12.5 13.5 L10.5 16 L12.5 18.5" />
  </>),
  ShortcutColumn: svg(<>
    <rect x="9" y="2" width="10" height="24" rx="4" />
    <path d="M10 8 h8 M10 13 h8 M10 18 h8" />
  </>),
  PipeSegment: svg(<>
    <path d="M2 11 L18 11 C21 11 21 14 18 14 L10 14 C7 14 7 17 10 17 L26 17" />
    <path d="M2 8 L2 14 M26 14 L26 20" />
  </>),
  ExtractionColumn: svg(<>
    <rect x="9" y="2" width="10" height="24" rx="4" />
    <circle cx="12.5" cy="8" r="1.2" /><circle cx="15.5" cy="12" r="1.2" />
    <circle cx="12.5" cy="16" r="1.2" /><circle cx="15.5" cy="20" r="1.2" />
    <path d="M5 5 L9 5 M19 23 L23 23" />
  </>),
  Absorber: svg(<>
    <rect x="9" y="2" width="10" height="24" rx="4" />
    <path d="M10 7.5 h8 M10 12 h8 M10 16.5 h8 M10 21 h8" strokeDasharray="2 1.5" />
    <path d="M5 5 L9 5 M19 23 L23 23" />
  </>),
  ReboiledAbsorber: svg(<>
    <rect x="9" y="2" width="10" height="20" rx="4" />
    <path d="M10 7 h8 M10 11.5 h8 M10 16 h8" strokeDasharray="2 1.5" />
    <path d="M9 24 c0 2 2.5 3 5 3 s5 -1 5 -3 M19 24 c2.5 0 2.5 -3 0 -3" />
  </>),
  RigorousColumn: svg(<>
    <rect x="9" y="2" width="10" height="24" rx="4" />
    <path d="M10 6.5 h8 M10 10 h8 M10 13.5 h8 M10 17 h8 M10 20.5 h8" />
    <path d="M19 5 c3 0 3 3.5 0 3.5 M9 23 c-3 0 -3 -3.5 0 -3.5" />
  </>),
  GibbsReactor: svg(<>
    <rect x="7" y="3" width="14" height="22" rx="6" />
    <path d="M11 17 c0 -2 1.5 -3 3 -3 s3 1 3 3" />
    <path d="M14 14 L14 9 M11.5 10.5 L14 8 L16.5 10.5" />
  </>),
  CSTR: svg(<>
    <rect x="7" y="5" width="14" height="20" rx="5" />
    <path d="M14 2 L14 14 M10 14 a4 2.2 0 0 0 8 0 M10 14 a4 2.2 0 0 1 8 0" />
  </>),
  PFR: svg(<>
    <rect x="2" y="9" width="24" height="10" rx="5" />
    <path d="M7 9 L7 19 M12 9 L12 19 M17 9 L17 19 M22 9 L22 19" strokeDasharray="1.5 2" />
  </>),
  Expander: svg(<>
    <path d="M5 11 L5 17 L23 22 L23 6 Z" />
    <path d="M2 14 L5 14 M23 14 L26 14" />
  </>),
  ComponentSplitter: svg(<>
    <rect x="8" y="4" width="12" height="20" rx="4" />
    <path d="M11 10 h6 M11 14 h6 M11 18 h6" strokeDasharray="1.5 1.5" />
    <path d="M20 8 L26 6 M20 20 L26 22" />
  </>),
  ThreePhaseSeparator: svg(<>
    <rect x="3" y="8" width="22" height="13" rx="6" />
    <path d="M5 15 h18" strokeDasharray="2.5 2" />
    <path d="M5 18 h18" strokeDasharray="1.2 1.8" />
    <path d="M9 8 L9 4 M14 21 L14 25 M20 21 L20 25" />
  </>),
  Balance: svg(<>
    <path d="M14 4 L14 24 M4 24 h20" />
    <path d="M5 9 L23 9 M8 9 l-3 5 h6 Z M20 9 l-3 5 h6 Z" />
  </>),
  Evaporator: svg(<>
    <rect x="5" y="7" width="18" height="16" rx="4" />
    <path d="M8 19 c2 -1.5 4 1.5 6 0 s4 1.5 6 0" />
    <path d="M10 14 L10 10 M14 15 L14 10 M18 14 L18 10" strokeDasharray="2 1.6" />
  </>),
  feed: svg(<>
    <path d="M3 14 L19 14 M13 8 L19 14 L13 20" />
    <path d="M23 6 L23 22" />
  </>),
  product: svg(<>
    <path d="M5 6 L5 22" />
    <path d="M9 14 L25 14 M19 8 L25 14 L19 20" />
  </>),
};

const FALLBACK = svg(<rect x="4" y="7" width="20" height="14" rx="3" />);

export function glyphFor(kind: string, unitType?: string): ReactNode {
  if (kind === "feed" || kind === "product") return GLYPHS[kind];
  return GLYPHS[unitType ?? ""] ?? FALLBACK;
}
