// Stream composition helpers. The engine stores per-stream composition `z`
// (mole fractions, SI) and the /solve payload carries per-component molar mass
// `molar_mass` (kg/mol). These format composition for display: normalize `z`
// defensively, derive per-component molar flow, and — when molar masses are
// available — mass fractions and mass flows.

export interface CompRow {
  comp: string;
  /** mole fraction, normalized to sum 1 */
  frac: number;
  /** per-component molar flow (mol/s, SI) or null if the stream flow is unknown */
  flow: number | null;
  /** mass fraction, or null if any component's molar mass is unavailable */
  massFrac: number | null;
  /** per-component mass flow (kg/s, SI), or null if unavailable */
  massFlow: number | null;
}

type MwMap = Record<string, number> | undefined | null;

/** Sorted (descending mole fraction) composition rows for a stream. Pass the
 *  molar-mass map (kg/mol) to also get mass fraction + mass flow per component. */
export function compositionRows(
  z: Record<string, number> | undefined | null,
  molarFlow: number | null | undefined,
  mw?: MwMap,
): CompRow[] {
  if (!z) return [];
  const entries = Object.entries(z).filter(([, v]) => Number.isFinite(v));
  const total = entries.reduce((a, [, v]) => a + v, 0);
  if (total <= 0) return [];

  // Mass basis needs a molar mass for EVERY component, else fractions are wrong.
  const haveMw = !!mw && entries.every(([c]) => Number.isFinite(mw[c]));
  // Σ x·M (kg/mol of the mixture) — used to split mass fractions.
  const mixM = haveMw
    ? entries.reduce((a, [c, v]) => a + (v / total) * (mw as Record<string, number>)[c], 0)
    : 0;

  return entries
    .map(([comp, v]) => {
      const frac = v / total;
      const massFrac = haveMw && mixM > 0
        ? (frac * (mw as Record<string, number>)[comp]) / mixM
        : null;
      return {
        comp,
        frac,
        flow: molarFlow != null ? molarFlow * frac : null,
        massFrac,
        massFlow: massFrac != null && molarFlow != null
          ? molarFlow * mixM * massFrac     // mol/s · kg/mol · (—) = kg/s
          : null,
      };
    })
    .sort((a, b) => b.frac - a.frac);
}

/** Total stream mass flow (kg/s), or null if molar flow / any molar mass is
 *  unavailable. = molar_flow · Σ(x·M). */
export function streamMassFlow(
  z: Record<string, number> | undefined | null,
  molarFlow: number | null | undefined,
  mw: MwMap,
): number | null {
  if (!z || molarFlow == null || !mw) return null;
  const entries = Object.entries(z).filter(([, v]) => Number.isFinite(v));
  const total = entries.reduce((a, [, v]) => a + v, 0);
  if (total <= 0) return null;
  if (!entries.every(([c]) => Number.isFinite(mw[c]))) return null;
  const mixM = entries.reduce((a, [c, v]) => a + (v / total) * mw[c], 0);
  return molarFlow * mixM;
}

/** Stable, union component order across many streams (for table columns). */
export function componentOrder(
  streams: { z?: Record<string, number> | null }[],
): string[] {
  const seen = new Set<string>();
  for (const s of streams) {
    if (!s.z) continue;
    for (const c of Object.keys(s.z)) seen.add(c);
  }
  return [...seen];
}

/** Format a mole/mass fraction for compact display (e.g. 0.618, 1.2e-4). */
export function fmtFrac(x: number): string {
  if (x === 0) return "0";
  if (x >= 0.001) return x.toFixed(4);
  return x.toExponential(1);
}
