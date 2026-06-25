// Stream composition helpers. The engine stores per-stream composition `z`
// (mole fractions, SI). These format it for display — normalize defensively
// (a solved product's z should already sum to 1, but a raw feed may not) and
// derive per-component molar flow. Mass fractions need molar masses, which the
// thin web client does not carry, so they are intentionally out of scope here.

export interface CompRow {
  comp: string;
  /** mole fraction, normalized to sum 1 */
  frac: number;
  /** per-component molar flow (mol/s, SI) or null if the stream flow is unknown */
  flow: number | null;
}

/** Sorted (descending mole fraction) composition rows for a stream. */
export function compositionRows(
  z: Record<string, number> | undefined | null,
  molarFlow: number | null | undefined,
): CompRow[] {
  if (!z) return [];
  const entries = Object.entries(z).filter(([, v]) => Number.isFinite(v));
  const total = entries.reduce((a, [, v]) => a + v, 0);
  if (total <= 0) return [];
  return entries
    .map(([comp, v]) => {
      const frac = v / total;
      return {
        comp,
        frac,
        flow: molarFlow != null ? molarFlow * frac : null,
      };
    })
    .sort((a, b) => b.frac - a.frac);
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
