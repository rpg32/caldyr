// Small display-formatting helpers for labels.

/** Turn an engine type id into spaced words for display, preserving acronyms.
 *  "EquilibriumReactor" -> "Equilibrium Reactor", "CSTR" -> "CSTR",
 *  "ThreePhaseSeparator" -> "Three Phase Separator", "PFR" -> "PFR". */
export function humanizeType(type: string): string {
  return type
    .replace(/([a-z0-9])([A-Z])/g, "$1 $2")   // camel boundary: xY -> x Y
    .replace(/([A-Z]+)([A-Z][a-z])/g, "$1 $2") // acronym then word: HTTPServer -> HTTP Server
    .trim();
}
