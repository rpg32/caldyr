// Plain-language definitions of the process-engineering + techno-economic jargon
// the app surfaces. Used by the inline <Term> tooltip and the Glossary dialog.

export interface GlossaryEntry {
  term: string;
  def: string;
  aka?: string[];   // alternate spellings/abbreviations matched by lookup
}

export const GLOSSARY: GlossaryEntry[] = [
  { term: "LCOP", def: "Levelized Cost of Product — the product price ($/kg) at which the project's net present value is exactly zero; the break-even cost per kg over the plant life." },
  { term: "TCI", def: "Total Capital Investment — the grassroots fixed capital plus working capital; the up-front money to build and start the plant." },
  { term: "ISBL", def: "Inside Battery Limits — the installed cost of the process equipment itself (the sum of the bare-module costs).", aka: ["inside battery limits"] },
  { term: "OSBL", def: "Outside Battery Limits — offsite / auxiliary facilities (utilities, storage, infrastructure); here ≈ 0.5 × Σ base bare-module cost.", aka: ["outside battery limits", "offsite"] },
  { term: "Cbm", def: "Bare-Module Cost — the installed cost of one equipment item: purchased cost × a bare-module factor accounting for installation, materials and pressure.", aka: ["CBM", "bare module", "bare-module"] },
  { term: "Fbm", def: "Bare-Module Factor — the multiplier applied to an item's purchased cost to get its installed cost (materials + pressure + labor).", aka: ["bare-module factor"] },
  { term: "COM", def: "Cost of Manufacturing — the total annual production cost. Turton's COM_d = 0.18·FCI + 2.73·labor + 1.23·(utilities + raw materials).", aka: ["COM_d", "manufacturing cost"] },
  { term: "FCI", def: "Fixed Capital Investment (grassroots) — the total fixed capital needed to build the plant on a new site.", aka: ["fixed capital"] },
  { term: "OPEX", def: "Operating Expenditure — annual operating cost: raw materials + utilities + fixed costs (labor, maintenance, overheads).", aka: ["opex", "operating cost"] },
  { term: "CAPEX", def: "Capital Expenditure — the one-time cost of building the plant (equipment, installation, offsites).", aka: ["capex", "capital cost"] },
  { term: "NPV", def: "Net Present Value — the sum of all project cash flows discounted to today; positive means value-creating at the chosen discount rate." },
  { term: "IRR", def: "Internal Rate of Return — the discount rate at which the project's NPV equals zero; compare it to your hurdle rate." },
  { term: "discount rate", def: "The annual rate used to convert future cash flows to present value (the time-value of money / cost of capital)." },
  { term: "grassroots", def: "The total cost to build the plant on a new, undeveloped (green-field) site, including offsites." },
  { term: "working capital", def: "Cash tied up in startup inventory and day-to-day operations; here ≈ 15% of the grassroots fixed capital." },
  { term: "CEPCI", def: "Chemical Engineering Plant Cost Index — an index that escalates historical equipment costs to a target year." },
  { term: "duty", def: "The heat or work rate (in watts) added to or removed from a unit — carried on its energy (duty) port." },
  { term: "vapor fraction", def: "The mole fraction of a stream that is vapor at its T and P (0 = all liquid, 1 = all vapor).", aka: ["vapor frac", "VF"] },
  { term: "tear stream", def: "In a recycle loop, the stream the sequential solver 'cuts' and iterates on until the loop converges." },
  { term: "reflux ratio", def: "In a distillation column, the ratio of liquid returned to the top stage (reflux L) to the distillate product taken off (D): R = L/D." },
  { term: "Murphree efficiency", def: "A tray's actual composition change divided by the equilibrium (ideal) change; 1.0 = an ideal equilibrium stage." },
  { term: "HETP", def: "Height Equivalent to a Theoretical Plate — the packed-column height that achieves one equilibrium stage." },
  { term: "LMTD", def: "Log-Mean Temperature Difference — the average driving-force ΔT used to size a heat exchanger: Q = U·A·LMTD." },
  { term: "U-value", def: "Overall heat-transfer coefficient (W/m²·K) — how readily a surface transfers heat; sets exchanger area for a given duty.", aka: ["overall U", "hx_overall_U"] },
  { term: "pinch", def: "In heat integration, the point of closest temperature approach between the hot and cold composite curves — it sets the minimum utility targets." },
  { term: "approach", def: "The minimum allowed temperature difference between two streams in an exchanger or heat-exchange network (ΔT_min).", aka: ["ΔT min", "dt_min", "min approach"] },
  { term: "residence time", def: "The average time material spends inside a vessel; a sizing basis for separators and reactors." },
  { term: "flooding", def: "The vapor velocity at which liquid can no longer drain down a column; trays are designed at a fraction (e.g. 80%) of it." },
  { term: "bubble point", def: "The temperature at which a liquid mixture first starts to boil (the first vapor bubble) at a given pressure." },
  { term: "dew point", def: "The temperature at which a vapor mixture first starts to condense (the first liquid drop) at a given pressure." },
  { term: "property package", def: "The thermodynamic model (PR, SRK, NRTL, UNIFAC…) used to compute phase equilibrium and properties." },
  { term: "bare-module", def: "See Cbm — the installed cost of an equipment item including installation, materials and pressure factors.", aka: ["bare module cost"] },
];

const INDEX: Map<string, GlossaryEntry> = (() => {
  const m = new Map<string, GlossaryEntry>();
  for (const e of GLOSSARY) {
    m.set(e.term.toLowerCase(), e);
    for (const a of e.aka ?? []) m.set(a.toLowerCase(), e);
  }
  return m;
})();

/** Definition for a term (matched case-insensitively against term + aka), or undefined. */
export const glossaryDef = (key: string): string | undefined =>
  INDEX.get(key.trim().toLowerCase())?.def;
