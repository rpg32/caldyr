# Thermodynamics

## Property packages

| id | Backend | Use for |
|---|---|---|
| `thermo:PR` | Peng-Robinson cubic EOS (`thermo` library) | Non-polar gases & hydrocarbons |
| `thermo:SRK` | Soave-Redlich-Kwong cubic EOS | Non-polar (alternative) |
| `thermo:NRTL` | NRTL gamma-phi, ChemSep parameters | Polar mixtures, azeotropes |
| `coolprop:Water` | IAPWS-95 steam tables (CoolProp) | Pure-water / steam systems |

All expose the same `PropertyPackage` protocol: enthalpy/entropy/volume,
PT/PH/PS flashes, bubble/dew (+ `bubble_point` for per-stage column use),
per-phase VLE compositions and enthalpies, three-phase `flash_pt_3p`
(PR/SRK only), and reaction `lnKeq`. Phases are identified by molar volume,
which stays robust at high pressure where label-based identification fails.

Notes per package:

- **`coolprop:Water`** is the HYSYS "NBS/ASME Steam" equivalent: a
  reference-quality pure-water Helmholtz EOS (IAPWS-95, the basis of the
  modern steam tables), far more accurate for water/steam than a cubic EOS —
  PR's water enthalpy departures are a few percent off the IAPWS values.
  **Single-component water only**; anything else raises a typed error. Its
  enthalpy is shifted by one constant so it sits on the same
  formation-inclusive basis as the `thermo` packages (anchored at the
  298.15 K ideal-gas state), so mixed-package comparisons of ΔH are honest.
- **`thermo:NRTL`** reads its binary interaction parameters from the ChemSep
  databank bundled with `thermo`. Pairs *without* ChemSep data silently fall
  back to an ideal-solution liquid (a warning is emitted only when **no**
  pair has data) — check coverage for your system before trusting an
  azeotrope prediction. Example: dimethyl ether has no ChemSep NRTL pairs
  with methanol or water, which is why the DME plant template runs on
  `thermo:PR`.

## ⚠ The enthalpy reference state (read this once)

**Caldyr enthalpies are formation-inclusive (absolute basis):** every stream's
`H` includes the composition-weighted ideal-gas enthalpy of formation,
`H = H_sensible + Σ zᵢ·ΔHf°ᵢ`.

Why: heats of reaction then emerge automatically from plain energy balances —
reactors need no separate ΔH_rxn bookkeeping, adiabatic flame/reactor
temperatures are right by construction, and every unit's balance closes on one
consistent basis.

Consequences to be aware of:

- **Do not compare `stream.H` against steam tables or simulator outputs that
  use a different reference** (e.g. HYSYS' elemental basis differs again). The
  *differences* between Caldyr enthalpies are physical; the absolute values
  carry the formation offset.
- For composition-conserved balances the offset cancels exactly, so nothing
  changes for mixers, heaters, exchangers, flashes.
- To recover a sensible-only enthalpy, subtract the mixture formation
  enthalpy: the offset is linear in composition,
  `H_sensible = H − Σ zᵢ·ΔHf°ᵢ`, with ΔHf° from the `chemicals` library
  (`chemicals.reaction.Hfg(CAS)`), the same source the engine uses.

This is documented in code at `engine/caldyr/thermo/_flasher.py`
(`formation_props`, `_hf_mix`) and verified by tests that close reactor energy
balances to machine precision.

## Bubble/dew points with permanent gases (non-condensables)

`thermo`'s native VF=0 / VF=1 ("PVF") flashes can fail outright when the
liquid carries dissolved supercritical components — H2 or N2 in a column
feed, HCl above its Tc on a hot stage. The engine layers documented fallbacks
behind the same `bubble_dew` / `bubble_point` API instead of crashing:

1. **Bubble point**: if the PVF flash fails, bisect the full PT flash's
   vapor fraction for the boiling edge (the PT flash carries stability
   analysis and works everywhere). VF(T) need not be monotone — dissolved
   permanent gases are *less* soluble cold — so the search targets the
   *upper* (boiling) transition.
2. **No all-liquid window at any T** (enough dissolved H2/N2/CH4 keeps
   VF > 0 everywhere): the classical non-condensable treatment (Seader 3e
   sec. 4.4) — take the bubble point of the *condensable submixture*
   (components with Tc ≥ 240 K), then fold the light gases into the
   incipient vapor by their K-values at that temperature.
3. **Dew point**: on PVF failure, bisect the all-vapor edge of the PT flash
   scanning down from high temperature.

Caveat: the non-condensable bubble *surrogate* can land above the true dew
edge; `bubble_dew` clamps the pair ordered (`min(bubble, dew), dew`) so
phase-envelope consumers see a degenerate-but-sane interval. These fallbacks
are what let bubble-point MESH columns tolerate traces of permanent gases —
but a feed with *percent*-level lights still deserves a degassing flash
first (see the cyclohexane plant example), and HCl-class near-critical
components may need extra reflux to keep them out of hot stage liquids (see
the VCM template).

## ⚠ Three-phase (VLLE) liquid labels can invert

`flash_pt_3p` (built on `thermo`'s `FlashVLN`; PR/SRK only) labels its two
liquids `light`/`heavy` by **mass density as the EOS computes it** — and
cubic-EOS liquid densities can be badly wrong for water (PR puts water near
toluene's density, a < 0.1% difference). When the two liquid densities are
near-degenerate, the light/heavy labels can come out **inverted** relative to
reality: the "aqueous" phase may land on `liquid_light` and the organic on
`liquid_heavy`.

Practical guidance:

- Identify phases by **composition**, not by which port they left from. The
  `ThreePhaseSeparator` stores the phase-split betas and both liquid
  densities on `unit.result` so a suspicious split is visible.
- The `ExtractionColumn` already defends itself: it tracks phases by
  composition continuity between stages instead of trusting the density
  labels.

## Validation anchors

- Normal boiling points vs CRC/NIST (water, ethanol, ±2 K).
- Ethanol/water azeotrope vs literature (NRTL captures it; PR cannot — by design).
- Haber-Bosch equilibrium (`lnKeq`) vs Smith, Van Ness & Abbott.
- Water/n-hexane VLLE structure vs Tsonopoulos (1999).
