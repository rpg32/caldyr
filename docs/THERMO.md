# Thermodynamics

## Property packages

| id | Backend | Use for |
|---|---|---|
| `thermo:PR` | Peng-Robinson cubic EOS (`thermo` library) | Non-polar gases & hydrocarbons |
| `thermo:SRK` | Soave-Redlich-Kwong cubic EOS | Non-polar (alternative) |
| `thermo:NRTL` | NRTL gamma-phi, ChemSep parameters | Polar mixtures, azeotropes |

All expose the same `PropertyPackage` protocol: enthalpy/entropy/volume,
PT/PH/PS flashes, bubble/dew (+ `bubble_point` for per-stage column use),
per-phase VLE compositions and enthalpies, three-phase `flash_pt_3p`
(PR/SRK only), and reaction `lnKeq`. Phases are identified by molar volume,
which stays robust at high pressure where label-based identification fails.

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

## Validation anchors

- Normal boiling points vs CRC/NIST (water, ethanol, ±2 K).
- Ethanol/water azeotrope vs literature (NRTL captures it; PR cannot — by design).
- Haber-Bosch equilibrium (`lnKeq`) vs Smith, Van Ness & Abbott.
- Water/n-hexane VLLE structure vs Tsonopoulos (1999).
