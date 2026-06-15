"""M15 demo: the "book chapter 12" solids operations batch (Hameed, *Chemical
Process Simulations using Aspen HYSYS*, Wiley 2025).

Three worked examples from the book, reproduced with Caldyr's v1 particle
model (the solid is an ordinary stream component; the particle size
distribution lives on the *unit* as the ``psd`` param — exactly how the book
enters it; PSDs do not propagate between units):

1. **Cyclone (sec. 12.1)** — 150 kmol/h of air at 45 C / 150 kPa carrying
   5 wt% carbon dust (rho_p = 1642 kg/m^3). The book's HYSYS design mode hits
   exactly its 95% efficiency spec with the High Output (Stairmand HT)
   geometry at D = 2.881 m and 3 parallel cyclones for the discrete PSD
   {1.5 mm: 40%, 2.0: 35%, 2.5: 15%, 3.0: 10%} (Fig. 12.9). Caldyr evaluates
   the *Lapple* cut-diameter model at that geometry: overall efficiency 98.4%
   (delta +3.4 points vs the book's Leith/Licht design point), with the full
   grade-efficiency table, the d50, and the Shepherd-Lapple pressure drop on
   ``unit.design``.

2. **Rotary vacuum filter (sec. 12.2)** — 4000 kg/h of slurry (30 wt%
   calcium in water) at 25 C / 160 kPa, dP = 10 kPa, 10-min cycle, 20%
   submergence, 0.4 m drum radius. HYSYS reports area 2.4368 m^2 and drum
   width 0.96957 m (Fig. 12.13); its internal cake resistance is not
   published, so alpha = 1.032e14 m/kg is backed out and the McCabe
   continuous-filtration design then reproduces the book's area and width.

3. **Baghouse filter + the cyclone-baghouse train (sec. 12.3)** — 60 kmol/h
   of air + 5 mol% sulfur at 35 C / 180 kPa, sized on the air-to-cloth
   ratio, with the filter-drag dP model giving the filtration time to a
   2 kPa dirty-bag limit (the book's case-study variable, Fig. 12.18). Then
   the book's Exercise-12.3 train (p. 408-409): 10,000 kg/h air + 2000 kg/h
   kaolin at 30 C / 5 bar through a 75% cyclone (gas keeps 500 kg/h, book
   value) and a baghouse that removes ~100% of the rest. Kaolin has no
   critical constants in the databank, so carbon stands in as the dust
   component with the book's kaolin PSD and density.

Everything is Turton-style costed at the end (the solids correlations are
order-of-magnitude — Couper/Walas, Peters-Timmerhaus-West and the EPA Cost
Manual; see caldyr/economics/data.py for the confidence notes).

Run from the repo root:

    python examples/21_solids.py
"""
import sys
from pathlib import Path

_ENGINE = Path(__file__).resolve().parent.parent / "engine"
if _ENGINE.is_dir() and str(_ENGINE) not in sys.path:
    sys.path.insert(0, str(_ENGINE))

from caldyr.core import Component, Flowsheet  # noqa: E402
from caldyr.economics.costing import cost_equipment  # noqa: E402
from caldyr.economics.sizing import size_flowsheet  # noqa: E402
from caldyr.thermo import make_package  # noqa: E402
from caldyr.unitops import BaghouseFilter, Cyclone, RotaryVacuumFilter  # noqa: E402

MW_C = 0.0120107
MW_CA = 0.040078
MW_W = 0.01801528
MW_S = 0.032065
MW_AIR = 0.028851

# -- 1. the book's sec. 12.1 cyclone (PSD case, High Output geometry) ---------
print("=" * 72)
print("1. Cyclone — Hameed sec. 12.1 (book Fig. 12.9: D=2.881 m, 3 cyclones)")
print("=" * 72)
_mc, _mair = 0.05 / MW_C, 0.95 / MW_AIR          # 5 wt% carbon -> mole fracs
xc = _mc / (_mc + _mair)
cyc_fs = Flowsheet(
    components=[Component("nitrogen"), Component("oxygen"), Component("carbon")],
    property_package="thermo:PR")
cyc_fs.add(Cyclone("CYC", {
    "solids": "carbon", "particle_density": 1642.0,
    "geometry": "Stairmand_HT", "body_diameter": 2.881, "n_cyclones": 3,
    "psd": [{"d_microns": 1500.0, "mass_frac": 0.40},
            {"d_microns": 2000.0, "mass_frac": 0.35},
            {"d_microns": 2500.0, "mass_frac": 0.15},
            {"d_microns": 3000.0, "mass_frac": 0.10}],
}))
cyc_fs.feed("F", "CYC:gas_in", T=318.15, P=150e3, molar_flow=150e3 / 3600.0,
            z={"nitrogen": (1 - xc) * 0.79, "oxygen": (1 - xc) * 0.21, "carbon": xc})
cyc_fs.connect("G", "CYC:gas_out", None)
cyc_fs.connect("S", "CYC:solids_out", None)
cyc_rep = cyc_fs.solve()
d = cyc_fs.units["CYC"].design
print(f"  geometry            {d['geometry']} (book 'High Output' ratios)")
print(f"  body diameter       {d['body_diameter_m']:.3f} m x {d['n_cyclones']} cyclones")
print(f"  inlet velocity      {d['inlet_velocity_m_s']:.3f} m/s")
print(f"  cut diameter d50    {d['d50_microns']:.1f} um")
print("  grade efficiencies:")
for g in d["grade"]:
    print(f"    {g['d_microns']:7.0f} um  ({g['mass_frac'] * 100:4.0f} wt%)  "
          f"eta = {g['efficiency'] * 100:.2f} %")
print(f"  overall efficiency  {d['overall_efficiency'] * 100:.2f} %   "
      f"(book design spec: 95.00 % -> delta +{(d['overall_efficiency'] - 0.95) * 100:.1f} pts)")
print(f"  pressure drop       {d['dP_Pa']:.1f} Pa (Shepherd-Lapple, "
      f"N_H = {d['NH_velocity_heads']:.1f})")
print(f"  dust captured       {d['solids_captured_kg_s'] * 3600:.1f} of "
      f"{d['solids_in_kg_s'] * 3600:.1f} kg/h")

# -- 2. the book's sec. 12.2 rotary vacuum filter ------------------------------
print()
print("=" * 72)
print("2. Rotary vacuum filter — Hameed sec. 12.2 (book Fig. 12.13)")
print("=" * 72)
n_ca = (1200.0 / 3600.0) / MW_CA                 # 30 wt% calcium of 4000 kg/h
n_w = (2800.0 / 3600.0) / MW_W
rvf_fs = Flowsheet(components=[Component("water"), Component("calcium")],
                   property_package="thermo:PR")
rvf_fs.add(RotaryVacuumFilter("RVF", {
    "solids": "calcium", "pressure_drop": 10e3, "cycle_time_s": 600.0,
    "submergence": 0.20, "cake_moisture": 0.5, "drum_radius_m": 0.4,
    # HYSYS computes alpha internally (particle size/sphericity); the book
    # does not publish it, so it is backed out to match the book's area.
    "alpha": 1.032e14,
}))
rvf_fs.feed("F", "RVF:slurry_in", T=298.15, P=160e3, molar_flow=n_ca + n_w,
            z={"water": n_w / (n_ca + n_w), "calcium": n_ca / (n_ca + n_w)})
rvf_fs.connect("L", "RVF:filtrate_out", None)
rvf_fs.connect("C", "RVF:cake_out", None)
rvf_rep = rvf_fs.solve()
d = rvf_fs.units["RVF"].design
print(f"  filtration area     {d['area_m2']:.4f} m^2   (book: 2.4368 m^2)")
print(f"  drum                R = {d['drum_radius_m']:.2f} m x "
      f"W = {d['drum_width_m']:.5f} m   (book width: 0.96957 m)")
print(f"  cycle-avg flux      {d['flux_m3_m2_s'] * 1e3:.3f} mm/s of filtrate")
print(f"  cake                {d['cake_solids_kg_s'] * 3600:.0f} kg/h solids + "
      f"{d['cake_liquid_kg_s'] * 3600:.0f} kg/h liquid "
      f"({d['cake_moisture'] * 100:.0f} % moisture)")
print(f"  filtrate            {d['filtrate_m3_s'] * 3600:.3f} m^3/h at "
      f"{rvf_fs.streams['L'].P / 1e3:.0f} kPa (vacuum side)")

# -- 3a. the book's sec. 12.3 baghouse ------------------------------------------
print()
print("=" * 72)
print("3a. Baghouse filter — Hameed sec. 12.3 (60 kmol/h air + 5 mol% sulfur)")
print("=" * 72)
bh_fs = Flowsheet(
    components=[Component("nitrogen"), Component("oxygen"), Component("sulfur")],
    property_package="thermo:PR")
bh_fs.add(BaghouseFilter("BH", {"solids": "sulfur", "dP_max": 2000.0}))
bh_fs.feed("F", "BH:gas_in", T=308.15, P=180e3, molar_flow=60e3 / 3600.0,
           z={"nitrogen": 0.95 * 0.79, "oxygen": 0.95 * 0.21, "sulfur": 0.05})
bh_fs.connect("G", "BH:gas_out", None)
bh_fs.connect("S", "BH:solids_out", None)
bh_rep = bh_fs.solve()
d = bh_fs.units["BH"].design
print(f"  gas flow            {d['Q_m3_s']:.3f} m^3/s")
print(f"  cloth area          {d['cloth_area_m2']:.1f} m^2 at "
      f"{d['face_velocity_m_s'] * 100:.0f} cm/s air-to-cloth (~{d['n_bags']} bags)")
print(f"  dust loading        {d['dust_loading_kg_m3'] * 1e3:.1f} g/m^3")
print(f"  collection          {d['overall_efficiency'] * 100:.1f} % -> "
      f"{d['solids_emitted_kg_s'] * 3600:.3f} kg/h escapes")
print(f"  filtration time to dP = {d['dP_max_Pa']:.0f} Pa: "
      f"{d['filtration_time_s']:.0f} s "
      f"({d['filtration_time_s'] / 3600:.2f} h between cleanings)")
print("  (the book's HYSYS number, 5.92e7 s, uses unpublished S_E/K2 "
      "defaults; ours are the cited Cooper & Alley mid-range values)")

# -- 3b. the Exercise-12.3 train: cyclone roughing + baghouse polishing ---------
print()
print("=" * 72)
print("3b. Cyclone + baghouse train — Hameed pp. 408-409 (kaolin from air)")
print("=" * 72)
PSD_KAOLIN = [{"d_microns": dm, "mass_frac": w / 100.0} for dm, w in
              [(22.23e-3, 0.03), (0.1209, 0.80), (0.6577, 7.05),
               (3.577, 24.20), (19.46, 35.83), (105.8, 24.20),
               (575.6, 7.05), (3131.0, 0.80), (17030.0, 0.03)]]
n_air = (10000.0 / 3600.0) / MW_AIR
n_k = (2000.0 / 3600.0) / MW_C                   # carbon stands in for kaolin
zk = n_k / (n_air + n_k)
train = Flowsheet(
    components=[Component("nitrogen"), Component("oxygen"), Component("carbon")],
    property_package="thermo:PR")
train.add(Cyclone("CYC", {"solids": "carbon", "particle_density": 2600.0,
                          "geometry": "Stairmand_HE", "body_diameter": 0.6559,
                          "psd": PSD_KAOLIN}))
train.add(BaghouseFilter("BH", {"solids": "carbon"}))
train.feed("F", "CYC:gas_in", T=303.15, P=5e5, molar_flow=n_air + n_k,
           z={"nitrogen": (1 - zk) * 0.79, "oxygen": (1 - zk) * 0.21, "carbon": zk})
train.connect("G", "CYC:gas_out", "BH:gas_in")
train.connect("S1", "CYC:solids_out", None)
train.connect("G2", "BH:gas_out", None)
train.connect("S2", "BH:solids_out", None)
train_rep = train.solve()


def dust_kg_h(stream) -> float:
    return stream.molar_flow * stream.z.get("carbon", 0.0) * MW_C * 3600.0


print(f"  cyclone efficiency  {train.units['CYC'].design['overall_efficiency'] * 100:.2f} % "
      f"(book spec: 75 %)")
print(f"  dust past cyclone   {dust_kg_h(train.streams['G']):.1f} kg/h   (book: 500 kg/h)")
print(f"  final emission      {dust_kg_h(train.streams['G2']):.2f} kg/h   "
      f"(book: 0.0 — HYSYS assumes 100 % baghouse capture; ours is 99.9 %)")
print("  (the PSD is a unit param and does not propagate: the baghouse uses "
      "its bulk efficiency, not the cyclone-depleted size distribution)")

# -- costing --------------------------------------------------------------------
print()
print("=" * 72)
print("Equipment costing (order-of-magnitude solids correlations; see data.py)")
print("=" * 72)
for fs, rep in ((cyc_fs, cyc_rep), (rvf_fs, rvf_rep), (bh_fs, bh_rep)):
    pp = make_package(fs.property_package, fs.component_ids)
    for size in size_flowsheet(fs, rep, pp):
        cost = cost_equipment(size)
        print(f"  {size.unit_id:5s} {size.equipment_type:22s} "
              f"{size.attribute:9.3f} {size.attribute_name:10s} "
              f"x{size.quantity}  Cbm = ${cost.bare_module:>10,.0f}  ({cost.year})")
        for note in size.notes:
            print(f"        - {note}")
