"""M13 tests: the PipeSegment unit op (Darcy-Weisbach + Churchill friction).

References:

* **Hameed, *Chemical Process Simulations using Aspen HYSYS* (Wiley 2025),
  sec. 7.1 worked example** — 700 kg/h of water at 60 F pumped through a piping
  network of 1-in and 0.5-in schedule-40 cast-iron pipe (Fig. 7.1). The PIPE-101
  pressure profile is tabulated in the book's Fig. 7.13 (HYSYS, ASME Steam):
  inlet (point C) 568.284 psia -> outlet (point D) 410.837 psia, with friction
  gradients of 0.675 psi per 150 ft in the 1-in pipe and 5.040 psi per 75 ft
  in the 0.5-in pipe. ``test_book_pipe101_*`` reproduce that table.
  Measured deltas of this engine (coolprop:Water + Churchill) vs the book,
  asserted with headroom below:
    - P at point D:               -0.04%   (asserted 0.2%)
    - total PIPE-101 dP:          +0.2%    (asserted 1%)
    - 1-in   friction gradient:   +0.3%    (asserted 2%)
    - 0.5-in friction gradient:   +0.2%    (asserted 2%)
* **Crane TP-410** (Flow of Fluids, 2009): fitting K factors (90-deg std elbow
  K = 30 fT, half-open gate valve K = 160 fT, with fT = 0.023 / 0.027 for 1-in /
  0.5-in pipe, A-26..A-29) and water properties at 60 F (rho = 62.371 lb/ft^3 =
  999.0 kg/m^3, mu = 1.12 cP, A-2/A-6).
* **Swamee & Jain**, J. Hydraulics Div. ASCE 102(5), 657-664 (1976) (also
  Perry's 8e sec. 6): explicit Colebrook approximation (+/-1%) used for the
  independent hand Darcy-Weisbach check.
"""
import math

import pytest

from caldyr.core import Component, Flowsheet
from caldyr.economics.costing import cost_equipment
from caldyr.economics.sizing import size_flowsheet
from caldyr.io import from_dict, to_dict
from caldyr.thermo import make_package
from caldyr.unitops import PipeFlowError, PipeSegment
from caldyr.unitops.pipe import churchill_friction_factor

FT = 0.3048               # m / ft
IN = 0.0254               # m / in
PSI = 6894.757            # Pa / psi
G = 9.80665               # m/s^2
T60F = (60.0 - 32.0) / 1.8 + 273.15          # 60 F = 288.706 K
M_H2O = 0.01801528                           # kg/mol
N_700KGH = (700.0 / 3600.0) / M_H2O          # 700 kg/h water, mol/s

# Book pipe geometry (sch-40 cast iron; Hameed Fig. 7.7: ID 1.049 in / 0.622 in,
# roughness 8.497e-4 ft).
D_1IN = 1.049 * IN
D_HALF = 0.622 * IN
ROUGH_CI = 8.497e-4 * FT


def water_pipe(P_in: float, n: float = N_700KGH, **params) -> Flowsheet:
    """One water-filled PipeSegment on the IAPWS steam-tables package."""
    fs = Flowsheet(components=[Component("water")],
                   property_package="coolprop:Water")
    fs.add(PipeSegment("PIPE", params))
    fs.feed("FEED", "PIPE:in1", T=T60F, P=P_in, molar_flow=n, z={"water": 1.0})
    fs.connect("OUT", "PIPE:out", None)
    return fs


# -- book sec. 7.1: the PIPE-101 network (Fig. 7.13) --------------------------
@pytest.fixture(scope="module")
def book_pipe101():
    """PIPE-101 of Hameed Fig. 7.1/7.13, anchored at the book's point-C
    pressure: 450 ft of 1-in pipe rising 300 ft (two std elbows + a union),
    then 375 ft of 0.5-in pipe ending in a half-open gate valve."""
    fs = Flowsheet(components=[Component("water")],
                   property_package="coolprop:Water")
    fs.add(PipeSegment("P1IN", {
        "length": 450 * FT, "diameter": D_1IN, "roughness": ROUGH_CI,
        "elevation_change": 300 * FT,
        # two 90-deg std elbows (K = 30*0.023, Crane) + a screwed union (~0.08)
        "fittings_K": 2 * 30 * 0.023 + 0.08, "segments": 15}))
    fs.add(PipeSegment("PHALF", {
        "length": 375 * FT, "diameter": D_HALF, "roughness": ROUGH_CI,
        # half-open gate valve, K = 160*fT = 160*0.027 (Crane A-27/A-29)
        "fittings_K": 160 * 0.027, "segments": 13}))
    fs.feed("C", "P1IN:in1", T=T60F, P=568.284 * PSI,    # book Fig. 7.13 row 1
            molar_flow=N_700KGH, z={"water": 1.0})
    fs.connect("MID", "P1IN:out", "PHALF:in1")
    fs.connect("D", "PHALF:out", None)
    rep = fs.solve()
    assert rep.converged
    return fs


def test_book_pipe101_point_d_pressure(book_pipe101):
    # Book Fig. 7.13 last row: 410.837 psia at point D (achieved: -0.04%).
    assert book_pipe101.streams["D"].P == pytest.approx(410.837 * PSI, rel=2e-3)
    # Total dP across PIPE-101: 568.284 - 410.837 = 157.447 psi (achieved +0.2%).
    d1 = book_pipe101.units["P1IN"].design
    d2 = book_pipe101.units["PHALF"].design
    assert d1["dP_total"] + d2["dP_total"] == pytest.approx(157.447 * PSI, rel=1e-2)


def test_book_pipe101_friction_gradients(book_pipe101):
    d1 = book_pipe101.units["P1IN"].design
    d2 = book_pipe101.units["PHALF"].design
    # 1-in pipe: book rows 0->150 ft drop 568.284 - 567.609 = 0.675 psi of pure
    # friction (no elevation there); the friction gradient is uniform along the
    # liquid-full pipe (achieved +0.3%).
    grad_1in = d1["dP_friction"] / d1["length"]
    assert grad_1in == pytest.approx(0.675 * PSI / (150 * FT), rel=2e-2)
    # 0.5-in pipe: book rows 450->525 ft drop 436.136 - 431.095 = 5.041 psi per
    # 75 ft (achieved +0.2%).
    grad_half = d2["dP_friction"] / d2["length"]
    assert grad_half == pytest.approx(5.041 * PSI / (75 * FT), rel=2e-2)
    # Elevation term of the riser: rho*g*dz, ~130 psi of the 157.
    assert d1["dP_elevation"] == pytest.approx(
        d1["density"] * G * 300 * FT, rel=5e-3)
    # Both lines are turbulent liquid-full flow (book: "Liquid Only").
    assert d1["flow_regime"] == d2["flow_regime"] == "turbulent"
    assert d1["phase"] == d2["phase"] == "liquid"


def test_book_pipe101_profile_shape(book_pipe101):
    d1 = book_pipe101.units["P1IN"].design
    assert len(d1["L_profile"]) == len(d1["P_profile"]) == 15 + 1
    assert d1["L_profile"][0] == 0.0
    assert d1["L_profile"][-1] == pytest.approx(450 * FT)
    assert d1["P_profile"][0] == pytest.approx(568.284 * PSI)
    # Friction + gravity both fight the flow up the riser: monotone decreasing.
    assert all(a > b for a, b in zip(d1["P_profile"], d1["P_profile"][1:]))


# -- independent Darcy-Weisbach hand check ------------------------------------
def test_hand_darcy_weisbach_commercial_steel():
    """100 m of 1-in commercial-steel pipe (default roughness 4.5e-5 m, Crane
    TP-410 / Perry's 8e Table 6-1), 700 kg/h water at 60 F. Hand calculation
    with Crane water properties (rho = 999.0 kg/m^3, mu = 1.12 cP) and the
    Swamee-Jain explicit friction factor."""
    L = 100.0
    fs = water_pipe(5e5, length=L, diameter=D_1IN)
    fs.solve()
    d = fs.units["PIPE"].design

    # The engine's transport properties match the Crane handbook values.
    assert d["density"] == pytest.approx(999.0, rel=5e-3)
    assert d["viscosity"] == pytest.approx(1.12e-3, rel=4e-2)

    # Hand numbers (fully independent of the engine):
    rho, mu = 999.0, 1.12e-3
    area = math.pi * D_1IN ** 2 / 4.0
    v = (700.0 / 3600.0) / (rho * area)                  # 0.3491 m/s
    re = rho * v * D_1IN / mu                            # ~8.29e3, turbulent
    rel_rough = 4.5e-5 / D_1IN
    f_sj = 0.25 / math.log10(rel_rough / 3.7 + 5.74 / re ** 0.9) ** 2
    dp_hand = f_sj * (L / D_1IN) * rho * v * v / 2.0     # ~8.04 kPa

    assert d["Re"] == pytest.approx(re, rel=4e-2)        # mu databank vs Crane
    assert d["friction_factor"] == pytest.approx(f_sj, rel=3e-2)  # Churchill vs S-J
    assert d["dP_friction"] == pytest.approx(dp_hand, rel=4e-2)
    assert d["dP_elevation"] == 0.0 and d["dP_fittings"] == 0.0
    assert d["dP_total"] == pytest.approx(d["dP_friction"], rel=1e-12)


# -- friction-factor correlation ----------------------------------------------
def test_churchill_laminar_limit_is_64_over_re():
    for re in (10.0, 100.0, 1000.0):
        assert churchill_friction_factor(re, 0.0) == pytest.approx(64.0 / re, rel=1e-4)
        assert churchill_friction_factor(re, 2e-3) == pytest.approx(64.0 / re, rel=1e-4)
    # Approaching transition the turbulent branch starts to contribute (the
    # correlation is smooth across regimes by design): ~0.1% high at Re 2000.
    assert churchill_friction_factor(2000.0, 0.0) == pytest.approx(64.0 / 2000.0,
                                                                   rel=2e-3)


def test_churchill_tracks_colebrook_turbulent():
    # Solve Colebrook iteratively and compare (Churchill is within ~2-3%).
    for re, rr in ((1e4, 0.0), (1e5, 1e-3), (1e6, 1e-4), (5e4, 5e-2)):
        f = 0.02
        for _ in range(60):
            f = (-2.0 * math.log10(rr / 3.7 + 2.51 / (re * math.sqrt(f)))) ** -2
        assert churchill_friction_factor(re, rr) == pytest.approx(f, rel=3e-2)


def test_churchill_rejects_nonpositive_re():
    with pytest.raises(ValueError, match="Reynolds"):
        churchill_friction_factor(0.0, 1e-3)
    with pytest.raises(ValueError, match="Reynolds"):
        churchill_friction_factor(-100.0, 1e-3)


# -- structural behavior -------------------------------------------------------
def test_friction_drop_scales_linearly_with_length():
    fs1 = water_pipe(10e5, length=50.0, diameter=D_1IN)
    fs2 = water_pipe(10e5, length=100.0, diameter=D_1IN)
    fs1.solve()
    fs2.solve()
    dp1 = fs1.units["PIPE"].design["dP_friction"]
    dp2 = fs2.units["PIPE"].design["dP_friction"]
    assert dp2 == pytest.approx(2.0 * dp1, rel=1e-3)     # liquid: ~incompressible


def test_turbulent_friction_scales_roughly_with_v_squared():
    # Doubling the flow doubles v; dP ~ f*v^2 with f falling slightly with Re
    # (and approaching constant in the fully-rough limit) -> ratio in (3.3, 4.0].
    fs1 = water_pipe(10e5, length=50.0, diameter=D_HALF, roughness=ROUGH_CI)
    fs2 = water_pipe(10e5, n=2 * N_700KGH, length=50.0, diameter=D_HALF,
                     roughness=ROUGH_CI)
    fs1.solve()
    fs2.solve()
    assert fs1.units["PIPE"].design["flow_regime"] == "turbulent"
    ratio = (fs2.units["PIPE"].design["dP_friction"]
             / fs1.units["PIPE"].design["dP_friction"])
    assert 3.3 < ratio <= 4.0


def test_laminar_regime_reports_f_64_over_re():
    # 36 kg/h of water in a 1-in pipe: Re = 4*mdot/(pi*D*mu) ~ 4e2, laminar.
    n = (36.0 / 3600.0) / M_H2O
    fs = water_pipe(5e5, n=n, length=20.0, diameter=D_1IN)
    fs.solve()
    d = fs.units["PIPE"].design
    assert d["flow_regime"] == "laminar"
    assert d["Re"] < 2300.0
    assert d["friction_factor"] == pytest.approx(64.0 / d["Re"], rel=1e-4)


def test_elevation_term_is_rho_g_dz_and_downhill_recovers_pressure():
    flat = water_pipe(5e5, length=10.0, diameter=D_1IN)
    up = water_pipe(5e5, length=10.0, diameter=D_1IN, elevation_change=10.0)
    down = water_pipe(5e5, length=10.0, diameter=D_1IN, elevation_change=-10.0)
    flat.solve()
    up.solve()
    down.solve()
    du = up.units["PIPE"].design
    assert du["dP_elevation"] == pytest.approx(du["density"] * G * 10.0, rel=1e-3)
    assert flat.units["PIPE"].design["dP_elevation"] == 0.0
    # A short downhill liquid line gains more head than friction eats.
    assert down.streams["OUT"].P > 5e5
    assert down.units["PIPE"].design["dP_total"] < 0.0


def test_fittings_term_is_k_velocity_heads():
    k = 7.5
    fs = water_pipe(5e5, length=1.0, diameter=D_1IN, fittings_K=k)
    fs.solve()
    d = fs.units["PIPE"].design
    assert d["dP_fittings"] == pytest.approx(
        k * d["density"] * d["velocity"] ** 2 / 2.0, rel=1e-3)


def test_rougher_pipe_drops_more_pressure():
    smooth = water_pipe(5e5, length=100.0, diameter=D_1IN, roughness=0.0)
    rough = water_pipe(5e5, length=100.0, diameter=D_1IN, roughness=ROUGH_CI)
    smooth.solve()
    rough.solve()
    assert (rough.units["PIPE"].design["dP_friction"]
            > smooth.units["PIPE"].design["dP_friction"])


def test_gas_line_accelerates_as_it_expands():
    """Methane at 3 bar down 2 km of 2-in line: the marching solution picks up
    the falling density, so the outlet velocity is ~P_in/P_out times the inlet
    velocity and the friction gradient steepens downstream."""
    fs = Flowsheet(components=[Component("methane")], property_package="thermo:PR")
    fs.add(PipeSegment("G", {"length": 2000.0, "diameter": 0.05, "segments": 30}))
    fs.feed("F", "G:in1", T=300.0, P=3e5, molar_flow=3.0, z={"methane": 1.0})
    fs.connect("OUT", "G:out", None)
    fs.solve()
    d = fs.units["G"].design
    out = fs.streams["OUT"]
    assert d["phase"] == "vapor"
    assert d["velocity_out"] == pytest.approx(
        d["velocity"] * 3e5 / out.P, rel=5e-2)            # ~ideal-gas expansion
    # The last marching step loses more pressure than the first.
    p = d["P_profile"]
    assert (p[-2] - p[-1]) > (p[0] - p[1])
    # Isothermal model: T preserved, enthalpy re-flashed at outlet P.
    assert out.T == pytest.approx(300.0)


# -- typed errors --------------------------------------------------------------
def test_two_phase_inlet_raises_typed_error():
    # 50/50 propane/n-decane at 350 K / 5 bar flashes to VF ~ 0.36 under PR.
    fs = Flowsheet(components=[Component("propane"), Component("n-decane")],
                   property_package="thermo:PR")
    fs.add(PipeSegment("P", {"length": 10.0, "diameter": 0.05}))
    fs.feed("F", "P:in1", T=350.0, P=5e5, molar_flow=10.0,
            z={"propane": 0.5, "n-decane": 0.5})
    fs.connect("OUT", "P:out", None)
    with pytest.raises(PipeFlowError, match="two-phase"):
        fs.solve()
    assert issubclass(PipeFlowError, ValueError)          # engine convention


def test_liquid_flashing_along_the_line_raises():
    # Pure liquid propane at 300 K just above its ~10-bar (PR) vapor pressure:
    # the friction drop pushes the line across the saturation curve, where a
    # pure component jumps straight from VF 0 to VF 1 — the model must refuse
    # rather than silently continue with vapor densities.
    fs = Flowsheet(components=[Component("propane")], property_package="thermo:PR")
    fs.add(PipeSegment("P", {"length": 15.0, "diameter": 0.025, "segments": 5}))
    fs.feed("F", "P:in1", T=300.0, P=10.2e5, molar_flow=20.0, z={"propane": 1.0})
    fs.connect("OUT", "P:out", None)
    with pytest.raises(PipeFlowError, match="changes phase"):
        fs.solve()


def test_pressure_drop_exhausting_inlet_raises():
    # 1 t/h of water through 1 km of 1-cm line at 1.5 bar: impossible.
    n = (1000.0 / 3600.0) / M_H2O
    fs = water_pipe(1.5e5, n=n, length=1000.0, diameter=0.01)
    with pytest.raises(PipeFlowError, match="exceeds the inlet"):
        fs.solve()


def test_parameter_validation():
    with pytest.raises(ValueError, match="length"):
        water_pipe(5e5, diameter=D_1IN).solve()           # missing length
    with pytest.raises(ValueError, match="positive"):
        water_pipe(5e5, length=-1.0, diameter=D_1IN).solve()
    with pytest.raises(ValueError, match="positive"):
        water_pipe(5e5, length=10.0, diameter=0.0).solve()
    with pytest.raises(ValueError, match="roughness"):
        water_pipe(5e5, length=10.0, diameter=D_1IN, roughness=-1e-5).solve()
    with pytest.raises(ValueError, match="segments"):
        water_pipe(5e5, length=10.0, diameter=D_1IN, segments=0).solve()


# -- conservation, .flow round-trip, economics ---------------------------------
def test_mass_conserved_and_isothermal():
    fs = water_pipe(5e5, length=100.0, diameter=D_1IN, elevation_change=5.0)
    fs.solve()
    feed, out = fs.streams["FEED"], fs.streams["OUT"]
    assert out.molar_flow == pytest.approx(feed.molar_flow, rel=1e-12)
    assert out.z == pytest.approx(feed.z)
    assert out.T == pytest.approx(feed.T)
    assert out.P < feed.P


def test_flow_roundtrip():
    fs = water_pipe(5e5, length=100.0, diameter=D_1IN,
                    roughness=ROUGH_CI, elevation_change=5.0,
                    fittings_K=1.38, segments=12)
    fs.solve()
    doc = to_dict(fs)
    assert doc == to_dict(from_dict(doc))
    fs2 = from_dict(doc)
    fs2.solve()
    assert fs2.streams["OUT"].P == pytest.approx(fs.streams["OUT"].P, rel=1e-9)


def test_sizer_and_costing_give_finite_installed_cost():
    L = 450 * FT
    fs = water_pipe(40e5, length=L, diameter=D_1IN, roughness=ROUGH_CI,
                    elevation_change=300 * FT)
    rep = fs.solve()
    pp = make_package(fs.property_package, fs.component_ids)
    sizes = size_flowsheet(fs, rep, pp)
    assert len(sizes) == 1
    size = sizes[0]
    assert size.equipment_type == "pipe"
    assert size.attribute_name == "L_m_x_D_m^0.74"
    assert size.attribute == pytest.approx(L * D_1IN ** 0.74, rel=1e-12)

    cost = cost_equipment(size)
    assert math.isfinite(cost.bare_module) and cost.bare_module > 0.0
    # The Sinnott correlation is an installed cost: Fbm = 1, so the bare-module
    # cost equals the (escalated) purchased cost, ~1416 USD2001 per m*D^0.74.
    assert cost.bare_module == pytest.approx(cost.purchased, rel=1e-12)
    cp0_2001 = 10 ** 3.1510 * size.attribute
    assert cost.purchased > cp0_2001                       # CEPCI-escalated up
