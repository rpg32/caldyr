"""Shared base for property packages backed by a `thermo` flasher object.

Both the cubic-EOS package (:mod:`caldyr.thermo.thermo_pkg`) and the
activity-coefficient package (:mod:`caldyr.thermo.activity_pkg`) wrap a
`thermo` flasher that answers the same `flash(...)` calls. Everything except how
that flasher is built is identical, and lives here.
"""
from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np

from .base import PhaseResult, ThreePhaseResult

_VF_EPS = 1e-9
_R = 8.314462618          # J/mol/K
_P_REF = 1e5             # Pa, standard state for ideal-gas activities


def _aligned(components: tuple[str, ...], values, label: str) -> list[float]:
    """Coerce a thermo per-component list to floats aligned to ``components``,
    defaulting a missing entry to 0 with a warning."""
    out: list[float] = []
    for comp, value in zip(components, values):
        if value is None:
            warnings.warn(
                f"no ideal-gas {label} for {comp!r}; using 0 J/mol. Reaction "
                f"calculations involving it will be wrong.",
                stacklevel=3,
            )
            out.append(0.0)
        else:
            out.append(float(value))
    return out


def formation_props(components: tuple[str, ...], hfgs, gfgs) -> tuple[list[float], list[float]]:
    """Ideal-gas formation enthalpy and Gibbs energy (J/mol) aligned to
    ``components``."""
    return _aligned(components, hfgs, "formation enthalpy"), \
        _aligned(components, gfgs, "formation Gibbs energy")


class FlasherPackage:
    """Implements the `PropertyPackage` protocol over a `thermo` flasher.

    Subclasses set ``self.components`` (ordered) and ``self._flasher`` before
    (or by calling) ``_init``.
    """

    components: list[str]
    _flasher: Any
    _hf: list[float]      # ideal-gas formation enthalpy per component, J/mol
    _gf: list[float]      # ideal-gas formation Gibbs energy per component, J/mol
    _flasher_3p: Any      # lazily built three-phase (VLN) flasher, or None
    _mws: list[float] | None   # molar masses (kg/mol) for liquid-density ordering

    def _init(self, components: list[str], flasher: Any,
              hf: list[float], gf: list[float]) -> None:
        if not components:
            raise ValueError(f"{type(self).__name__} requires at least one component")
        self.components = list(components)
        self._flasher = flasher
        self._hf = list(hf)
        self._gf = list(gf)
        self._flasher_3p = None
        self._mws = None

    # -- enthalpy basis ----------------------------------------------------
    # `thermo`'s flash H() is sensible-only (zero for each pure ideal gas at
    # 298.15 K), so heats of reaction would vanish. We shift every enthalpy to a
    # formation-inclusive absolute basis by adding the composition-weighted
    # ideal-gas formation enthalpy. Because the shift is linear in composition,
    # it cancels in any balance where composition is conserved (mixers, heaters,
    # pumps...) yet correctly carries the heat of reaction across a reactor.
    def _hf_mix(self, zs: list[float]) -> float:
        return sum(zi * hi for zi, hi in zip(zs, self._hf))

    # -- internals ---------------------------------------------------------
    def _zs(self, z: dict[str, float]) -> list[float]:
        """Project a composition dict onto this package's ordered components,
        normalizing to sum 1. Raises on unknown components or empty totals."""
        unknown = set(z) - set(self.components)
        if unknown:
            raise ValueError(f"composition has components not in package: {sorted(unknown)}")
        zs = [float(z.get(c, 0.0)) for c in self.components]
        total = sum(zs)
        if total <= 0.0:
            raise ValueError(f"composition sums to {total}; expected > 0")
        return [zi / total for zi in zs]

    @staticmethod
    def _phase_of(vf: float) -> str:
        if vf <= _VF_EPS:
            return "liquid"
        if vf >= 1.0 - _VF_EPS:
            return "vapor"
        return "VLE"

    def _comp(self, zs) -> dict[str, float]:
        return {c: float(v) for c, v in zip(self.components, zs)}

    def _result(self, r) -> PhaseResult:
        bulk_zs = list(r.zs)
        vf, y, hv, x, hl = self._split_phases(r)
        return PhaseResult(
            T=float(r.T),
            P=float(r.P),
            H=float(r.H()) + self._hf_mix(bulk_zs),
            phase=self._phase_of(vf),
            vapor_fraction=vf,
            x=x,
            y=y,
            H_liquid=hl,
            H_vapor=hv,
        )

    def _split_phases(self, r):
        """Classify a flash result into a vapor and a (possibly lumped) liquid by
        molar volume — the *least dense* phase is the vapor. This is robust where
        thermo's gas/liquid labelling is not: at high pressure it can tag a
        light, vapor-like phase as a "liquid", which would otherwise hide the
        overhead vapor of a separator. Returns
        ``(vapor_fraction, y, H_vapor, x, H_liquid)`` with absent phases as None.
        """
        phases = list(r.phases)
        betas = [float(b) for b in r.betas]
        if len(phases) == 1:
            ph = phases[0]
            zs = list(ph.zs)
            comp = self._comp(zs)
            h = float(ph.H()) + self._hf_mix(zs)
            if getattr(r, "gas", None) is not None and not getattr(r, "liquids", None):
                return 1.0, comp, h, None, None
            return 0.0, None, None, comp, h

        # Multiphase: vapor = single least-dense phase; lump the rest as liquid.
        vols = [float(ph.V()) for ph in phases]
        vidx = max(range(len(phases)), key=lambda i: vols[i])
        vapor = phases[vidx]
        vf = betas[vidx]
        y = self._comp(list(vapor.zs))
        hv = float(vapor.H()) + self._hf_mix(list(vapor.zs))

        liq_idx = [i for i in range(len(phases)) if i != vidx]
        beta_liq = sum(betas[i] for i in liq_idx) or 1.0
        x_zs = [
            sum(betas[i] * float(phases[i].zs[k]) for i in liq_idx) / beta_liq
            for k in range(len(self.components))
        ]
        h_liq = sum(betas[i] * float(phases[i].H()) for i in liq_idx) / beta_liq
        x = self._comp(x_zs)
        hl = h_liq + self._hf_mix(x_zs)
        return vf, y, hv, x, hl

    # -- PropertyPackage protocol -----------------------------------------
    def enthalpy(self, T: float, P: float, z: dict[str, float]) -> float:
        zs = self._zs(z)
        return float(self._flasher.flash(T=T, P=P, zs=zs).H()) + self._hf_mix(zs)

    def entropy(self, T: float, P: float, z: dict[str, float]) -> float:
        return float(self._flasher.flash(T=T, P=P, zs=self._zs(z)).S())

    def volume(self, T: float, P: float, z: dict[str, float]) -> float:
        """Bulk molar volume, m^3/mol."""
        return float(self._flasher.flash(T=T, P=P, zs=self._zs(z)).V())

    def flash_pt(self, T: float, P: float, z: dict[str, float]) -> PhaseResult:
        return self._result(self._flasher.flash(T=T, P=P, zs=self._zs(z)))

    def flash_ph(self, P: float, H: float, z: dict[str, float]) -> PhaseResult:
        # H arrives on the absolute (formation-inclusive) basis; convert back to
        # thermo's sensible reference before flashing.
        zs = self._zs(z)
        return self._result(self._flasher.flash(P=P, H=H - self._hf_mix(zs), zs=zs))

    def flash_ps(self, P: float, S: float, z: dict[str, float]) -> PhaseResult:
        return self._result(self._flasher.flash(P=P, S=S, zs=self._zs(z)))

    def bubble_dew(self, P: float, z: dict[str, float]) -> tuple[float, float]:
        zs = self._zs(z)
        bubble = float(self._flasher.flash(P=P, VF=0.0, zs=zs).T)
        dew = float(self._flasher.flash(P=P, VF=1.0, zs=zs).T)
        return bubble, dew

    def bubble_point(self, P: float, x: dict[str, float]) -> PhaseResult:
        """Saturated-liquid state of ``x`` at ``P``: a single VF=0 flash whose
        result carries the bubble temperature, the incipient vapor composition
        (``y``, so K_i = y_i/x_i) and both saturated phase enthalpies. See
        :meth:`caldyr.thermo.base.PropertyPackage.bubble_point`."""
        return self._result(self._flasher.flash(P=P, VF=0.0, zs=self._zs(x)))

    # -- per-phase properties at an arbitrary (T, P) -------------------------
    # Evaluated directly on the flasher's phase models (no stability analysis,
    # no flash iteration), so a liquid can be evaluated above its bubble point
    # and a vapor below its dew point. This is exactly what energy-balance-
    # driven MESH methods (sum-rates absorbers; Seader, Henley & Roper 3e
    # ch. 10.4) and tray hydraulic sizing need: stage temperatures there are
    # set by heat balances, not by saturation.
    def _phase_models(self) -> tuple[Any, Any]:
        """``(liquid_phase, gas_phase)`` model objects of the wrapped flasher."""
        gas = self._flasher.gas
        liq = getattr(self._flasher, "liquid", None)
        if liq is None:
            liq = self._flasher.liquids[0]
        return liq, gas

    def k_values(self, T: float, P: float, x: dict[str, float],
                 y: dict[str, float]) -> dict[str, float]:
        """Phi-phi K-values K_i = phi_i^L(T,P,x)/phi_i^V(T,P,y) — the liquid
        fugacity coefficient at composition ``x`` over the vapor's at ``y``.
        For the gamma-phi (activity) package the liquid lnphi already folds in
        gamma_i * Psat_i / P, so the same expression is correct there too."""
        liq, gas = self._phase_models()
        lnphi_l = liq.to(T=T, P=P, zs=self._zs(x)).lnphis()
        lnphi_v = gas.to(T=T, P=P, zs=self._zs(y)).lnphis()
        return {c: math.exp(float(ll) - float(lv))
                for c, ll, lv in zip(self.components, lnphi_l, lnphi_v)}

    def enthalpy_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Molar enthalpy (J/mol, formation-inclusive) of ``x`` as a liquid."""
        liq, _ = self._phase_models()
        xs = self._zs(x)
        return float(liq.to(T=T, P=P, zs=xs).H()) + self._hf_mix(xs)

    def enthalpy_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Molar enthalpy (J/mol, formation-inclusive) of ``y`` as a vapor."""
        _, gas = self._phase_models()
        ys = self._zs(y)
        return float(gas.to(T=T, P=P, zs=ys).H()) + self._hf_mix(ys)

    def volume_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Molar volume (m^3/mol) of ``x`` as a liquid at (T, P)."""
        liq, _ = self._phase_models()
        return float(liq.to(T=T, P=P, zs=self._zs(x)).V())

    def volume_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Molar volume (m^3/mol) of ``y`` as a vapor at (T, P)."""
        _, gas = self._phase_models()
        return float(gas.to(T=T, P=P, zs=self._zs(y)).V())

    # -- three-phase (VLLE) flashes -----------------------------------------
    # Built on `thermo`'s FlashVLN (two trial liquids + a gas). Only the
    # cubic-EOS backends implement `_build_3p`; the base raises a clear
    # NotImplementedError so e.g. the NRTL activity package fails loudly
    # instead of silently mislabeling phases (PR/SRK only for now).
    def _build_3p(self) -> tuple[Any, list[float]]:
        """Return ``(three_phase_flasher, molar_masses_kg_per_mol)``."""
        raise NotImplementedError(
            f"{type(self).__name__} does not support three-phase (VLLE) flashes; "
            f"only the cubic-EOS packages do for now — select 'thermo:PR' or "
            f"'thermo:SRK' as the flowsheet property package."
        )

    def _get_3p(self) -> Any:
        if self._flasher_3p is None:
            self._flasher_3p, self._mws = self._build_3p()
        return self._flasher_3p

    def flash_pt_3p(self, T: float, P: float, z: dict[str, float]) -> ThreePhaseResult:
        return self._result_3p(self._get_3p().flash(T=T, P=P, zs=self._zs(z)))

    def flash_ph_3p(self, P: float, H: float, z: dict[str, float]) -> ThreePhaseResult:
        # H arrives on the absolute (formation-inclusive) basis, as in flash_ph.
        zs = self._zs(z)
        return self._result_3p(self._get_3p().flash(P=P, H=H - self._hf_mix(zs), zs=zs))

    def _rho_mass(self, ph: Any) -> float:
        """Mass density of a thermo phase, kg/m^3."""
        assert self._mws is not None
        mw = sum(zi * mi for zi, mi in zip(ph.zs, self._mws))
        return mw / float(ph.V())

    def _result_3p(self, r: Any) -> ThreePhaseResult:
        """Classify a (up to) three-phase flash result into vapor + light liquid
        + heavy liquid. The vapor is thermo's gas phase; the liquids are ordered
        by *mass density* (light = less dense, e.g. organic; heavy = denser,
        e.g. aqueous). Absent phases get beta 0 and None fields, so a fully
        miscible system degrades to two (or one) phases instead of failing.
        """
        betas = [float(b) for b in r.betas]
        beta_v, y, hv = 0.0, None, None
        liqs: list[tuple[Any, float]] = []
        for ph, beta in zip(r.phases, betas):
            if r.gas is not None and ph is r.gas:
                beta_v = beta
                y = self._comp(list(ph.zs))
                hv = float(ph.H()) + self._hf_mix(list(ph.zs))
            else:
                liqs.append((ph, beta))
        liqs.sort(key=lambda t: self._rho_mass(t[0]))

        def liq_fields(entry: tuple[Any, float] | None):
            if entry is None:
                return 0.0, None, None, None
            ph, beta = entry
            zs = list(ph.zs)
            return (beta, self._comp(zs), float(ph.H()) + self._hf_mix(zs),
                    self._rho_mass(ph))

        b_l, x_l, h_l, rho_l = liq_fields(liqs[0] if liqs else None)
        b_h, x_h, h_h, rho_h = liq_fields(liqs[1] if len(liqs) > 1 else None)
        return ThreePhaseResult(
            T=float(r.T), P=float(r.P),
            H=float(r.H()) + self._hf_mix(list(r.zs)),
            beta_vapor=beta_v, beta_light=b_l, beta_heavy=b_h,
            y=y, x_light=x_l, x_heavy=x_h,
            H_vapor=hv, H_light=h_l, H_heavy=h_h,
            rho_light=rho_l, rho_heavy=rho_h,
        )

    # -- reaction equilibrium ---------------------------------------------
    def lnKeq(self, stoich: dict[str, float], T: float) -> float:
        """ln of the ideal-gas equilibrium constant (standard state 1 bar) for a
        reaction ``{component: nu}`` at temperature ``T``.

        Anchored at 298.15 K from formation Gibbs energies, then propagated with
        the van't Hoff / Gibbs-Helmholtz relation using the *temperature-dependent*
        heat of reaction (from the formation-inclusive enthalpies), so the result
        is far better than a constant-ΔH fit:

            ln K(T) = ln K(298) + (1/R) ∫_298^T ΔH_rxn(T') / T'^2 dT'

        (ΔH_rxn < 0 for an exothermic reaction makes the integral negative, so K
        falls with temperature — Le Chatelier.)
        """
        idx = {c: self.components.index(c) for c in stoich}
        dG298 = sum(nu * self._gf[idx[c]] for c, nu in stoich.items())
        ln_k298 = -dG298 / (_R * 298.15)

        def dH_rxn(temp: float) -> float:
            return sum(nu * self.enthalpy(temp, _P_REF, {c: 1.0}) for c, nu in stoich.items())

        temps = np.linspace(298.15, T, 11)
        integrand = np.array([dH_rxn(float(t)) / (t * t) for t in temps])
        # numpy >= 2 renamed trapz -> trapezoid; getattr keeps both stub
        # generations happy (the missing name is a dead branch either way).
        trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")
        integral = float(trapezoid(integrand, temps))
        return ln_k298 + integral / _R
