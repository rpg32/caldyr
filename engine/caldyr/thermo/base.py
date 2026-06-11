from dataclasses import dataclass
from typing import Protocol, runtime_checkable


@dataclass
class PhaseResult:
    """Resolved single-flash state. SI units (K, Pa, J/mol).

    Carries the full resolved (T, P) so a flash_ph result alone is enough to
    fully define a stream — the caller never needs a second flash to recover T.
    For a two-phase result, the per-phase compositions and molar enthalpies are
    populated too (``x``/``y`` and ``H_liquid``/``H_vapor``) so a flash drum can
    split the phases; for a single phase they describe the one phase present.
    """
    T: float                   # K
    P: float                   # Pa
    H: float                   # J/mol (bulk)
    phase: str                 # "vapor" | "liquid" | "VLE"
    vapor_fraction: float
    x: dict[str, float] | None = None        # liquid-phase mole fractions
    y: dict[str, float] | None = None        # vapor-phase mole fractions
    H_liquid: float | None = None            # J/mol of the liquid phase
    H_vapor: float | None = None             # J/mol of the vapor phase


@dataclass
class ThreePhaseResult:
    """Resolved three-phase (vapor / light liquid / heavy liquid) flash state.
    SI units (K, Pa, J/mol, kg/m^3).

    The two liquids are identified by *mass density*: the light liquid is the
    less dense one (e.g. the organic phase of a water/hydrocarbon split), the
    heavy liquid the denser (e.g. the aqueous phase). Absent phases carry a
    zero ``beta`` and ``None`` composition/enthalpy/density, so a fully
    miscible system degrades gracefully to a two-phase (or single-phase)
    result rather than failing.
    """
    T: float                   # K
    P: float                   # Pa
    H: float                   # J/mol (bulk, formation-inclusive)
    beta_vapor: float          # molar phase fractions; sum to 1
    beta_light: float
    beta_heavy: float
    y: dict[str, float] | None = None         # vapor mole fractions
    x_light: dict[str, float] | None = None   # light-liquid mole fractions
    x_heavy: dict[str, float] | None = None   # heavy-liquid mole fractions
    H_vapor: float | None = None              # J/mol of each phase
    H_light: float | None = None
    H_heavy: float | None = None
    rho_light: float | None = None            # kg/m^3 (mass density orders the liquids)
    rho_heavy: float | None = None


@runtime_checkable
class PropertyPackage(Protocol):
    """The boundary that keeps the engine out of re-deriving thermodynamics.
    Implementations wrap `thermo`/`CoolProp`/`Cantera`. Selected per-flowsheet.

    Compositions are passed as ``{component_id: mole_fraction}`` dicts; the
    implementation owns its component ordering and normalizes internally.
    """
    def enthalpy(self, T: float, P: float, z: dict[str, float]) -> float: ...
    def entropy(self, T: float, P: float, z: dict[str, float]) -> float: ...
    def volume(self, T: float, P: float, z: dict[str, float]) -> float: ...
    def flash_pt(self, T: float, P: float, z: dict[str, float]) -> PhaseResult: ...
    def flash_ph(self, P: float, H: float, z: dict[str, float]) -> PhaseResult: ...
    def flash_ps(self, P: float, S: float, z: dict[str, float]) -> PhaseResult: ...
    def bubble_dew(self, P: float, z: dict[str, float]) -> tuple[float, float]: ...
    def bubble_point(self, P: float, x: dict[str, float]) -> PhaseResult:
        """Saturated-liquid state of liquid composition ``x`` at pressure ``P``.

        Returns a PhaseResult at the bubble temperature with ``x`` the liquid,
        ``y`` the *incipient* vapor composition (so ``K_i = y_i / x_i`` at the
        bubble point) and ``H_liquid``/``H_vapor`` the saturated per-phase molar
        enthalpies. One call gives a stage's temperature, K-values and both
        phase enthalpies — the workhorse of tray-by-tray (MESH) columns.
        """
        ...
    # Per-phase properties at an *arbitrary* (T, P) — not restricted to the
    # saturation locus like bubble_point. These are what energy-balance-driven
    # MESH methods (sum-rates absorbers, Seader 3e ch. 10.4) and tray hydraulic
    # sizing need: the stage temperature there is set by the heat balance, so
    # K-values and phase enthalpies must be evaluable off saturation.
    def k_values(self, T: float, P: float, x: dict[str, float],
                 y: dict[str, float]) -> dict[str, float]:
        """K_i = phi_i^L(T, P, x) / phi_i^V(T, P, y) — phi-phi K-values of a
        liquid of composition ``x`` against a vapor of composition ``y``."""
        ...
    def enthalpy_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Molar enthalpy (J/mol) of composition ``x`` forced to the liquid
        phase at (T, P), on the formation-inclusive basis."""
        ...
    def enthalpy_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Molar enthalpy (J/mol) of composition ``y`` forced to the vapor
        phase at (T, P), on the formation-inclusive basis."""
        ...
    def volume_liquid(self, T: float, P: float, x: dict[str, float]) -> float:
        """Molar volume (m^3/mol) of composition ``x`` forced to the liquid
        phase at (T, P)."""
        ...
    def volume_vapor(self, T: float, P: float, y: dict[str, float]) -> float:
        """Molar volume (m^3/mol) of composition ``y`` forced to the vapor
        phase at (T, P)."""
        ...
    # Three-phase (VLLE) flashes. Only the cubic-EOS backends implement these;
    # others raise NotImplementedError with a clear message (PR/SRK only for now).
    def flash_pt_3p(self, T: float, P: float, z: dict[str, float]) -> ThreePhaseResult: ...
    def flash_ph_3p(self, P: float, H: float, z: dict[str, float]) -> ThreePhaseResult: ...
