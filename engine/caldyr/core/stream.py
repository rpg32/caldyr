from dataclasses import dataclass, field, replace

_Z_TOL = 1e-6


@dataclass
class Stream:
    """A material stream. SI units (K, Pa, mol/s, J/mol). State is defined by
    (T or H, P, molar_flow, z); the property package resolves the rest."""
    id: str
    components: list[str]
    T: float | None = None
    P: float | None = None
    molar_flow: float | None = None
    z: dict[str, float] = field(default_factory=dict)
    # resolved by the property package after solve:
    H: float | None = None
    phase: str | None = None
    vapor_fraction: float | None = None

    def require_state(self) -> tuple[float, float, float]:
        """Return (T, P, molar_flow) asserting all are specified. Unit ops call
        this at solve time, where an unspecified stream is a wiring error."""
        if self.T is None or self.P is None or self.molar_flow is None:
            raise ValueError(
                f"stream {self.id!r} is underspecified at solve "
                f"(T={self.T}, P={self.P}, molar_flow={self.molar_flow})"
            )
        return self.T, self.P, self.molar_flow

    def normalized_z(self) -> dict[str, float]:
        """Composition rescaled to sum exactly 1. Raises if it sums to <= 0."""
        total = sum(self.z.values())
        if total <= 0.0:
            raise ValueError(f"stream {self.id!r} composition sums to {total}")
        return {k: v / total for k, v in self.z.items()}

    def validate(self) -> None:
        """Sanity-check a feed spec: composition present and ~normalized,
        non-negative flow. Raises ValueError with a diagnostic (no silent fail)."""
        if not self.z:
            raise ValueError(f"stream {self.id!r} has no composition")
        total = sum(self.z.values())
        if abs(total - 1.0) > _Z_TOL:
            raise ValueError(
                f"stream {self.id!r} composition sums to {total:.6f}, expected 1.0"
            )
        if self.molar_flow is not None and self.molar_flow < 0:
            raise ValueError(f"stream {self.id!r} has negative molar_flow {self.molar_flow}")

    def with_(self, **overrides) -> "Stream":
        """Return a copy with fields overridden (e.g. ``s.with_(T=350.0)``)."""
        return replace(self, **overrides)


@dataclass
class EnergyStream:
    """A pure energy stream carrying a duty in watts."""
    id: str
    duty: float | None = None     # W
