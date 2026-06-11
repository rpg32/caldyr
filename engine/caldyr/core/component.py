from dataclasses import dataclass


@dataclass(frozen=True)
class Component:
    """A chemical species. Property backends resolve by id/cas; no property
    data is stored on the component itself."""
    id: str                 # canonical key, e.g. "water"
    name: str = ""
    formula: str | None = None
    cas: str | None = None
