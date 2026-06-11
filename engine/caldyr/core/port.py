from dataclasses import dataclass


@dataclass(frozen=True)
class Port:
    name: str                     # e.g. "in1", "out", "duty"
    direction: str                # "inlet" | "outlet"
    kind: str = "material"        # "material" | "energy"
