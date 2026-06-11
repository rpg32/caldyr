from abc import ABC, abstractmethod

from .port import Port
from .stream import EnergyStream, Stream

# What a unit op may emit on a port: a material stream or an energy duty.
PortStream = Stream | EnergyStream


class UnitOp(ABC):
    """Base class for all unit operations. The single solve() contract is what
    lets the sequential-modular solver and the AI tool layer stay simple."""
    id: str
    ports: list[Port]
    params: dict

    def __init__(self, id: str, params: dict | None = None) -> None:
        self.id = id
        self.params = params or {}
        self.ports = self.define_ports()

    @abstractmethod
    def define_ports(self) -> list[Port]:
        """Declare inlet/outlet ports for this unit op."""

    @abstractmethod
    def solve(self, inlets: dict[str, Stream], pp) -> dict[str, PortStream]:
        """Given inlet material streams keyed by port name and a PropertyPackage,
        compute and return outlet streams (material or energy) keyed by port name."""
