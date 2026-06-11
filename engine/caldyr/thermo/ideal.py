"""Back-compat shim.

The M0 property package is implemented in :mod:`caldyr.thermo.thermo_pkg` as a
real cubic-EOS wrapper around `thermo` (the stub here originally recommended
exactly that). Import :class:`ThermoPackage` from there or from the package root.
"""
from .thermo_pkg import ThermoPackage  # noqa: F401
