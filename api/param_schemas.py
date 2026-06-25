"""Per-unit parameter schemas served by ``GET /unit-types`` so the GUI can build
typed parameter forms (number / boolean / select) with units, defaults, ranges,
help text, and conditional applicability.

The ENGINE remains the single source of truth for validation — every unit op
validates its own params in ``solve``/``_read_params`` and raises a typed error.
These schemas are a UI convenience (labels, widgets, defaults); a unit not listed
here simply advertises no param metadata (the GUI falls back to a generic editor).

A param entry:
  ``name``        the engine param key (round-trips in ``.flow``)
  ``label``       human label for the form field
  ``type``        "number" | "int" | "boolean" | "select" | "string"
  ``default``     default value (omitted => no default; usually means required)
  ``required``    bool (default False)
  ``unit``        SI unit suffix shown in the form (engine is SI throughout)
  ``min`` / ``max``  numeric bounds (advisory; the engine enforces the hard ones)
  ``options``     list of allowed values for ``type == "select"``
  ``requires``    ``{param: value}`` — this field only APPLIES when another param
                  holds the given value (e.g. ``condenser_T`` requires
                  ``decant_condenser == True``); the GUI hides it otherwise and
                  treats ``required`` as conditional on the same predicate
  ``help``        one-line description
"""
from __future__ import annotations

_COLUMN_METHODS = ["bubble_point", "sum_rates", "inside_out", "naphtali_sandholm"]

RIGOROUS_COLUMN_PARAMS: list[dict] = [
    {"name": "n_stages", "label": "Stages", "type": "int", "required": True,
     "min": 3, "help": "Equilibrium stages including condenser + reboiler."},
    {"name": "feed_stage", "label": "Feed stage", "type": "int", "min": 1,
     "help": "Feed stage (1 = top); or use 'feeds' for multiple feeds."},
    {"name": "reflux_ratio", "label": "Reflux ratio", "type": "number", "min": 0.0,
     "help": "Molar reflux ratio L/D (the organic-layer ratio in decant mode)."},
    {"name": "distillate_rate", "label": "Distillate rate", "type": "number",
     "unit": "mol/s", "min": 0.0, "help": "Distillate molar flow spec."},
    {"name": "P", "label": "Pressure", "type": "number", "unit": "Pa", "min": 0.0,
     "help": "Top-stage pressure (default: the feed pressure)."},
    {"name": "method", "label": "Solver method", "type": "select",
     "options": _COLUMN_METHODS, "default": "naphtali_sandholm",
     "help": "MESH solution method; decant mode requires naphtali_sandholm."},
    {"name": "reboiled", "label": "Reboiled", "type": "boolean", "default": True,
     "help": "Column has a reboiler (required for decant mode)."},
    {"name": "partial_condenser", "label": "Partial condenser", "type": "boolean",
     "default": False, "help": "Vapor distillate (mutually exclusive with decant)."},
    {"name": "max_iter", "label": "Max iterations", "type": "int", "min": 1,
     "help": "Solver iteration cap."},
    # -- integrated decanting condenser (heteroazeotropic entrainer column) -----
    {"name": "decant_condenser", "label": "Decant condenser", "type": "boolean",
     "default": False,
     "help": "Integrate a DECANTING condenser: the overhead settles into an "
             "organic + aqueous layer; the organic layer is refluxed in full and "
             "the aqueous layer is the distillate (anhydrous-ethanol entrainer "
             "columns, Hameed §9.5.6)."},
    {"name": "condenser_T", "label": "Condenser T", "type": "number", "unit": "K",
     "min": 150.0, "max": 1500.0, "required": True,
     "requires": {"decant_condenser": True},
     "help": "Temperature the overhead is condensed + decanted at. Required when "
             "the decanting condenser is on."},
    {"name": "reflux_layer", "label": "Reflux layer", "type": "select",
     "options": ["organic", "aqueous"], "default": "organic",
     "requires": {"decant_condenser": True},
     "help": "Which settled layer is refluxed in full (organic = the "
             "entrainer-rich layer for a cyclohexane entrainer)."},
]

PARAM_SCHEMAS: dict[str, list[dict]] = {
    "RigorousColumn": RIGOROUS_COLUMN_PARAMS,
}
