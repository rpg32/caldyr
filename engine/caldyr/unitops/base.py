"""Unit-op registry so `.flow` "type" strings map to classes, and the AI layer
can enumerate available units."""
REGISTRY: dict[str, type] = {}


def register(name: str):
    def deco(cls):
        REGISTRY[name] = cls
        return cls
    return deco
