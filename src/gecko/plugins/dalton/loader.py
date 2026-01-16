# src/gecko/plugins/dalton/loader.py
from __future__ import annotations
from pathlib import Path
from gecko.core.model import Calculation
from gecko.plugins.dalton.detect import can_load

def load(path: Path) -> Calculation:
    root = path.expanduser().resolve()
    if not can_load(root):
        raise ValueError(f"Not a DALTON run directory: {root}")
    return Calculation(code="dalton", root=root, artifacts={}, data={}, meta={})
