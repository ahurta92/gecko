from __future__ import annotations

from pathlib import Path

from gecko.core.model import Calculation
from gecko.plugins.madness.detect import can_load as madness_can_load
from gecko.plugins.dalton.detect import can_load as dalton_can_load

from gecko.plugins.madness.loader import load as load_madness
from gecko.plugins.dalton.loader import load as load_dalton


def load_calc(path: str | Path) -> Calculation:
    """
    Load and parse a calculation directory.

    Step 2 behavior:
    - Detect whether the directory is MADNESS or DALTON
    - Delegate to the appropriate plugin loader, which parses artifacts
    """
    root = Path(path).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Path does not exist: {root}")

    if root.is_file():
        # For now, enforce directory-only to keep things simple.
        # (We can add file-path support later if you want.)
        raise ValueError(f"Expected a directory, got a file: {root}")

    if madness_can_load(root):
        return load_madness(root)

    if dalton_can_load(root):
        return load_dalton(root)

    raise ValueError(
        "Could not detect calculation type (madness/dalton) from directory. "
        f"Path: {root}"
    )
