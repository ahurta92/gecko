from __future__ import annotations

from pathlib import Path


def can_load(path: Path) -> bool:
    """
    Heuristic detection for DALTON outputs (migration fixtures first).
    """
    if not path.exists() or not path.is_dir():
        return False

    # Fixture marker:
    if any(path.glob("*.out")):
        return True

    # Some DALTON runs use DALTON.OUT
    if (path / "DALTON.OUT").exists():
        return True

    return False
