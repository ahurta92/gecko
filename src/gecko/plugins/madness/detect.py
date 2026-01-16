from __future__ import annotations

from pathlib import Path


def can_load(path: Path) -> bool:
    """
    Heuristic detection for MADNESS outputs (migration fixtures first).

    We will make this more robust once we inspect real production run layouts.
    """
    if not path.exists() or not path.is_dir():
        return False

    # Fixture markers you currently have:
    # - n12_mad_output.json
    # - ... or any *_mad_output.json
    if any(path.glob("*.calc_info.json")):
        return True

    # Common future marker candidates you might have in real runs:
    # - responses/metadata.json
    if (path / "responses" / "metadata.json").exists():
        return True

    return False
