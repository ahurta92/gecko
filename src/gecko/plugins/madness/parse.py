from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from gecko.core.model import Calculation


def parse_run(calc: Calculation) -> None:
    """
    Populate calc.data/meta from MADQC-style MADNESS artifacts.
    Migration-first: keep dicts as dicts; normalize later.
    """
    root = calc.root

    # 1) calc_info is the anchor for MADQC-style runs
    ci_path = calc.artifacts.get("calc_info_json")
    if ci_path and ci_path.exists():
        calc_info = _read_json(ci_path)
        calc.data["calc_info"] = calc_info
        calc.meta["style"] = "madqc"

        # light metadata extraction (best-effort)
        _populate_meta_from_calc_info(calc, calc_info)

    # 2) optional: mad_output_json (older/newer auxiliary)
    mo_path = calc.artifacts.get("mad_output_json")
    if mo_path and mo_path.exists():
        calc.data["mad_output"] = _read_json(mo_path)

    # 3) optional: responses metadata (if present)
    rm_path = calc.artifacts.get("responses_metadata_json")
    if rm_path and rm_path.exists():
        calc.data["responses_metadata"] = _read_json(rm_path)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _populate_meta_from_calc_info(calc: Calculation, ci: dict[str, Any]) -> None:
    """
    Best-effort extraction. Don't overfit yet—just grab what’s commonly useful.
    """
    # You can adjust these keys once we look at a real calc_info schema
    # Typical candidates: molecule name, basis, method, protocol, frequencies, etc.

    # Example patterns (safe lookups):
    for key in ("molecule", "mol", "name"):
        if key in ci and "molecule" not in calc.meta:
            calc.meta["molecule"] = ci[key]

    for key in ("basis", "basis_set"):
        if key in ci and "basis" not in calc.meta:
            calc.meta["basis"] = ci[key]

    for key in ("protocol", "threshold", "dconv", "econv"):
        if key in ci and key not in calc.meta:
            calc.meta[key] = ci[key]
