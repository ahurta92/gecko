from __future__ import annotations

import json
from pathlib import Path

from gecko.core.model import Calculation
from gecko.plugins.madness.detect import can_load
from gecko.plugins.madness.parse import parse_run


def load(path: Path) -> Calculation:
    root = path.expanduser().resolve()

    if not can_load(root):
        raise ValueError(f"Not a MADNESS run directory: {root}")

    artifacts = _discover_artifacts(root)
    calc = Calculation(code="madness", root=root, artifacts=artifacts, data={}, meta={})

    # Fill data/meta using parse_run
    parse_run(calc)
    return calc


def _discover_artifacts(root: Path) -> dict[str, Path]:
    artifacts: dict[str, Path] = {}

    # MADQC "marker": {stem}.calc_info.json with paired {stem}.in (required by contract).
    for p in sorted(root.glob("*.calc_info.json")):
        stem = p.name[: -len(".calc_info.json")]
        input_in = root / f"{stem}.in"
        if input_in.exists():
            artifacts["calc_info_json"] = p
            artifacts["input_in"] = input_in
            raw_out = root / f"{stem}.out"
            if raw_out.exists():
                artifacts["raw_out"] = raw_out
            break

    # If you still produce something like "n12_mad_output.json" or "*_mad_output.json"
    # keep this for compatibility.
    mad_out = next(iter(root.glob("*_mad_output.json")), None)
    if mad_out:
        artifacts["mad_output_json"] = mad_out

    # Future / common location
    resp_meta = root / "responses" / "metadata.json"
    if resp_meta.exists():
        artifacts["responses_metadata_json"] = resp_meta
    else:
        legacy_resp_meta = root / "response_metadata.json"
        if legacy_resp_meta.exists():
            artifacts["responses_metadata_json"] = legacy_resp_meta

    # Legacy molresponse markers
    # Prefer output.json, but accept outputs.json (older databases).
    legacy = (root / "output.json")
    if legacy.exists():
        artifacts["output_json"] = legacy
    else:
        legacy_plural = (root / "outputs.json")
        if legacy_plural.exists():
            artifacts["output_json"] = legacy_plural

    input_json = root / "input.json"
    if input_json.exists():
        artifacts["input_json"] = input_json

    return artifacts
