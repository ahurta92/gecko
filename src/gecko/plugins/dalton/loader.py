# src/gecko/plugins/dalton/loader.py
from __future__ import annotations
from pathlib import Path
from typing import Any

from gecko.core.model import Calculation
from gecko.plugins.dalton.detect import can_load
from gecko.plugins.dalton.parse import parse_run

def load(
    path: Path,
    *,
    output_file: Path | None = None,
    run_id: str | None = None,
    meta: dict | None = None,
) -> Calculation:
    root = path.expanduser().resolve()
    if output_file is None and not can_load(root):
        raise ValueError(f"Not a DALTON run directory: {root}")

    artifacts = _discover_artifacts(root, output_file=output_file)
    calc = Calculation(code="dalton", root=root, artifacts=artifacts, data={}, meta={})
    if meta:
        calc.meta.update(meta)
    if run_id:
        calc.meta["run_id"] = run_id
        calc.meta["calc_key"] = f"{root}:{run_id}"
    parse_run(calc)
    return calc


def _discover_artifacts(root: Path, *, output_file: Path | None = None) -> dict[str, Any]:
    artifacts: dict[str, Any] = {}

    dal_files = sorted(root.glob("*.dal"))
    mol_files = sorted(root.glob("*.mol"))
    out_files = sorted(root.glob("*.out"))
    dalton_upper = root / "DALTON.OUT"
    if dalton_upper.exists() and dalton_upper not in out_files:
        out_files.insert(0, dalton_upper)

    artifacts["dalton_dal_files"] = dal_files
    artifacts["dalton_mol_files"] = mol_files
    artifacts["dalton_out_files"] = out_files

    pairs: list[dict[str, Path]] = []
    for dal in dal_files:
        for mol in mol_files:
            out_path = root / f"{dal.stem}_{mol.stem}.out"
            if out_path.exists():
                pairs.append({"dal": dal, "mol": mol, "out": out_path})
    artifacts["dalton_pairs"] = pairs

    if output_file is not None:
        out_path = Path(output_file).expanduser().resolve()
        artifacts["out"] = out_path
        artifacts["dalton_out"] = out_path
    else:
        preferred = root / "DALTON.OUT"
        if preferred.exists():
            artifacts["out"] = preferred
            artifacts["dalton_out"] = preferred
        else:
            if out_files:
                artifacts["out"] = out_files[0]
                artifacts["dalton_out"] = out_files[0]

    if "out" in artifacts:
        out_path: Path = artifacts["out"]
        quad_candidates = [
            p
            for p in root.glob("*.out")
            if p != out_path
            and any(tok in p.name.lower() for tok in ("quad", "qr", "response"))
        ]
        if quad_candidates:
            artifacts["dalton_quad_out"] = quad_candidates[0]

    return artifacts
