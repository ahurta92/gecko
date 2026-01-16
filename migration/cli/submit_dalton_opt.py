#!/usr/bin/env python3
from __future__ import annotations
import argparse
import subprocess
from pathlib import Path
from typing import List


def submit_dalton_opt_jobs(
    db_root: Path,
    basis: str,
    molecules: List[str],
    run_script: Path,
) -> None:
    """
    For each molecule, submit a Dalton optimization job using:

      RamanBenchDB/<mol>/dalton/optimize.dal
      RamanBenchDB/<mol>/dalton/<mol>_<basis>.mol

    Jobs are submitted via the provided run_script (run_dalton.sh),
    and the working directory is set to the molecule's dalton directory.
    """
    db_root = db_root.resolve()
    run_script = run_script.resolve()

    if not run_script.is_file():
        raise FileNotFoundError(f"run_dalton.sh not found at: {run_script}")

    for mol in molecules:
        mol_dir = db_root / mol
        dalton_dir = mol_dir / "dalton"

        if not dalton_dir.is_dir():
            print(f"[WARN] Dalton directory not found for molecule '{mol}': {dalton_dir}")
            continue

        input_file = dalton_dir / "optimize.dal"
        mol_file = dalton_dir / f"{mol}_{basis}.mol"

        if not input_file.is_file():
            print(f"[WARN] optimize.dal not found for '{mol}': {input_file}")
            continue

        if not mol_file.is_file():
            print(f"[WARN] .mol file not found for '{mol}' with basis '{basis}': {mol_file}")
            continue

        print(f"[INFO] Submitting optimization for {mol}")
        print(f"       cwd: {dalton_dir}")
        print(f"       dalton input: {input_file.name}")
        print(f"       mol file:     {mol_file.name}")

        # Call run_dalton.sh from within the dalton directory
        cmd = [
            str(run_script),
            input_file.name,   # note: basename, since cwd = dalton_dir
            mol_file.name,
        ]

        # This will invoke sbatch inside the dalton_dir
        subprocess.run(cmd, cwd=dalton_dir, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit Dalton optimization runs for a list of molecules."
    )
    parser.add_argument(
        "--db-root",
        type=Path,
        default=Path("../RamanBenchDB"),
        help="Path to RamanBenchDB root directory (default: ../RamanBenchDB)",
    )
    parser.add_argument(
        "--basis",
        type=str,
        default="d-aug-cc-pVTZ",
        help="Basis set tag used in the Dalton .mol file names (default: d-aug-cc-pVTZ)",
    )
    parser.add_argument(
        "--run-script",
        type=Path,
        default=Path("./run_dalton.sh"),
        help="Path to the run_dalton.sh wrapper script (default: ./run_dalton.sh)",
    )
    parser.add_argument(
        "molecules",
        nargs="+",
        help="List of molecule names, e.g. H2 H2O BH3 C2H4",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    submit_dalton_opt_jobs(
        db_root=args.db_root,
        basis=args.basis,
        molecules=args.molecules,
        run_script=args.run_script,
    )

