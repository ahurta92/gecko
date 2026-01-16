#!/usr/bin/env python3
from __future__ import annotations
import argparse
import subprocess
from pathlib import Path
from typing import List


def submit_madness_jobs(
    db_root: Path,
    molecules: List[str],
    run_script: Path,
    input_file_name: str = "mad.raman.in",
):
    """
    Submit MADNESS runs for each molecule.

    Expected directory structure:

      RamanBenchDB/<mol>/madness/
         mad.raman.in
         <mol>.madmol

    The script will submit:
        run_madness.sh mad.raman.in
    with cwd set to the madness directory.
    """
    db_root = db_root.resolve()
    run_script = run_script.resolve()

    if not run_script.is_file():
        raise FileNotFoundError(f"run_madness.sh not found at: {run_script}")

    for mol in molecules:
        madness_dir = db_root / mol / "madness"

        if not madness_dir.is_dir():
            print(f"[WARN] No madness directory for {mol}: {madness_dir}")
            continue

        input_file = madness_dir / input_file_name

        if not input_file.is_file():
            print(f"[WARN] Missing MADNESS input file for {mol}: {input_file}")
            continue

        print(f"[INFO] Submitting MADNESS run for {mol}")
        print(f"       cwd: {madness_dir}")
        print(f"       input: {input_file_name}")

        cmd = [
            str(run_script),
            input_file_name,   # pass just filename, since cwd=madness_dir
        ]

        # submit job
        subprocess.run(cmd, cwd=madness_dir, check=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit MADNESS Raman jobs for a list of molecules."
    )
    parser.add_argument(
        "--db-root",
        type=Path,
        default=Path("../RamanBenchDB"),
        help="Path to RamanBenchDB root directory",
    )
    parser.add_argument(
        "--run-script",
        type=Path,
        default=Path("./run_madness.sh"),
        help="Path to the MADNESS run script (default: ./run_madness.sh)",
    )
    parser.add_argument(
        "molecules",
        nargs="+",
        help="List of molecule names, e.g. H2 H2O BH3 C2H4",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    submit_madness_jobs(
        db_root=args.db_root,
        molecules=args.molecules,
        run_script=args.run_script,
    )

