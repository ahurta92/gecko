from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Optional, Sequence

import qcelemental as qcel

from ramanbench.db.dalton_write_inputs import build_dalton_mol_string
from ramanbench.db.layout import RamanBenchDB
from ramanbench.parsers import dalton as dalton_parser


DaltonActionType = Literal["done", "submit_raman", "submit_optimize", "blocked"]


_ATOMTYPE_CHARGE_RE = re.compile(r"^\s*Atomtype=\d+\s+\S+\s+Charge=([+-]?\d+)\b")


@dataclass(frozen=True)
class DaltonCheckpointAction:
    mol: str
    basis: str
    action: DaltonActionType
    reason: str
    optimize_out: Optional[Path] = None
    raman_out: Optional[Path] = None
    opt_mol: Optional[Path] = None
    script: Optional[Path] = None


def _read_lines(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


def can_parse_dalton_raman_output(path: Path) -> bool:
    """
    "Done" means we can parse at least one Raman table with at least one row.
    """
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        lines = _read_lines(path)
        tables, _ = dalton_parser.parse_all_raman_tables(lines)
        return any(len(rows) > 0 for rows in tables.values())
    except Exception:
        return False


def parse_final_optimized_geometry_from_opt_output(path: Path) -> qcel.models.Molecule:
    """
    Extract the final optimized geometry (Å) from a Dalton optimization output.
    """
    lines = _read_lines(path)
    return dalton_parser.parse_optimized_geometry(lines)


def _read_charge_from_dalton_mol(mol_file: Path) -> Optional[int]:
    if not mol_file.exists():
        return None
    for line in mol_file.read_text(errors="replace").splitlines():
        m = _ATOMTYPE_CHARGE_RE.match(line)
        if m:
            return int(m.group(1))
    return None


def write_dalton_opt_mol(
    db: RamanBenchDB,
    mol: str,
    basis: str,
    optimized_geometry: qcel.models.Molecule,
) -> Path:
    """
    Write `<mol>_<basis>_opt.mol` in the Dalton directory, preserving molecular charge
    from the original `<mol>_<basis>.mol` when available.
    """
    dalton_dir = db.dalton_dir(mol)
    dalton_dir.mkdir(parents=True, exist_ok=True)

    original_mol = db.dalton_mol(mol, basis)
    charge = _read_charge_from_dalton_mol(original_mol)
    mol_text = build_dalton_mol_string(optimized_geometry, basis, charge=charge)

    opt_mol = db.dalton_opt_mol(mol, basis)
    opt_mol.write_text(mol_text)
    return opt_mol


def _default_job_script_body(
    *,
    job: Literal["optimize", "raman"],
    input_dal: str,
    mol_file: str,
    base_name: str,
    out_file: str,
    err_file: str,
    nodes: int,
    ntasks_per_node: int,
    walltime: str,
    partition: str,
    dalton_gb: int,
    load_env_command: str,
) -> str:
    """
    Generates a bash wrapper that submits a Dalton job via Slurm (sbatch heredoc),
    matching the `ramanbench/cli/run_dalton.sh` style on this system.
    """
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            'cd "$(dirname "$0")"',
            f'input_file="{input_dal}"',
            f'mol_file="{mol_file}"',
            "",
            "sbatch <<EOF",
            "#!/bin/bash",
            f"#SBATCH --job-name=dalton_{base_name}",
            f"#SBATCH --nodes={nodes}",
            f"#SBATCH --ntasks-per-node={ntasks_per_node}",
            f"#SBATCH --time={walltime}",
            f"#SBATCH --output={out_file}",
            f"#SBATCH --error={err_file}",
            f"#SBATCH --partition={partition}",
            "",
            "module purge",
            load_env_command,
            "",
            'echo "Running Dalton job: ' + base_name + '"',
            'echo "Input file: ${input_file}"',
            'echo "Molecule file: ${mol_file}"',
            'echo "--------------------------------------"',
            "",
            f"dalton -N \\$SLURM_NTASKS -gb {dalton_gb} -dal $input_file -mol $mol_file",
            "",
            'echo "Job completed successfully."',
            "EOF",
            "",
        ]
    )


def ensure_job_script(
    db: RamanBenchDB,
    mol: str,
    basis: str,
    job: Literal["optimize", "raman"],
    *,
    nodes: int = 6,
    ntasks_per_node: int = 96,
    walltime: str = "48:00:00",
    partition: str = "hbm-long-96core",
    dalton_gb: int = 10,
    load_env_command: str = "source ~/load_xeonmax.sh",
    overwrite: bool = False,
) -> Path:
    dalton_dir = db.dalton_dir(mol)
    dalton_dir.mkdir(parents=True, exist_ok=True)

    input_dal = db.dalton_opt_input(mol).name if job == "optimize" else db.dalton_raman_input(mol).name
    mol_file = db.dalton_mol(mol, basis).name if job == "optimize" else db.dalton_opt_mol(mol, basis).name

    # Name outputs in a basis-aware way; layout.py will find them.
    base_name = f"dalton_{job}_{mol}_{basis}" if job == "optimize" else f"dalton_{job}_{mol}_opt_{basis}"
    out_file = f"{base_name}.out"
    err_file = f"{base_name}.err"

    script = dalton_dir / f"run_{job}_{basis}.sh"
    if script.exists() and not overwrite:
        return script

    script.write_text(
        _default_job_script_body(
            job=job,
            input_dal=input_dal,
            mol_file=mol_file,
            base_name=base_name,
            out_file=out_file,
            err_file=err_file,
            nodes=nodes,
            ntasks_per_node=ntasks_per_node,
            walltime=walltime,
            partition=partition,
            dalton_gb=dalton_gb,
            load_env_command=load_env_command,
        )
    )
    script.chmod(0o755)
    return script


def plan_dalton_checkpoints(
    db_root: Path | str,
    molecules: Iterable[str],
    basis_list: Sequence[str],
    *,
    write_opt_mol: bool = False,
) -> list[DaltonCheckpointAction]:
    db = RamanBenchDB(Path(db_root))
    actions: list[DaltonCheckpointAction] = []

    for mol in molecules:
        for basis in basis_list:
            raman_out = db.dalton_output(mol, "raman", basis)
            if raman_out and can_parse_dalton_raman_output(raman_out):
                actions.append(
                    DaltonCheckpointAction(
                        mol=mol,
                        basis=basis,
                        action="done",
                        reason=f"Parsed Raman tables from `{raman_out.name}`.",
                        raman_out=raman_out,
                    )
                )
                continue

            optimize_out = db.dalton_output(mol, "optimize", basis)
            if not optimize_out:
                actions.append(
                    DaltonCheckpointAction(
                        mol=mol,
                        basis=basis,
                        action="submit_optimize",
                        reason="No optimization output found; run optimization first.",
                    )
                )
                continue

            try:
                optimized = parse_final_optimized_geometry_from_opt_output(optimize_out)
            except Exception as e:
                actions.append(
                    DaltonCheckpointAction(
                        mol=mol,
                        basis=basis,
                        action="blocked",
                        reason=f"Optimization output exists but optimized geometry not parseable ({e}).",
                        optimize_out=optimize_out,
                    )
                )
                continue

            opt_mol = db.dalton_opt_mol(mol, basis)
            if write_opt_mol:
                opt_mol = write_dalton_opt_mol(db, mol, basis, optimized)

            actions.append(
                DaltonCheckpointAction(
                    mol=mol,
                    basis=basis,
                    action="submit_raman",
                    reason=f"Optimized geometry parsed from `{optimize_out.name}`; Raman not parsed yet.",
                    optimize_out=optimize_out,
                    raman_out=raman_out,
                    opt_mol=opt_mol,
                )
            )

    return actions


def execute_actions(
    actions: Sequence[DaltonCheckpointAction],
    db_root: Path | str,
    *,
    submit_command_template: Optional[str] = None,
    overwrite_scripts: bool = False,
    dry_run: bool = True,
    nodes: int = 6,
    ntasks_per_node: int = 96,
    walltime: str = "48:00:00",
    partition: str = "hbm-long-96core",
    dalton_gb: int = 10,
    load_env_command: str = "source ~/load_xeonmax.sh",
) -> None:
    """
    Generate job scripts and optionally submit them.

    `submit_command_template` can be e.g.:
      - "bash {script}"

    Available fields:
      {script} {mol} {basis} {action}
    """
    db = RamanBenchDB(Path(db_root))

    for a in actions:
        if a.action not in ("submit_optimize", "submit_raman"):
            continue

        job: Literal["optimize", "raman"] = "optimize" if a.action == "submit_optimize" else "raman"
        script = ensure_job_script(
            db,
            a.mol,
            a.basis,
            job,
            nodes=nodes,
            ntasks_per_node=ntasks_per_node,
            walltime=walltime,
            partition=partition,
            dalton_gb=dalton_gb,
            load_env_command=load_env_command,
            overwrite=overwrite_scripts,
        )

        if not submit_command_template:
            continue

        cmd = submit_command_template.format(
            script=str(script),
            mol=a.mol,
            basis=a.basis,
            action=a.action,
        )
        if dry_run:
            print(cmd)
        else:
            subprocess.run(cmd, shell=True, check=True)


def _cli() -> int:
    p = argparse.ArgumentParser(
        description="Checkpoint workflow for Dalton optimize->raman runs in RamanBenchDB."
    )
    p.add_argument("--db-root", type=Path, default=Path("RamanBenchDB"))
    p.add_argument("--basis", action="append", required=True, help="Repeatable; e.g. --basis d-aug-cc-pVTZ")
    p.add_argument("--molecule", action="append", help="Repeatable; defaults to all in db-root")
    p.add_argument("--write-opt-mol", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--submit", default=None, help="Submit command template, e.g. 'sbatch {script}' or 'bash {script}'")
    p.add_argument("--execute", action="store_true", help="Actually submit (otherwise dry-run prints commands).")
    p.add_argument("--overwrite-scripts", action="store_true")
    p.add_argument("--nodes", type=int, default=6)
    p.add_argument("--ntasks-per-node", type=int, default=96)
    p.add_argument("--time", dest="walltime", default="48:00:00")
    p.add_argument("--partition", default="hbm-long-96core")
    p.add_argument("--dalton-gb", type=int, default=10)
    p.add_argument(
        "--load-env",
        dest="load_env_command",
        default="source ~/load_xeonmax.sh",
        help="Command to load the cluster environment (default: 'source ~/load_xeonmax.sh')",
    )
    args = p.parse_args()

    db = RamanBenchDB(args.db_root)
    molecules = args.molecule or db.molecules()

    actions = plan_dalton_checkpoints(
        args.db_root, molecules=molecules, basis_list=args.basis, write_opt_mol=args.write_opt_mol
    )
    for a in actions:
        print(f"{a.mol:10s} {a.basis:18s} {a.action:14s} {a.reason}")

    if args.submit:
        execute_actions(
            actions,
            args.db_root,
            submit_command_template=args.submit,
            overwrite_scripts=args.overwrite_scripts,
            dry_run=not args.execute,
            nodes=args.nodes,
            ntasks_per_node=args.ntasks_per_node,
            walltime=args.walltime,
            partition=args.partition,
            dalton_gb=args.dalton_gb,
            load_env_command=args.load_env_command,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
