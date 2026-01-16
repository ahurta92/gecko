from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Optional, Sequence

from ramanbench.db.layout import RamanBenchDB


MadnessActionType = Literal["done", "submit", "blocked"]


@dataclass(frozen=True)
class MadnessCheckpointAction:
    mol: str
    action: MadnessActionType
    reason: str
    calc_info: Optional[Path] = None
    stdout: Optional[Path] = None
    script: Optional[Path] = None


def _default_job_script_body(
    *,
    input_file: str,
    base_name: str,
    out_file: str,
    err_file: str,
    nodes: int,
    ntasks_per_node: int,
    walltime: str,
    partition: str,
    load_env_command: str,
    mad_num_threads: int,
) -> str:
    """
    Generates a bash wrapper that submits a MADNESS job via Slurm (sbatch heredoc).
    """
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            'cd "$(dirname "$0")"',
            f'input_file="{input_file}"',
            "",
            "sbatch <<EOF",
            "#!/bin/bash",
            f"#SBATCH --job-name={base_name}",
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
            f"export MAD_NUM_THREADS={mad_num_threads}",
            "",
            'echo "Running MADNESS job: ${input_file}"',
            'echo "--------------------------------------"',
            "",
            "mpirun --map-by numa madqc --wf=response --input=${input_file}",
            "",
            'echo "Job completed successfully."',
            "EOF",
            "",
        ]
    )


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def can_parse_madness_calc_info(path: Path) -> tuple[bool, str]:
    """
    Check that calc_info exists, is readable, and contains a Raman table we can parse.
    """
    if not path.exists():
        return False, f"{path.name} not found"
    if path.stat().st_size == 0:
        return False, f"{path.name} is empty"

    try:
        data = _read_json(path)
    except Exception as exc:
        return False, f"calc_info unreadable ({exc.__class__.__name__})"

    if not isinstance(data, dict):
        return False, "calc_info JSON root is not an object"
    tasks = data.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        return False, "calc_info missing tasks array"

    try:
        # Reuse the production parser to ensure Raman data are actually readable.
        from ramanbench.parsers.madness import madqc_parser

        parsed = madqc_parser(path)
        raman_table = getattr(parsed, "raman_pivot", None)
        if raman_table is not None:
            try:
                row_count = len(raman_table)
            except Exception:
                row_count = None
            if not getattr(raman_table, "empty", False) and (row_count is None or row_count > 0):
                detail = f"{row_count} rows" if row_count is not None else "non-empty"
                return True, f"Raman table parsed ({detail})"
        return False, "calc_info missing Raman table (parse yielded empty result)"
    except Exception as exc:
        return False, f"Raman table parse failed ({exc.__class__.__name__})"


def ensure_madness_job_script(
    db: RamanBenchDB,
    mol: str,
    *,
    input_file_name: str = "mad.raman.in",
    nodes: int = 16,
    ntasks_per_node: int = 8,
    walltime: str = "08:00:00",
    partition: str = "hbm-large-96core",
    load_env_command: str = "source ~/load_xeonmax.sh",
    mad_num_threads: int = 20,
    overwrite: bool = False,
) -> Path:
    madness_dir = db.madness_dir(mol)
    madness_dir.mkdir(parents=True, exist_ok=True)

    base_name = Path(input_file_name).stem
    out_file = f"{base_name}.out"
    err_file = f"{base_name}.err"
    # if out_file and err_file exist, mv them to next available .old.N
    for f in (out_file, err_file):
        f_path = madness_dir / f
        if f_path.exists() and not overwrite:
            n = 1
            while True:
                old_path = madness_dir / f"{f}.old.{n}"
                if not old_path.exists():
                    f_path.rename(old_path)
                    break
                n += 1

    script = madness_dir / "run_madness.sh"
    if script.exists() and not overwrite:
        return script

    script.write_text(
        _default_job_script_body(
            input_file=input_file_name,
            base_name=base_name,
            out_file=out_file,
            err_file=err_file,
            nodes=nodes,
            ntasks_per_node=ntasks_per_node,
            walltime=walltime,
            partition=partition,
            load_env_command=load_env_command,
            mad_num_threads=mad_num_threads,
        )
    )
    script.chmod(0o755)
    return script


def plan_madness_checkpoints(
    db_root: Path | str,
    molecules: Iterable[str],
    *,
    input_file_name: str = "mad.raman.in",
) -> list[MadnessCheckpointAction]:
    db = RamanBenchDB(Path(db_root))
    actions: list[MadnessCheckpointAction] = []

    for mol in molecules:
        calc_info = db.madness_calc_info_json(mol)
        print(calc_info)
        stdout = db.madness_stdout(mol)

        parse_ok = False
        parse_reason = "calc_info JSON not found."
        if calc_info:
            parse_ok, parse_reason = can_parse_madness_calc_info(calc_info)

        if parse_ok:
            actions.append(
                MadnessCheckpointAction(
                    mol=mol,
                    action="done",
                    reason=f"calc_info OK ({parse_reason})",
                    calc_info=calc_info,
                    stdout=stdout,
                )
            )
            continue

        madness_dir = db.madness_dir(mol)
        input_path = madness_dir / input_file_name
        geom_path = madness_dir / f"{mol}.madmol"
        if not madness_dir.exists():
            actions.append(
                MadnessCheckpointAction(
                    mol=mol,
                    action="blocked",
                    reason=f"No MADNESS directory found: {madness_dir}",
                )
            )
            continue
        if not input_path.exists():
            actions.append(
                MadnessCheckpointAction(
                    mol=mol,
                    action="blocked",
                    reason=f"Missing input file: {input_path.name}",
                )
            )
            continue
        if not geom_path.exists():
            actions.append(
                MadnessCheckpointAction(
                    mol=mol,
                    action="blocked",
                    reason=f"Missing geometry file: {geom_path.name}",
                )
            )
            continue

        actions.append(
            MadnessCheckpointAction(
                mol=mol,
                action="submit",
                reason=f"{parse_reason}; resubmit MADNESS run.",
                calc_info=calc_info,
                stdout=stdout,
            )
        )

    return actions


def execute_actions(
    actions: Sequence[MadnessCheckpointAction],
    db_root: Path | str,
    *,
    submit_command_template: Optional[str] = None,
    dry_run: bool = True,
    input_file_name: str = "mad.raman.in",
    nodes: int = 16,
    ntasks_per_node: int = 8,
    walltime: str = "08:00:00",
    partition: str = "hbm-large-96core",
    load_env_command: str = "source ~/load_xeonmax.sh",
    mad_num_threads: int = 20,
    overwrite_scripts: bool = False,
) -> None:
    db = RamanBenchDB(Path(db_root))

    for a in actions:
        if a.action != "submit":
            continue

        script = ensure_madness_job_script(
            db,
            a.mol,
            input_file_name=input_file_name,
            nodes=nodes,
            ntasks_per_node=ntasks_per_node,
            walltime=walltime,
            partition=partition,
            load_env_command=load_env_command,
            mad_num_threads=mad_num_threads,
            overwrite=overwrite_scripts,
        )

        if not submit_command_template:
            print(f"[DRY-RUN] Would submit {script}")
            continue

        cmd = submit_command_template.format(
            script=str(script),
            mol=a.mol,
            action=a.action,
        )
        if dry_run:
            print(cmd)
        else:
            subprocess.run(cmd, shell=True, check=True)


def _cli() -> int:
    p = argparse.ArgumentParser(
        description="Checkpoint workflow for MADNESS Raman runs in RamanBenchDB."
    )
    p.add_argument("--db-root", type=Path, default=Path("RamanBenchDB"))
    p.add_argument("--molecule", action="append", help="Repeatable; defaults to all in db-root")
    p.add_argument(
        "--input-file",
        default="mad.raman.in",
        help="MADNESS input file name (default: mad.raman.in)",
    )
    p.add_argument(
        "--submit",
        default=None,
        help="Submit command template, e.g. 'bash {script}' or 'sbatch {script}'",
    )
    p.add_argument(
        "--execute",
        action="store_true",
        help="Actually submit (otherwise dry-run prints commands).",
    )
    p.add_argument("--overwrite-scripts", action="store_true")
    p.add_argument("--nodes", type=int, default=16)
    p.add_argument("--ntasks-per-node", type=int, default=8)
    p.add_argument("--time", dest="walltime", default="08:00:00")
    p.add_argument("--partition", default="hbm-large-96core")
    p.add_argument(
        "--load-env",
        dest="load_env_command",
        default="source ~/load_xeonmax.sh",
        help="Command to load the cluster environment (default: 'source ~/load_xeonmax.sh')",
    )
    p.add_argument(
        "--threads", type=int, default=20, help="MAD_NUM_THREADS to export inside the job."
    )
    args = p.parse_args()

    db = RamanBenchDB(args.db_root)
    molecules = args.molecule or db.molecules()

    actions = plan_madness_checkpoints(
        args.db_root, molecules=molecules, input_file_name=args.input_file
    )
    for a in actions:
        print(f"{a.mol:10s} {a.action:10s} {a.reason}")

    if args.submit:
        execute_actions(
            actions,
            args.db_root,
            submit_command_template=args.submit,
            dry_run=not args.execute,
            input_file_name=args.input_file,
            nodes=args.nodes,
            ntasks_per_node=args.ntasks_per_node,
            walltime=args.walltime,
            partition=args.partition,
            load_env_command=args.load_env_command,
            mad_num_threads=args.threads,
            overwrite_scripts=args.overwrite_scripts,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
