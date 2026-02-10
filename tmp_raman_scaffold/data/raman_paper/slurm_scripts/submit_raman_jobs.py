#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

import qcelemental as qcel


@dataclass(frozen=True)
class Job:
    code: str
    stage: str
    molecule: str
    basis: str
    run_dir: Path
    cmd: list[str]
    reason: str


def _candidate_gecko_src_paths(script_path: Path) -> list[Path]:
    out: list[Path] = []

    env_src = os.environ.get("GECKO_SRC")
    if env_src:
        out.append(Path(env_src).expanduser())

    for parent in script_path.resolve().parents:
        out.append(parent / "gecko" / "src")

    out.append(Path.cwd() / "src")
    out.append(Path.cwd() / "gecko" / "src")

    unique: list[Path] = []
    seen: set[Path] = set()
    for p in out:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        unique.append(rp)
    return unique


def _bootstrap_gecko(gecko_src: Path | None, script_path: Path) -> Callable:
    try:
        from gecko.core.load import load_calc  # type: ignore

        return load_calc
    except Exception:
        pass

    candidates: list[Path] = []
    if gecko_src is not None:
        candidates.append(gecko_src.expanduser().resolve())
    candidates.extend(_candidate_gecko_src_paths(script_path))

    for candidate in candidates:
        if not candidate.exists():
            continue
        if str(candidate) not in sys.path:
            sys.path.insert(0, str(candidate))
        try:
            from gecko.core.load import load_calc  # type: ignore

            return load_calc
        except Exception:
            continue

    raise SystemExit(
        "Could not import gecko. Install gecko in this environment or pass --gecko-src /path/to/gecko/src"
    )


def _sorted_molecule_dirs(data_root: Path, selected: set[str] | None, mra_dir: str) -> list[Path]:
    def _is_molecule_dir(path: Path) -> bool:
        if not path.is_dir():
            return False
        if (path / mra_dir).is_dir():
            return True

        for child in path.iterdir():
            if not child.is_dir() or child.name == mra_dir:
                continue
            if any(child.glob("*.mol")) and any(child.glob("*.dal")):
                return True
        return False

    def _key(name: str):
        m = re.fullmatch(r"n(\d+)", name)
        if m:
            return (0, int(m.group(1)))
        return (1, name)

    dirs = [p for p in data_root.iterdir() if _is_molecule_dir(p)]
    if selected is not None:
        dirs = [p for p in dirs if p.name in selected]
    return sorted(dirs, key=lambda p: _key(p.name))


def _find_madness_input(mra_dir: Path, molecule: str, explicit: str | None) -> Path | None:
    if not mra_dir.exists():
        return None

    if explicit:
        p = mra_dir / explicit
        if p.exists():
            return p
        return None

    preferred = mra_dir / f"{molecule}_raman.in"
    if preferred.exists():
        return preferred

    ins = sorted(mra_dir.glob("*.in"))
    if len(ins) == 1:
        return ins[0]
    if len(ins) > 1:
        for p in ins:
            if "raman" in p.name.lower():
                return p
        return ins[0]

    return None


def _select_base_mol_file(molecule: str, basis_dir: Path) -> Path | None:
    expected = basis_dir / f"{molecule}_{basis_dir.name}.mol"
    if expected.exists():
        return expected

    mols = sorted([p for p in basis_dir.glob("*.mol") if not p.name.startswith("opt_")])
    if len(mols) == 1:
        return mols[0]
    if len(mols) > 1:
        for p in mols:
            if p.stem.startswith(f"{molecule}_"):
                return p
        return mols[0]

    return None


def _basis_from_mol_file(mol_path: Path, fallback: str) -> str:
    try:
        lines = mol_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return fallback

    for i, line in enumerate(lines):
        if line.strip().upper() != "BASIS":
            continue
        for j in range(i + 1, len(lines)):
            text = lines[j].strip()
            if text:
                return text
    return fallback


def _has_nosymmetry(mol_path: Path) -> bool:
    try:
        text = mol_path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return True
    return "nosymmetry" in text.lower()


def _optimized_mol_path(basis_dir: Path, molecule: str, template: str) -> Path:
    return basis_dir / template.format(mol=molecule, basis=basis_dir.name)


def _optimization_out_candidates(basis_dir: Path) -> list[Path]:
    cands = sorted(basis_dir.glob("optimize*.out"), key=lambda p: p.stat().st_mtime, reverse=True)
    if cands:
        return cands

    generic = [
        p
        for p in basis_dir.glob("*.out")
        if "raman" not in p.name.lower() and "quad" not in p.name.lower()
    ]
    generic = sorted(generic, key=lambda p: p.stat().st_mtime, reverse=True)
    return generic


def _parse_optimized_molecule(lines: Sequence[str]):
    from gecko.plugins.dalton.legacy.dalton import parse_optimized_geometry
    from gecko.plugins.dalton.parse import parse_last_molecular_geometry

    errors: list[str] = []
    for parser in (parse_optimized_geometry, parse_last_molecular_geometry):
        try:
            return parser(lines), None
        except Exception as exc:
            errors.append(f"{parser.__name__}: {type(exc).__name__}: {exc}")
    return None, "; ".join(errors)


def _write_optimized_mol(
    out_path: Path,
    angstrom_mol,
    *,
    basis: str,
    no_orient: bool,
    dry_run: bool,
) -> None:
    from gecko.plugins.dalton.legacy.dalton_write_inputs import to_string as dalton_mol_to_string

    ang_to_bohr = qcel.constants.conversion_factor("angstrom", "bohr")
    bohr_geom = [
        [float(coord * ang_to_bohr) for coord in atom]
        for atom in angstrom_mol.geometry.tolist()
    ]

    qmol = qcel.models.Molecule(
        symbols=[str(s) for s in angstrom_mol.symbols],
        geometry=bohr_geom,
        molecular_charge=float(getattr(angstrom_mol, "molecular_charge", 0.0)),
        molecular_multiplicity=int(round(float(getattr(angstrom_mol, "molecular_multiplicity", 1)))),
        fix_orientation=bool(no_orient),
    )
    qmol.extras["no_orient"] = bool(no_orient)

    text = dalton_mol_to_string(qmol, basis, units="angstrom")
    if dry_run:
        return

    out_path.write_text(text.rstrip() + "\n", encoding="utf-8")


def _ensure_optimized_mol(
    *,
    basis_dir: Path,
    molecule: str,
    base_mol: Path,
    opt_mol: Path,
    dry_run: bool,
) -> tuple[bool, str]:
    if opt_mol.exists():
        return True, f"optimized molecule present: {opt_mol.name}"

    cands = _optimization_out_candidates(basis_dir)
    if not cands:
        return False, "optimized molecule missing and no optimize output found"

    basis = _basis_from_mol_file(base_mol, fallback=basis_dir.name)
    no_orient = _has_nosymmetry(base_mol)

    parse_errors: list[str] = []
    for out_file in cands:
        try:
            lines = out_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        except Exception as exc:
            parse_errors.append(f"{out_file.name}: read failed ({type(exc).__name__}: {exc})")
            continue

        mol, err = _parse_optimized_molecule(lines)
        if mol is None:
            if err:
                parse_errors.append(f"{out_file.name}: {err}")
            continue

        _write_optimized_mol(
            opt_mol,
            mol,
            basis=basis,
            no_orient=no_orient,
            dry_run=dry_run,
        )
        action = "would write" if dry_run else "wrote"
        return True, f"{action} {opt_mol.name} from {out_file.name}"

    preview = "; ".join(parse_errors[:2])
    if len(parse_errors) > 2:
        preview += "; ..."
    return False, f"optimize output found but final geometry parse failed ({preview})"


def _has_raman(calc: object) -> bool:
    data = getattr(calc, "data", None)
    if not isinstance(data, dict) or "raman" not in data:
        return False

    raman = data.get("raman")
    if raman is None:
        return False

    if isinstance(raman, dict):
        raman_by_freq = raman.get("raman_by_freq")
        if isinstance(raman_by_freq, dict) and raman_by_freq:
            return True

        pol = raman.get("polarization_frequencies")
        if pol is not None:
            try:
                return len(pol) > 0
            except Exception:
                return True

        vib = raman.get("vibrational_frequencies")
        if vib is not None:
            try:
                return len(vib) > 0
            except Exception:
                return True

        return bool(raman)

    try:
        return len(raman) > 0
    except Exception:
        return True


def _run_complete_raman(run_dir: Path, load_calc: Callable) -> tuple[bool, str]:
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            calc = load_calc(run_dir)
    except Exception as exc:
        return False, f"gecko.load_calc failed: {type(exc).__name__}: {exc}"

    if _has_raman(calc):
        return True, "raman present"

    return False, "calculation loaded but raman is missing"


def _running_workdirs() -> set[Path]:
    user = os.environ.get("USER")
    if not user:
        return set()

    try:
        proc = subprocess.run(
            ["squeue", "-h", "-u", user, "-o", "%Z"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return set()

    if proc.returncode != 0:
        return set()

    out: set[Path] = set()
    for line in proc.stdout.splitlines():
        wd = line.strip()
        if not wd:
            continue
        try:
            out.add(Path(wd).expanduser().resolve())
        except Exception:
            continue
    return out


def _submit(job: Job, *, dry_run: bool) -> tuple[bool, str]:
    cmd_str = " ".join(shlex.quote(x) for x in job.cmd)
    if dry_run:
        return True, f"DRY-RUN: (cd {job.run_dir} && {cmd_str})"

    proc = subprocess.run(
        job.cmd,
        cwd=job.run_dir,
        check=False,
        capture_output=True,
        text=True,
    )
    out = (proc.stdout or "").strip()
    err = (proc.stderr or "").strip()

    if proc.returncode != 0:
        msg = err or out or f"command failed with exit code {proc.returncode}"
        return False, msg

    return True, out or "submitted"


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    raman_root = script_dir.parent

    parser = argparse.ArgumentParser(
        description=(
            "Submit missing Raman MADNESS/DALTON jobs. "
            "DALTON runs are two-stage: optimize -> generate opt_<mol>_<basis>.mol -> raman."
        )
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=raman_root / "data",
        help="Root containing molecule folders (default: <raman_root>/data)",
    )
    parser.add_argument(
        "--mra-dir",
        default="mra-p07",
        help="MADNESS subdirectory name per molecule (default: mra-p07)",
    )
    parser.add_argument(
        "--madness-input",
        default=None,
        help=(
            "MADNESS input filename in each MRA dir. "
            "If omitted, script auto-detects <mol>_raman.in or *.in"
        ),
    )
    parser.add_argument(
        "--optimize-input",
        default="optimize.dal",
        help="DALTON optimize input filename in each basis dir (default: optimize.dal)",
    )
    parser.add_argument(
        "--raman-input",
        default="raman.dal",
        help="DALTON Raman input filename in each basis dir (default: raman.dal)",
    )
    parser.add_argument(
        "--optimized-mol-template",
        default="opt_{mol}_{basis}.mol",
        help=(
            "Filename template for optimized DALTON molecule file. "
            "Available fields: {mol}, {basis}"
        ),
    )
    parser.add_argument(
        "--run-madness-script",
        type=Path,
        default=script_dir / "run_madness.sh",
        help="Submission helper script for MADNESS jobs",
    )
    parser.add_argument(
        "--run-dalton-script",
        type=Path,
        default=script_dir / "run_dalton.sh",
        help="Submission helper script for DALTON jobs",
    )
    parser.add_argument(
        "--molecules",
        nargs="*",
        default=None,
        help="Optional subset of molecules (e.g. H2O CH4 NH3)",
    )
    parser.add_argument(
        "--bases",
        nargs="*",
        default=None,
        help="Optional subset of DALTON basis dirs (e.g. aug-cc-pVDZ)",
    )
    parser.add_argument(
        "--skip-madness",
        action="store_true",
        help="Do not check/submit MADNESS runs",
    )
    parser.add_argument(
        "--skip-dalton",
        action="store_true",
        help="Do not check/submit DALTON runs",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Submit response jobs even when gecko says they are complete",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Actually submit with sbatch. Default mode is dry-run.",
    )
    parser.add_argument(
        "--max-submit",
        type=int,
        default=None,
        help="Maximum number of jobs to submit (after filtering)",
    )
    parser.add_argument(
        "--gecko-src",
        type=Path,
        default=None,
        help="Optional path to gecko/src if gecko is not importable",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    raman_root = script_dir.parent

    data_root = args.data_root.expanduser().resolve()
    if not data_root.exists() and data_root == (raman_root / "data").resolve():
        # Backward-compatible fallback for legacy layout: <raman_root>/<mol>/...
        legacy = raman_root.resolve()
        if legacy.exists():
            data_root = legacy

    run_mad_script = args.run_madness_script.expanduser().resolve()
    run_dal_script = args.run_dalton_script.expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"Data root not found: {data_root}")
    if not args.skip_madness and not run_mad_script.exists():
        raise SystemExit(f"MADNESS submission script not found: {run_mad_script}")
    if not args.skip_dalton and not run_dal_script.exists():
        raise SystemExit(f"DALTON submission script not found: {run_dal_script}")

    load_calc = _bootstrap_gecko(args.gecko_src, Path(__file__))

    selected = set(args.molecules) if args.molecules else None
    selected_bases = set(args.bases) if args.bases else None

    molecule_dirs = _sorted_molecule_dirs(data_root, selected, args.mra_dir)
    if not molecule_dirs:
        raise SystemExit(f"No molecule directories found in {data_root}")

    running_dirs = _running_workdirs()

    jobs: list[Job] = []
    complete = 0
    skipped_running = 0
    skipped_missing_inputs = 0

    dry_run = not args.submit

    for mol_dir in molecule_dirs:
        mol = mol_dir.name

        if not args.skip_madness:
            mra_dir = mol_dir / args.mra_dir
            mad_input = _find_madness_input(mra_dir, mol, args.madness_input)

            if not mra_dir.exists() or mad_input is None:
                print(f"[skip][madness][{mol}] missing run dir or input")
                skipped_missing_inputs += 1
            else:
                is_done, reason = _run_complete_raman(mra_dir, load_calc)
                if args.force:
                    is_done = False
                    reason = "forced resubmission"

                if is_done:
                    print(f"[ok][madness][{mol}] {reason}")
                    complete += 1
                elif mra_dir.resolve() in running_dirs:
                    print(f"[skip-running][madness][{mol}] {reason}")
                    skipped_running += 1
                else:
                    jobs.append(
                        Job(
                            code="madness",
                            stage="raman",
                            molecule=mol,
                            basis=args.mra_dir,
                            run_dir=mra_dir,
                            cmd=["bash", str(run_mad_script), mad_input.name],
                            reason=reason,
                        )
                    )

        if args.skip_dalton:
            continue

        basis_dirs = [p for p in mol_dir.iterdir() if p.is_dir() and p.name != args.mra_dir]
        if selected_bases is not None:
            basis_dirs = [p for p in basis_dirs if p.name in selected_bases]

        for basis_dir in sorted(basis_dirs, key=lambda p: p.name):
            basis = basis_dir.name
            optimize_input = basis_dir / args.optimize_input
            raman_input = basis_dir / args.raman_input
            base_mol = _select_base_mol_file(mol, basis_dir)
            opt_mol = _optimized_mol_path(basis_dir, mol, args.optimized_mol_template)

            # If Raman is already complete, both DALTON stages are effectively complete.
            if base_mol is not None and raman_input.exists() and not args.force:
                raman_done, raman_reason = _run_complete_raman(basis_dir, load_calc)
                if raman_done:
                    print(f"[ok][dalton-optimize][{mol}/{basis}] implied by raman completion")
                    print(f"[ok][dalton-raman][{mol}/{basis}] {raman_reason}")
                    complete += 2
                    continue

            # Stage 1: optimize
            if base_mol is None or not optimize_input.exists():
                print(
                    f"[skip][dalton-optimize][{mol}/{basis}] "
                    f"missing input or base mol (input={optimize_input.exists()}, mol={base_mol is not None})"
                )
                skipped_missing_inputs += 1
                optimize_ready = False
                optimize_reason = "missing optimize input or base molecule"
            else:
                optimize_ready, optimize_reason = _ensure_optimized_mol(
                    basis_dir=basis_dir,
                    molecule=mol,
                    base_mol=base_mol,
                    opt_mol=opt_mol,
                    dry_run=dry_run,
                )

                if optimize_ready:
                    print(f"[ok][dalton-optimize][{mol}/{basis}] {optimize_reason}")
                    complete += 1
                elif basis_dir.resolve() in running_dirs:
                    print(f"[skip-running][dalton-optimize][{mol}/{basis}] {optimize_reason}")
                    skipped_running += 1
                else:
                    jobs.append(
                        Job(
                            code="dalton",
                            stage="optimize",
                            molecule=mol,
                            basis=basis,
                            run_dir=basis_dir,
                            cmd=["bash", str(run_dal_script), args.optimize_input, base_mol.name],
                            reason=optimize_reason,
                        )
                    )

            # Stage 2: raman
            if not raman_input.exists():
                print(f"[skip][dalton-raman][{mol}/{basis}] missing input: {raman_input.name}")
                skipped_missing_inputs += 1
                continue

            opt_exists_or_planned = opt_mol.exists() or (optimize_ready and dry_run)
            if not opt_exists_or_planned:
                print(f"[wait][dalton-raman][{mol}/{basis}] waiting on optimized molecule")
                continue

            raman_done, raman_reason = _run_complete_raman(basis_dir, load_calc)
            if args.force:
                raman_done = False
                raman_reason = "forced resubmission"

            if raman_done:
                print(f"[ok][dalton-raman][{mol}/{basis}] {raman_reason}")
                complete += 1
            elif basis_dir.resolve() in running_dirs:
                print(f"[skip-running][dalton-raman][{mol}/{basis}] {raman_reason}")
                skipped_running += 1
            else:
                jobs.append(
                    Job(
                        code="dalton",
                        stage="raman",
                        molecule=mol,
                        basis=basis,
                        run_dir=basis_dir,
                        cmd=["bash", str(run_dal_script), args.raman_input, opt_mol.name],
                        reason=raman_reason,
                    )
                )

    print("\n=== Submission Plan ===")
    print(f"molecule dirs scanned : {len(molecule_dirs)}")
    print(f"already complete      : {complete}")
    print(f"missing inputs        : {skipped_missing_inputs}")
    print(f"skipped (running)     : {skipped_running}")
    print(f"to submit             : {len(jobs)}")

    if args.max_submit is not None:
        jobs = jobs[: args.max_submit]
        print(f"capped submissions    : {len(jobs)} (via --max-submit)")

    success = 0
    failed = 0

    for job in jobs:
        ok, msg = _submit(job, dry_run=dry_run)
        tag = "submit" if ok else "error"
        print(f"[{tag}][{job.code}:{job.stage}][{job.molecule}/{job.basis}] {msg}")
        if ok:
            success += 1
        else:
            failed += 1

    print("\n=== Result ===")
    print(f"mode                 : {'dry-run' if dry_run else 'submit'}")
    print(f"jobs attempted       : {len(jobs)}")
    print(f"jobs succeeded       : {success}")
    print(f"jobs failed          : {failed}")

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
