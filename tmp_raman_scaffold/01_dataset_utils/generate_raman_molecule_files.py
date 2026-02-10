#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import qcelemental as qcel

try:
    from gecko.plugins.dalton.legacy.dalton_write_inputs import to_string as dalton_mol_to_string
    from gecko.plugins.madness.legacy.madness_molecule import MADMolecule
except ImportError:
    gecko_src = Path(__file__).resolve().parents[2] / "gecko" / "src"
    if gecko_src.exists():
        sys.path.insert(0, str(gecko_src))
        from gecko.plugins.dalton.legacy.dalton_write_inputs import to_string as dalton_mol_to_string
        from gecko.plugins.madness.legacy.madness_molecule import MADMolecule
    else:
        raise

DEFAULT_BASIS_SETS = (
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "aug-cc-pVQZ",
)


def _parse_bool(value: str, default: bool = True) -> bool:
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return default


def parse_madness_mol_as_qcel(
    mol_path: Path,
    *,
    charge: int,
    multiplicity: int,
    no_orient_default: bool,
) -> tuple[qcel.models.Molecule, bool]:
    symbols: list[str] = []
    coords: list[list[float]] = []
    units = "angstrom"
    no_orient = no_orient_default

    for raw in mol_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        lower = line.lower()
        if lower in {"geometry", "molecule"}:
            continue
        if lower == "end":
            break

        parts = line.split()
        key = parts[0].lower()
        if key == "units" and len(parts) >= 2:
            units = parts[1].lower()
            continue
        if key == "no_orient" and len(parts) >= 2:
            no_orient = _parse_bool(parts[1], default=no_orient_default)
            continue
        if key in {"eprec", "field", "psp_calc", "pure_ae", "symtol", "core_type"}:
            continue

        if len(parts) != 4:
            continue

        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if not symbols:
        raise ValueError(f"No atom coordinates parsed from {mol_path}")

    if units in {"angstrom", "angs"}:
        factor = qcel.constants.conversion_factor("angstrom", "bohr")
        coords = [[float(c * factor) for c in xyz] for xyz in coords]
    elif units in {"bohr", "atomic", "au"}:
        pass
    else:
        raise ValueError(f"Unsupported units in {mol_path}: {units}")

    qmol = qcel.models.Molecule(
        symbols=symbols,
        geometry=coords,
        molecular_charge=charge,
        molecular_multiplicity=multiplicity,
        fix_orientation=no_orient,
        name=mol_path.stem,
    )
    qmol.extras["no_orient"] = bool(no_orient)
    return qmol, bool(no_orient)


def to_madness_molecule(
    qmol: qcel.models.Molecule,
    *,
    eprec: float,
    no_orient: bool,
) -> MADMolecule:
    bohr_to_ang = qcel.constants.conversion_factor("bohr", "angstrom")
    geometry = [[float(coord * bohr_to_ang) for coord in atom] for atom in qmol.geometry.tolist()]

    return MADMolecule(
        name=str(qmol.name),
        symbols=[str(symbol) for symbol in qmol.symbols],
        geometry=geometry,
        parameters={
            "eprec": float(eprec),
            "field": [0.0, 0.0, 0.0],
            "no_orient": bool(no_orient),
            "psp_calc": False,
            "pure_ae": True,
            "symtol": -1e-2,
            "core_type": "none",
            "units": "angstrom",
        },
    )


def parse_madness_frequencies(template_text: str) -> list[float]:
    m = re.search(r"dipole\.frequencies\s*\[([^\]]+)\]", template_text, flags=re.IGNORECASE)
    if m is None:
        raise ValueError("Could not find dipole.frequencies [...] block in MADNESS template")

    payload = m.group(1)
    tokens = [tok.strip() for tok in payload.replace("\n", " ").split(",")]
    freqs = [float(tok) for tok in tokens if tok]
    if not freqs:
        raise ValueError("MADNESS dipole.frequencies block is empty")
    return freqs


def parse_dalton_frequency_blocks(template_text: str) -> list[list[float]]:
    lines = template_text.splitlines()
    blocks: list[list[float]] = []

    i = 0
    while i < len(lines):
        if lines[i].strip().upper() != ".FREQUE":
            i += 1
            continue

        j = i + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines):
            raise ValueError("Malformed DALTON template: missing frequency count after .FREQUE")

        try:
            nfreq = int(lines[j].strip())
        except ValueError as exc:
            raise ValueError("Malformed DALTON template: invalid frequency count after .FREQUE") from exc

        values: list[float] = []
        k = j + 1
        while k < len(lines) and len(values) < nfreq:
            text = lines[k].strip()
            if not text:
                k += 1
                continue
            if text.startswith("*") or text.startswith("."):
                break
            values.extend(float(tok) for tok in text.split())
            k += 1

        if len(values) != nfreq:
            raise ValueError(
                f"Malformed DALTON template: .FREQUE expected {nfreq} values, parsed {len(values)}"
            )

        blocks.append(values)
        i = k

    if not blocks:
        raise ValueError("Could not find any .FREQUE blocks in DALTON template")

    return blocks


def validate_frequency_consistency(mad_freqs: list[float], dalton_blocks: list[list[float]]) -> None:
    tol = 1e-12
    for idx, block in enumerate(dalton_blocks, start=1):
        if len(block) != len(mad_freqs):
            raise ValueError(
                f"Frequency mismatch: MADNESS has {len(mad_freqs)} values but DALTON block {idx} has {len(block)}"
            )
        for a, b in zip(mad_freqs, block):
            if abs(a - b) > tol:
                raise ValueError(
                    "Frequency mismatch between MADNESS and DALTON templates: "
                    f"MADNESS={mad_freqs}, DALTON block {idx}={block}"
                )


def _render_mad_input(template_text: str, mol_text: str) -> str:
    body = template_text.rstrip()
    if re.search(r"(?mi)^\s*(geometry|molecule)\s*$", body) is not None:
        return body + "\n"
    return body + "\n\n" + mol_text.strip() + "\n"


def parse_args() -> argparse.Namespace:
    default_project_root = Path(__file__).resolve().parents[1] / "data" / "raman_paper"

    parser = argparse.ArgumentParser(
        description=(
            "Generate Raman project run directories from source molecules and templates. "
            "Also validates that MADNESS and DALTON frequency lists are identical."
        )
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=default_project_root,
        help=f"Raman project root (default: {default_project_root})",
    )
    parser.add_argument(
        "--molecules-dir",
        type=Path,
        default=None,
        help="Directory containing source MADNESS molecule files (default: <project-root>/molecules)",
    )
    parser.add_argument(
        "--inputs-dir",
        type=Path,
        default=None,
        help="Directory containing shared input templates (default: <project-root>/inputs)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Root where molecule run directories are written (default: <project-root>/data)",
    )
    parser.add_argument(
        "--mra-name",
        default="mra-p07",
        help="MADNESS run subdirectory name (default: mra-p07)",
    )
    parser.add_argument(
        "--basis-set",
        dest="basis_sets",
        nargs="+",
        default=list(DEFAULT_BASIS_SETS),
        help="DALTON basis sets used to generate runs",
    )
    parser.add_argument(
        "--mad-template",
        default="mad.raman.in",
        help="MADNESS template filename in inputs dir (default: mad.raman.in)",
    )
    parser.add_argument(
        "--raman-template",
        default="raman.dal",
        help="DALTON Raman template filename in inputs dir (default: raman.dal)",
    )
    parser.add_argument(
        "--optimize-template",
        default="optimize.dal",
        help="DALTON optimize template filename in inputs dir (default: optimize.dal)",
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Molecular charge for parsed source molecules (default: 0)",
    )
    parser.add_argument(
        "--multiplicity",
        type=int,
        default=1,
        help="Molecular multiplicity for parsed source molecules (default: 1)",
    )
    parser.add_argument(
        "--eprec",
        type=float,
        default=1e-6,
        help="MADNESS eprec written to generated molecule files (default: 1e-6)",
    )
    parser.add_argument(
        "--no-orient",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write no_orient/fix_orientation in generated molecules (default: True)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing generated files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview actions without writing files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    project_root = args.project_root.expanduser().resolve()
    molecules_dir = (args.molecules_dir or (project_root / "molecules")).expanduser().resolve()
    inputs_dir = (args.inputs_dir or (project_root / "inputs")).expanduser().resolve()
    output_root = (args.output_root or (project_root / "data")).expanduser().resolve()

    if not molecules_dir.exists():
        raise SystemExit(f"Molecules directory not found: {molecules_dir}")
    if not inputs_dir.exists():
        raise SystemExit(f"Inputs directory not found: {inputs_dir}")

    mad_template_path = inputs_dir / args.mad_template
    raman_template_path = inputs_dir / args.raman_template
    optimize_template_path = inputs_dir / args.optimize_template

    for p in (mad_template_path, raman_template_path, optimize_template_path):
        if not p.exists():
            raise SystemExit(f"Template file not found: {p}")

    mad_template_text = mad_template_path.read_text(encoding="utf-8")
    raman_template_text = raman_template_path.read_text(encoding="utf-8")
    optimize_template_text = optimize_template_path.read_text(encoding="utf-8")

    mad_freqs = parse_madness_frequencies(mad_template_text)
    dalton_blocks = parse_dalton_frequency_blocks(raman_template_text)
    validate_frequency_consistency(mad_freqs, dalton_blocks)

    print(f"[freq] MADNESS frequencies : {mad_freqs}")
    print(f"[freq] DALTON blocks found : {len(dalton_blocks)}")

    mol_files = sorted(molecules_dir.glob("*.mol"))
    if not mol_files:
        raise SystemExit(f"No .mol files found in {molecules_dir}")

    written = 0
    skipped = 0

    for mol_path in mol_files:
        molecule_name = mol_path.stem
        qmol, _source_no_orient = parse_madness_mol_as_qcel(
            mol_path,
            charge=args.charge,
            multiplicity=args.multiplicity,
            no_orient_default=args.no_orient,
        )
        no_orient = bool(args.no_orient)
        extras = dict(getattr(qmol, "extras", {}) or {})
        extras["no_orient"] = bool(no_orient)
        qmol = qmol.copy(
            update={
                "name": molecule_name,
                "fix_orientation": bool(no_orient),
                "extras": extras,
            }
        )

        mra_dir = output_root / molecule_name / args.mra_name
        mad_mol_path = mra_dir / f"{molecule_name}.mol"
        mad_in_path = mra_dir / f"{molecule_name}_raman.in"

        mad_mol = to_madness_molecule(qmol, eprec=args.eprec, no_orient=no_orient)

        if not args.dry_run:
            mra_dir.mkdir(parents=True, exist_ok=True)

        if mad_mol_path.exists() and not args.overwrite:
            skipped += 1
            print(f"[skip] {mad_mol_path}")
        else:
            if not args.dry_run:
                mad_mol.to_molfile(mad_mol_path)
            written += 1
            print(f"[write] {mad_mol_path}")

        if mad_in_path.exists() and not args.overwrite:
            skipped += 1
            print(f"[skip] {mad_in_path}")
        else:
            if not args.dry_run:
                if not mad_mol_path.exists():
                    mad_mol.to_molfile(mad_mol_path)
                mol_text = mad_mol_path.read_text(encoding="utf-8")
                mad_input_text = _render_mad_input(mad_template_text, mol_text)
                mad_in_path.write_text(mad_input_text, encoding="utf-8")
            written += 1
            print(f"[write] {mad_in_path}")

        for basis in args.basis_sets:
            basis_dir = output_root / molecule_name / basis
            mol_out = basis_dir / f"{molecule_name}_{basis}.mol"
            optimize_out = basis_dir / "optimize.dal"
            raman_out = basis_dir / "raman.dal"

            if not args.dry_run:
                basis_dir.mkdir(parents=True, exist_ok=True)

            if mol_out.exists() and not args.overwrite:
                skipped += 1
                print(f"[skip] {mol_out}")
            else:
                if not args.dry_run:
                    mol_text = dalton_mol_to_string(qmol, basis, units="angstrom")
                    mol_out.write_text(f"{mol_text}\n", encoding="utf-8")
                written += 1
                print(f"[write] {mol_out}")

            for src_text, out_path in (
                (optimize_template_text, optimize_out),
                (raman_template_text, raman_out),
            ):
                if out_path.exists() and not args.overwrite:
                    skipped += 1
                    print(f"[skip] {out_path}")
                    continue
                if not args.dry_run:
                    out_path.write_text(src_text.rstrip() + "\n", encoding="utf-8")
                written += 1
                print(f"[write] {out_path}")

    mode = "Planned" if args.dry_run else "Wrote"
    print(f"{mode} {written} files (skipped {skipped}).")


if __name__ == "__main__":
    main()
