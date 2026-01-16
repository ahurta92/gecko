from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import qcelemental as qcel

from .dalton_write_inputs import build_dalton_mol_string


def _write_text(path: Path, content: str, *, overwrite: bool) -> bool:
    """
    Write `content` to `path` unless it already exists and overwrite=False.

    Returns True if the file was written (either created or overwritten).
    """
    if path.exists() and not overwrite:
        return False
    path.write_text(content)
    return True


def _parse_geometry_file(path: Path) -> qcel.models.Molecule:
    """
    Parse a simple geometry block (geometry ... end) into a qcelemental Molecule.
    """
    text_lines = path.read_text().splitlines()
    in_geometry = False
    units = "angstrom"
    eprec: float | None = None
    symbols: list[str] = []
    coords: list[list[float]] = []

    for raw in text_lines:
        line = raw.strip()
        if not line:
            continue

        lower = line.lower()
        if lower == "geometry":
            in_geometry = True
            continue
        if not in_geometry:
            continue
        if lower == "end":
            break
        if lower.startswith("eprec"):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    eprec = float(parts[1])
                except ValueError:
                    pass
            continue
        if lower.startswith("units"):
            parts = line.split()
            if len(parts) >= 2:
                units = parts[1]
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
            except ValueError:
                continue
            symbols.append(parts[0])
            coords.append(xyz)

    if not symbols:
        raise ValueError(f"No atomic coordinates found in {path}")

    unit_lower = units.lower()
    if unit_lower.startswith("ang"):
        factor = 1 / qcel.constants.bohr2angstroms
    elif unit_lower in {"au", "bohr"}:
        factor = 1.0
    else:
        raise ValueError(f"Unsupported geometry units '{units}' in {path}")

    geometry = [[c * factor for c in xyz] for xyz in coords]
    extras = {"madness_parameters": {"units": units}}
    if eprec is not None:
        extras["madness_parameters"]["eprec"] = eprec

    return qcel.models.Molecule(
        symbols=symbols,
        geometry=geometry,
        extras=extras,
        fix_com=False,
        fix_orientation=False,
    )


def load_molecules_from_dir(
    molecule_dir: Path, names: Iterable[str] | None = None
) -> dict[str, qcel.models.Molecule]:
    """
    Load a mapping of molecule name -> qcelemental Molecule from a directory of .mol files.
    """
    molecule_dir = molecule_dir.expanduser().resolve()
    if not molecule_dir.is_dir():
        raise FileNotFoundError(f"Molecule directory not found: {molecule_dir}")

    selected = list(names) if names else [p.stem for p in sorted(molecule_dir.glob("*.mol"))]
    if not selected:
        raise FileNotFoundError(f"No .mol files found in {molecule_dir}")

    mol_dict: dict[str, qcel.models.Molecule] = {}
    for name in selected:
        path = molecule_dir / f"{name}.mol"
        if not path.is_file():
            raise FileNotFoundError(f"Molecule file not found: {path}")
        mol_dict[name] = _parse_geometry_file(path)
    return mol_dict


def _extract_madness_eprec(mol: qcel.models.Molecule) -> float | None:
    extras = getattr(mol, "extras", None) or {}
    for key in ("madness_parameters", "parameters"):
        params = extras.get(key, {})
        if isinstance(params, dict) and params.get("eprec") is not None:
            try:
                return float(params["eprec"])
            except (TypeError, ValueError):
                return None
    if extras.get("eprec") is not None:
        try:
            return float(extras["eprec"])
        except (TypeError, ValueError):
            return None
    return None


def build_madness_geometry_string(mol: qcel.models.Molecule, *, eprec: float | None) -> str:
    """
    Render a MADNESS molecule string and inject an `eprec` line when provided.
    """
    geom = mol.to_string("madness")
    target_eprec = _extract_madness_eprec(mol) if eprec is None else eprec
    if target_eprec is None:
        return geom

    lines = geom.splitlines()
    out: list[str] = []
    inserted = False
    replaced = False
    for line in lines:
        stripped = line.strip().lower()
        if stripped.startswith("eprec") and not replaced:
            out.append(f"eprec {target_eprec}")
            replaced = True
            continue
        out.append(line)
        if stripped == "molecule" and not inserted:
            out.append(f"eprec {target_eprec}")
            inserted = True
    if replaced:
        return "\n".join(out)
    if not inserted:
        out.insert(0, f"eprec {target_eprec}")
    return "\n".join(out)

# ----------------------------------------------------------------------
# Database setup helper
# ----------------------------------------------------------------------
def setup_ramanbench_db(
    db_root: Path | str,
    mol_dict: dict[str, qcel.models.Molecule],
    madness_raman_input: str,
    dalton_raman_input: str,
    dalton_opt_input: str,
    dalton_basis: str = "d-aug-cc-pVTZ",
    madness_eprec: float | None = None,
    overwrite: bool = False,
) -> None:
    """
    Create a simple RamanBenchDB-style layout:

    RamanBenchDB/
      <mol_name>/
        madness/
          mad.raman.in          # per-molecule MADNESS Raman input
          <mol_name>.madmol     # MADNESS geometry/molecule
        dalton/
          <mol_name>_<basis>.mol  # Dalton .mol (geometry + basis)
          optimize.dal
          raman.dal
    """
    db_root = Path(db_root).resolve()
    db_root.mkdir(parents=True, exist_ok=True)

    for mol_name, mol in mol_dict.items():
        # Top-level molecule directory
        mol_dir = db_root / mol_name
        madness_dir = mol_dir / "madness"
        dalton_dir = mol_dir / "dalton"

        madness_dir.mkdir(parents=True, exist_ok=True)
        dalton_dir.mkdir(parents=True, exist_ok=True)

        # ---------------- MADNESS setup ----------------
        # Geometry in MADNESS format
        madness_geom_str = build_madness_geometry_string(mol, eprec=madness_eprec)
        full_madness_input = f"{madness_raman_input.rstrip()}\n\n{madness_geom_str}"
        madness_input_path = madness_dir / "mad.raman.in"
        mad_input_written = _write_text(madness_input_path, full_madness_input, overwrite=overwrite)
        mad_geom_path = madness_dir / f"{mol_name}.madmol"
        mad_geom_written = _write_text(mad_geom_path, madness_geom_str, overwrite=overwrite)

        # ---------------- DALTON setup ----------------
        # Build .mol file with the chosen basis
        dalton_mol_str = build_dalton_mol_string(mol, dalton_basis)
        dalton_mol_path = dalton_dir / f"{mol_name}_{dalton_basis}.mol"
        dalton_mol_written = _write_text(dalton_mol_path, dalton_mol_str, overwrite=overwrite)

        # Copy optimization and Raman input templates
        dalton_opt_path = dalton_dir / "optimize.dal"
        dalton_raman_path = dalton_dir / "raman.dal"

        dalton_opt_written = _write_text(dalton_opt_path, dalton_opt_input, overwrite=overwrite)
        dalton_raman_written = _write_text(dalton_raman_path, dalton_raman_input, overwrite=overwrite)

        # Later:
        #  - Run Dalton geometry optimization using optimize.dal + .mol
        #  - Parse final optimized geometry and write a new
        #    <mol_name>_<basis>_opt.mol (or similar) into dalton_dir
        #  - Then use that opt geometry + raman.dal for final Raman runs

        status_bits = [
            f"mad.raman.in {'wrote' if mad_input_written else 'kept'}",
            f"{mad_geom_path.name} {'wrote' if mad_geom_written else 'kept'}",
            f"{dalton_mol_path.name} {'wrote' if dalton_mol_written else 'kept'}",
            f"{dalton_opt_path.name} {'wrote' if dalton_opt_written else 'kept'}",
            f"{dalton_raman_path.name} {'wrote' if dalton_raman_written else 'kept'}",
        ]
        print(f"Set up molecule {mol_name} in {mol_dir} ({', '.join(status_bits)})")


def _cli() -> int:
    parser = argparse.ArgumentParser(description="Initialize RamanBenchDB inputs.")
    parser.add_argument("--db-root", type=Path, default=Path("RamanBenchDB"), help="Database root directory.")
    parser.add_argument(
        "--molecule-dir",
        type=Path,
        default=Path("molecules"),
        help="Directory containing input geometries (*.mol).",
    )
    parser.add_argument(
        "--molecule",
        action="append",
        dest="molecules",
        help="Molecule name to initialize (repeatable). Defaults to all .mol in molecule-dir.",
    )
    parser.add_argument("--basis", default="d-aug-cc-pVTZ", help="Basis set label for Dalton inputs.")
    parser.add_argument(
        "--madness-template",
        type=Path,
        default=Path("templates/mad.raman.in"),
        help="MADNESS Raman input template.",
    )
    parser.add_argument(
        "--dalton-opt-template",
        type=Path,
        default=Path("templates/optimize.dal"),
        help="Dalton optimization input template.",
    )
    parser.add_argument(
        "--dalton-raman-template",
        type=Path,
        default=Path("templates/raman.dal"),
        help="Dalton Raman input template.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite any existing inputs instead of keeping them.",
    )
    parser.add_argument(
        "--madness-eprec",
        type=float,
        default=None,
        help="Override eprec (nuclear smoothing) added to MADNESS geometry blocks.",
    )
    args = parser.parse_args()

    mol_dict = load_molecules_from_dir(args.molecule_dir, names=args.molecules)
    setup_ramanbench_db(
        db_root=args.db_root,
        mol_dict=mol_dict,
        madness_raman_input=Path(args.madness_template).read_text(),
        dalton_raman_input=Path(args.dalton_raman_template).read_text(),
        dalton_opt_input=Path(args.dalton_opt_template).read_text(),
        dalton_basis=args.basis,
        madness_eprec=args.madness_eprec,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
