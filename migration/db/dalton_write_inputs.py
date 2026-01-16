from __future__ import annotations

from pathlib import Path
from typing import Optional

import qcelemental as qcel


def _extract_units(extras: Optional[dict]) -> str:
    if not isinstance(extras, dict):
        return "Angstrom"

    for key in ("units", "madness_units"):
        if extras.get(key):
            return str(extras[key]).capitalize()

    for key in ("parameters", "madness_parameters"):
        params = extras.get(key)
        if isinstance(params, dict) and params.get("units"):
            return str(params["units"]).capitalize()

    return "Angstrom"


def _has_no_orient(extras: Optional[dict], default: bool) -> bool:
    if not isinstance(extras, dict):
        return default

    if extras.get("no_orient") is not None:
        value = extras["no_orient"]
        return str(value).lower() == "true" if isinstance(value, str) else bool(value)

    for key in ("parameters", "madness_parameters"):
        params = extras.get(key)
        if isinstance(params, dict) and params.get("no_orient") is not None:
            value = params["no_orient"]
            return str(value).lower() == "true" if isinstance(value, str) else bool(value)

    return default


def build_dalton_mol_string(
    molecule: qcel.models.Molecule, basis: str, *, charge: Optional[int] = None
) -> str:
    """
    Return the Dalton .mol string for a qcelemental Molecule.
    """
    if not isinstance(molecule, qcel.models.Molecule):
        raise TypeError("molecule must be an instance of qcelemental.models.Molecule")

    dalton_lines = ["BASIS", str(basis), "blah", "blah"]

    extras = getattr(molecule, "extras", None) or {}
    units = _extract_units(extras)

    charge_value = molecule.molecular_charge if charge is None else charge
    charge_value = int(round(charge_value))

    atom_coords: dict[str, list[str]] = {}
    coordinates = molecule.geometry.reshape((-1, 3))
    for symbol, coord in zip(molecule.symbols, coordinates):
        formatted = " ".join(f"{float(c):.12f}" for c in coord)
        atom_coords.setdefault(symbol, []).append(formatted)

    general_line = f"Atomtype={len(atom_coords)} {units} Charge={charge_value}"
    #if _has_no_orient(extras, bool(getattr(molecule, "fix_orientation", True))):
    general_line += " Nosymmetry"
    dalton_lines.append(general_line)

    for atom, coords in atom_coords.items():
        atomic_number = qcel.periodictable.to_atomic_number(atom)
        dalton_lines.append(f"Charge={atomic_number} Atoms={len(coords)}")
        suffix = ord("a")
        for coord in coords:
            dalton_lines.append(f"{atom}_{chr(suffix)} {coord}")
            suffix += 1

    return "\n".join(dalton_lines)


def write_dalton_mol_file(
    molecule: qcel.models.Molecule,
    basis: str,
    destination: str | Path,
    *,
    charge: Optional[int] = None,
) -> str:
    """
    Write the Dalton .mol representation to `destination`.

    Returns the string that was written for convenience.
    """
    mol_string = build_dalton_mol_string(molecule, basis, charge=charge)
    Path(destination).write_text(mol_string)
    return mol_string
