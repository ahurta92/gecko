from __future__ import annotations

from hashlib import sha1
from typing import Optional

import qcelemental as qcel

from gecko.core.model import Calculation


def geom_id_from_molecule(mol: qcel.models.Molecule | None) -> Optional[str]:
    if mol is None:
        return None
    return str(mol.name)


def geom_id(molecule: qcel.models.Molecule) -> str:
    return str(molecule.name)


def calc_id(calc: Calculation) -> str:
    root = calc.root.resolve()
    base = f"{calc.code}|{root}"
    digest = sha1(base.encode("utf-8")).hexdigest()[:16]
    return f"calc_{digest}"


def mol_id_from_molecule(mol: qcel.models.Molecule | None) -> Optional[str]:
    if mol is None:
        return None
    formula = str(mol.get_molecular_formula())
    if formula is not None:
        return formula 
    return None


def mol_id(calc: Calculation) -> Optional[str]:
    return mol_id_from_molecule(calc.molecule)
