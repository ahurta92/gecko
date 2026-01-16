import numpy as np
import qcelemental as qcel


def to_qcel(mol):
    """
    Convert a molecule to a QCEngine schema.
    """
    units = mol.get("parameters", {}).get("units", "Angstrom")
    # convert to bohr if needed
    if units.lower() == "angstrom":
        geometry = (1 / qcel.constants.bohr2angstroms) * np.array(mol["geometry"])
    else:
        geometry = np.array(mol["geometry"])

    return qcel.models.Molecule(**{"symbols": mol["symbols"], "geometry": geometry})
