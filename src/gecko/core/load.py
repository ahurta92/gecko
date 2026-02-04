from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from gecko.core.model import Calculation
from gecko.plugins.madness.detect import can_load as madness_can_load
from gecko.plugins.dalton.detect import can_load as dalton_can_load

from gecko.plugins.madness.loader import load as load_madness
from gecko.plugins.dalton.loader import load as load_dalton


def _maybe_attach_calc_info(calc: Calculation) -> None:
    """Attach optional calc_info.json payload without affecting parsing."""
    info_path = calc.root / "calc_info.json"
    if not info_path.exists():
        return
    try:
        import json
        calc.data.setdefault("calc_info", json.loads(info_path.read_text(encoding="utf-8")))
    except Exception as exc:
        calc.meta.setdefault("warnings", []).append(
            f"Failed to read calc_info.json: {info_path} ({type(exc).__name__}: {exc})"
        )


def _maybe_load_molecule_from_calc_dir(calc: Calculation) -> Optional[Any]:
    if calc.molecule is not None:
        return calc.molecule

    root = calc.root
    preferred = root / "molecule.mol"
    candidates: list[Path] = []
    if preferred.exists():
        candidates.append(preferred)
    candidates.extend(sorted(p for p in root.glob("*.mol") if p != preferred))

    if not candidates:
        return None

    from gecko.mol.io import read_mol

    for path in candidates:
        try:
            mol = read_mol(path)
            calc.meta.setdefault("molecule_source", "calc_dir")
            calc.meta.setdefault("molecule_path", str(path))
            return mol
        except Exception as exc:
            calc.meta.setdefault("warnings", []).append(
                f"Failed to read molecule .mol file: {path} ({type(exc).__name__}: {exc})"
            )
    return None


def _maybe_attach_input_output_molecules(calc: Calculation) -> None:
    """
    Attach standardized molecules for downstream analysis:
      - calc.data["input_molecule"]
      - calc.data["output_molecule"]

    Conventions:
      - input_molecule: best-effort from input artifacts
      - output_molecule: best-effort from outputs; falls back to input_molecule
      - calc.molecule: "current" geometry, prefer output_molecule
    """
    import numpy as np

    input_mol = calc.data.get("input_molecule")
    output_mol = calc.data.get("output_molecule")

    def _canon(symbols, geometry):
        from gecko.molecule.canonical import canonicalize_atom_order

        symbols_sorted, geometry_sorted = canonicalize_atom_order(
            list(symbols), np.asarray(geometry, dtype=float).reshape(-1, 3), decimals=10
        )
        return symbols_sorted, np.asarray(geometry_sorted, dtype=float).reshape(-1, 3)

    # -------------------------
    # input_molecule
    # -------------------------
    if input_mol is None:
        if calc.code == "madness":
            input_json = calc.artifacts.get("input_json")
            if isinstance(input_json, Path) and input_json.exists():
                try:
                    import json
                    import qcelemental as qcel

                    data = json.loads(input_json.read_text(encoding="utf-8"))
                    mol = data.get("molecule")
                    if isinstance(mol, dict):
                        symbols = mol.get("symbols")
                        geometry = mol.get("geometry")
                        if symbols is not None and geometry is not None:
                            units = mol.get("units") or mol.get("parameters", {}).get("units")
                            coords = np.asarray(geometry, dtype=float).reshape(-1, 3)
                            if isinstance(units, str) and units.lower() in ("bohr", "atomic", "au"):
                                coords = coords * qcel.constants.bohr2angstroms
                            sym, coords = _canon(symbols, coords)
                            input_mol = qcel.models.Molecule(symbols=sym, geometry=coords)
                            calc.meta.setdefault("input_molecule_source", "input.json")
                            calc.meta.setdefault("input_molecule_path", str(input_json))
                except Exception as exc:
                    calc.meta.setdefault("warnings", []).append(
                        f"Failed to read input molecule from input.json: {input_json} ({type(exc).__name__}: {exc})"
                    )

            if input_mol is None:
                from gecko.mol.io import read_mol

                preferred = calc.root / "molecule.mol"
                candidates: list[Path] = []
                if preferred.exists():
                    candidates.append(preferred)
                candidates.extend(sorted(p for p in calc.root.glob("*.mol") if p != preferred))
                for p in candidates:
                    try:
                        input_mol = read_mol(p)
                        calc.meta.setdefault("input_molecule_source", "mol")
                        calc.meta.setdefault("input_molecule_path", str(p))
                        break
                    except Exception:
                        continue

            # Fallback: if the calc payload includes an embedded molecule block (common in MADQC calc_info),
            # treat that as the input molecule when no dedicated input artifact exists.
            if input_mol is None:
                raw = calc.data.get("raw_json")
                if isinstance(raw, dict):
                    try:
                        import qcelemental as qcel

                        def _find_first_molecule_dict(obj):
                            if isinstance(obj, dict):
                                if isinstance(obj.get("symbols"), list) and obj.get("geometry") is not None:
                                    return obj
                                for v in obj.values():
                                    r = _find_first_molecule_dict(v)
                                    if r is not None:
                                        return r
                            elif isinstance(obj, list):
                                for it in obj:
                                    r = _find_first_molecule_dict(it)
                                    if r is not None:
                                        return r
                            return None

                        m = _find_first_molecule_dict(raw.get("tasks") or raw)
                        if isinstance(m, dict):
                            units = m.get("units") or m.get("parameters", {}).get("units")
                            coords = np.asarray(m.get("geometry"), dtype=float).reshape(-1, 3)
                            if isinstance(units, str) and units.lower() in ("bohr", "atomic", "au"):
                                coords = coords * qcel.constants.bohr2angstroms
                            sym, coords = _canon(m.get("symbols"), coords)
                            input_mol = qcel.models.Molecule(symbols=sym, geometry=coords)
                            calc.meta.setdefault("input_molecule_source", "raw_json")
                    except Exception:
                        pass

        elif calc.code == "dalton":
            from gecko.plugins.dalton.parse import read_dalton_mol

            preferred_mol: Path | None = None
            out_path = calc.artifacts.get("out")
            pairs = calc.artifacts.get("dalton_pairs")

            if isinstance(out_path, Path) and isinstance(pairs, list):
                for entry in pairs:
                    if not isinstance(entry, dict):
                        continue
                    if entry.get("out") == out_path and isinstance(entry.get("mol"), Path):
                        preferred_mol = entry["mol"]
                        calc.meta.setdefault("input_molecule_source", "dalton_pairs")
                        calc.meta.setdefault("input_molecule_path", str(preferred_mol))
                        break

            if preferred_mol is None and isinstance(pairs, list) and pairs:
                first = pairs[0]
                if isinstance(first, dict) and isinstance(first.get("mol"), Path):
                    preferred_mol = first["mol"]
                    calc.meta.setdefault("input_molecule_source", "dalton_pairs_first")
                    calc.meta.setdefault("input_molecule_path", str(preferred_mol))

            if preferred_mol is None:
                mols = calc.artifacts.get("dalton_mol_files")
                if isinstance(mols, list) and mols:
                    candidate = mols[0]
                    if isinstance(candidate, Path):
                        preferred_mol = candidate
                        calc.meta.setdefault("input_molecule_source", "dalton_mol_files")
                        calc.meta.setdefault("input_molecule_path", str(preferred_mol))

            if preferred_mol is not None and preferred_mol.exists():
                try:
                    input_mol = read_dalton_mol(preferred_mol)
                except Exception as exc:
                    calc.meta.setdefault("warnings", []).append(
                        f"Failed to read Dalton input molecule: {preferred_mol} ({type(exc).__name__}: {exc})"
                    )

    # --------------------------
    # output_molecule
    # --------------------------
    if output_mol is None:
        mol_from_data = calc.data.get("molecule")
        if mol_from_data is not None:
            output_mol = mol_from_data
            calc.meta.setdefault("output_molecule_source", "data.molecule")
        elif calc.molecule is not None:
            output_mol = calc.molecule
            calc.meta.setdefault("output_molecule_source", "calc.molecule")

    if output_mol is None and input_mol is not None:
        output_mol = input_mol
        calc.meta.setdefault("output_molecule_source", "input_fallback")

    if input_mol is not None:
        calc.data["input_molecule"] = input_mol
    if output_mol is not None:
        calc.data["output_molecule"] = output_mol

    calc.molecule = output_mol or input_mol or calc.molecule


def _finalize_calc(calc: Calculation) -> Calculation:
    """Finalize a Calculation after loading and parsing.
    Attach optional artifacts, ensure molecule is set, and populate meta.
    """
    _maybe_attach_calc_info(calc)

    if calc.molecule is None and calc.data.get("molecule") is not None:
        calc.molecule = calc.data.get("molecule")

    _maybe_attach_input_output_molecules(calc)

    if calc.molecule is None:
        calc.molecule = _maybe_load_molecule_from_calc_dir(calc)

    _maybe_attach_input_output_molecules(calc)

    if calc.molecule is None:
        calc.meta.setdefault("mol_source", "missing")
    else:
        calc.meta.setdefault("mol_source", "embedded")

    if calc.molecule is not None and calc.meta.get("molecule_id") is None:
        from gecko.ids import geom_id_from_molecule

        calc.meta["molecule_id"] = geom_id_from_molecule(calc.molecule)

    return calc


def load_calc(
    path: str | Path,
) -> Calculation:
    """
    Load and parse a calculation directory.

    Step 2 behavior:
    - Detect whether the directory is MADNESS or DALTON
    - Delegate to the appropriate plugin loader, which parses artifacts
    """
    root = Path(path).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Path does not exist: {root}")

    if madness_can_load(root):
        calc = load_madness(root)
        return _finalize_calc(calc)

    if dalton_can_load(root):
        calc = load_dalton(root)
        return _finalize_calc(calc)

    raise ValueError(
        "Could not detect calculation type (madness/dalton) from directory. "
        f"Path: {root}"
    )

