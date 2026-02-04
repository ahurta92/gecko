from __future__ import annotations

from pathlib import Path

import pytest
import logging

from gecko.core.load import load_calc


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"

logger=logging.getLogger(__name__)


def _assert_common(calc) -> None:
    assert calc.code in {"madness", "dalton"}
    assert isinstance(calc.root, Path)
    assert calc.root.exists()
    assert isinstance(calc.artifacts, dict)
    assert isinstance(calc.data, dict)
    assert isinstance(calc.meta, dict)

    if calc.molecule is None:
        assert calc.meta.get("mol_source") == "missing"
    else:
        assert calc.meta.get("molecule_id") is not None


@pytest.mark.parametrize(
    "relpath, expected_code",
    [
        ("00_mra-d06_bh2cl", "madness"),
        ("01_mra-d04_n2", "madness"),
        ("02_aug-cc-pVDZ_n2", "dalton"),
        ("03_mra-raman_h2o", "madness"),
        ("05_dalton_raman_h2o", "dalton"),
    ],
)
def test_load_calc_directory_contract(relpath: str, expected_code: str) -> None:
    root = FIXTURES / relpath
    calc = load_calc(root)
    _assert_common(calc)
    assert calc.code == expected_code


def test_load_calc_file_contract_dalton_out() -> None:
    out_path = FIXTURES / "02_aug-cc-pVDZ_n2"
    calc = load_calc(out_path)
    _assert_common(calc)
    assert calc.code == "dalton"
    assert calc.root == out_path


def test_load_calc_file_contract_madness_outputs_json() -> None:
    out_path = FIXTURES / "00_mra-d06_bh2cl" 
    calc = load_calc(out_path)
    _assert_common(calc)
    assert calc.code == "madness"
    assert calc.root == out_path

def test_load_calc_dalton_pairs_discovered() -> None:
    logger.info("Testing DALTON pairs discovery")

    root = FIXTURES / "05_dalton_raman_h2o"
    calc = load_calc(root)
    _assert_common(calc)
    assert calc.code == "dalton"

    dal_files = calc.artifacts.get("dalton_dal_files")
    mol_files = calc.artifacts.get("dalton_mol_files")
    out_files = calc.artifacts.get("dalton_out_files")
    pairs = calc.artifacts.get("dalton_pairs")


    logger.info(f"DALTON pairs discovered: ")
    if pairs is not None:
        logger.info(f" Total pairs found: {len(pairs)}")
        for p in pairs:
            logger.info(f" {p['dal'].name}, {p['mol'].name}, {p['out'].name}")

    assert isinstance(dal_files, list) and len(dal_files) >= 2
    assert isinstance(mol_files, list) and len(mol_files) >= 2
    assert isinstance(out_files, list) and len(out_files) >= 2
    assert isinstance(pairs, list) and len(pairs) >= 2

    pair_names = {(p["dal"].name, p["mol"].name, p["out"].name) for p in pairs}
    assert (
        "optimize.dal",
        "H2O_d-aug-cc-pV6Z.mol",
        "optimize_H2O_d-aug-cc-pV6Z.out",
    ) in pair_names
    assert (
        "raman.dal",
        "opt_H2O_d-aug-cc-pV6Z.mol",
        "raman_opt_H2O_d-aug-cc-pV6Z.out",
    ) in pair_names


def test_load_calc_madqc_requires_paired_input_in_artifact() -> None:
    root = FIXTURES / "03_mra-raman_h2o"
    calc = load_calc(root)
    _assert_common(calc)
    assert calc.code == "madness"
    assert calc.artifacts.get("calc_info_json") is not None
    assert calc.artifacts.get("input_in") is not None
    assert Path(calc.artifacts["input_in"]).name.endswith(".in")
