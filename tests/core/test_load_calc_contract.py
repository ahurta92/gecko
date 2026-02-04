from __future__ import annotations

from pathlib import Path

import pytest

from gecko.core.load import load_calc


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


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
    out_path = FIXTURES / "02_aug-cc-pVDZ_n2" / "quad.out"
    calc = load_calc(out_path)
    _assert_common(calc)
    assert calc.code == "dalton"
    assert calc.root == out_path.parent.resolve()


def test_load_calc_file_contract_madness_outputs_json() -> None:
    json_path = FIXTURES / "00_mra-d06_bh2cl" / "outputs.json"
    calc = load_calc(json_path)
    _assert_common(calc)
    assert calc.code == "madness"
    assert calc.root == json_path.parent.resolve()

