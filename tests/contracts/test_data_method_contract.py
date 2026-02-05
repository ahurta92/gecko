from __future__ import annotations

from pathlib import Path

import pytest

from gecko.core.load import load_calc

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


@pytest.mark.parametrize(
    "relpath, code, expected_method",
    [
        ("03_mra-raman_h2o", "madness", "HF"),
        ("02_aug-cc-pVDZ_n2", "dalton", "HF"),
        ("05_dalton_raman_h2o", "dalton", "HF"),
    ],
)
def test_method_contract(relpath: str, code: str, expected_method: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    assert calc.code == code
    method = calc.meta.get("method")
    assert method, "Expected method to be present"
    assert str(method).upper() == expected_method
