from __future__ import annotations

from pathlib import Path

import pytest

from gecko.core.load import load_calc
import logging 
logger = logging.getLogger(__name__)


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


@pytest.mark.parametrize(
    "relpath, expected_basis",
    [
        ("00_mra-d06_bh2cl", "mra-d06"),
        ("01_mra-d04_n2", "mra-d04"),
        ("02_aug-cc-pVDZ_n2", "aug-cc-pVDZ"),
        ("03_mra-raman_h2o", "mra-p07"),
        ("05_dalton_raman_h2o", "d-aug-cc-pV6Z"),
    ],
)
def test_basis_contract(relpath: str, expected_basis: str) -> None:
    calc = load_calc(FIXTURES / relpath)

    # Basis must exist and match the contract expectation for these fixtures.
    assert isinstance(calc.basis, str)
    assert calc.basis == expected_basis

    logger.info(f"Tested basis for {relpath}: {calc.basis}")
    

    # Calculations should be identifiable by input molecule + basis.
    assert calc.data.get("input_molecule") is not None
