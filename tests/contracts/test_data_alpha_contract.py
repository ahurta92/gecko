from __future__ import annotations

from pathlib import Path

import pytest
import numpy as np
import logging 
logger = logging.getLogger(__name__)

from gecko.core.load import load_calc


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"

@pytest.mark.parametrize(
    "relpath, code, expected_basis",
    [
        ("03_mra-raman_h2o", "madness", "mra-p07"),
        ("05_dalton_raman_h2o","dalton", "d-aug-cc-pV6Z"),
    ]
)



def test_alpha_contract_dalton(relpath:str,code:str, expected_basis:str) -> None:

        calc = load_calc(FIXTURES / relpath)

        # Alpha data must exist and match the contract expectation for these fixtures.
        assert calc.code == code 
        assert calc.basis == expected_basis

        alpha = calc.data.get("alpha")
        assert isinstance(alpha, dict)
        assert set(alpha.keys()) >= {"omega", "components", "values", "shape"}
        assert alpha["shape"] == ("freq", "component")

        omega = np.asarray(alpha["omega"], dtype=float).reshape(-1)
        values = np.asarray(alpha["values"], dtype=float)
        components = list(alpha["components"])

        assert omega.ndim == 1
        assert values.ndim == 2
        assert values.shape[0] == omega.shape[0]
        assert values.shape[1] == 9

        assert components == ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]

    
        [logger.info(f"\n omega: {omi} values shape:\n {alpha_i.reshape(3,3)}\n") for omi, alpha_i in zip(omega, values)]












