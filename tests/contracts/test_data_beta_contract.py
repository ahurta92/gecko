from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import logging


from gecko.core.load import load_calc


logger=logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"
_COMP_RE = re.compile(r"^[xyz]{3}$")


def _assert_beta_schema(beta: dict) -> None:
    assert isinstance(beta, dict)
    assert set(beta.keys()) >= {"omega", "components", "values", "shape"}
    assert beta["shape"] == ("freq", "component")

    omega = np.asarray(beta["omega"], dtype=float)
    values = np.asarray(beta["values"], dtype=float)
    components = list(beta["components"])

    assert omega.ndim == 2
    assert omega.shape[1] == 3
    assert values.ndim == 2
    assert values.shape[0] == omega.shape[0]
    assert values.shape[1] == len(components)

    assert components == sorted(components)
    assert len(components) > 0
    for comp in components:
        assert isinstance(comp, str)
        assert _COMP_RE.match(comp) is not None
        


def test_beta_contract_madness() -> None:
    calc = load_calc(FIXTURES / "01_mra-d04_n2")
    assert calc.code == "madness"
    beta = calc.data.get("beta")
    assert beta, "Expected beta data to be present for this fixture"
    _assert_beta_schema(beta)

    logger.info(f"Test beta data for madness fixture: omega shape {beta['omega'].shape}, values shape {beta['values'].shape}, components {beta['components']}")
    logger.info(f"Beta omega:\n{beta['omega']}\nBeta components:\n{beta['components']}\nBeta values:\n{beta['values']}")


def test_beta_contract_dalton() -> None:
    calc = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")
    assert calc.code == "dalton"
    beta = calc.data.get("beta")
    assert beta, "Expected beta data to be present for this fixture"
    _assert_beta_schema(beta)
    logger.info(f"Test beta data for madness fixture: omega shape {beta['omega'].shape}, values shape {beta['values'].shape}, components {beta['components']}")
    logger.info(f"Beta omega:\n{beta['omega']}\nBeta components:\n{beta['components']}\nBeta values:\n{beta['values']}")


