from __future__ import annotations

from pathlib import Path

import pytest

from gecko.core.load import load_calc


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


@pytest.mark.parametrize(
    "relpath",
    [
        "01_mra-d04_n2",
        "02_aug-cc-pVDZ_n2",
        "03_mra-raman_h2o",
        "05_dalton_raman_h2o",
    ],
)
def test_ground_state_energy_present(relpath: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    e = calc.meta.get("ground_state_energy")
    assert isinstance(e, (int, float))


def test_ground_state_energy_optional(relpath: str = "00_mra-d06_bh2cl") -> None:
    calc = load_calc(FIXTURES / relpath)
    e = calc.meta.get("ground_state_energy")
    assert e is None or isinstance(e, (int, float))

