from __future__ import annotations

from pathlib import Path
import logging

import numpy as np
import pytest

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


@pytest.mark.parametrize(
    "relpath, code",
    [
        ("01_mra-d04_n2", "madness"),
        ("02_aug-cc-pVDZ_n2", "dalton"),
        ("03_mra-raman_h2o", "madness"),
        ("05_dalton_raman_h2o", "dalton"),
    ],
)
def test_table_energy_contract(relpath: str, code: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    assert calc.code == code

    tb = TableBuilder([calc])
    df = tb.build_energy()
    assert not df.empty

    required_cols = {
        "calc_id",
        "geom_id",
        "mol_id",
        "molecule_id",
        "label",
        "code",
        "root",
        "basis",
        "method",
        "energy",
    }
    assert set(df.columns) >= required_cols

    assert (df["code"] == code).all()

    vals = df["energy"].astype(float)
    assert np.isfinite(vals).all()

    view_cols = ["code", "basis", "energy"]
    df_view = df[view_cols].reset_index(drop=True)
    logger.info("Energy table %s rows=%d", relpath, len(df))
    logger.info("Energy table (key columns):\n%s", df_view.to_string(index=False))
