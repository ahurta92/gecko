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
    ],
)
def test_table_beta_contract(relpath: str, code: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    assert calc.code == code

    tb = TableBuilder([calc])
    df = tb.build_beta()
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
        "omegaA",
        "omegaB",
        "omegaC",
        "ijk",
        "value",
    }
    assert set(df.columns) >= required_cols

    assert (df["code"] == code).all()

    vals = df["value"].astype(float)
    assert np.isfinite(vals).all()

    view_cols = ["code", "basis", "omegaA", "omegaB", "omegaC", "ijk", "value"]
    df_view = df[view_cols].sort_values(["omegaA", "omegaB", "omegaC", "ijk"]).reset_index(drop=True)
    logger.info("Beta table %s rows=%d", relpath, len(df))
    logger.info("Beta table (key columns):\n%s", df_view.to_string(index=False))
