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
        ("03_mra-raman_h2o", "madness"),
        ("05_dalton_raman_h2o", "dalton"),
    ],
)
def test_table_raman_contract(relpath: str, code: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    assert calc.code == code

    tb = TableBuilder([calc])
    df = tb.build_raman()
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
        "omega_pol",
        "mode",
        "freq_cm1",
        "alpha2",
        "beta2",
        "pol_int",
        "depol_int",
        "dep_ratio",
        "alpha_iso",
        "dalpha_iso",
    }
    assert set(df.columns) >= required_cols

    assert (df["code"] == code).all()

    for col in ["omega_pol", "freq_cm1", "alpha2", "beta2", "pol_int", "depol_int", "dep_ratio"]:
        vals = df[col].astype(float)
        assert np.isfinite(vals).all()

    assert (df["mode"].astype(int) >= 1).all()

    logger.info("Raman table %s rows=%d", relpath, len(df))
    view_cols = [
        "code",
        "basis",
        "omega_pol",
        "mode",
        "freq_cm1",
        "alpha2",
        "beta2",
        "pol_int",
        "depol_int",
        "dep_ratio",
        "alpha_iso",
        "dalpha_iso",
    ]
    df_view = df[view_cols].sort_values(["mode", "freq_cm1", "omega_pol"])\
        .reset_index(drop=True)
    logger.info("Raman table (key columns):\n%s", df_view.to_string(index=False))
