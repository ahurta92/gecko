from __future__ import annotations

from pathlib import Path
import logging

import numpy as np

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def test_table_beta_compare_long_contract() -> None:
    calc_a = load_calc(FIXTURES / "01_mra-d04_n2")
    calc_b = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")

    tb = TableBuilder([calc_a, calc_b])
    df = tb.compare_beta_long(ref_basis="mra-d04", keys=["mol_id", "omegaA", "omegaB", "omegaC", "ijk"])
    assert not df.empty

    required_cols = {"mol_id", "omegaA", "omegaB", "omegaC", "ijk", "ref_basis", "basis", "value", "ref_value", "delta", "rel"}
    assert set(df.columns) >= required_cols

    vals = df["value"].astype(float)
    assert np.isfinite(vals).all()

    view_cols = ["mol_id", "omegaA", "omegaB", "omegaC", "ijk", "ref_basis", "basis", "value", "ref_value", "delta", "rel"]
    df_view = df[view_cols].sort_values(["omegaA", "omegaB", "omegaC", "ijk", "basis"]).reset_index(drop=True)
    logger.info("Beta comparison long table rows=%d", len(df_view))
    logger.info("Beta comparison long table:\n%s", df_view.to_string(index=False))
