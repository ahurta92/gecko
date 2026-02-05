from __future__ import annotations

from pathlib import Path
import logging

import numpy as np

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def test_table_alpha_compare_long_contract() -> None:
    calc_a = load_calc(FIXTURES / "03_mra-raman_h2o")
    calc_b = load_calc(FIXTURES / "05_dalton_raman_h2o")

    tb = TableBuilder([calc_a, calc_b])
    df = tb.compare_alpha_long(ref_basis="mra-p07", keys=["mol_id", "omega", "ij"])
    assert not df.empty

    required_cols = {"mol_id", "omega", "ij", "ref_basis", "basis", "value", "ref_value", "delta", "rel"}
    assert set(df.columns) >= required_cols

    for col in ["value", "ref_value", "delta"]:
        vals = df[col].astype(float)
        assert np.isfinite(vals).all()

    view_cols = ["mol_id", "omega", "ij", "ref_basis", "basis", "value", "ref_value", "delta", "rel"]
    df_view = df[view_cols].sort_values(["omega", "ij", "basis"]).reset_index(drop=True)
    logger.info("Alpha comparison long table rows=%d", len(df_view))
    logger.info("Alpha comparison long table:\n%s", df_view.to_string(index=False))
