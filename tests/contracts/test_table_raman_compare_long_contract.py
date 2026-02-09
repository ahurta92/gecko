from __future__ import annotations

from pathlib import Path
import logging

import numpy as np

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def _run_property(property_name: str) -> None:
    calc_a = load_calc(FIXTURES / "03_mra-raman_h2o")
    calc_b = load_calc(FIXTURES / "05_dalton_raman_h2o")

    tb = TableBuilder([calc_a, calc_b])
    df = tb.compare_raman_long(ref_basis="mra-p07", property_name=property_name, keys=["mol_id", "omega_pol", "mode"])
    assert not df.empty

    required_cols = {"mol_id", "omega_pol", "mode", "freq_cm1", "ref_freq_cm1", "ref_basis", "basis", "property", "value", "ref_value", "delta", "rel"}
    assert set(df.columns) >= required_cols

    vals = df["value"].astype(float)
    assert np.isfinite(vals.dropna()).all()

    view_cols = ["mol_id", "omega_pol", "mode", "freq_cm1", "ref_freq_cm1", "ref_basis", "basis", "property", "value", "ref_value", "delta", "rel"]
    df_view = df[view_cols].sort_values(["omega_pol", "mode", "basis"]).reset_index(drop=True)
    logger.info("Raman comparison long table (%s) rows=%d", property_name, len(df_view))
    logger.info("Raman comparison long table (%s):\n%s", property_name, df_view.to_string(index=False))


def test_table_raman_compare_long_alpha2() -> None:
    _run_property("alpha2")


def test_table_raman_compare_long_depol() -> None:
    _run_property("dep_ratio")
