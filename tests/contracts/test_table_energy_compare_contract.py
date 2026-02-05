from __future__ import annotations

from pathlib import Path
import logging

import numpy as np

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def test_table_energy_compare_contract() -> None:
    calc_a = load_calc(FIXTURES / "01_mra-d04_n2")
    calc_b = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")

    tb = TableBuilder([calc_a, calc_b])
    df = tb.compare_energy(ref_basis="mra-d04", keys=["mol_id"])
    assert not df.empty

    assert "mra-d04" in df.columns
    assert "aug-cc-pVDZ" in df.columns
    assert "delta_aug-cc-pVDZ" in df.columns
    assert "rel_aug-cc-pVDZ" in df.columns

    for col in ["mra-d04", "aug-cc-pVDZ", "delta_aug-cc-pVDZ"]:
        vals = df[col].astype(float)
        assert np.isfinite(vals).all()

    logger.info("Energy comparison table rows=%d", len(df))
    logger.info("Energy comparison table:\n%s", df.to_string(index=False))

    calc_c = load_calc(FIXTURES / "03_mra-raman_h2o")
    calc_d = load_calc(FIXTURES / "05_dalton_raman_h2o")

    tb_h2o = TableBuilder([calc_c, calc_d])
    df_h2o = tb_h2o.compare_energy(ref_basis="mra-p07", keys=["mol_id"])
    assert not df_h2o.empty

    assert "mra-p07" in df_h2o.columns
    assert "d-aug-cc-pV6Z" in df_h2o.columns
    assert "delta_d-aug-cc-pV6Z" in df_h2o.columns
    assert "rel_d-aug-cc-pV6Z" in df_h2o.columns

    for col in ["mra-p07", "d-aug-cc-pV6Z", "delta_d-aug-cc-pV6Z"]:
        vals = df_h2o[col].astype(float)
        assert np.isfinite(vals).all()

    logger.info("Energy comparison table (H2O) rows=%d", len(df_h2o))
    logger.info("Energy comparison table (H2O):\n%s", df_h2o.to_string(index=False))


def test_table_energy_compare_multibasis_beta_data() -> None:
    beta_root = FIXTURES.parent / "beta_data"
    mol_dirs = [p for p in sorted(beta_root.iterdir()) if p.is_dir()]

    for mol_dir in mol_dirs:
        calc_dirs = [p for p in sorted(mol_dir.iterdir()) if p.is_dir()]
        calcs = [load_calc(p) for p in calc_dirs]
        tb = TableBuilder(calcs)
        df = tb.compare_energy(ref_basis="mra-d06", keys=["mol_id"])
        assert not df.empty

        basis_cols = [
            c
            for c in df.columns
            if c not in {"mol_id"} and not c.startswith("delta_") and not c.startswith("rel_")
        ]
        if "mra-d06" not in df.columns and basis_cols:
            df = tb.compare_energy(ref_basis=basis_cols[0], keys=["mol_id"])
            basis_cols = [
                c
                for c in df.columns
                if c not in {"mol_id"} and not c.startswith("delta_") and not c.startswith("rel_")
            ]
        assert basis_cols
        assert len(basis_cols) >= 3

        logger.info("Energy comparison table (beta_data/%s) rows=%d", mol_dir.name, len(df))
        logger.info("Energy comparison table (beta_data/%s):\n%s", mol_dir.name, df.to_string(index=False))
