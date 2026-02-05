from __future__ import annotations

from pathlib import Path
import logging

import numpy as np

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def test_table_energy_compare_long_contract() -> None:
    calc_a = load_calc(FIXTURES / "01_mra-d04_n2")
    calc_b = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")

    tb = TableBuilder([calc_a, calc_b])
    df = tb.compare_energy_long(ref_basis="mra-d04", keys=["mol_id"])
    assert not df.empty

    required_cols = {"mol_id", "ref_basis", "basis", "energy", "ref_energy", "delta", "rel"}
    assert set(df.columns) >= required_cols

    for col in ["energy", "ref_energy", "delta"]:
        vals = df[col].astype(float)
        assert np.isfinite(vals).all()

    logger.info("Energy comparison long table rows=%d", len(df))
    view_cols = ["mol_id", "ref_basis", "basis", "energy", "ref_energy", "delta", "rel"]
    logger.info("Energy comparison long table:\n%s", df[view_cols].to_string(index=False))


def test_table_energy_compare_long_multibasis_beta_data() -> None:
    beta_root = FIXTURES.parent / "beta_data"
    mol_dirs = [p for p in sorted(beta_root.iterdir()) if p.is_dir()]

    for mol_dir in mol_dirs:
        calc_dirs = [p for p in sorted(mol_dir.iterdir()) if p.is_dir()]
        calcs = [load_calc(p) for p in calc_dirs]
        tb = TableBuilder(calcs)

        ref_basis = "mra-d06"
        df = tb.compare_energy_long(ref_basis=ref_basis, keys=["mol_id"])
        if df.empty:
            # fallback to first available basis in this molecule set
            energies = tb.build_energy()
            if energies.empty:
                continue
            basis_cols = sorted(set(energies["basis"].dropna().tolist()))
            if not basis_cols:
                continue
            ref_basis = basis_cols[0]
            df = tb.compare_energy_long(ref_basis=ref_basis, keys=["mol_id"])

        assert not df.empty
        assert (df["ref_basis"] == ref_basis).all()

        logger.info("Energy comparison long table (beta_data/%s) rows=%d", mol_dir.name, len(df))
        view_cols = ["mol_id", "ref_basis", "basis", "energy", "ref_energy", "delta", "rel"]
        logger.info(
            "Energy comparison long table (beta_data/%s):\n%s",
            mol_dir.name,
            df[view_cols].to_string(index=False),
        )
