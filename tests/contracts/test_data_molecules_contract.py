from __future__ import annotations

from pathlib import Path

import numpy as np

from gecko.core.load import load_calc


FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def _sig(mol) -> tuple[list[str], np.ndarray]:
    symbols = [str(s) for s in mol.symbols]
    geom = np.asarray(mol.geometry, dtype=float).reshape(-1, 3)
    return symbols, geom


def _assert_same_geometry(m1, m2, *, atol: float = 1e-12) -> None:
    s1, g1 = _sig(m1)
    s2, g2 = _sig(m2)
    assert s1 == s2
    assert g1.shape == g2.shape
    assert np.allclose(g1, g2, atol=atol, rtol=0.0)


def test_input_molecule_same_for_n2_fixtures() -> None:
    mra = load_calc(FIXTURES / "01_mra-d04_n2")
    bs = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")

    in1 = mra.data.get("input_molecule")
    in2 = bs.data.get("input_molecule")
    assert in1 is not None
    assert in2 is not None
    _assert_same_geometry(in1, in2)


def test_output_molecule_exists_and_matches_input_for_n2_fixtures() -> None:
    mra = load_calc(FIXTURES / "01_mra-d04_n2")
    bs = load_calc(FIXTURES / "02_aug-cc-pVDZ_n2")

    for calc in (mra, bs):
        inp = calc.data.get("input_molecule")
        out = calc.data.get("output_molecule")
        assert inp is not None
        assert out is not None
        _assert_same_geometry(inp, out)

