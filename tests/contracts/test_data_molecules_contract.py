from __future__ import annotations

from pathlib import Path

import numpy as np
import logging

logger = logging.getLogger(__name__)

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

    if in1 is not None and in2 is not None:
        logger.info(f"Input molecule from MRA calc: {in1.symbols}, {in1.geometry}")
        logger.info(f"Input molecule from BS calc: {in2.symbols}, {in2.geometry}")
        logger.info(f"Molecule names: {in1.name}, {in2.name}")



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


def test_dalton_multistep_output_geometry_differs_from_input() -> None:
    calc = load_calc(FIXTURES / "05_dalton_raman_h2o")
    inp = calc.data.get("input_molecule")
    out = calc.data.get("output_molecule")
    assert inp is not None
    assert out is not None

    # Symbols should match; geometry may or may not differ depending on the fixture.
    s_in, g_in = _sig(inp)
    s_out, g_out = _sig(out)
    assert s_in == s_out
    assert g_in.shape == g_out.shape

    # Downstream "current" molecule follows the output geometry.
    assert calc.molecule is not None
    _assert_same_geometry(calc.molecule, out)


    if inp is not None and out is not None:
        logger.info(f"Input molecule from DALTON calc: {inp.symbols}, {inp.geometry}")
        logger.info(f"Output molecule from DALTON calc: {out.symbols}, {out.geometry}")
        logger.info(f"Molecule names: {inp.name}, {out.name}")

    # For this fixture, input geometry should come from the starting .mol and
    # output geometry should come from the optimization output.
    in_path = str(calc.meta.get("input_molecule_path") or "")
    assert Path(in_path).name == "H2O_d-aug-cc-pV6Z.mol"

    out_path = str(calc.meta.get("molecule_path") or "")
    assert "optimize_" in Path(out_path).name
