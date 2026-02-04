from __future__ import annotations

from pathlib import Path
import logging

import numpy as np
import pytest

from gecko.core.load import load_calc

logger = logging.getLogger(__name__)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "load_calc"


def _assert_alpha_present(alpha: dict) -> None:
    assert isinstance(alpha, dict)
    assert set(alpha.keys()) >= {"omega", "components", "values", "shape"}
    assert alpha["shape"] == ("freq", "component")

    omega = np.asarray(alpha["omega"], dtype=float).reshape(-1)
    values = np.asarray(alpha["values"], dtype=float)
    components = list(alpha["components"])

    assert omega.ndim == 1
    assert values.ndim == 2
    assert values.shape[0] == omega.shape[0]
    assert values.shape[1] == 9
    assert components == ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]


def _assert_raman_schema(raman: dict) -> None:
    assert isinstance(raman, dict)
    assert set(raman.keys()) >= {
        "polarization_frequencies",
        "vibrational_frequencies",
        "polarizability_derivatives",
        "polarizability_derivatives_by_mode",
        "raman_by_freq",
    }

    pol_freqs = np.asarray(raman["polarization_frequencies"], dtype=float).reshape(-1)
    vib_freqs = np.asarray(raman["vibrational_frequencies"], dtype=float).reshape(-1)
    assert pol_freqs.ndim == 1
    assert vib_freqs.ndim == 1
    assert len(pol_freqs) > 0
    assert len(vib_freqs) > 0

    pder = raman["polarizability_derivatives"]
    if isinstance(pder, list):
        assert len(pder) > 0
        for item in pder:
            np.asarray(item, dtype=float)
    else:
        np.asarray(pder, dtype=float)

    raman_by_freq = raman["raman_by_freq"]
    assert isinstance(raman_by_freq, dict)
    assert len(raman_by_freq) > 0

    required_keys = {
        "mode",
        "freq_cm1",
        "alpha2",
        "beta2",
        "pol_int",
        "depol_int",
        "dep_ratio",
    }

    freq_keys = []
    for freq_key, rows in raman_by_freq.items():
        freq_keys.append(float(freq_key))
        assert isinstance(rows, list)
        assert len(rows) > 0
        for row in rows:
            assert isinstance(row, dict)
            assert set(row.keys()) >= required_keys
            assert int(row["mode"]) >= 1
            float(row["freq_cm1"])
            float(row["alpha2"])
            float(row["beta2"])
            float(row["pol_int"])
            float(row["depol_int"])
            float(row["dep_ratio"])

    if len(freq_keys) == len(pol_freqs):
        assert np.allclose(sorted(freq_keys), sorted(pol_freqs))


@pytest.mark.parametrize(
    "relpath, code",
    [
        ("03_mra-raman_h2o", "madness"),
        ("05_dalton_raman_h2o", "dalton"),
    ],
)
def test_raman_contract(relpath: str, code: str) -> None:
    calc = load_calc(FIXTURES / relpath)
    assert calc.code == code

    raman = calc.data.get("raman")
    assert raman, "Expected raman data to be present for this fixture"
    _assert_raman_schema(raman)

    alpha = calc.data.get("alpha")
    assert alpha, "Expected polarizability (alpha) data for Raman fixtures"
    _assert_alpha_present(alpha)

    logger.info(
        "Raman contract fixture %s: pol_freqs=%s vib_freqs=%s rows=%s",
        relpath,
        np.asarray(raman["polarization_frequencies"]).shape,
        np.asarray(raman["vibrational_frequencies"]).shape,
        len(raman["raman_by_freq"]),
    )

    logger.info(
        f"Raman data for {relpath}:\n"
        f"Polarization frequencies:\n{raman['polarization_frequencies']}\n"
        f"Vibrational frequencies:\n{raman['vibrational_frequencies']}\n"
        f"Polarizability derivatives:\n{raman['polarizability_derivatives']}\n"
        f"Raman by frequency:\n{raman['raman_by_freq']}"
    )





