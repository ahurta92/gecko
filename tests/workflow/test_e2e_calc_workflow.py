"""End-to-end CLI workflow tests: calc init → submit → results.

Covers three molecules not present in the existing fixture set:
SO2, NH3, CO2.

- TestCalcInit    : verifies input files and SLURM scripts are generated
- TestCalcSubmit  : verifies job submission is recorded in the jobstore
- TestCalcResults : verifies property extraction using the BH2Cl alpha fixture
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import qcelemental as qcel

from gecko.cli import main

# ---------------------------------------------------------------------------
# Molecule fixtures  (not in existing fixture set)
# ---------------------------------------------------------------------------

_ANG2BOHR = qcel.constants.conversion_factor("angstrom", "bohr")


def _mol(symbols: list[str], coords_ang: list[list[float]]) -> qcel.models.Molecule:
    return qcel.models.Molecule(
        symbols=symbols,
        geometry=(np.array(coords_ang) * _ANG2BOHR).flatten().tolist(),
    )


_SO2 = _mol(
    ["S", "O", "O"],
    [[0.000,  0.000,  0.000],
     [1.432,  0.000, -0.572],
     [-1.432, 0.000, -0.572]],
)
_NH3 = _mol(
    ["N", "H", "H", "H"],
    [[0.000,  0.000,  0.116],
     [0.000,  0.940, -0.271],
     [0.814, -0.470, -0.271],
     [-0.814, -0.470, -0.271]],
)
_CO2 = _mol(
    ["C", "O", "O"],
    [[0.000, 0.000,  0.000],
     [0.000, 0.000,  1.162],
     [0.000, 0.000, -1.162]],
)

MOLECULES = [
    ("SO2", _SO2),
    ("NH3", _NH3),
    ("CO2", _CO2),
]

# Path to the BH2Cl alpha fixture (used for calc results tests)
_BHLCL_FIXTURE = (
    Path(__file__).parent.parent / "fixtures" / "load_calc" / "00_mra-d06_bh2cl"
)


# ---------------------------------------------------------------------------
# TestCalcInit
# ---------------------------------------------------------------------------


class TestCalcInit:
    @pytest.mark.parametrize("mol_name, mol", MOLECULES)
    def test_generates_madness_input_file(self, tmp_path, mol_name, mol):
        with patch("gecko.workflow.geometry.fetch_geometry", return_value=mol):
            rc = main([
                "calc", "init",
                "--molecule", mol_name,
                "--property", "alpha",
                "--tier", "medium",
                "--out", str(tmp_path / mol_name),
            ])
        assert rc == 0
        in_files = list((tmp_path / mol_name).rglob("*.in"))
        assert len(in_files) >= 1

    @pytest.mark.parametrize("mol_name, mol", MOLECULES)
    def test_input_has_dft_response_molecule_sections(self, tmp_path, mol_name, mol):
        with patch("gecko.workflow.geometry.fetch_geometry", return_value=mol):
            main([
                "calc", "init",
                "--molecule", mol_name,
                "--property", "alpha",
                "--tier", "medium",
                "--out", str(tmp_path / mol_name),
            ])
        in_file = next((tmp_path / mol_name).rglob("*.in"))
        content = in_file.read_text()
        assert "dft" in content
        assert "response" in content
        assert "molecule" in content

    @pytest.mark.parametrize("mol_name, mol", MOLECULES)
    def test_slurm_script_generated_with_manual_flags(self, tmp_path, mol_name, mol):
        with patch("gecko.workflow.geometry.fetch_geometry", return_value=mol):
            rc = main([
                "calc", "init",
                "--molecule", mol_name,
                "--property", "alpha",
                "--tier", "medium",
                "--out", str(tmp_path / mol_name),
                "--slurm",
                "--partition", "hbm-short-96core",
                "--nodes", "1",
                "--tasks-per-node", "2",
                "--walltime", "00:30:00",
            ])
        assert rc == 0
        scripts = list((tmp_path / mol_name).rglob("run_*.sh"))
        assert len(scripts) >= 1


# ---------------------------------------------------------------------------
# TestCalcSubmit
# ---------------------------------------------------------------------------


class TestCalcSubmit:
    def _init(self, tmp_path, mol_name, mol):
        """Helper: run calc init with SLURM scripts for a molecule."""
        with patch("gecko.workflow.geometry.fetch_geometry", return_value=mol):
            main([
                "calc", "init",
                "--molecule", mol_name,
                "--property", "alpha",
                "--tier", "medium",
                "--out", str(tmp_path / mol_name),
                "--slurm",
                "--partition", "hbm-short-96core",
                "--nodes", "1",
                "--tasks-per-node", "2",
                "--walltime", "00:30:00",
            ])
        return tmp_path / mol_name

    @pytest.mark.parametrize("mol_name, mol", MOLECULES)
    def test_submit_records_job_in_store(self, tmp_path, mol_name, mol):
        out_dir = self._init(tmp_path, mol_name, mol)

        sbatch = MagicMock(returncode=0, stdout="Submitted batch job 99999\n", stderr="")
        with patch("subprocess.run", return_value=sbatch):
            rc = main(["calc", "submit", str(out_dir)])
        assert rc == 0

        from gecko.workflow.jobstore import load_store
        records = load_store(out_dir).all()
        assert len(records) >= 1
        assert records[0].job_id == "99999"

    @pytest.mark.parametrize("mol_name, mol", MOLECULES)
    def test_submit_fails_gracefully_on_no_scripts(self, tmp_path, mol_name, mol):
        """If there are no run_*.sh scripts, submit should return non-zero."""
        out_dir = tmp_path / mol_name
        out_dir.mkdir()
        rc = main(["calc", "submit", str(out_dir)])
        assert rc != 0


# ---------------------------------------------------------------------------
# TestCalcResults
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _BHLCL_FIXTURE.exists(),
    reason="BH2Cl load_calc fixture not present",
)
class TestCalcResults:
    def test_alpha_table_printed_to_stdout(self, capsys):
        rc = main(["calc", "results", str(_BHLCL_FIXTURE), "--property", "alpha"])
        assert rc == 0
        out = capsys.readouterr().out
        assert len(out.strip()) > 0

    def test_alpha_csv_written_to_file(self, tmp_path):
        out_csv = tmp_path / "alpha.csv"
        rc = main([
            "calc", "results", str(_BHLCL_FIXTURE),
            "--property", "alpha",
            "--format", "csv",
            "--out", str(out_csv),
        ])
        assert rc == 0
        assert out_csv.exists()
        assert len(out_csv.read_text().strip()) > 0

    def test_csv_has_header_row(self, tmp_path):
        out_csv = tmp_path / "alpha.csv"
        main([
            "calc", "results", str(_BHLCL_FIXTURE),
            "--property", "alpha",
            "--format", "csv",
            "--out", str(out_csv),
        ])
        lines = out_csv.read_text().strip().splitlines()
        assert len(lines) >= 2  # header + at least one data row
