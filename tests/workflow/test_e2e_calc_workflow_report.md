# BTS-66 End-to-End Workflow Test Report

**Date:** 2026-04-02
**Platform:** Linux, Python 3.12.13, pytest 9.0.2
**Test file:** `gecko/tests/workflow/test_e2e_calc_workflow.py`

---

## Summary

| Result | Count |
|--------|-------|
| Passed | 18    |
| Failed | 0     |
| Total  | 18    |
| Duration | 6.12s |

---

## Molecules Tested

All three molecules are new — not present in the existing fixture set.

| Molecule | Formula | Shape | Atoms |
|----------|---------|-------|-------|
| SO2 | Sulfur dioxide | Bent | S, O, O |
| NH3 | Ammonia | Pyramidal | N, H, H, H |
| CO2 | Carbon dioxide | Linear | C, O, O |

---

## Test Classes & Results

### `TestCalcInit` — `gecko calc init`

**Source:** `gecko/src/gecko/cli.py:_calc_init_command`
**Geometry mock:** `gecko.workflow.geometry.fetch_geometry`

| Test | SO2 | NH3 | CO2 |
|------|-----|-----|-----|
| `test_generates_madness_input_file` | PASSED | PASSED | PASSED |
| `test_input_has_dft_response_molecule_sections` | PASSED | PASSED | PASSED |
| `test_slurm_script_generated_with_manual_flags` | PASSED | PASSED | PASSED |

Verified that `calc init --tier medium` produces a valid `.in` file with `dft`, `response`, and `molecule` sections, and that `--slurm` produces a `run_*.sh` script.

---

### `TestCalcSubmit` — `gecko calc submit`

**Source:** `gecko/src/gecko/cli.py:_calc_submit_command`
**Jobstore:** `gecko/src/gecko/workflow/jobstore.py`
**sbatch mock:** `subprocess.run` → returns job ID `99999`

| Test | SO2 | NH3 | CO2 |
|------|-----|-----|-----|
| `test_submit_records_job_in_store` | PASSED | PASSED | PASSED |
| `test_submit_fails_gracefully_on_no_scripts` | PASSED | PASSED | PASSED |

Verified that a submitted job is recorded in the jobstore with the correct ID, and that the command exits non-zero when no `run_*.sh` scripts are found.

---

### `TestCalcResults` — `gecko calc results`

**Source:** `gecko/src/gecko/cli.py:_calc_results_command`
**Fixture used:** `gecko/tests/fixtures/load_calc/00_mra-d06_bh2cl/` (BH2Cl alpha data)

| Test | Result |
|------|--------|
| `test_alpha_table_printed_to_stdout` | PASSED |
| `test_alpha_csv_written_to_file` | PASSED |
| `test_csv_has_header_row` | PASSED |

Verified that `--format table` prints to stdout and `--format csv --out alpha.csv` writes a valid CSV with at least a header and one data row.

---

## Key Files Modified / Created

| File | Change |
|------|--------|
| `gecko/src/gecko/workflow/params.py` | Added `_FIXTURES_DIR`, `_load_tier()` |
| `gecko/src/gecko/mcp_server.py` | Removed duplicate `_load_tier`, now imports from `params` |
| `gecko/src/gecko/workflow/hpc.py` | Added `load_slurm_profile()` |
| `gecko/src/gecko/cli.py` | Added `--tier`, `--cluster` to `calc init`; added `calc results` subcommand |
| `gecko/tests/workflow/test_e2e_calc_workflow.py` | New — this test file |
