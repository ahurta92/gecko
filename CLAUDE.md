# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install (uv recommended)
uv venv && uv pip install -e .

# Run tests
pytest
pytest --cov                          # with coverage
pytest tests/contracts/test_table_beta_contract.py  # single test file

# Lint & format
ruff format src tests
ruff check src tests
mypy src

# CLI
gecko shg build --db /path/to/calcs --out /path/to/output

# Visualization apps
python -m gecko.viz.apps.beta_viewer --shg-csv data/csv_data/shg_ijk.csv
python -m gecko.viz.apps.raman_dashboard --db-dir /path/to/calcs
```

## Architecture

Gecko is a plugin-based Python toolkit for loading, comparing, and visualizing quantum-chemical response data from MADNESS and DALTON electronic structure codes.

### Plugin Dispatch

`load_calc(path)` (exported from `src/gecko/__init__.py`) is the main entry point. It auto-detects the calculation type (MADNESS vs DALTON) via `plugins/*/detect.py` and delegates loading to the appropriate `plugins/*/loader.py`. After loading, a finalization pipeline in `core/load.py` handles molecule resolution, artifact attachment, and metadata enrichment.

### Data Model

`Calculation` (in `core/model.py`) is the central container. `Calculation.data` holds raw legacy parser outputs; `.artifacts` holds discovered file paths. The design is migration-first — legacy parsers live in `plugins/madness/legacy/` and `plugins/dalton/legacy/`.

### Tables

`tables/builder.py` contains `TableBuilder`, which uses extractors from `tables/extractors.py` to pull response properties (beta, alpha, Raman, energy, timing) from loaded calculations into pandas DataFrames. `recipes/shg_csv.py` is a high-level workflow that builds beta tables from a directory of calculations.

### Visualization

Interactive web dashboards in `viz/apps/` use [Trame](https://trame.readthedocs.io/) with VTK rendering. `viz/state.py` manages app state; `viz/vtk_scene.py` handles VTK setup.

### Testing

Tests are contract-based in `tests/contracts/`. Fixtures are real calculation directories under `tests/fixtures/` — `01_mra-d04_n2/` for MADNESS and `02_aug-cc-pVDZ_n2/` for DALTON.
