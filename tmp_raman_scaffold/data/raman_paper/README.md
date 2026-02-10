# Raman Project Workflow

This directory contains the setup and automated Slurm submission workflow for Raman calculations.

## Layout

- `molecules/`: source molecule files (`*.mol`) in MADNESS geometry format.
- `inputs/`: shared templates.
  - `mad.raman.in`
  - `optimize.dal`
  - `raman.dal`
- `data/`: generated per-molecule run directories.
  - `data/<mol>/mra-p07/` for MADNESS
  - `data/<mol>/<basis>/` for DALTON
- `slurm_scripts/`
  - `run_madness.sh`
  - `run_dalton.sh`
  - `submit_raman_jobs.py`
  - `use_gecko_env.sh`
- `progress_tracker.ipynb`: status dashboard + plots.

## 1. Generate Run Directories

```bash
python /home/ahurta92/Projects/project_data/01_dataset_utils/generate_raman_molecule_files.py
```

What this does:

- Reads molecules from `molecules/*.mol`.
- Builds `data/<mol>/...` run trees.
- Writes DALTON molecule/input files for configured basis sets.
- Writes MADNESS input files (`<mol>_raman.in`) with molecule block appended.
- Validates that frequencies in `inputs/mad.raman.in` and all `.FREQUE` blocks in `inputs/raman.dal` are identical.

Defaults:

- MRA dir: `mra-p07`
- Basis sets: `aug-cc-pVDZ`, `aug-cc-pVTZ`, `aug-cc-pVQZ`
- MADNESS `eprec`: `1e-6`

Useful flags:

- `--dry-run`
- `--overwrite`
- `--basis-set aug-cc-pVDZ aug-cc-pVTZ`
- `--output-root /path/to/custom/root`

## 2. Load Gecko Environment

```bash
cd /home/ahurta92/Projects/project_data/data/raman_paper
source slurm_scripts/use_gecko_env.sh /path/to/gecko
```

If gecko is already importable in your active Python environment, this step is optional.

## 3. Plan and Submit Jobs

Dry-run (recommended first):

```bash
python slurm_scripts/submit_raman_jobs.py --max-submit 20
```

Submit real jobs:

```bash
python slurm_scripts/submit_raman_jobs.py --submit
```

The DALTON flow is stage-aware:

1. `optimize.dal` is submitted first if needed.
2. Once optimize output exists, final geometry is parsed and `opt_<mol>_<basis>.mol` is written automatically.
3. `raman.dal` is then submitted with the optimized molecule.

Useful flags:

- `--molecules H2O CH4 NH3`
- `--bases aug-cc-pVDZ`
- `--skip-madness`
- `--skip-dalton`
- `--force`
- `--gecko-src /path/to/gecko/src`

## 4. Track Progress

Open:

- `progress_tracker.ipynb`

It reports statuses for:

- `madness_raman`
- `dalton_optimize`
- `dalton_raman`

and includes summary plots plus an optional submit cell.

## Legacy Layout Note

If you still use the legacy layout (`raman_paper/<mol>/...` without `data/`), the submit script can still run by pointing `--data-root` to the legacy root:

```bash
python slurm_scripts/submit_raman_jobs.py --data-root .
```
