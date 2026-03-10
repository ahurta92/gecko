# Gecko Analysis Workspace

Canonical location for analysis notebooks, submission scripts, and Slurm logs.

## Layout

- `notebooks/`
  - Analysis notebooks (H2O, frequency development, timing plots, etc.)
- `scripts/`
  - Reusable submit helpers (e.g. `run_madness_hbm96_medium.sh`)
- `slurm_logs/`
  - Centralized Slurm stdout/stderr
  - Organized as `slurm_logs/<partition>/<jobname>-<jobid>.out|err`

## Standard log policy

All updated project-data submit scripts now default to:

`/gpfs/projects/rjh/adrian/development/gecko/analysis/slurm_logs`

Override at submit time if needed:

```bash
export SLURM_LOG_ROOT=/path/to/custom/log/root
```

## Updated submit scripts

- `/gpfs/scratch/ahurtado/project_data/freq_dev/slurm_scripts/run_madness.sh`
- `/gpfs/scratch/ahurtado/project_data/response_scaling_dev/slurm_scripts/run_madness.sh`
- `/gpfs/scratch/ahurtado/project_data/data/raman_paper/slurm_scripts/run_madness.sh`
- `/gpfs/scratch/ahurtado/project_data/data/NLO/slurm_scripts/run_madness.sh`
- `/gpfs/scratch/ahurtado/project_data/freq_dev/slurm_scripts/run_dalton.sh`
- `/gpfs/scratch/ahurtado/project_data/response_scaling_dev/slurm_scripts/run_dalton.sh`
- `/gpfs/scratch/ahurtado/project_data/data/raman_paper/slurm_scripts/run_dalton.sh`
- `/gpfs/scratch/ahurtado/project_data/data/NLO/slurm_scripts/run_dalton.sh`
- `/gpfs/scratch/ahurtado/project_data/excited_state_dev/slurm_scripts/run_dalton.sh`
