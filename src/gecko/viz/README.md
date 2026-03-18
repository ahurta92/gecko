# gecko.viz

Reusable visualization helpers and the Trame-based Gecko apps.

## Run with an existing CSV

If running from the repo root, ensure the package is installed in your environment:

```bash
python -m pip install -e .
```

```bash
python -m gecko.viz.apps.beta_viewer --shg-csv data/csv_data/shg_ijk.csv
```

## Build from a database directory

```bash
python -m gecko.viz.apps.beta_viewer --db-dir /path/to/calcs
```

## Build a reusable bundle

```bash
python -m gecko.viz.apps.beta_viewer --db-dir /path/to/calcs --write-bundle ./beta_bundle
python -m gecko.viz.apps.beta_viewer --bundle-dir ./beta_bundle
```

## Raman dashboard

`gecko.viz.apps.raman_dashboard` provides a data dashboard for frequency-dependent
polarizability, polarizability-derivative trends, Raman intensity trends, and
basis-set error summaries relative to a selectable reference basis.

```bash
python -m gecko.viz.apps.raman_dashboard \
  --db-dir /gpfs/scratch/ahurtado/project_data/data/raman_paper/data
```

## Notes
- SHG omega indexing uses $\omega_B = \omega_C$ and starts at 0.
- MADNESS basis is labeled as "MRA".
