# SHG Pipeline Example (Reuse-First Recipe)

This tutorial shows the **canonical, reuse-first** SHG workflow in Gecko. It uses the existing recipe entrypoint in `gecko.recipes.shg_csv` and the iterator utilities in `gecko.core.iterators`.

## 1) Input: Mixed MADNESS + Dalton calculations

Assume a directory named `calc_nlo_beta/` that contains a mix of MADNESS and Dalton calculations in nested subdirectories.

## 2) Build `shg_ijk.csv` using the recipe entrypoint

```python
from pathlib import Path
from gecko.core.iterators import iter_calc_dirs
from gecko.recipes.shg_csv import build_beta_table

root = Path("calc_nlo_beta")
calc_dirs = list(iter_calc_dirs(root))

shg_df = build_beta_table(
    calc_dirs,
    shg_only=True,
    add_shg_omega=True,
    shg_start_at=0,
    shg_tol=1e-12,
    include_geometry=True,
    app_compat=True,
    verbose=True,
)

print(shg_df.head())
```

What this does:
- Uses the iterator to discover MADNESS and Dalton calculations
- Builds a long-form SHG table from parsed beta data
- Filters to SHG rows where $|\omega_B - \omega_C| \leq 10^{-12}$
- Adds per-molecule `omega` indices starting at 0
- Adds legacy app columns (`molecule`, `Beta`) for visualization compatibility

## 3) Verify omega indexing (per molecule)

```python
shg_df[["mol_id", "omegaB", "omega"]].drop_duplicates().sort_values(
    ["mol_id", "omega"]
)
```

## 4) Write the CSV

```python
out = Path("data/csv_data")
out.mkdir(parents=True, exist_ok=True)

shg_df.to_csv(out / "shg_ijk.csv", index=False)
```

## 5) Use with the visualization app

Point the trame / PyVista app to the CSV:
- `data/csv_data/shg_ijk.csv`

If you are using the migration app, it reads that path by default. From the repo root:

```
python migration/viz/application.py
```

To use a custom location, pass the new CLI argument:

```
python migration/viz/application.py --shg-csv /path/to/shg_ijk.csv
```

The app expects these columns:
- `molecule` (formula) and `Beta`
- `basis`, `omega`, `ijk`

The recipe keeps the scientific columns as well:
- `geom_id`, `mol_id`, `code`
- `omegaA`, `omegaB`, `omegaC`
- `value`

No additional parsing is required once the CSV is generated.
