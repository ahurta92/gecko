# Gecko

gecko is a modular Python toolkit for loading, comparing, and visualizing quantum-chemical response data across electronic structure codes.

## `gecko` src

```graphsql
src/gecko/
  __init__.py                 # exposes load_calc + a couple convenience imports
  core/
    model.py                  # Calculation (minimal), maybe Molecule later
    load.py                   # load_calc + registry
    iterators.py              # directory scanning helpers (read-only mode)
  plugins/
    madness/
      __init__.py
      detect.py               # can_load rules
      parse.py                # thin wrapper around migrated legacy parser
      legacy/                 # (copy of migration/parsers/madness.py + helpers)
    dalton/
      __init__.py
      detect.py
      parse.py
      legacy/                 # (copy of migration/parsers/dalton*.py)
  recipes/
    shg_csv.py                # replaces the notebook workflow gradually
  viz/
    unit_sphere.py            # core plotting / mapping helpers
    metrics.py                # field_error equivalents
    io.py                     # load shg_ijk.csv / dataframe helpers
  apps/
    trame_shg_viewer.py       # migration/viz/application.py
  workflows/
    templates/                # packaged templates
      dalton/
      madness/
    legacy/                   # quarantine: db/ and cli/ initially
```
