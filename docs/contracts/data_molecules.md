# Contract: `calc.data["input_molecule"]` / `calc.data["output_molecule"]`

This contract defines how Gecko represents **starting** and **final** geometries on a loaded `Calculation`.

Motivation:
- Basis-set comparisons should be performed between calculations that share the **same input geometry**.
- Some calculations are single-step (no geometry change): then `output_molecule == input_molecule`.
- Some calculations are multi-step (e.g. optimize → raman): the final property should be associated with the **initial geometry**.

## Scope

Applies to calculations loaded via:
- `gecko.core.load.load_calc`

Non-goals (not part of this contract yet):
- Full workflow step graph representation
- Storing *all* intermediate geometries

## Representation

If Gecko can determine the geometry from inputs/outputs, it stores:

- `calc.data["input_molecule"]`: a `qcelemental.models.Molecule` for the starting geometry
- `calc.data["output_molecule"]`: a `qcelemental.models.Molecule` for the final geometry

Selection rules (best-effort):

1) `input_molecule` is derived from input artifacts:
   - MADNESS: `input.json` molecule block if present, else a `.mol` file in the directory
   - Dalton: the `.mol` paired to the chosen output (via `{dal_stem}_{mol_stem}.out`) if possible, else a `.mol` file in the directory

2) `output_molecule` is derived from parsed outputs if present; otherwise it falls back to `input_molecule`.

3) `calc.molecule` is the “current geometry” for downstream analysis, and should match `output_molecule` when available.

## Fixture-backed examples (must pass)

Starting geometry consistency (basis comparisons):
- `tests/fixtures/load_calc/01_mra-d04_n2` and `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`
  must produce equal `input_molecule` (same symbols and coordinates within tolerance).

Single-step behavior:
- In the same fixtures above, `output_molecule` must exist and match `input_molecule`.

Multi-step behavior (optimize → property):
- `tests/fixtures/load_calc/05_dalton_raman_h2o`:
  - `input_molecule` must exist and come from the starting `.mol`
  - `output_molecule` must exist and represent the optimized geometry (it may differ from `input_molecule`)
  - `calc.molecule` should match `output_molecule` (final geometry for downstream analysis)

## How to change this contract

If you want to add support for intermediate steps (optimize geometry distinct from raman geometry, etc.):
1) Update this document.
2) Update tests.
3) Update loader/parser logic accordingly.
