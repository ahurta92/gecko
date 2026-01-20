# Developer Notes

## Canonical molecule ordering

Gecko canonicalizes atom ordering before creating `qcelemental.models.Molecule` so that
molecule hashes are stable across codes (Dalton/MADNESS) when only atom order differs.

Ordering rule: sort by `(atomic_number, x, y, z)` where coordinates are rounded to
`decimals=10` for stable ordering under tiny floating-point noise. The original (unrounded)
coordinates are preserved in the final molecule.

This does **not** address:
- Translation or rotation differences between geometry representations
- Symmetry-related degeneracy issues

If matching must tolerate rigid transforms, implement a geometry key (e.g., distance-matrix
based) and do not rely solely on QCElemental hashes.

## Molecule identifiers

- `geom_id` now comes directly from `qcelemental.models.Molecule.get_hash()`.
- `mol_id` now comes from `qcelemental.models.Molecule.formula` (composition only).

This means:
- Molecule parsing + canonicalization is required for stable IDs across codes.
- `mol_id` does not distinguish isomers; it is composition-only.
- `geom_id` is geometry-sensitive and reflects exact coordinates after canonicalization.

## MADNESS basis metadata

For MADNESS calculations, `calc.meta["basis"]` is currently forced to "MRA".
This is a temporary label until we standardize MRA accuracy metadata (e.g.,
orbital residual + density thresholds).

## SHG pipeline findings (2026-01-20)

Findings:
- New SHG logic was introduced under `gecko.tables.shg` (now deprecated).
- A CLI path routed SHG generation through the tables module.
- Tutorial docs referenced the tables pipeline.

Introduced modules/files (now deprecated or redirected):
- src/gecko/tables/shg.py
- src/gecko/tables/__init__.py (exports)
- src/gecko/cli.py (SHG command path)
- docs/dev/shg_pipeline.md
- docs/tutorials/shg_pipeline.md

Resolution:
- SHG CSV generation is routed exclusively through `gecko.recipes.shg_csv.build_beta_table`.
- Calculation discovery uses `gecko.core.iterators.iter_calc_dirs`.
- SHG omega indexing starts at 0 with deterministic sorted frequency mapping.
- No molecule labels are used in SHG outputs (use `geom_id` + `mol_id` instead).
