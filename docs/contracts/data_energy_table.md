# Contract: Energy table (TableBuilder)

This contract defines the standard table emitted by `TableBuilder.build_energy()`.

## Representation

`TableBuilder.build_energy()` returns a Pandas DataFrame in long format with at least the following columns:

Envelope columns (from `make_envelope`):
- `calc_id`, `geom_id`, `mol_id`, `molecule_id`, `label`, `code`, `root`, `basis`, `method`

Energy columns:
- `energy`: ground-state energy (float, native units)

## Types

- `energy` must be convertible to float.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/01_mra-d04_n2`
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`
- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update table extraction accordingly.
