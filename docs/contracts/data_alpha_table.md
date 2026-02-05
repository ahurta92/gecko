# Contract: Alpha table (TableBuilder)

This contract defines the standard, long-form table emitted by `TableBuilder.build_alpha()`.

## Representation

`TableBuilder.build_alpha()` returns a Pandas DataFrame in long format with at least the following columns:

Envelope columns (from `make_envelope`):
- `calc_id`, `geom_id`, `mol_id`, `molecule_id`, `label`, `code`, `root`, `basis`, `method`

Alpha columns:
- `omega`: frequency (au)
- `ij`: tensor component (`xx`, `xy`, ...)
- `value`: component value

## Types

- `omega` and `value` must be convertible to float.
- `ij` must be a lowercase 2-letter component.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update table extraction accordingly.
