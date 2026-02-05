# Contract: Beta table (TableBuilder)

This contract defines the standard, long-form table emitted by `TableBuilder.build_beta()`.

## Representation

`TableBuilder.build_beta()` returns a Pandas DataFrame in long format with at least the following columns:

Envelope columns (from `make_envelope`):
- `calc_id`, `geom_id`, `mol_id`, `molecule_id`, `label`, `code`, `root`, `basis`, `method`

Beta columns:
- `omegaA`, `omegaB`, `omegaC`: frequency triplet (au)
- `ijk`: tensor component (`xxx`, `xxy`, ...)
- `value`: component value

## Types

- `omegaA`, `omegaB`, `omegaC`, and `value` must be convertible to float.
- `ijk` must be a lowercase 3-letter component.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/01_mra-d04_n2`
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update table extraction accordingly.
