# Contract: Beta comparison table (long form)

This contract defines the long-form table emitted by `TableBuilder.compare_beta_long()`.

## Representation

`TableBuilder.compare_beta_long(ref_basis=..., keys=...)` returns a long table with:
- Key columns (default: `mol_id`, `method`, `omegaA`, `omegaB`, `omegaC`, `ijk`)
- `ref_basis`: the reference basis label
- `basis`: basis label for the compared value
- `value`: beta component value for `basis`
- `ref_value`: beta component value for `ref_basis`
- `delta`: `value - ref_value`
- `rel`: relative delta (`delta / |ref_value|`)

## Types

- `value`, `ref_value`, and `delta` must be convertible to float.
- `rel` may be `NaN` when reference value is small.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/01_mra-d04_n2`
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update comparison logic accordingly.
