# Contract: Alpha comparison table (long form)

This contract defines the long-form table emitted by `TableBuilder.compare_alpha_long()`.

## Representation

`TableBuilder.compare_alpha_long(ref_basis=..., keys=...)` returns a long table with:
- Key columns (default: `mol_id`, `method`, `omega`, `ij`)
- `ref_basis`: the reference basis label
- `basis`: basis label for the compared value
- `value`: alpha component value for `basis`
- `ref_value`: alpha component value for `ref_basis`
- `delta`: `value - ref_value`
- `rel`: relative delta (`delta / |ref_value|`)

## Types

- `value`, `ref_value`, and `delta` must be convertible to float.
- `rel` may be `NaN` when reference value is near zero.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update comparison logic accordingly.
