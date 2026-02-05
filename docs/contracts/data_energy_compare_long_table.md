# Contract: Energy comparison table (long form)

This contract defines the long-form table emitted by `TableBuilder.compare_energy_long()`.

## Representation

`TableBuilder.compare_energy_long(ref_basis=..., keys=...)` returns a long table with:
- Key columns (default: `mol_id`, `method`)
- `ref_basis`: the reference basis label
- `basis`: basis label for the compared value
- `energy`: energy for `basis`
- `ref_energy`: energy for `ref_basis`
- `delta`: `energy - ref_energy`
- `rel`: relative delta (`delta / |ref_energy|`)

## Types

- Energy and delta columns must be convertible to float.
- `rel` may be `NaN` when reference energy is near zero.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/01_mra-d04_n2` and `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`
- `tests/fixtures/beta_data/*`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update comparison logic accordingly.
