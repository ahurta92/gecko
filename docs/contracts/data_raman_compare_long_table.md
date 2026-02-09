# Contract: Raman comparison table (long form)

This contract defines the long-form table emitted by `TableBuilder.compare_raman_long()`.

## Representation

`TableBuilder.compare_raman_long(ref_basis=..., property_name=..., keys=...)` returns a long table with:
- Key columns (default: `mol_id`, `method`, `omega_pol`, `mode`)
- `ref_basis`: the reference basis label
- `basis`: basis label for the compared value
- `property`: Raman property name (e.g., `alpha2`, `beta2`, `pol_int`, `depol_int`, `dep_ratio`, `alpha_iso`, `dalpha_iso`)
- `value`: property value for `basis`
- `ref_value`: property value for `ref_basis`
- `freq_cm1`: basis vibrational frequency (cm⁻¹)
- `ref_freq_cm1`: reference vibrational frequency (cm⁻¹)
- `delta`: `value - ref_value`
- `rel`: relative delta (`delta / |ref_value|`)

## Types

- `value`, `ref_value`, and `delta` must be convertible to float.
- `rel` may be `NaN` when reference values are small.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update comparison logic accordingly.
