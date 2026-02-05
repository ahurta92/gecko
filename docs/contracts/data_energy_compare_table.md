# Contract: Energy comparison table (TableBuilder)

This contract defines the standard table emitted by `TableBuilder.compare_energy()`.

## Representation

`TableBuilder.compare_energy(ref_basis=..., keys=...)` returns a wide table with:
- Key columns (default: `mol_id`, `method`)
- One column per basis containing the energy value
- Delta columns `delta_<basis>` for each basis vs `ref_basis`
- Relative delta columns `rel_<basis>` for each basis vs `ref_basis`

## Types

- Energy and delta columns must be convertible to float.
- Relative deltas may be `NaN` when reference energy is near zero.

## Fixture-backed example (must pass)

- `tests/fixtures/load_calc/01_mra-d04_n2` and `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update comparison logic accordingly.
