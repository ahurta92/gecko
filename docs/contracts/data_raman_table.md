# Contract: Raman table (TableBuilder)

This contract defines the standard, long-form table emitted by `TableBuilder.build_raman()`.

## Scope

Applies to:
- MADNESS Raman fixtures
- Dalton Raman fixtures

## Representation

`TableBuilder.build_raman()` returns a Pandas DataFrame in long format with at least the following columns:

Envelope columns (from `make_envelope`):
- `calc_id`, `geom_id`, `mol_id`, `molecule_id`, `label`, `code`, `root`, `basis`, `method`

Raman columns:
- `omega_pol`: polarization frequency (au)
- `mode`: vibrational mode index (1-based)
- `freq_cm1`: vibrational frequency (cm⁻¹)
- `alpha2`
- `beta2`
- `pol_int`
- `depol_int`
- `dep_ratio`
- `alpha_iso`: isotropic polarizability at `omega_pol` (trace/3)
- `dalpha_iso`: isotropic polarizability derivative for `mode` at `omega_pol`

## Types

- Raman numeric columns must be convertible to float.
- `mode` must be convertible to int.
- `alpha_iso`/`dalpha_iso` may be `None` if not available.

## Fixture-backed examples (must pass)

- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update table extraction accordingly.
