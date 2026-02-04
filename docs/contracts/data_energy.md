# Contract: `calc.meta["ground_state_energy"]`

This contract defines the minimal, standardized representation for a calculation’s ground-state energy.

## Scope

Applies to calculations loaded via:
- `gecko.core.load.load_calc`

Non-goals (not part of this contract yet):
- Units (assumed “native” to the backend; typically Hartree)
- Guaranteed availability (some legacy fixtures may not include energy)

## Representation

If ground-state energy is available, it must be stored as:

- `calc.meta["ground_state_energy"]`: `float`

The key may be absent or set to `None` when unavailable.

## Fixture-backed examples (must pass)

Fixtures that must provide energy:
- `tests/fixtures/load_calc/01_mra-d04_n2`
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`
- `tests/fixtures/load_calc/03_mra-raman_h2o`
- `tests/fixtures/load_calc/05_dalton_raman_h2o`

Fixture where energy is allowed to be missing:
- `tests/fixtures/load_calc/00_mra-d06_bh2cl`

## How to change this contract

If you want a richer energy schema (method, reference, components, etc.):
1) Update this document.
2) Update tests.
3) Update parsers accordingly.

