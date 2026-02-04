# Contract: `calc.data["alpha"]` (Polarizability)

This contract defines the standard representation for polarizability data (`α`) when it is available on a loaded `Calculation`.

## Scope

Applies to calculations loaded via:
- `gecko.core.load.load_calc`

Non-goals (not part of this contract yet):
- Units and unit conversions
- Guarantee that every calculation has polarizability data

## Representation

If polarizability data is available, it must be stored as:

`calc.data["alpha"] = { ... }` with required keys:
- `omega`: array-like of shape `(n_freq,)` containing frequencies
- `components`: list of length `9` containing `["xx","xy","xz","yx","yy","yz","zx","zy","zz"]` (lowercase)
- `values`: array-like of shape `(n_freq, 9)` containing numeric values
- `shape`: the literal tuple `("freq", "component")`

### Types

- `omega` and `values` must be convertible to `numpy.ndarray` with `dtype=float`.
- Missing components may be represented as `NaN` in `values`.

## Fixture-backed examples (must pass)

- Dalton: `tests/fixtures/load_calc/05_dalton_raman_h2o`
- MRA: `tests/fixtures/load_calc/03_mra-raman_h2o`

## How to change this contract

If you want to include tensor symmetry, isotropic values, or complex polarizabilities:
1) Update this document.
2) Update tests.
3) Update parsers accordingly.

