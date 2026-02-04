# Contract: `calc.data["beta"]` (Hyperpolarizability)

This contract defines the **standard, code-agnostic shape** for hyperpolarizability data (`β`) stored on a loaded `Calculation`.

It is enforced by tests using fixtures in `tests/fixtures/load_calc/`.

## Scope

Applies to:
- MADNESS runs that provide hyperpolarizability (legacy `output.json`/`outputs.json`, or MADQC `*.calc_info.json`)
- Dalton quadratic-response runs that print `beta(...)` lines

Non-goals (not part of this contract yet):
- Units and unit conversions
- Complex-valued tensors (currently stored as real floats)
- Guarantee that every run has beta data (absence is allowed and should be handled upstream)

## Representation

If beta data is available, it must be stored as:

`calc.data["beta"] = { ... }` with required keys:
- `omega`: array-like of shape `(n_freq, 3)` containing `[omegaA, omegaB, omegaC]` per row
- `components`: list of length `n_comp` of tensor components (strings)
- `values`: array-like of shape `(n_freq, n_comp)` containing numeric values
- `shape`: the literal tuple `("freq", "component")`

### Component naming

- Components are lowercase strings matching `^[xyz]{3}$` (e.g. `"xxx"`, `"xxy"`, ...).
- `components` is sorted (lexicographic) and `values[:, j]` corresponds to `components[j]`.

### Types

- `omega` and `values` must be convertible to `numpy.ndarray` with `dtype=float`.
- Missing components for a frequency may be represented as `NaN` in `values`.

## Fixture-backed examples (must pass)

- MADNESS: `tests/fixtures/load_calc/01_mra-d04_n2`
- Dalton: `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2`

## How to change this contract

If you want to add richer semantics (units, complex values, symmetry, etc.):
1) Update this document.
2) Update tests.
3) Update parsers accordingly.

