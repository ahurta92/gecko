# Contract: `calc.data["raman"]` (Raman spectra)

This contract defines the **standard, code-agnostic shape** for Raman spectra data stored on a loaded `Calculation`.

It is enforced by tests using fixtures in `tests/fixtures/load_calc/`.

## Scope

Applies to:
- MADNESS runs that provide Raman spectra (MADQC `*.calc_info.json` or legacy response JSON)
- Dalton Raman runs that print Raman tables and polarizability gradients

Non-goals (not part of this contract yet):
- Units and unit conversions
- Complex-valued tensors
- Guarantee that every run has Raman data (absence is allowed for non-Raman calculations)

## Representation

If Raman data is available, it must be stored as:

`calc.data["raman"] = { ... }` with required keys:
- `polarization_frequencies`: array-like of shape `(n_pol,)` with polarization frequencies (au)
- `vibrational_frequencies`: array-like of shape `(n_modes,)` with vibrational frequencies (cm⁻¹)
- `polarizability_derivatives`: list/array containing polarizability derivatives at each polarization frequency
- `polarizability_derivatives_by_mode`: list/array containing polarizability derivatives resolved by normal modes
- `raman_by_freq`: dict mapping polarization frequency → list of Raman mode rows

### Raman row schema

Each row in `raman_by_freq[freq]` must include:
- `mode`: integer mode index (1-based)
- `freq_cm1`: vibrational frequency in cm⁻¹
- `alpha2`: Raman $\alpha^2$
- `beta2`: Raman $\beta^2$
- `pol_int`: polarizability intensity
- `depol_int`: depolarization intensity
- `dep_ratio`: depolarization ratio

### Types

- `polarization_frequencies` and `vibrational_frequencies` must be convertible to `numpy.ndarray` with `dtype=float`.
- `polarizability_derivatives` entries must be convertible to `numpy.ndarray` with `dtype=float`.
- `raman_by_freq` keys must be convertible to `float`.

### Related requirement: Polarizability

For Raman fixtures, polarizability must also be available and conform to [docs/contracts/data_alpha.md](docs/contracts/data_alpha.md).

## Fixture-backed examples (must pass)

- MADNESS: `tests/fixtures/load_calc/03_mra-raman_h2o`
- Dalton: `tests/fixtures/load_calc/05_dalton_raman_h2o`

## How to change this contract

If you want to add richer semantics (units, complex values, symmetry, etc.):
1) Update this document.
2) Update tests.
3) Update parsers accordingly.
