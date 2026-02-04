# Contract: `calc.meta["basis"]` (Basis / MRA label)

This contract defines how Gecko determines the **basis label** for a loaded `Calculation`.

The goal is a single, human-readable string that is stable enough to group calculations for analysis.

## Scope

Applies to calculations loaded via:
- `gecko.core.load.load_calc`

Non-goals (not part of this contract yet):
- Full method labeling (HF/DFT functional, etc.)
- Per-step basis labeling for multi-step workflows (may be added later)

## Representation

After `load_calc(...)`, the returned `calc` must have:
- `calc.meta["basis"]`: a non-empty string

### Dalton rule (basis-set calculations)

For Dalton calculations, the basis label must be read from the Dalton-format `.mol` file, using:

1) If the directory contains discovered `{dal}_{mol}.out` pairs (`calc.artifacts["dalton_pairs"]`):
   - Choose the `.mol` associated with the primary output (`calc.artifacts["out"]`) when possible.
   - Otherwise, choose the first paired `.mol`.
2) Parse the `.mol` contents:
   - Find a line equal to `BASIS` (case-insensitive), and take the next non-empty line as the basis name.

Examples: `aug-cc-pVDZ`, `d-aug-cc-pV6Z`.

### MADNESS rule (MRA calculations)

For MADNESS calculations, the basis label must be determined by this priority order:

0) **paired input `.in` file** (MADQC):
   - If a paired `{stem}.in` exists, it is the preferred source of `dconv` / `protocol`

1) **dconv** (preferred):
   - If a `dconv` value is available (e.g. `1e-6`, `1e-4`), label as `mra-dXX`
   - Example: `dconv=1e-4` → `mra-d04`

2) **protocol final threshold** (fallback):
   - If a protocol list exists, take the final element as the threshold and label as `mra-pXX`
   - Example: final protocol value `1e-8` → `mra-p08`

3) Otherwise:
   - Use `mra`

Notes:
- `dconv` may come from an input file (`input.json`) or from modern calc-info payloads
  (e.g. `convergence.converged_for_dconv`).

## Fixture-backed examples (must pass)

Expected basis labels:
- `tests/fixtures/load_calc/00_mra-d06_bh2cl` → `mra-d06`
- `tests/fixtures/load_calc/01_mra-d04_n2` → `mra-d04`
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2` → `aug-cc-pVDZ`
- `tests/fixtures/load_calc/03_mra-raman_h2o` → `mra-p07`
- `tests/fixtures/load_calc/05_dalton_raman_h2o` → `d-aug-cc-pV6Z`

## How to change this contract

If you want to encode additional semantics (method, grid, eprec, etc.):
1) Update this document.
2) Update tests.
3) Update the parser/loader logic accordingly.
