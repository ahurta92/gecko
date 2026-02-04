# Contract: `gecko.core.load.load_calc`

This contract defines the **minimum stable behavior** of Gecko’s primary entrypoint for loading finished calculations.

It is enforced by tests (mirrored under `tests/`) using fixtures in `tests/fixtures/load_calc/`.

## Scope

Public surface:
- `gecko.core.load.load_calc(path)`
- `gecko.core.model.Calculation` (as returned by `load_calc`)

Non-goals (not part of this contract yet):
- Exact scientific values (energies/tensors)
- Full schema coverage of every legacy JSON variant
- Job submission / workflow execution

## Inputs

`load_calc(path)` accepts:
- A **directory** containing either a MADNESS run or a Dalton run (auto-detected).
- A **single file**:
  - Dalton: `*.out` (loads its parent directory, preferring that output)
  - MADNESS: `output.json` / `outputs.json`, `*.calc_info.json`, `*_mad_output.json`

## Outputs

`load_calc(...)` returns a `Calculation` with these invariants:

- `calc.code` is one of: `"madness"`, `"dalton"`
- `calc.root` is the resolved calculation directory (`Path`)
- `calc.artifacts` contains the discovered primary input artifacts (best-effort). It may also include structured groupings.
- `calc.data` contains parser outputs (best-effort), and may be empty if parsing fails
- `calc.molecule` is:
  - a `qcelemental.models.Molecule` when geometry can be determined, otherwise `None`
  - when `None`, `calc.meta["mol_source"] == "missing"`
- `calc.meta["molecule_id"]` is set when `calc.molecule` is set

Missing/optional data should be represented as:
- Absent keys in `calc.data`/`calc.meta`, or `None` values
- Non-fatal issues collected in `calc.meta["warnings"]` (list of strings)

## Fixture-backed examples (must pass)

Directory loads:
- `tests/fixtures/load_calc/00_mra-d06_bh2cl` → MADNESS legacy (`outputs.json`) + `.mol` fallback
- `tests/fixtures/load_calc/01_mra-d04_n2` → MADNESS legacy (`output.json`)
- `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2` → Dalton (single `quad.out`)
- `tests/fixtures/load_calc/03_mra-raman_h2o` → MADNESS MADQC (`*.calc_info.json`)
- `tests/fixtures/load_calc/05_dalton_raman_h2o` → Dalton (multi-`.out` directory)

File loads:
- `.../02_aug-cc-pVDZ_n2/quad.out` loads as Dalton with `calc.root == parent_dir`
- `.../00_mra-d06_bh2cl/outputs.json` loads as MADNESS with `calc.root == parent_dir`

## Dalton artifact discovery (multi-step)

For Dalton directories with multiple steps (e.g. optimize + raman), `load_calc(dir)` should discover **all available inputs and outputs**.

Minimum expectations (best-effort; presence depends on what files exist):
- `calc.artifacts["dalton_dal_files"]`: list of `Path` for `*.dal` in the directory
- `calc.artifacts["dalton_mol_files"]`: list of `Path` for `*.mol` in the directory
- `calc.artifacts["dalton_out_files"]`: list of `Path` for `*.out` in the directory (including `DALTON.OUT` if present)
- `calc.artifacts["dalton_pairs"]`: list of dicts of the form `{"dal": Path, "mol": Path, "out": Path}`

Pairing rule (initial contract):
- For each `X.dal` and `Y.mol`, if an output exists named `X_Y.out`, then that triple appears in `dalton_pairs`.

## How to change this contract

If you want to change any of the invariants above:
1) Update this contract doc.
2) Update the corresponding tests.
3) Then refactor implementation freely until tests pass.
