# Contract: `calc.meta["method"]`

This contract defines the minimal representation for the electronic structure method (level of theory).

## Representation

If available, `calc.meta["method"]` must be a non-empty string such as:
- `"HF"`
- `"MP2"`
- `"CCSD"`

For DFT calculations, `method` should be the **XC functional name** (e.g., `"PBE"`, `"B3LYP"`) rather than a generic `"DFT"`.

## Fixture-backed examples (must pass)

- MADNESS: `tests/fixtures/load_calc/03_mra-raman_h2o` → `"HF"`
- Dalton: `tests/fixtures/load_calc/02_aug-cc-pVDZ_n2` → `"HF"`
- Dalton: `tests/fixtures/load_calc/05_dalton_raman_h2o` → `"HF"`

## How to change this contract

1) Update this document.
2) Update tests.
3) Update method inference accordingly.
