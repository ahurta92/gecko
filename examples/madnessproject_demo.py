"""Demo: madnessproject API — mirrors daltonproject pattern for MADNESS.

Run from the gecko repo root:

    .venv/bin/python examples/madnessproject_demo.py

Dalton (for reference):

    molecule = dp.Molecule(input_file='water.xyz')
    basis    = dp.Basis(basis='cc-pVDZ')
    method   = dp.QCMethod('HF')
    props    = dp.Property(polarizabilities=True)
    settings = dp.ComputeSettings(work_dir='...', mpi_num_procs=90)
    result   = Dalton.compute(molecule, basis, method, props, compute_settings=settings)
    print('alpha (static):', result.polarizabilities.values[0])

MADNESS equivalent (this demo):

    molecule = mad.Molecule(input_file='water.xyz')
    calc     = mad.CalculationParameters(dipole=True, econv=1e-7, ...)
    response = mad.ResponseParameters()
    settings = mad.ComputeSettings(work_dir='...', MAD_NUM_THREADS=5, mpi_num_procs=8)
    result   = mad.Madness.compute(molecule, calc, response, compute_settings=settings)
    print('alpha:', result.alpha)
"""

from pathlib import Path

from gecko import madnessproject as mad

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"

# ============================================================
# 1. Object construction (one object per .in block)
# ============================================================

# Molecule -> molecule ... end block
molecule = mad.Molecule(
    atoms="O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692",
)
print("Molecule:          ", molecule)

# CalculationParameters -> dft ... end block
calc = mad.CalculationParameters(
    dipole=True,
    econv=1e-7,
    maxiter=20,
    protocol=[1e-4, 1e-6, 1e-7],
    gopt=True,
)
print("CalculationParams: ", calc)

# Protocol helper (can pass to CalculationParameters(protocol=...))
protocol = mad.Protocol([1e-4, 1e-6])
print("Protocol:          ", protocol)

# ResponseParameters -> response ... end block
response = mad.ResponseParameters(
    dipole=True,
    frequencies=[0.0, 0.05],
)
print("ResponseParams:    ", response)

# ComputeSettings (runtime config, not part of .in file)
settings = mad.ComputeSettings(
    work_dir="/tmp/gecko_madness_demo",
    MAD_NUM_THREADS=5,
    mpi_num_procs=8,
)
print("ComputeSettings:   ", settings)

# ============================================================
# 2. Generate .in file
# ============================================================

print("\n===== Generated MADNESS .in file =====\n")
text = mad.madness_input(calc, molecule, response)
print(text)

# ============================================================
# 3. Write to disk
# ============================================================

path = mad.write_input(calc, molecule, response,
                       output_path="/tmp/gecko_madness_demo/mad.in")
print(f"Written to: {path}\n")

# ============================================================
# 4. Dry-run compute (writes input, doesn't execute MADNESS)
# ============================================================

result = mad.Madness.compute(
    molecule, calc, response,
    compute_settings=settings,
    dry_run=True,
)
print("Dry-run result:", result)

# ============================================================
# 5. Read existing JSON output
# ============================================================

print("\n===== Reading existing MADNESS outputs =====\n")

outputs_json = FIXTURES / "beta_data" / "H2O" / "mra-d06" / "outputs.json"
if outputs_json.exists():
    r = mad.load_output(outputs_json)
    print(f"Source:  {outputs_json}")
    print(f"Energy:  {r.energy}")
    print(f"Alpha:   {list(r.alpha.keys()) if r.alpha else None}")
    print(f"Beta:    {list(r.beta.keys()) if r.beta else None}")
else:
    print(f"(fixture not found: {outputs_json})")

calc_info = FIXTURES / "load_calc" / "03_mra-raman_h2o" / "mad.h2o_gopt.calc_info.json"
if calc_info.exists():
    r2 = mad.load_output(calc_info)
    print(f"\nSource:  {calc_info}")
    print(f"Energy:  {r2.energy}")
    print(f"Molecule: {r2.molecule}")
else:
    print(f"(fixture not found: {calc_info})")
