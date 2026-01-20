import gecko

calc = gecko.load_calc("tests/fixtures/calc_nlo_beta/NLO/madness/n2")

beta = calc.data.get("beta", {})
print("beta keys:", beta.keys())
print("omega shape:", None if not beta else beta["omega"].shape)
print("values shape:", None if not beta else beta["values"].shape)
print("first components:", None if not beta else beta["components"][:5])

from pathlib import Path

from gecko.core.iterators import iter_calc_dirs
from gecko.recipes.shg_csv import build_beta_table

root = Path("tests/fixtures/calc_nlo_beta/NLO")
calc_dirs = list(iter_calc_dirs(root))

shg_df = build_beta_table(
    calc_dirs,
    shg_only=True,
    add_shg_omega=True,
    shg_start_at=0,
    shg_tol=1e-12,
    include_geometry=True,
    app_compat=True,
    verbose=True,
)

print(shg_df.head())

out = Path("data/csv_data")
out.mkdir(parents=True, exist_ok=True)

shg_df.to_csv(out / "shg_ijk.csv", index=False)
# print("beta['omega']:", beta["omega"])
# print("beta['values']", beta["values"])
