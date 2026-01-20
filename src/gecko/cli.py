from __future__ import annotations

import argparse
from pathlib import Path

from gecko.core.iterators import iter_calc_dirs
from gecko.recipes.shg_csv import build_beta_table


def _build_shg_command(args: argparse.Namespace) -> int:
    db_dir = Path(args.db)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    calc_dirs = list(iter_calc_dirs(db_dir))
    shg_df = build_beta_table(
        calc_dirs,
        shg_only=True,
        add_shg_omega=True,
        shg_start_at=args.start_at,
        shg_tol=args.tol,
        app_compat=True,
        verbose=args.verbose,
        fail_fast=args.fail_fast,
    )

    out_path = out_dir / "shg_ijk.csv"
    shg_df.to_csv(out_path, index=False)
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gecko")
    subparsers = parser.add_subparsers(dest="command", required=True)

    shg_parser = subparsers.add_parser("shg", help="SHG table utilities")
    shg_subparsers = shg_parser.add_subparsers(dest="shg_command", required=True)

    build_parser = shg_subparsers.add_parser("build", help="Build SHG ijk table")
    build_parser.add_argument("--db", required=True, help="Database root directory")
    build_parser.add_argument("--out", required=True, help="Output directory")
    build_parser.add_argument("--start-at", type=int, default=0)
    build_parser.add_argument("--tol", type=float, default=1e-12)
    build_parser.add_argument("--fail-fast", action="store_true", default=False)
    build_parser.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=True)
    build_parser.set_defaults(func=_build_shg_command)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
