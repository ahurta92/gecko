"""Trame dashboard for Raman/polarizability convergence analysis.

This app is intentionally lightweight and built on existing Gecko patterns:
- data loading through `iter_calc_dirs` + `load_calc` + `TableBuilder`
- Vuetify layout via trame's `SinglePageWithDrawerLayout`
- matplotlib-rendered plots served as in-memory data URLs

Run:
    python -m gecko.viz.apps.raman_dashboard --db-dir /path/to/project_data
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import errno
import io
import socket
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from trame.app import get_server
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import vuetify

from gecko.core.iterators import iter_calc_dirs
from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder
from gecko.tables.extractors import make_envelope


_DEFAULT_DB_DIR = Path("/gpfs/scratch/ahurtado/project_data/data/raman_paper/data")
_DB_DIR: Path | None = None

_COMPONENTS = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"]
_SMALL_THRESH = 1e-10


@dataclass
class DashboardData:
    root: Path
    n_calc_dirs: int
    n_loaded_calcs: int
    load_errors: list[str]
    alpha_df: pd.DataFrame
    raman_df: pd.DataFrame
    dalpha_df: pd.DataFrame
    alpha_fit_df: pd.DataFrame
    dalpha_fit_df: pd.DataFrame
    molecules: list[str]
    bases: list[str]


def _repo_root() -> Path:
    # src/gecko/viz/apps/raman_dashboard.py -> apps -> viz -> gecko -> src -> repo
    return Path(__file__).resolve().parents[4]


def _resolve_db_dir() -> Path:
    if _DB_DIR is not None:
        return _DB_DIR
    if _DEFAULT_DB_DIR.exists():
        return _DEFAULT_DB_DIR
    return _repo_root() / "data"


def _safe_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _clean_string(value: Any) -> str | None:
    if value is None:
        return None
    out = str(value).strip()
    return out or None


def _is_mra_basis(basis: str) -> bool:
    return "mra" in str(basis).lower()


def _infer_molecule_from_root(root: str | Path) -> str:
    p = Path(root)
    parts = p.parts
    if "data" in parts:
        idx = parts.index("data")
        if idx + 1 < len(parts):
            return str(parts[idx + 1])
    if p.parent != p:
        return str(p.parent.name)
    return str(p.name)


def _attach_molecule_column(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        out = df.copy()
        out["molecule"] = []
        return out

    def _pick(row: pd.Series) -> str:
        for key in ("mol_id", "molecule_id", "label"):
            val = _clean_string(row.get(key))
            if val:
                return val
        root_val = row.get("root")
        if root_val is not None:
            return _infer_molecule_from_root(root_val)
        return "unknown"

    out = df.copy()
    out["molecule"] = out.apply(_pick, axis=1)
    if "basis" in out.columns:
        out["basis"] = out["basis"].astype(str).str.strip()
    return out


def _iter_derivative_tensors(
    pol_freqs: np.ndarray,
    derivs_by_mode: Any,
) -> list[tuple[float, np.ndarray]]:
    rows: list[tuple[float, np.ndarray]] = []

    if isinstance(derivs_by_mode, dict):
        for key, value in derivs_by_mode.items():
            om = _safe_float(key)
            if om is None:
                continue
            rows.append((om, np.asarray(value, dtype=float)))
        rows.sort(key=lambda item: item[0])
        return rows

    if isinstance(derivs_by_mode, list):
        for idx, value in enumerate(derivs_by_mode):
            omega = float(idx)
            if idx < len(pol_freqs):
                parsed = _safe_float(pol_freqs[idx])
                if parsed is not None:
                    omega = parsed
            rows.append((omega, np.asarray(value, dtype=float)))
    return rows


def _coerce_derivative_matrix(arr: np.ndarray) -> np.ndarray | None:
    if arr.ndim == 3 and arr.shape[0] == 3 and arr.shape[1] == 3:
        # (3,3,n_modes) -> (9,n_modes)
        return arr.reshape(9, arr.shape[2])

    if arr.ndim == 2:
        if arr.shape[0] == 9:
            return arr
        if arr.shape[1] == 9:
            return arr.T

    return None


def _mode_frequency_map(raman_by_freq: Any) -> dict[float, dict[int, float]]:
    out: dict[float, dict[int, float]] = {}
    if not isinstance(raman_by_freq, dict):
        return out

    for omega_key, rows in raman_by_freq.items():
        omega = _safe_float(omega_key)
        if omega is None or not isinstance(rows, list):
            continue
        mode_map: dict[int, float] = {}
        for row in rows:
            if not isinstance(row, dict):
                continue
            mode_raw = row.get("mode")
            freq_raw = row.get("freq_cm1")
            try:
                mode = int(mode_raw)
            except (TypeError, ValueError):
                continue
            freq = _safe_float(freq_raw)
            if freq is None:
                continue
            mode_map[mode] = freq
        out[omega] = mode_map
    return out


def _extract_polarizability_derivative_df(calcs: Iterable[Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for calc in calcs:
        raman = calc.data.get("raman") or {}
        derivs_by_mode = raman.get("polarizability_derivatives_by_mode")
        if derivs_by_mode is None:
            continue

        pol_freqs = np.asarray(raman.get("polarization_frequencies", []), dtype=float).reshape(-1)
        mode_freq_map = _mode_frequency_map(raman.get("raman_by_freq") or {})
        env = make_envelope(calc)

        for omega, raw_arr in _iter_derivative_tensors(pol_freqs, derivs_by_mode):
            matrix = _coerce_derivative_matrix(raw_arr)
            if matrix is None:
                continue

            n_modes = int(matrix.shape[1])
            for mode_idx in range(n_modes):
                mode = mode_idx + 1
                freq_cm1 = mode_freq_map.get(float(omega), {}).get(mode)
                for comp_idx, ij in enumerate(_COMPONENTS):
                    value = _safe_float(matrix[comp_idx, mode_idx])
                    if value is None:
                        continue
                    rows.append(
                        {
                            **env,
                            "omega": float(omega),
                            "mode": int(mode),
                            "freq_cm1": float(freq_cm1) if freq_cm1 is not None else np.nan,
                            "ij": str(ij),
                            "value": float(value),
                        }
                    )

    if not rows:
        return pd.DataFrame(
            columns=[
                "calc_id",
                "geom_id",
                "mol_id",
                "molecule_id",
                "label",
                "code",
                "root",
                "basis",
                "method",
                "omega",
                "mode",
                "freq_cm1",
                "ij",
                "value",
            ]
        )

    return pd.DataFrame(rows)


def _power_series_fit(omega: np.ndarray, values: np.ndarray) -> dict[str, float] | None:
    mask = np.isfinite(omega) & np.isfinite(values)
    if int(mask.sum()) < 3:
        return None

    om = omega[mask]
    val = values[mask]

    x2 = om**2
    x4 = om**4
    X = np.column_stack([np.ones_like(om), x2, x4])

    try:
        coeff, _, rank, _ = np.linalg.lstsq(X, val, rcond=None)
    except np.linalg.LinAlgError:
        return None

    if int(rank) < 3:
        return None

    c0, c2, c4 = (float(coeff[0]), float(coeff[1]), float(coeff[2]))
    fit = X @ coeff

    ss_res = float(np.sum((val - fit) ** 2))
    ss_tot = float(np.sum((val - np.mean(val)) ** 2))
    r2 = np.nan if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean((val - fit) ** 2)))

    A = np.nan
    B = np.nan
    if abs(c0) > _SMALL_THRESH:
        A = c2 / c0
        B = c4 / c0

    return {
        "f0": c0,
        "A": float(A),
        "B": float(B),
        "coef2": c2,
        "coef4": c4,
        "rmse": rmse,
        "r2": float(r2),
        "n_points": float(mask.sum()),
    }


def _fit_df(
    df: pd.DataFrame,
    *,
    group_cols: list[str],
    omega_col: str,
    value_col: str,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=[*group_cols, "f0", "A", "B", "coef2", "coef4", "rmse", "r2", "n_points"])

    rows: list[dict[str, Any]] = []
    grouped = df.groupby(group_cols, dropna=False, sort=True)
    for key, chunk in grouped:
        omega = np.asarray(chunk[omega_col], dtype=float)
        value = np.asarray(chunk[value_col], dtype=float)
        fit = _power_series_fit(omega, value)
        if fit is None:
            continue

        row: dict[str, Any] = {}
        if isinstance(key, tuple):
            for idx, col in enumerate(group_cols):
                row[col] = key[idx]
        else:
            row[group_cols[0]] = key
        row.update(fit)
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=[*group_cols, "f0", "A", "B", "coef2", "coef4", "rmse", "r2", "n_points"])

    return pd.DataFrame(rows)


def _collect_calcs(root: Path) -> tuple[list[Any], list[str], int]:
    calc_dirs = list(iter_calc_dirs(root))
    calcs: list[Any] = []
    errors: list[str] = []

    for calc_dir in calc_dirs:
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                calc = load_calc(calc_dir)
            calcs.append(calc)
        except Exception as exc:  # pragma: no cover - best-effort cataloging
            errors.append(f"{calc_dir}: {type(exc).__name__}: {exc}")

    return calcs, errors, len(calc_dirs)


@lru_cache(maxsize=1)
def _data() -> DashboardData:
    root = _resolve_db_dir().expanduser().resolve()
    if not root.exists():
        empty = pd.DataFrame()
        return DashboardData(
            root=root,
            n_calc_dirs=0,
            n_loaded_calcs=0,
            load_errors=[f"Dataset root does not exist: {root}"],
            alpha_df=empty,
            raman_df=empty,
            dalpha_df=empty,
            alpha_fit_df=empty,
            dalpha_fit_df=empty,
            molecules=[],
            bases=[],
        )

    calcs, load_errors, n_calc_dirs = _collect_calcs(root)
    tb = TableBuilder(calcs)

    alpha_df = _attach_molecule_column(tb.build_alpha())
    raman_df = _attach_molecule_column(tb.build_raman())
    dalpha_df = _attach_molecule_column(_extract_polarizability_derivative_df(calcs))

    alpha_fit_df = _fit_df(
        alpha_df,
        group_cols=["molecule", "basis", "ij"],
        omega_col="omega",
        value_col="value",
    )
    dalpha_fit_df = _fit_df(
        dalpha_df,
        group_cols=["molecule", "basis", "mode", "ij"],
        omega_col="omega",
        value_col="value",
    )

    molecules = sorted(
        set(alpha_df.get("molecule", []))
        | set(raman_df.get("molecule", []))
        | set(dalpha_df.get("molecule", []))
    )
    bases = sorted(
        set(alpha_df.get("basis", []))
        | set(raman_df.get("basis", []))
        | set(dalpha_df.get("basis", []))
    )
    bases = [b for b in bases if _clean_string(b)]

    return DashboardData(
        root=root,
        n_calc_dirs=n_calc_dirs,
        n_loaded_calcs=len(calcs),
        load_errors=load_errors,
        alpha_df=alpha_df,
        raman_df=raman_df,
        dalpha_df=dalpha_df,
        alpha_fit_df=alpha_fit_df,
        dalpha_fit_df=dalpha_fit_df,
        molecules=[str(m) for m in molecules if _clean_string(m)],
        bases=[str(b) for b in bases],
    )


def _error_long(
    df: pd.DataFrame,
    *,
    key_cols: list[str],
    value_col: str,
    ref_basis: str,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    if "basis" not in df.columns or value_col not in df.columns:
        return pd.DataFrame()

    work = df.copy()
    work = work[work["basis"].astype(str).str.len() > 0]
    if work.empty:
        return pd.DataFrame()

    pivot = work.pivot_table(index=key_cols, columns="basis", values=value_col, aggfunc="first")
    if ref_basis not in pivot.columns:
        return pd.DataFrame()

    pivot = pivot.reset_index()
    basis_cols = [col for col in pivot.columns if col not in key_cols]

    rows: list[dict[str, Any]] = []
    for _, row in pivot.iterrows():
        ref_value = _safe_float(row.get(ref_basis))
        if ref_value is None:
            continue
        for basis in basis_cols:
            if basis == ref_basis:
                continue
            value = _safe_float(row.get(basis))
            if value is None:
                continue

            rel = np.nan
            if abs(ref_value) > _SMALL_THRESH:
                rel = (value - ref_value) / abs(ref_value)

            out = {col: row[col] for col in key_cols}
            out.update(
                {
                    "basis": str(basis),
                    "ref_basis": str(ref_basis),
                    "value": float(value),
                    "ref_value": float(ref_value),
                    "rel": float(rel) if np.isfinite(rel) else np.nan,
                    "abs_rel": float(abs(rel)) if np.isfinite(rel) else np.nan,
                }
            )
            rows.append(out)

    return pd.DataFrame(rows)


def _coverage_frame(data: DashboardData) -> pd.DataFrame:
    molecules = list(data.molecules)
    bases = list(data.bases)
    if not molecules or not bases:
        return pd.DataFrame(
            columns=[
                "molecule",
                "basis",
                "alpha_points",
                "alpha_omega_count",
                "dalpha_points",
                "dalpha_omega_count",
                "dalpha_mode_count",
                "raman_points",
                "raman_omega_count",
                "raman_mode_count",
                "has_alpha",
                "has_dalpha",
                "has_raman",
                "complete_row",
            ]
        )

    grid = pd.MultiIndex.from_product([molecules, bases], names=["molecule", "basis"]).to_frame(index=False)

    def _count(df: pd.DataFrame, value_col: str, omega_col: str, mode_col: str | None = None) -> pd.DataFrame:
        if df.empty:
            cols = ["molecule", "basis", f"{value_col}_points", f"{value_col}_omega_count"]
            if mode_col is not None:
                cols.append(f"{value_col}_mode_count")
            return pd.DataFrame(columns=cols)

        agg_map: dict[str, tuple[str, str]] = {
            f"{value_col}_points": (value_col, "count"),
            f"{value_col}_omega_count": (omega_col, pd.Series.nunique),
        }
        if mode_col is not None:
            agg_map[f"{value_col}_mode_count"] = (mode_col, pd.Series.nunique)

        return df.groupby(["molecule", "basis"], as_index=False).agg(**agg_map)

    alpha_cov = _count(data.alpha_df, value_col="value", omega_col="omega")
    alpha_cov = alpha_cov.rename(columns={"value_points": "alpha_points", "value_omega_count": "alpha_omega_count"})

    dalpha_cov = _count(data.dalpha_df, value_col="value", omega_col="omega", mode_col="mode")
    dalpha_cov = dalpha_cov.rename(
        columns={
            "value_points": "dalpha_points",
            "value_omega_count": "dalpha_omega_count",
            "value_mode_count": "dalpha_mode_count",
        }
    )

    raman_cov = _count(data.raman_df, value_col="pol_int", omega_col="omega_pol", mode_col="mode")
    raman_cov = raman_cov.rename(
        columns={
            "pol_int_points": "raman_points",
            "pol_int_omega_count": "raman_omega_count",
            "pol_int_mode_count": "raman_mode_count",
        }
    )

    out = grid.merge(alpha_cov, on=["molecule", "basis"], how="left")
    out = out.merge(dalpha_cov, on=["molecule", "basis"], how="left")
    out = out.merge(raman_cov, on=["molecule", "basis"], how="left")

    for col in (
        "alpha_points",
        "alpha_omega_count",
        "dalpha_points",
        "dalpha_omega_count",
        "dalpha_mode_count",
        "raman_points",
        "raman_omega_count",
        "raman_mode_count",
    ):
        out[col] = out[col].fillna(0).astype(int)

    out["has_alpha"] = out["alpha_points"] > 0
    out["has_dalpha"] = out["dalpha_points"] > 0
    out["has_raman"] = out["raman_points"] > 0
    out["complete_row"] = out["has_alpha"] & out["has_dalpha"] & out["has_raman"]

    return out


def _summary_frame(data: DashboardData, ref_basis: str) -> pd.DataFrame:
    coverage = _coverage_frame(data)
    if coverage.empty:
        return coverage

    alpha_err = _error_long(
        data.alpha_df,
        key_cols=["molecule", "omega", "ij"],
        value_col="value",
        ref_basis=ref_basis,
    )
    dalpha_err = _error_long(
        data.dalpha_df,
        key_cols=["molecule", "omega", "mode", "ij"],
        value_col="value",
        ref_basis=ref_basis,
    )
    raman_err = _error_long(
        data.raman_df,
        key_cols=["molecule", "omega_pol", "mode"],
        value_col="pol_int",
        ref_basis=ref_basis,
    )

    def _agg(err_df: pd.DataFrame, prefix: str) -> pd.DataFrame:
        if err_df.empty:
            return pd.DataFrame(columns=["molecule", "basis", f"{prefix}_mae", f"{prefix}_max"])
        return (
            err_df.groupby(["molecule", "basis"], as_index=False)
            .agg(
                **{
                    f"{prefix}_mae": ("abs_rel", "mean"),
                    f"{prefix}_max": ("abs_rel", "max"),
                }
            )
        )

    summary = coverage.copy()
    summary = summary.merge(_agg(alpha_err, "alpha"), on=["molecule", "basis"], how="left")
    summary = summary.merge(_agg(dalpha_err, "dalpha"), on=["molecule", "basis"], how="left")
    summary = summary.merge(_agg(raman_err, "pol_int"), on=["molecule", "basis"], how="left")

    for col in [
        "alpha_mae",
        "alpha_max",
        "dalpha_mae",
        "dalpha_max",
        "pol_int_mae",
        "pol_int_max",
    ]:
        summary[col] = summary[col].astype(float)
        summary.loc[summary["basis"] == ref_basis, col] = 0.0

    summary["is_mra"] = summary["basis"].map(_is_mra_basis)
    return summary.sort_values(["molecule", "is_mra", "basis"], ascending=[True, False, True]).reset_index(drop=True)


def _data_url_from_figure(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _global_error_plot_data_url(
    summary: pd.DataFrame,
    *,
    ref_basis: str,
    metric_key: str,
    yscale: str,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    if summary.empty:
        return ""

    metric_label = {
        "alpha_mae": "Alpha mean |rel error| [%]",
        "dalpha_mae": "dAlpha mean |rel error| [%]",
        "pol_int_mae": "Raman pol_int mean |rel error| [%]",
    }.get(metric_key, metric_key)

    nonref = summary[summary["basis"] != ref_basis].copy()
    if nonref.empty:
        return ""

    molecules = sorted(nonref["molecule"].astype(str).unique())
    bases = sorted(nonref["basis"].astype(str).unique())

    fig, ax = plt.subplots(figsize=(9.4, 4.5), dpi=160)

    xs = np.arange(len(molecules), dtype=float)
    for basis in bases:
        y_values: list[float] = []
        for mol in molecules:
            row = nonref[(nonref["molecule"] == mol) & (nonref["basis"] == basis)]
            if row.empty:
                y_values.append(np.nan)
                continue
            val = _safe_float(row.iloc[0].get(metric_key))
            y_values.append(100.0 * val if val is not None else np.nan)

        arr = np.asarray(y_values, dtype=float)
        if not np.isfinite(arr).any():
            continue
        ax.plot(xs, arr, marker="o", linewidth=1.6, label=str(basis))

    ax.set_xticks(xs)
    ax.set_xticklabels(molecules, rotation=25, ha="right")
    ax.set_ylabel(metric_label)
    ax.set_title(f"Global basis-set error vs molecule (reference: {ref_basis})")
    ax.grid(True, axis="y", alpha=0.25)

    if yscale == "symlog":
        ax.set_yscale("symlog", linthresh=1e-3)

    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = _data_url_from_figure(fig)
    plt.close(fig)
    return out


def _evaluate_power_fit(omega: np.ndarray, *, f0: float, A: float, B: float) -> np.ndarray:
    return f0 * (1.0 + A * (omega**2) + B * (omega**4))


def _basis_color_map(bases: list[str]) -> dict[str, Any]:
    try:
        import matplotlib.pyplot as plt

        cmap = plt.get_cmap("tab10")
        return {basis: cmap(i % 10) for i, basis in enumerate(bases)}
    except Exception:
        return {basis: None for basis in bases}


def _molecule_plot_data_url(
    data: DashboardData,
    *,
    molecule: str,
    ref_basis: str,
    alpha_component: str,
    dalpha_mode: int,
    dalpha_component: str,
    yscale: str,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    alpha_mol = data.alpha_df[(data.alpha_df["molecule"] == molecule) & (data.alpha_df["ij"] == alpha_component)]
    dalpha_mol = data.dalpha_df[
        (data.dalpha_df["molecule"] == molecule)
        & (data.dalpha_df["ij"] == dalpha_component)
        & (data.dalpha_df["mode"].astype(int) == int(dalpha_mode))
    ]
    raman_mol = data.raman_df[
        (data.raman_df["molecule"] == molecule)
        & (data.raman_df["mode"].astype(int) == int(dalpha_mode))
    ]

    if alpha_mol.empty and dalpha_mol.empty and raman_mol.empty:
        return ""

    basis_order = sorted(
        set(alpha_mol.get("basis", []))
        | set(dalpha_mol.get("basis", []))
        | set(raman_mol.get("basis", []))
    )
    if ref_basis in basis_order:
        basis_order.remove(ref_basis)
        basis_order = [ref_basis] + basis_order

    color_map = _basis_color_map(basis_order)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9.6, 9.6), dpi=160)

    # Panel 1: alpha(omega) + fit
    for basis in basis_order:
        chunk = alpha_mol[alpha_mol["basis"] == basis].copy()
        if chunk.empty:
            continue
        chunk = chunk.sort_values("omega")
        om = np.asarray(chunk["omega"], dtype=float)
        val = np.asarray(chunk["value"], dtype=float)
        ax1.plot(om, val, marker="o", linewidth=1.5, label=str(basis), color=color_map.get(basis))

        fit_row = data.alpha_fit_df[
            (data.alpha_fit_df["molecule"] == molecule)
            & (data.alpha_fit_df["basis"] == basis)
            & (data.alpha_fit_df["ij"] == alpha_component)
        ]
        if not fit_row.empty:
            f0 = _safe_float(fit_row.iloc[0].get("f0"))
            A = _safe_float(fit_row.iloc[0].get("A"))
            B = _safe_float(fit_row.iloc[0].get("B"))
            if f0 is not None and A is not None and B is not None:
                xfit = np.linspace(float(np.nanmin(om)), float(np.nanmax(om)), 200)
                yfit = _evaluate_power_fit(xfit, f0=f0, A=A, B=B)
                ax1.plot(xfit, yfit, linestyle="--", linewidth=1.0, color=color_map.get(basis))

    ax1.set_title(f"Molecule detail: {molecule}  (reference basis: {ref_basis})")
    ax1.set_ylabel(f"alpha[{alpha_component}]")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best", fontsize=8)

    # Panel 2: d(alpha)/dQ(omega) + fit
    for basis in basis_order:
        chunk = dalpha_mol[dalpha_mol["basis"] == basis].copy()
        if chunk.empty:
            continue
        chunk = chunk.sort_values("omega")
        om = np.asarray(chunk["omega"], dtype=float)
        val = np.asarray(chunk["value"], dtype=float)
        ax2.plot(om, val, marker="o", linewidth=1.5, label=str(basis), color=color_map.get(basis))

        fit_row = data.dalpha_fit_df[
            (data.dalpha_fit_df["molecule"] == molecule)
            & (data.dalpha_fit_df["basis"] == basis)
            & (data.dalpha_fit_df["mode"].astype(int) == int(dalpha_mode))
            & (data.dalpha_fit_df["ij"] == dalpha_component)
        ]
        if not fit_row.empty:
            f0 = _safe_float(fit_row.iloc[0].get("f0"))
            A = _safe_float(fit_row.iloc[0].get("A"))
            B = _safe_float(fit_row.iloc[0].get("B"))
            if f0 is not None and A is not None and B is not None:
                xfit = np.linspace(float(np.nanmin(om)), float(np.nanmax(om)), 200)
                yfit = _evaluate_power_fit(xfit, f0=f0, A=A, B=B)
                ax2.plot(xfit, yfit, linestyle="--", linewidth=1.0, color=color_map.get(basis))

    ax2.set_ylabel(f"dalpha[{dalpha_component}]  (mode {int(dalpha_mode)})")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="best", fontsize=8)

    # Panel 3: Raman intensity proxy from derivatives
    for basis in basis_order:
        chunk = raman_mol[raman_mol["basis"] == basis].copy()
        if chunk.empty:
            continue
        chunk = chunk.sort_values("omega_pol")
        om = np.asarray(chunk["omega_pol"], dtype=float)
        val = np.asarray(chunk["pol_int"], dtype=float)
        ax3.plot(om, val, marker="o", linewidth=1.5, label=str(basis), color=color_map.get(basis))

    ax3.set_xlabel("omega")
    ax3.set_ylabel(f"pol_int  (mode {int(dalpha_mode)})")
    ax3.grid(True, alpha=0.25)
    ax3.legend(loc="best", fontsize=8)

    if yscale == "symlog":
        ax1.set_yscale("symlog", linthresh=1e-6)
        ax2.set_yscale("symlog", linthresh=1e-6)
        ax3.set_yscale("symlog", linthresh=1e-6)

    fig.tight_layout()
    out = _data_url_from_figure(fig)
    plt.close(fig)
    return out


def _default_reference_basis(data: DashboardData) -> str:
    if not data.bases:
        return ""
    for basis in data.bases:
        if _is_mra_basis(basis):
            return basis
    return data.bases[0]


def _mode_items_for_molecule(data: DashboardData, molecule: str) -> list[int]:
    if not molecule or data.dalpha_df.empty:
        return []
    chunk = data.dalpha_df[data.dalpha_df["molecule"] == molecule]
    if chunk.empty:
        return []
    return sorted(set(int(m) for m in chunk["mode"].dropna().astype(int).tolist()))


def _summary_table_items(summary: pd.DataFrame) -> list[dict[str, Any]]:
    if summary.empty:
        return []

    items: list[dict[str, Any]] = []
    for _, row in summary.iterrows():
        def _pct(col: str) -> float | None:
            v = _safe_float(row.get(col))
            return None if v is None else 100.0 * v

        items.append(
            {
                "molecule": str(row.get("molecule")),
                "basis": str(row.get("basis")),
                "complete": "yes" if bool(row.get("complete_row")) else "no",
                "alpha_omega": int(row.get("alpha_omega_count", 0)),
                "dalpha_modes": int(row.get("dalpha_mode_count", 0)),
                "raman_modes": int(row.get("raman_mode_count", 0)),
                "alpha_mae_pct": _pct("alpha_mae"),
                "dalpha_mae_pct": _pct("dalpha_mae"),
                "pol_int_mae_pct": _pct("pol_int_mae"),
            }
        )

    return items


def _refresh_dashboard() -> None:
    data = _data()

    state.dataset_root = str(data.root)
    state.molecule_items = list(data.molecules)
    state.basis_items = list(data.bases)

    try:
        state.status = "loading"
        state.status_color = "orange"
        state.last_error = ""

        if not data.molecules or not data.bases:
            state.status = "no data"
            state.status_color = "red"
            state.last_error = "No alpha/raman data found in dataset root."
            state.mode_items = []
            state.summary_items = []
            state.global_plot_src = ""
            state.molecule_plot_src = ""
            return

        if not state.selected_molecule or state.selected_molecule not in data.molecules:
            state.selected_molecule = data.molecules[0]

        if not state.reference_basis or state.reference_basis not in data.bases:
            state.reference_basis = _default_reference_basis(data)

        mode_items = _mode_items_for_molecule(data, str(state.selected_molecule))
        state.mode_items = mode_items
        if mode_items:
            if int(state.selected_mode) not in mode_items:
                state.selected_mode = int(mode_items[0])
        else:
            state.selected_mode = 1

        summary = _summary_frame(data, str(state.reference_basis))
        state.summary_items = _summary_table_items(summary)

        # KPIs
        state.kpi_total_molecules = int(len(data.molecules))

        complete_molecules = 0
        if not summary.empty:
            for mol in sorted(summary["molecule"].astype(str).unique()):
                chunk = summary[summary["molecule"] == mol]
                has_mra = bool((chunk["complete_row"] & chunk["is_mra"]).any())
                has_bs = bool((chunk["complete_row"] & (~chunk["is_mra"])).any())
                if has_mra and has_bs:
                    complete_molecules += 1
        state.kpi_complete_molecules = int(complete_molecules)

        nonref = summary[summary["basis"] != str(state.reference_basis)] if not summary.empty else summary
        alpha_mean = float(np.nanmean(nonref["alpha_mae"].to_numpy(dtype=float))) if not nonref.empty else np.nan
        pol_mean = float(np.nanmean(nonref["pol_int_mae"].to_numpy(dtype=float))) if not nonref.empty else np.nan
        state.kpi_alpha_mae_pct = None if not np.isfinite(alpha_mean) else 100.0 * alpha_mean
        state.kpi_pol_mae_pct = None if not np.isfinite(pol_mean) else 100.0 * pol_mean

        state.global_plot_src = _global_error_plot_data_url(
            summary,
            ref_basis=str(state.reference_basis),
            metric_key=str(state.global_metric),
            yscale=str(state.global_scale),
        )
        state.molecule_plot_src = _molecule_plot_data_url(
            data,
            molecule=str(state.selected_molecule),
            ref_basis=str(state.reference_basis),
            alpha_component=str(state.alpha_component),
            dalpha_mode=int(state.selected_mode),
            dalpha_component=str(state.dalpha_component),
            yscale=str(state.fit_scale),
        )

        state.status = "ready"
        state.status_color = "green"
    except Exception as exc:  # pragma: no cover - defensive UI path
        state.status = "error"
        state.status_color = "red"
        state.last_error = f"{type(exc).__name__}: {exc}"


# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

state.setdefault("status", "starting")
state.setdefault("status_color", "orange")
state.setdefault("last_error", "")

state.setdefault("dataset_root", str(_resolve_db_dir()))
state.setdefault("molecule_items", [])
state.setdefault("basis_items", [])
state.setdefault("mode_items", [])

state.setdefault("selected_molecule", "")
state.setdefault("reference_basis", "")
state.setdefault("alpha_component", "xx")
state.setdefault("dalpha_component", "xx")
state.setdefault("selected_mode", 1)
state.setdefault("global_metric", "alpha_mae")
state.setdefault("global_scale", "linear")
state.setdefault("fit_scale", "linear")

state.setdefault("summary_items", [])
state.setdefault("global_plot_src", "")
state.setdefault("molecule_plot_src", "")

state.setdefault("kpi_total_molecules", 0)
state.setdefault("kpi_complete_molecules", 0)
state.setdefault("kpi_alpha_mae_pct", None)
state.setdefault("kpi_pol_mae_pct", None)

state.setdefault(
    "summary_headers",
    [
        {"text": "Molecule", "value": "molecule"},
        {"text": "Basis", "value": "basis"},
        {"text": "Complete", "value": "complete"},
        {"text": "alpha omegas", "value": "alpha_omega"},
        {"text": "dalpha modes", "value": "dalpha_modes"},
        {"text": "Raman modes", "value": "raman_modes"},
        {"text": "alpha |rel| %", "value": "alpha_mae_pct"},
        {"text": "dalpha |rel| %", "value": "dalpha_mae_pct"},
        {"text": "pol_int |rel| %", "value": "pol_int_mae_pct"},
    ],
)


@state.change("selected_molecule")
def _on_molecule_change(**_):
    _refresh_dashboard()


@state.change(
    "reference_basis",
    "alpha_component",
    "dalpha_component",
    "selected_mode",
    "global_metric",
    "global_scale",
    "fit_scale",
)
def _on_controls_change(**_):
    _refresh_dashboard()


def _fmt_pct_expr(key: str) -> str:
    return f"{{{{ {key} === null || {key} === undefined ? 'n/a' : {key}.toFixed(2) + '%' }}}}"


with SinglePageWithDrawerLayout(server) as layout:
    layout.title.set_text("Raman Convergence Dashboard")

    with layout.toolbar:
        vuetify.VSpacer()
        vuetify.VChip(
            "{{ status }}",
            small=True,
            outlined=True,
            color=("status_color", "orange"),
            classes="mr-2",
        )

    with layout.drawer as drawer:
        drawer.width = 380
        vuetify.VContainer(fluid=True, classes="pa-2")

        vuetify.VAlert(
            "Dataset root: {{ dataset_root }}",
            type="info",
            dense=True,
            outlined=True,
            classes="mb-2",
        )

        vuetify.VAlert(
            "Assumptions: Raman intensity error uses pol_int; derivative fits use component tensors from polarizability_derivatives_by_mode when available.",
            type="warning",
            dense=True,
            outlined=True,
            classes="mb-2",
        )

        vuetify.VAlert(
            "{{ last_error }}",
            type="error",
            dense=True,
            outlined=True,
            v_if="last_error",
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Molecule",
            v_model=("selected_molecule", state.selected_molecule),
            items=("molecule_items", []),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Reference basis",
            v_model=("reference_basis", state.reference_basis),
            items=("basis_items", []),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Global error metric",
            v_model=("global_metric", state.global_metric),
            items=(
                "metric_items",
                [
                    {"text": "alpha mean |rel|", "value": "alpha_mae"},
                    {"text": "dalpha mean |rel|", "value": "dalpha_mae"},
                    {"text": "pol_int mean |rel|", "value": "pol_int_mae"},
                ],
            ),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Alpha component",
            v_model=("alpha_component", state.alpha_component),
            items=("alpha_component_items", _COMPONENTS),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Derivative mode",
            v_model=("selected_mode", state.selected_mode),
            items=("mode_items", []),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Derivative component",
            v_model=("dalpha_component", state.dalpha_component),
            items=("dalpha_component_items", _COMPONENTS),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Global plot y-scale",
            v_model=("global_scale", state.global_scale),
            items=("global_scale_items", ["linear", "symlog"]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Molecule plot y-scale",
            v_model=("fit_scale", state.fit_scale),
            items=("fit_scale_items", ["linear", "symlog"]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

    with layout.content:
        with vuetify.VContainer(fluid=True, classes="pa-2"):
            with vuetify.VRow(dense=True, classes="mb-2"):
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Molecules", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_total_molecules }}", classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Complete Molecules", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_complete_molecules }}", classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Global alpha Error", classes="text-subtitle-2")
                        vuetify.VCardText(_fmt_pct_expr("kpi_alpha_mae_pct"), classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Global pol_int Error", classes="text-subtitle-2")
                        vuetify.VCardText(_fmt_pct_expr("kpi_pol_mae_pct"), classes="text-h5")

            with vuetify.VRow(dense=True, classes="mb-2"):
                with vuetify.VCol(cols=6):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Global Basis-Error Trend")
                        vuetify.VImg(
                            src=("global_plot_src", ""),
                            contain=True,
                            v_if="global_plot_src",
                            style="width: 100%;",
                        )
                with vuetify.VCol(cols=6):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Molecule Trends + Power-Series Fits")
                        vuetify.VImg(
                            src=("molecule_plot_src", ""),
                            contain=True,
                            v_if="molecule_plot_src",
                            style="width: 100%;",
                        )

            with vuetify.VRow(dense=True):
                with vuetify.VCol(cols=12):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Breakdown Table (coverage + basis errors)")
                        vuetify.VDataTable(
                            dense=True,
                            items=("summary_items", []),
                            headers=("summary_headers", []),
                            hide_default_footer=True,
                            disable_pagination=True,
                        )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def _parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db-dir",
        default=None,
        help=(
            "Calculation database root. Example: "
            "/gpfs/scratch/ahurtado/project_data/data/raman_paper/data"
        ),
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host interface to bind (default: 127.0.0.1)",
    )
    parser.add_argument(
        "--port",
        default=9020,
        type=int,
        help="Port to bind (default: 9020; use 0 for auto)",
    )
    return parser.parse_args(argv)


def _select_free_port(host: str, preferred_port: int, *, max_tries: int = 50) -> int:
    bind_host = host or "127.0.0.1"
    if preferred_port == 0:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind((bind_host, 0))
            return int(sock.getsockname()[1])

    port = int(preferred_port)
    for _ in range(max_tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                sock.bind((bind_host, port))
                return port
            except OSError as exc:
                if exc.errno == errno.EADDRINUSE:
                    port += 1
                    continue
                raise

    raise RuntimeError(
        f"Could not find a free port starting at {preferred_port} after {max_tries} attempts"
    )


def main(argv: list[str] | None = None) -> int:
    global _DB_DIR

    args = _parse_args(argv)
    if args.db_dir:
        _DB_DIR = Path(args.db_dir).expanduser().resolve()

    # Ensure cache and UI state are rebuilt using the final CLI dataset path.
    _data.cache_clear()
    _refresh_dashboard()

    host = str(args.host)
    port = _select_free_port(host, int(args.port))

    url_host = host
    if url_host in ("0.0.0.0", "::"):
        url_host = "127.0.0.1"

    print(f"Raman dashboard running at http://{url_host}:{port}/")
    server.start(host=host, port=port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
