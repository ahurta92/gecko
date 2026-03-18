"""Dataset dashboard for frequency-dependent polarizability and Raman trends.

Purpose
-------
Provide a single launch point for paper-oriented exploration of project_data
Raman calculations, reusing Gecko's existing data-loading and Trame patterns.

This dashboard intentionally mirrors the H2O notebook plotting style and
extends it to all systems/bases in the dataset:
- alpha_xx, alpha_yy, alpha_zz, alpha_avg vs frequency (+ optional fit)
- Raman pol_int vs frequency for up to 3 modes (+ optional fit)
- KPI cards + trend charts + breakdown table

Run
---
python -m gecko.viz.apps.polar_raman_dashboard \
  --db-dir /gpfs/scratch/ahurtado/project_data/data/raman_paper/data
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import errno
import io
import json
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


_DEFAULT_DB_DIR = Path("/gpfs/scratch/ahurtado/project_data/data/raman_paper/data")
_DB_DIR: Path | None = None

_ALPHA_DIAG = ["xx", "yy", "zz"]
_SMALL_THRESH = 1e-12


@dataclass
class DashboardData:
    root: Path
    n_calc_dirs: int
    n_loaded_calcs: int
    load_errors: list[str]
    alpha_df: pd.DataFrame
    raman_df: pd.DataFrame
    molecules: list[str]
    bases: list[str]


def _repo_root() -> Path:
    # src/gecko/viz/apps/polar_raman_dashboard.py -> apps -> viz -> gecko -> src -> repo
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
    if isinstance(value, float) and np.isnan(value):
        return None
    out = str(value).strip()
    if not out:
        return None
    if out.lower() in {"nan", "none", "null"}:
        return None
    return out


def _infer_molecule_from_root(root: str | Path) -> str:
    p = Path(root)
    parts = p.parts
    data_indices = [i for i, part in enumerate(parts) if part == "data"]
    if data_indices:
        idx = data_indices[-1]
        if idx + 1 < len(parts):
            return str(parts[idx + 1])
    if p.parent != p:
        return str(p.parent.name)
    return str(p.name)


def _infer_run_key_from_root(root: str | Path) -> str:
    p = Path(root)
    parts = p.parts
    data_indices = [i for i, part in enumerate(parts) if part == "data"]
    if data_indices:
        idx = data_indices[-1]
        # .../data/<molecule>/<run_key>/...
        if idx + 2 < len(parts):
            return str(parts[idx + 2])
    # Fallback: parent folder name if no canonical data layout is detected.
    if p.parent != p:
        return str(p.parent.name)
    return str(p.name)


def _decorate_basis_label(code: Any, basis: Any, run_key: Any) -> str:
    basis_label = _clean_string(basis) or "unknown"
    code_label = _clean_string(code) or ""
    run_label = _clean_string(run_key) or ""

    # Keep conventional Dalton basis labels unchanged.
    if code_label.lower() != "madness":
        return basis_label

    # If parser only reports a generic MRA basis, infer a more specific base
    # label from the run key when available (e.g., mra-p07-hbm... -> mra-p07).
    if basis_label.lower() in {"mra", "unknown"} and run_label.startswith("mra-"):
        parts = run_label.split("-")
        if len(parts) >= 2:
            basis_label = "-".join(parts[:2])

    # For MADNESS, expose each run as a distinct effective basis variant.
    if run_label and run_label != basis_label:
        return f"{basis_label} | {run_label}"
    return basis_label


def _attach_molecule_column(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        out = df.copy()
        out["molecule"] = []
        out["run_key"] = []
        return out

    def _pick_molecule(row: pd.Series) -> str:
        for key in ("mol_id", "molecule_id", "label"):
            val = _clean_string(row.get(key))
            if val:
                return val
        root_val = row.get("root")
        if root_val is not None:
            return _infer_molecule_from_root(root_val)
        return "unknown"

    def _pick_run_key(row: pd.Series) -> str:
        root_val = row.get("root")
        if root_val is None:
            return ""
        return _infer_run_key_from_root(root_val)

    out = df.copy()
    out["molecule"] = out.apply(_pick_molecule, axis=1)
    out["run_key"] = out.apply(_pick_run_key, axis=1)

    if "basis" in out.columns:
        out["basis"] = out.apply(
            lambda row: _decorate_basis_label(row.get("code"), row.get("basis"), row.get("run_key")),
            axis=1,
        )

    return out


def _alpha_from_properties_rows(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    # Expected flat-row schema example:
    # {"property": "polarizability", "freqB": 0.02, "component": ["x", "x"], "value": ...}
    entries: list[tuple[float, str, float]] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        if str(row.get("property", "")).strip().lower() != "polarizability":
            continue

        comp_raw = row.get("component")
        if not isinstance(comp_raw, list) or len(comp_raw) != 2:
            continue

        try:
            i = str(comp_raw[0]).strip().lower().replace("dipole_", "")
            j = str(comp_raw[1]).strip().lower().replace("dipole_", "")
            ij = f"{i}{j}"
            freq = float(row.get("freqB"))
            value = float(row.get("value"))
        except Exception:
            continue

        if len(ij) != 2:
            continue
        entries.append((freq, ij, value))

    if not entries:
        return None

    freqs = sorted(set(freq for freq, _, _ in entries))
    comps = sorted(set(comp for _, comp, _ in entries))
    if not freqs or not comps:
        return None

    vals = np.full((len(freqs), len(comps)), np.nan, dtype=float)
    fi = {f: idx for idx, f in enumerate(freqs)}
    cj = {c: idx for idx, c in enumerate(comps)}

    for freq, comp, value in entries:
        vals[fi[freq], cj[comp]] = value

    return {
        "omega": np.asarray(freqs, dtype=float),
        "components": comps,
        "values": vals,
        "shape": ("freq", "component"),
    }


def _hydrate_calc_from_properties(calc: Any) -> bool:
    """Best-effort fallback for MADNESS run fragments that only have properties.json."""
    if not isinstance(getattr(calc, "data", None), dict):
        return False

    alpha = calc.data.get("alpha")
    has_alpha = False
    if isinstance(alpha, dict):
        try:
            omega = np.asarray(alpha.get("omega", []), dtype=float).reshape(-1)
        except Exception:
            omega = np.asarray([], dtype=float)
        components = alpha.get("components", [])
        has_alpha = omega.size > 0 and len(list(components)) > 0
    if has_alpha:
        return False

    root = Path(getattr(calc, "root", ""))
    props_path = root / "properties.json"
    if not props_path.exists():
        return False

    try:
        payload = json.loads(props_path.read_text(encoding="utf-8"))
    except Exception:
        return False

    if not isinstance(payload, list):
        return False

    alpha_tensor = _alpha_from_properties_rows(payload)
    if alpha_tensor is None:
        return False

    calc.data["alpha"] = alpha_tensor
    calc.meta.setdefault("alpha_source", "properties.json")
    return True


def _collect_calcs(root: Path) -> tuple[list[Any], list[str], int]:
    calc_dirs = list(iter_calc_dirs(root))
    calcs: list[Any] = []
    errors: list[str] = []

    for calc_dir in calc_dirs:
        try:
            # Some parsers print diagnostics to stdout; keep dashboard logs clean.
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                calc = load_calc(calc_dir)
            _hydrate_calc_from_properties(calc)
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
            molecules=[],
            bases=[],
        )

    calcs, load_errors, n_calc_dirs = _collect_calcs(root)
    tb = TableBuilder(calcs)

    alpha_df = _attach_molecule_column(tb.build_alpha())
    raman_df = _attach_molecule_column(tb.build_raman())

    molecules = sorted(set(alpha_df.get("molecule", [])) | set(raman_df.get("molecule", [])))
    bases = sorted(set(alpha_df.get("basis", [])) | set(raman_df.get("basis", [])))
    bases = [b for b in bases if _clean_string(b)]

    return DashboardData(
        root=root,
        n_calc_dirs=n_calc_dirs,
        n_loaded_calcs=len(calcs),
        load_errors=load_errors,
        alpha_df=alpha_df,
        raman_df=raman_df,
        molecules=[str(m) for m in molecules if _clean_string(m)],
        bases=[str(b) for b in bases],
    )


def _as_bool(value: Any, *, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        out = [_clean_string(v) for v in value]
        return [str(v) for v in out if v]
    item = _clean_string(value)
    return [item] if item else []


def _power_series_fit(omega: np.ndarray, values: np.ndarray) -> dict[str, float] | None:
    mask = np.isfinite(omega) & np.isfinite(values)
    if int(mask.sum()) < 3:
        return None

    om = omega[mask]
    val = values[mask]

    X = np.column_stack([np.ones_like(om), om**2, om**4])
    try:
        coeff, _, rank, _ = np.linalg.lstsq(X, val, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if int(rank) < 3:
        return None

    c0 = float(coeff[0])
    c2 = float(coeff[1])
    c4 = float(coeff[2])
    fit = X @ coeff

    ss_res = float(np.sum((val - fit) ** 2))
    ss_tot = float(np.sum((val - np.mean(val)) ** 2))
    r2 = np.nan if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot

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
        "r2": float(r2),
    }


def _evaluate_fit(omega: np.ndarray, fit: dict[str, float]) -> np.ndarray:
    f0 = _safe_float(fit.get("f0"))
    A = _safe_float(fit.get("A"))
    B = _safe_float(fit.get("B"))
    if f0 is None or A is None or B is None:
        return np.full_like(omega, np.nan, dtype=float)
    return f0 * (1.0 + A * (omega**2) + B * (omega**4))


def _coverage_frame(data: DashboardData) -> pd.DataFrame:
    molecules = list(data.molecules)
    bases = list(data.bases)
    if not molecules or not bases:
        return pd.DataFrame(
            columns=[
                "molecule",
                "basis",
                "alpha_omega_count",
                "alpha_points",
                "raman_omega_count",
                "raman_mode_count",
                "raman_points",
                "has_alpha",
                "has_raman",
                "complete_row",
            ]
        )

    grid = pd.MultiIndex.from_product([molecules, bases], names=["molecule", "basis"]).to_frame(index=False)

    alpha_diag = data.alpha_df[data.alpha_df.get("ij", "").isin(_ALPHA_DIAG)].copy()
    if alpha_diag.empty:
        alpha_cov = pd.DataFrame(columns=["molecule", "basis", "alpha_omega_count", "alpha_points"])
    else:
        alpha_cov = (
            alpha_diag.groupby(["molecule", "basis"], as_index=False)
            .agg(
                alpha_omega_count=("omega", pd.Series.nunique),
                alpha_points=("value", "count"),
            )
        )

    if data.raman_df.empty:
        raman_cov = pd.DataFrame(columns=["molecule", "basis", "raman_omega_count", "raman_mode_count", "raman_points"])
    else:
        raman_cov = (
            data.raman_df.groupby(["molecule", "basis"], as_index=False)
            .agg(
                raman_omega_count=("omega_pol", pd.Series.nunique),
                raman_mode_count=("mode", pd.Series.nunique),
                raman_points=("pol_int", "count"),
            )
        )

    out = grid.merge(alpha_cov, on=["molecule", "basis"], how="left")
    out = out.merge(raman_cov, on=["molecule", "basis"], how="left")

    for col in [
        "alpha_omega_count",
        "alpha_points",
        "raman_omega_count",
        "raman_mode_count",
        "raman_points",
    ]:
        out[col] = out[col].fillna(0).astype(int)

    out["has_alpha"] = out["alpha_points"] > 0
    out["has_raman"] = out["raman_points"] > 0
    out["complete_row"] = out["has_alpha"] & out["has_raman"]

    return out.sort_values(["molecule", "basis"]).reset_index(drop=True)


def _summary_items(cov: pd.DataFrame) -> list[dict[str, Any]]:
    if cov.empty:
        return []

    items: list[dict[str, Any]] = []
    for _, row in cov.iterrows():
        items.append(
            {
                "molecule": str(row.get("molecule", "")),
                "basis": str(row.get("basis", "")),
                "complete": "yes" if bool(row.get("complete_row", False)) else "no",
                "alpha_omega": int(row.get("alpha_omega_count", 0)),
                "alpha_points": int(row.get("alpha_points", 0)),
                "raman_omega": int(row.get("raman_omega_count", 0)),
                "raman_modes": int(row.get("raman_mode_count", 0)),
                "raman_points": int(row.get("raman_points", 0)),
            }
        )
    return items


def _data_url_from_figure(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _basis_items_for_molecule(data: DashboardData, molecule: str) -> list[str]:
    if not molecule:
        return []
    cov = _coverage_frame(data)
    if cov.empty:
        return []
    chunk = cov[cov["molecule"] == molecule]
    if chunk.empty:
        return []
    return sorted(chunk["basis"].astype(str).unique().tolist())


def _alpha_compare_plot_data_url(
    data: DashboardData,
    *,
    molecule: str,
    bases: list[str],
    reference_basis: str,
    yscale: str,
    show_fit: bool,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    if not bases:
        return ""

    chunk = data.alpha_df[
        (data.alpha_df["molecule"] == molecule)
        & (data.alpha_df["basis"].isin(bases))
        & (data.alpha_df["ij"].isin(_ALPHA_DIAG))
    ].copy()
    if chunk.empty:
        return ""

    avg_by_basis: dict[str, pd.Series] = {}
    for basis in bases:
        bchunk = chunk[chunk["basis"] == basis].copy()
        if bchunk.empty:
            continue
        pivot = bchunk.pivot_table(index="omega", columns="ij", values="value", aggfunc="first").sort_index()
        if pivot.empty:
            continue
        available = [c for c in _ALPHA_DIAG if c in pivot.columns]
        if not available:
            continue
        avg_values = np.nanmean(np.column_stack([np.asarray(pivot[c], dtype=float) for c in available]), axis=1)
        avg_by_basis[str(basis)] = pd.Series(avg_values, index=np.asarray(pivot.index, dtype=float)).sort_index()

    if not avg_by_basis:
        return ""

    has_reference = reference_basis in avg_by_basis and len(avg_by_basis) > 1
    nrows = 2 if has_reference else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(9.8, 4.5 + (2.8 if has_reference else 0.0)), dpi=160, sharex=has_reference)
    if isinstance(axes, np.ndarray):
        ax_main = axes[0]
        ax_delta = axes[1] if has_reference else None
    else:
        ax_main = axes
        ax_delta = None

    cmap = plt.get_cmap("tab10")
    basis_order = [b for b in bases if b in avg_by_basis]

    for idx, basis in enumerate(basis_order):
        series = avg_by_basis[basis]
        om = np.asarray(series.index, dtype=float)
        val = np.asarray(series.values, dtype=float)
        color = cmap(idx % 10)
        is_ref = basis == reference_basis
        ax_main.plot(
            om,
            val,
            marker="o",
            linewidth=2.2 if is_ref else 1.5,
            alpha=0.95 if is_ref else 0.85,
            label=basis,
            color=color,
        )

        if show_fit:
            fit = _power_series_fit(om, val)
            if fit is not None:
                xfit = np.linspace(float(np.nanmin(om)), float(np.nanmax(om)), 220)
                yfit = _evaluate_fit(xfit, fit)
                ax_main.plot(xfit, yfit, linestyle="--", linewidth=1.0, color=color, alpha=0.7)

    ax_main.set_title(f"alpha_avg comparison: {molecule}")
    ax_main.set_ylabel("alpha_avg (a.u.)")
    ax_main.grid(True, alpha=0.25)
    if yscale == "symlog":
        ax_main.set_yscale("symlog", linthresh=1e-6)
    ax_main.legend(loc="best", fontsize=8)

    if has_reference and ax_delta is not None:
        ref_series = avg_by_basis[reference_basis]
        ref_index = set(float(x) for x in ref_series.index)
        for idx, basis in enumerate(basis_order):
            if basis == reference_basis:
                continue
            series = avg_by_basis[basis]
            common = sorted(ref_index.intersection(float(x) for x in series.index))
            if not common:
                continue
            om = np.asarray(common, dtype=float)
            delta = np.asarray([series[o] - ref_series[o] for o in om], dtype=float)
            color = cmap(idx % 10)
            ax_delta.plot(om, delta, marker="o", linewidth=1.4, alpha=0.90, label=f"{basis} - {reference_basis}", color=color)

        ax_delta.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
        ax_delta.set_ylabel("Delta alpha_avg")
        ax_delta.set_xlabel("Frequency (a.u.)")
        ax_delta.grid(True, alpha=0.25)
        if yscale == "symlog":
            ax_delta.set_yscale("symlog", linthresh=1e-6)
        ax_delta.legend(loc="best", fontsize=8)
    else:
        ax_main.set_xlabel("Frequency (a.u.)")

    fig.tight_layout()
    out = _data_url_from_figure(fig)
    plt.close(fig)
    return out


def _raman_compare_plot_data_url(
    data: DashboardData,
    *,
    molecule: str,
    bases: list[str],
    reference_basis: str,
    yscale: str,
    show_fit: bool,
    max_modes: int = 3,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    if not bases:
        return ""

    chunk = data.raman_df[
        (data.raman_df["molecule"] == molecule)
        & (data.raman_df["basis"].isin(bases))
    ].copy()
    if chunk.empty:
        return ""

    modes = sorted(set(int(m) for m in chunk["mode"].dropna().astype(int).tolist()))
    if not modes:
        return ""
    modes = modes[: max(1, int(max_modes))]

    fig, axes = plt.subplots(len(modes), 1, figsize=(9.8, 2.8 * len(modes) + 0.8), dpi=160, sharex=True)
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])

    cmap = plt.get_cmap("tab10")
    basis_order = [b for b in bases if b in set(chunk["basis"].astype(str))]

    for ridx, mode in enumerate(modes):
        ax = axes[ridx]
        mchunk = chunk[chunk["mode"].astype(int) == int(mode)].copy()
        if mchunk.empty:
            continue

        for bidx, basis in enumerate(basis_order):
            bchunk = mchunk[mchunk["basis"] == basis].copy()
            if bchunk.empty:
                continue
            bchunk = bchunk.sort_values("omega_pol")

            om = np.asarray(bchunk["omega_pol"], dtype=float)
            pol = np.asarray(bchunk["pol_int"], dtype=float)
            color = cmap(bidx % 10)
            is_ref = basis == reference_basis

            ax.plot(
                om,
                pol,
                marker="o",
                linewidth=2.2 if is_ref else 1.4,
                alpha=0.95 if is_ref else 0.85,
                label=basis,
                color=color,
            )

            if show_fit:
                fit = _power_series_fit(om, pol)
                if fit is not None:
                    xfit = np.linspace(float(np.nanmin(om)), float(np.nanmax(om)), 220)
                    yfit = _evaluate_fit(xfit, fit)
                    ax.plot(xfit, yfit, linestyle="--", linewidth=1.0, color=color, alpha=0.7)

        freq_cm1 = np.nan
        if "freq_cm1" in mchunk.columns:
            vals = np.asarray(mchunk["freq_cm1"], dtype=float)
            finite = vals[np.isfinite(vals)]
            if finite.size > 0:
                freq_cm1 = float(np.nanmean(finite))

        title = f"Mode {int(mode)}"
        if np.isfinite(freq_cm1):
            title += f" (~{freq_cm1:.1f} cm^-1)"

        ax.set_title(title, loc="left", fontsize=11)
        ax.set_ylabel("pol_int")
        ax.grid(True, alpha=0.25)
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=1e-6)
        ax.legend(loc="best", fontsize=8)

    axes[-1].set_xlabel("Frequency (a.u.)")
    fig.suptitle(f"Raman intensity comparison: {molecule} (reference: {reference_basis})", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out = _data_url_from_figure(fig)
    plt.close(fig)
    return out


def _mode_frequency_plot_data_url(
    data: DashboardData,
    *,
    molecule: str,
    bases: list[str],
    reference_basis: str,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    if not bases:
        return ""

    chunk = data.raman_df[
        (data.raman_df["molecule"] == molecule)
        & (data.raman_df["basis"].isin(bases))
    ].copy()
    if chunk.empty or "freq_cm1" not in chunk.columns:
        return ""

    chunk = chunk[np.isfinite(chunk["freq_cm1"].astype(float))].copy()
    if chunk.empty:
        return ""

    summary = (
        chunk.groupby(["basis", "mode"], as_index=False)
        .agg(freq_cm1=("freq_cm1", "mean"))
        .sort_values(["mode", "basis"])
    )
    if summary.empty:
        return ""

    modes = sorted(set(int(m) for m in summary["mode"].dropna().astype(int).tolist()))
    basis_order = [b for b in bases if b in set(summary["basis"].astype(str))]

    has_reference = reference_basis in basis_order and len(basis_order) > 1
    nrows = 2 if has_reference else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(9.8, 4.2 + (2.4 if has_reference else 0.0)), dpi=160, sharex=True)
    if isinstance(axes, np.ndarray):
        ax_main = axes[0]
        ax_delta = axes[1] if has_reference else None
    else:
        ax_main = axes
        ax_delta = None

    cmap = plt.get_cmap("tab10")

    for idx, basis in enumerate(basis_order):
        bsum = summary[summary["basis"] == basis].copy().sort_values("mode")
        if bsum.empty:
            continue
        x = np.asarray(bsum["mode"], dtype=float)
        y = np.asarray(bsum["freq_cm1"], dtype=float)
        color = cmap(idx % 10)
        is_ref = basis == reference_basis
        ax_main.plot(
            x,
            y,
            marker="o",
            linewidth=2.2 if is_ref else 1.5,
            alpha=0.95 if is_ref else 0.85,
            label=basis,
            color=color,
        )

    ax_main.set_title(f"Mode frequency by basis: {molecule}")
    ax_main.set_ylabel("Mode frequency (cm^-1)")
    ax_main.grid(True, alpha=0.25)
    ax_main.legend(loc="best", fontsize=8)

    if has_reference and ax_delta is not None:
        ref = summary[summary["basis"] == reference_basis].copy().set_index("mode")
        ref_modes = set(int(m) for m in ref.index.tolist())
        for idx, basis in enumerate(basis_order):
            if basis == reference_basis:
                continue
            bsum = summary[summary["basis"] == basis].copy().set_index("mode")
            common = sorted(ref_modes.intersection(int(m) for m in bsum.index.tolist()))
            if not common:
                continue
            x = np.asarray(common, dtype=float)
            y = np.asarray([float(bsum.loc[m, "freq_cm1"] - ref.loc[m, "freq_cm1"]) for m in common], dtype=float)
            color = cmap(idx % 10)
            ax_delta.plot(x, y, marker="o", linewidth=1.4, alpha=0.9, label=f"{basis} - {reference_basis}", color=color)

        ax_delta.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
        ax_delta.set_ylabel("Delta freq (cm^-1)")
        ax_delta.set_xlabel("Mode index")
        ax_delta.grid(True, alpha=0.25)
        ax_delta.legend(loc="best", fontsize=8)
    else:
        ax_main.set_xlabel("Mode index")

    fig.tight_layout()
    out = _data_url_from_figure(fig)
    plt.close(fig)
    return out


def _refresh_dashboard() -> None:

    data = _data()

    try:
        state.status = "loading"
        state.status_color = "orange"
        state.last_error = ""

        if not data.molecules or not data.bases:
            state.status = "no data"
            state.status_color = "red"
            state.last_error = "No alpha/raman records found in dataset root."
            state.summary_items = []
            state.alpha_plot_src = ""
            state.raman_plot_src = ""
            return

        if not state.selected_molecule or state.selected_molecule not in data.molecules:
            state.selected_molecule = data.molecules[0]

        basis_items = _basis_items_for_molecule(data, str(state.selected_molecule))
        if not basis_items:
            basis_items = data.bases
        state.basis_items = basis_items

        selected_bases = [b for b in _as_list(state.selected_bases) if b in basis_items]
        if not selected_bases:
            selected_bases = basis_items[: min(4, len(basis_items))]
        state.selected_bases = selected_bases

        if not state.reference_basis or state.reference_basis not in selected_bases:
            state.reference_basis = selected_bases[0]

        cov = _coverage_frame(data)
        state.summary_items = _summary_items(cov)

        state.kpi_total_molecules = int(len(data.molecules))
        state.kpi_total_basis_rows = int(len(cov))
        state.kpi_complete_rows = int(cov["complete_row"].sum()) if not cov.empty else 0
        state.kpi_load_errors = int(len(data.load_errors))

        state.alpha_plot_src = _alpha_compare_plot_data_url(
            data,
            molecule=str(state.selected_molecule),
            bases=_as_list(state.selected_bases),
            reference_basis=str(state.reference_basis),
            yscale=str(state.plot_scale),
            show_fit=_as_bool(state.show_fit, default=True),
        )
        state.raman_plot_src = _raman_compare_plot_data_url(
            data,
            molecule=str(state.selected_molecule),
            bases=_as_list(state.selected_bases),
            reference_basis=str(state.reference_basis),
            yscale=str(state.plot_scale),
            show_fit=_as_bool(state.show_fit, default=True),
            max_modes=3,
        )
        state.mode_freq_plot_src = _mode_frequency_plot_data_url(
            data,
            molecule=str(state.selected_molecule),
            bases=_as_list(state.selected_bases),
            reference_basis=str(state.reference_basis),
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

state.setdefault("selected_molecule", "")
state.setdefault("selected_bases", [])
state.setdefault("reference_basis", "")
state.setdefault("basis_items", [])
state.setdefault("show_fit", True)
state.setdefault("plot_scale", "linear")

state.setdefault("summary_items", [])
state.setdefault("alpha_plot_src", "")
state.setdefault("raman_plot_src", "")
state.setdefault("mode_freq_plot_src", "")

state.setdefault("kpi_total_molecules", 0)
state.setdefault("kpi_total_basis_rows", 0)
state.setdefault("kpi_complete_rows", 0)
state.setdefault("kpi_load_errors", 0)

state.setdefault(
    "summary_headers",
    [
        {"text": "Molecule", "value": "molecule"},
        {"text": "Basis", "value": "basis"},
        {"text": "Complete", "value": "complete"},
        {"text": "alpha omegas", "value": "alpha_omega"},
        {"text": "alpha points", "value": "alpha_points"},
        {"text": "raman omegas", "value": "raman_omega"},
        {"text": "raman modes", "value": "raman_modes"},
        {"text": "raman points", "value": "raman_points"},
    ],
)


@state.change("selected_molecule")
def _on_molecule_change(**_):
    _refresh_dashboard()


@state.change("selected_bases", "reference_basis", "show_fit", "plot_scale")
def _on_controls_change(**_):
    _refresh_dashboard()


with SinglePageWithDrawerLayout(server) as layout:
    layout.title.set_text("Polarizability + Raman Dataset Dashboard")

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
        drawer.width = 390
        vuetify.VContainer(fluid=True, classes="pa-2")

        vuetify.VAlert(
            f"Dataset root: {_resolve_db_dir()}",
            type="info",
            dense=True,
            outlined=True,
            classes="mb-2",
        )

        vuetify.VAlert(
            "Assumptions: simplest paper-launch view. Trends mirror H2O notebook style; Raman panel shows up to first 3 modes; fits use f(omega)=c0+c2*w^2+c4*w^4.",
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
            items=("molecule_items", _data().molecules),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Basis sets (compare)",
            v_model=("selected_bases", state.selected_bases),
            items=("basis_items", state.basis_items),
            multiple=True,
            chips=True,
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Reference basis",
            v_model=("reference_basis", state.reference_basis),
            items=("selected_bases", state.selected_bases),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Plot y-scale",
            v_model=("plot_scale", state.plot_scale),
            items=("plot_scale_items", ["linear", "symlog"]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSwitch(
            label="Overlay power-series fit",
            v_model=("show_fit", state.show_fit),
            dense=True,
            hide_details=True,
            classes="mb-2",
        )

    with layout.content:
        with vuetify.VContainer(fluid=True, classes="pa-2"):
            # KPI row
            with vuetify.VRow(dense=True, classes="mb-2"):
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Molecules", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_total_molecules }}", classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Molecule/Basis Rows", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_total_basis_rows }}", classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Complete Rows", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_complete_rows }}", classes="text-h5")
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VCardTitle("Load Errors", classes="text-subtitle-2")
                        vuetify.VCardText("{{ kpi_load_errors }}", classes="text-h5")

            # Trends
            with vuetify.VRow(dense=True, classes="mb-2"):
                with vuetify.VCol(cols=6):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Polarizability Comparison (alpha_avg)")
                        vuetify.VImg(
                            src=("alpha_plot_src", ""),
                            contain=True,
                            v_if="alpha_plot_src",
                            style="width: 100%;",
                        )
                with vuetify.VCol(cols=6):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Raman Comparison (pol_int, Modes 1-3)")
                        vuetify.VImg(
                            src=("raman_plot_src", ""),
                            contain=True,
                            v_if="raman_plot_src",
                            style="width: 100%;",
                        )


            with vuetify.VRow(dense=True, classes="mb-2"):
                with vuetify.VCol(cols=12):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Mode Frequency By Basis")
                        vuetify.VImg(
                            src=("mode_freq_plot_src", ""),
                            contain=True,
                            v_if="mode_freq_plot_src",
                            style="width: 100%;",
                        )

            # Breakdown table
            with vuetify.VRow(dense=True):
                with vuetify.VCol(cols=12):
                    with vuetify.VCard(outlined=True, classes="pa-2"):
                        vuetify.VSubheader("Breakdown Table (coverage by molecule/basis)")
                        vuetify.VDataTable(
                            dense=True,
                            items=("summary_items", []),
                            headers=("summary_headers", []),
                            hide_default_footer=True,
                            disable_pagination=True,
                        )

            _refresh_dashboard()


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
        default=9030,
        type=int,
        help="Port to bind (default: 9030; use 0 for auto)",
    )
    parser.add_argument(
        "--server",
        action="store_true",
        help="Compatibility flag for headless runs; browser auto-open is disabled by default.",
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

    host = str(args.host)
    port = _select_free_port(host, int(args.port))

    url_host = host
    if url_host in ("0.0.0.0", "::"):
        url_host = "127.0.0.1"

    print(f"Polar/Raman dashboard running at http://{url_host}:{port}/")
    try:
        server.start(host=host, port=port, open_browser=False)
    except TypeError:
        # Backward compatibility with older Trame versions.
        server.start(host=host, port=port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
