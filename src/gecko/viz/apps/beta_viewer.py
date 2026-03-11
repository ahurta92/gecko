"""Trame server app for viewing SHG unit-sphere fields and error metrics.

This replaces the old demo that read a static `.vtu` file.

Run:
    python -m gecko.viz.apps.beta_viewer

Then open the printed URL in a browser.

"""

from __future__ import annotations

import argparse
import base64
import json
import io
import os
import re
import time
from types import SimpleNamespace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Tuple


_DEBUG_STARTUP = os.environ.get("BETA_TRAME_DEBUG", "0") == "1"
_MOL_RESOLVER = None

_GLOBAL_CLIM_PERCENTILE = 99.0


def _dprint(msg: str) -> None:
    if _DEBUG_STARTUP:
        print(msg, flush=True)

try:
    import numpy as np
except ModuleNotFoundError as e:  # pragma: no cover
    raise SystemExit(
        "Missing Python dependencies (numpy). You're likely running with the wrong interpreter.\n"
        "Try one of:\n"
        "  - conda activate madqc && python notebooks/scripts/application.py\n"
        "  - conda run -n madqc python notebooks/scripts/application.py\n"
        "Or select the same environment you use in the notebooks."
    ) from e

_dprint("[startup] numpy imported")

from trame.app import get_server
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import vtk, vuetify

_dprint("[startup] trame imported")

from trame_vtk.modules.vtk.serializers import configure_serializer

from vtkmodules.vtkCommonCore import vtkLookupTable
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkCamera,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)

_dprint("[startup] vtk modules imported")

# Rendering backend import (safe to include)
import vtkmodules.vtkRenderingOpenGL2  # noqa


configure_serializer(encode_lut=True, skip_light=True)

_dprint("[startup] serializer configured")


def _repo_root() -> Path:
    # src/gecko/viz/apps/beta_viewer.py -> apps -> viz -> gecko -> src -> repo
    return Path(__file__).resolve().parents[4]

from gecko.viz import vtk_scene as _vtk_scene
from gecko.viz.fields import (
    ErrorSettings,
    compute_error_fields,
    evaluate_field,
    load_lebedev_grid,
    tensor_from_long as _tensor_from_long_impl,
)
from gecko.viz.state import (
    FIELD_CHOICES as _FIELD_CHOICES,
    METRIC_CHOICES as _METRIC_CHOICES,
    auto_clim as _auto_clim_impl,
    default_state as _default_state,
    metric_style as _metric_style_impl,
)


_SHG_CSV_PATH: Path | None = None
_SHG_DB_DIR: Path | None = None
_GEOM_MAP_PATH: Path | None = None
_BUNDLE_DIR: Path | None = None
_WRITE_BUNDLE_DIR: Path | None = None

_GLOBAL_BREAKDOWN_MAX_ROWS = 300
_GLOBAL_TREND_PERCENTILE = 95.0
_GLOBAL_CLUSTER_MAX = 8

GLOBAL_DASH_METRIC_CHOICES: List[Tuple[str, str]] = [
    ("L2 relative error (rel_L2)", "rel_L2"),
    ("Mean signed parallel-relative error", "mean_par_rel_signed"),
    ("P95 signed parallel-relative error", "p95_par_rel_signed"),
    ("Mean angle [deg]", "mean_ang"),
    ("P95 angle [deg]", "p95_ang"),
]


@lru_cache(maxsize=1)
def _shg_long():
    from gecko.viz.io import load_shg_df_from_csv

    csv_path = _SHG_CSV_PATH or (_repo_root() / "data" / "csv_data" / "shg_ijk.csv")
    print(f"Loading SHG tensors from {csv_path} ...", flush=True)
    try:
        return load_shg_df_from_csv(csv_path)
    except FileNotFoundError:
        _dprint(f"[shg] CSV not found: {csv_path}")
        import pandas as pd

        return pd.DataFrame(columns=["molecule", "basis", "omega", "ijk", "Beta"])


def _normalize_label(label: str | None) -> str | None:
    if label is None:
        return None
    cleaned = str(label).strip()
    cleaned = re.sub(r"\s+", "", cleaned)
    return cleaned or None


def _load_geometry_map_from_json(path: Path) -> Dict[str, Dict[str, Any]]:
    try:
        raw = path.read_text(encoding="utf-8")
        data = json.loads(raw)
    except Exception as exc:
        _dprint(f"[geom] failed to read geometry map: {path}: {exc}")
        return {}
    if not isinstance(data, dict):
        _dprint(f"[geom] geometry map must be a JSON object: {path}")
        return {}
    mapped: Dict[str, Dict[str, Any]] = {}
    for key, value in data.items():
        norm = _normalize_label(str(key))
        if norm is None:
            continue
        if isinstance(value, dict) and "symbols" in value and "geometry" in value:
            mapped[norm] = value
    return mapped


@lru_cache(maxsize=1)
def _geometry_map() -> Dict[str, Dict[str, Any]]:
    if _GEOM_MAP_PATH is not None and _GEOM_MAP_PATH.exists():
        return _load_geometry_map_from_json(_GEOM_MAP_PATH)
    if _BUNDLE_DIR is not None:
        geom_path = _BUNDLE_DIR / "geometries.json"
        if geom_path.exists():
            return _load_geometry_map_from_json(geom_path)

    df = _data()
    if "geometry" not in df.columns:
        return {}
    from gecko.viz.io import geometry_map_from_df

    return geometry_map_from_df(df, key="molecule")


def _beta_df_to_np(beta_df) -> np.ndarray:
    ijk_map = {"X": 0, "Y": 1, "Z": 2}
    beta_np = np.zeros((3, 3, 3), dtype=float)
    for ijk, value in beta_df.items():
        i, j, k = ijk_map[ijk[0]], ijk_map[ijk[1]], ijk_map[ijk[2]]
        try:
            beta_np[i, j, k] = float(value)
        except Exception:
            beta_np[i, j, k] = 0.0
    return beta_np


def _tensor_from_long(shg_long, mol: str, basis: str, omega: float | int) -> np.ndarray:
    key = (str(mol), str(basis), float(omega))
    tensor = _tensor_lookup().get(key)
    if tensor is None:
        tensor = _tensor_from_long_impl(shg_long, mol, basis, omega)
    if not np.any(tensor):
        _dprint(f"[shg] Missing tensor for (mol={mol}, basis={basis}, omega={omega}); using zeros")
    return tensor


@lru_cache(maxsize=1)
def _data():
    # Intentionally only load the SHG tensor table here.
    if _SHG_DB_DIR is not None:
        from gecko.viz.io import build_shg_df_from_db, write_beta_viewer_bundle

        df = build_shg_df_from_db(_SHG_DB_DIR, include_geometry=True)
        if _WRITE_BUNDLE_DIR is not None:
            try:
                write_beta_viewer_bundle(df, _WRITE_BUNDLE_DIR)
                _dprint(f"[bundle] wrote beta viewer bundle to {_WRITE_BUNDLE_DIR}")
            except Exception as exc:
                _dprint(f"[bundle] failed to write bundle: {exc}")
        return df

    return _shg_long()


@lru_cache(maxsize=1)
def _tensor_lookup() -> Dict[Tuple[str, str, float], np.ndarray]:
    df = _data()
    required = {"molecule", "basis", "omega", "ijk", "Beta"}
    if not required.issubset(set(df.columns)):
        return {}

    cols = df[["molecule", "basis", "omega", "ijk", "Beta"]].dropna(
        subset=["molecule", "basis", "omega", "ijk"]
    )

    lookup: Dict[Tuple[str, str, float], np.ndarray] = {}
    for (mol, basis, omega), group in cols.groupby(
        ["molecule", "basis", "omega"], sort=False
    ):
        series = {
            str(ijk): beta
            for ijk, beta in zip(group["ijk"], group["Beta"], strict=False)
        }
        lookup[(str(mol), str(basis), float(omega))] = _beta_df_to_np(series)
    return lookup


def _active_data_source_label() -> str:
    if _SHG_DB_DIR is not None:
        return f"db:{_SHG_DB_DIR}"
    if _BUNDLE_DIR is not None:
        return f"bundle:{_BUNDLE_DIR}"
    if _SHG_CSV_PATH is not None:
        return f"csv:{_SHG_CSV_PATH}"
    return "csv:data/csv_data/shg_ijk.csv"


@lru_cache(maxsize=256)
def _load_molecule(mol_name: str):
    from gecko.mol.io import read_mol

    norm_label = _normalize_label(mol_name)
    geom_map = _geometry_map()
    if norm_label is not None and norm_label in geom_map:
        try:
            import qcelemental as _qcel

            payload = geom_map[norm_label]
            symbols = payload.get("symbols")
            geometry = payload.get("geometry")
            if symbols is not None and geometry is not None:
                return _qcel.models.Molecule(symbols=symbols, geometry=geometry)
        except Exception as exc:
            _dprint(f"[molecule] failed to build molecule from geometry map: {exc}")

    if _MOL_RESOLVER is not None:
        _dprint("[molecule] MoleculeResolver support removed; ignoring _MOL_RESOLVER")

    mol_path = _repo_root() / "data" / "molecules" / f"{mol_name}.mol"
    if mol_path.exists():
        try:
            return read_mol(mol_path)
        except Exception as exc:
            _dprint(f"[molecule] failed to load fallback .mol ({mol_path}): {exc}")

    return None


# -----------------------------------------------------------------------------
# VTK scene
# -----------------------------------------------------------------------------


def _new_scene(
    *,
    renderWindow: vtkRenderWindow,
    viewport: Tuple[float, float, float, float],
    background: Tuple[float, float, float] = (0.08, 0.08, 0.10),
) -> Dict[str, Any]:
    return _vtk_scene.new_scene(renderWindow=renderWindow, viewport=viewport, background=background)


def _create_render_window() -> vtkRenderWindow:
    return _vtk_scene.create_render_window(debug_fn=_dprint)


_dprint("[startup] creating VTK render window")
renderWindow = _create_render_window()

# Keep an interactor attached for compatibility with some trame-vtk helpers.
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

_shared_camera = vtkCamera()
_camera_initialized = False


def _rotate_beta(beta: np.ndarray, R: np.ndarray) -> np.ndarray:
        """Rotate beta tensor into a new coordinate system.

        Conventions:
        - Vectors are treated as column vectors in the math, but in this code we
            often store coordinates as row vectors.
        - If we use `coords_new = coords_old @ R` for row vectors, then the column
            transform is `v_new = R.T @ v_old`.
        - Therefore beta transforms as:
                beta_new[i,j,k] = Rt[i,a] Rt[j,b] Rt[k,c] beta_old[a,b,c]
            where Rt = R.T.
        """
        beta = np.asarray(beta, dtype=float)
        R = np.asarray(R, dtype=float)
        Rt = R.T
        return np.einsum("ia,jb,kc,abc->ijk", Rt, Rt, Rt, beta, optimize=True)


def _principal_axes_rotation(mol) -> np.ndarray:
        """Compute a deterministic principal-axes frame with a phase convention.

        - Uses an inertia-tensor principal axes definition.
        - Enforces a right-handed basis.
        - Fixes axis signs by aligning axes with a reference direction derived from
            the heaviest atom.
        """

        import numpy as _np
        import qcelemental as _qcel

        if mol is None:
            return _np.eye(3, dtype=float)

        symbols = list(mol.symbols)
        coords = _np.asarray(mol.geometry, dtype=float)
        coords = coords - _np.mean(coords, axis=0, keepdims=True)

        masses = _np.array([
                float(_qcel.periodictable.to_mass(sym)) if sym is not None else 1.0 for sym in symbols
        ])
        if not _np.all(_np.isfinite(masses)):
                masses = _np.ones(len(symbols), dtype=float)

        # Inertia tensor about the origin: I = Σ m (r^2 I - r r^T)
        I = _np.zeros((3, 3), dtype=float)
        for r, m in zip(coords, masses, strict=False):
                r2 = float(_np.dot(r, r))
                I += m * (r2 * _np.eye(3) - _np.outer(r, r))

        evals, evecs = _np.linalg.eigh(I)
        order = _np.argsort(evals)
        R = evecs[:, order]

        # Phase convention: use the direction to the heaviest atom.
        ref_idx = int(_np.argmax(masses)) if masses.size else 0
        ref = coords[ref_idx] if coords.size else _np.array([1.0, 0.0, 0.0])
        if _np.linalg.norm(ref) < 1e-12:
                ref = _np.array([1.0, 0.0, 0.0])
        dots = R.T @ ref
        signs = _np.where(dots >= 0, 1.0, -1.0)
        R = R * signs

        # Ensure right-handed
        if _np.linalg.det(R) < 0:
                R[:, 2] *= -1.0
        return R

_scene_left = _new_scene(renderWindow=renderWindow, viewport=(0.0, 0.0, 0.5, 1.0))
_scene_mid = _new_scene(renderWindow=renderWindow, viewport=(0.5, 0.0, 1.0, 1.0))
_scene_right = _new_scene(renderWindow=renderWindow, viewport=(0.0, 0.0, 0.0, 0.0))
_scene_left["renderer"].SetActiveCamera(_shared_camera)
_scene_mid["renderer"].SetActiveCamera(_shared_camera)
_scene_right["renderer"].SetActiveCamera(_shared_camera)

_dprint("[startup] VTK render window created")


def _polydata_from_points(points: np.ndarray) -> vtkPolyData:
    return _vtk_scene.polydata_from_points(points)


def _set_vectors(pd: vtkPolyData, name: str, vectors: np.ndarray) -> None:
    _vtk_scene.set_vectors(pd, name, vectors)


def _set_scalars(pd: vtkPolyData, name: str, scalars: np.ndarray) -> None:
    _vtk_scene.set_scalars(pd, name, scalars)


def _add_array(pd: vtkPolyData, name: str, arr: np.ndarray) -> None:
    _vtk_scene.add_array(pd, name, arr)


def _element_rgb(symbol: str) -> Tuple[float, float, float]:
    return _vtk_scene.element_rgb(symbol)


def _default_clim(values: np.ndarray, *, symmetric: bool) -> Tuple[float, float]:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (0.0, 1.0)
    vmin, vmax = float(np.min(finite)), float(np.max(finite))
    if symmetric:
        bound = max(abs(vmin), abs(vmax))
        if bound == 0:
            bound = 1.0
        return (-bound, bound)
    if vmin == vmax:
        vmax = vmin + 1.0
    return (vmin, vmax)


def _auto_clim(values: np.ndarray, *, mode: str) -> Tuple[float, float]:
    return _auto_clim_impl(values, mode=mode)


def _metric_style(mode: str) -> Dict[str, Any]:
    """Return styling metadata for a scalar metric."""
    return _metric_style_impl(mode)


def _build_lut(*, kind: str, clim: Tuple[float, float], scalar_mode: str) -> vtkLookupTable:
    """Create a LUT consistent with our metric semantics."""
    return _vtk_scene.build_lut(kind=kind, clim=clim, scalar_mode=scalar_mode)


def _build_glyph_actor(
    *,
    n_hat: np.ndarray,
    vectors: np.ndarray,
    scalars: np.ndarray,
    scale: np.ndarray,
    glyph_factor: float,
    clim: Tuple[float, float],
    lut_kind: str,
    scalar_mode: str,
) -> vtkActor:
    return _vtk_scene.build_glyph_actor(
        points=n_hat,
        vectors=vectors,
        scalars=scalars,
        scale=scale,
        scale_factor=float(glyph_factor),
        clim=clim,
        lut_kind=str(lut_kind),
        scalar_mode=str(scalar_mode),
    )


def _build_atom_actor(points: np.ndarray, *, radius: float, rgb: Tuple[float, float, float]) -> vtkActor:
    return _vtk_scene.build_atom_actor(points, radius=radius, rgb=rgb)


# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------


server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# During module import, callbacks may run before the VTK view is created.
# Provide safe defaults until `ctrl.view_update` is set to the actual view's update.
ctrl.view_update = lambda **_: None


def _reset_camera() -> None:
    global _camera_initialized
    _scene_left["renderer"].ResetCamera()
    _camera_initialized = True
    ctrl.view_update()


ctrl.reset_camera = _reset_camera

state.setdefault("status", "starting")
state.setdefault("status_color", "orange")
state.setdefault("last_error", "")


def _basis_list() -> List[str]:
    df = _data()
    if "basis" not in df.columns:
        return []
    return sorted(set(str(x) for x in df["basis"].unique()))


def _molecule_list() -> List[str]:
    df = _data()
    if "molecule" not in df.columns:
        return []
    return sorted(set(str(x) for x in df["molecule"].unique()))


def _omega_list(mol: str) -> List[float]:
    df = _data()
    if "molecule" in df.columns and "omega" in df.columns:
        sub = df[df["molecule"] == str(mol)]
        if not sub.empty:
            return sorted(set(float(x) for x in sub["omega"].unique()))
    if "omega" not in df.columns:
        return []
    return sorted(set(float(x) for x in df["omega"].unique()))


METRIC_CHOICES = _METRIC_CHOICES
FIELD_CHOICES = _FIELD_CHOICES


def _scalar_mode_label(mode: str) -> str:
    for label, value in METRIC_CHOICES:
        if value == mode:
            return str(label)
    return str(mode)


def _rel_norm_label() -> str:
    kind = str(state.rel_norm)
    if kind == "pointwise":
        return "|v_ref(n)|"
    if kind == "global_max":
        return "max(|v_ref|)"
    if kind == "global_mean":
        return "mean(|v_ref|)"
    if kind == "global_rms":
        return "rms(|v_ref|)"
    return kind


def _finite_concat(chunks: List[np.ndarray]) -> np.ndarray:
    if not chunks:
        return np.asarray([], dtype=float)
    arr = np.concatenate([np.asarray(x, dtype=float).reshape(-1) for x in chunks], axis=0)
    return arr[np.isfinite(arr)]


@lru_cache(maxsize=64)
def _global_metric_clims(
    mol: str,
    omega: float,
    ref_basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
) -> Dict[str, Tuple[float, float]]:
    """Compute robust clims for each metric across all basis sets.

    Returned clims are intended for visualization only.
    """
    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()

    beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))

    R = None
    if coord_system == "principal_axes":
        qmol = _load_molecule(mol)
        R = _principal_axes_rotation(qmol)
        beta_ref = _rotate_beta(beta_ref, R)

    if isinstance(grid, tuple) and len(grid) == 2:
        n_hat = np.asarray(grid[0], dtype=float)
        w = np.asarray(grid[1], dtype=float)
    else:
        n_hat = np.asarray(grid.n_hat, dtype=float)
        w = np.asarray(grid.w, dtype=float)

    v_ref = evaluate_field(beta_ref, n_hat)
    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )

    # Collect values for each metric across all basis sets.
    metric_values: Dict[str, List[np.ndarray]] = {}
    for basis in _basis_list():
        if str(basis) == str(ref_basis):
            continue
        try:
            beta_bas = _tensor_from_long(df, mol, str(basis), float(omega))
        except Exception:
            continue
        if R is not None:
            beta_bas = _rotate_beta(beta_bas, R)
        v_bas = evaluate_field(beta_bas, n_hat)
        arrays, _ = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)

        for key in (
            "ref_mag",
            "bas_mag",
            "abs_err",
            "rel_err",
            "ang_err",
            "par_rel",
            "perp_rel",
        ):
            if key in arrays:
                metric_values.setdefault(key, []).append(np.asarray(arrays[key]))

    clims: Dict[str, Tuple[float, float]] = {}
    for key, chunks in metric_values.items():
        vals = _finite_concat(chunks)
        style = _metric_style(key)
        if vals.size == 0:
            clims[key] = (-1.0, 1.0) if style["signed"] else (0.0, 1.0)
            continue
        if style["signed"]:
            bound = float(np.percentile(np.abs(vals), _GLOBAL_CLIM_PERCENTILE))
            if not np.isfinite(bound) or bound <= 0:
                bound = 1.0
            clims[key] = (-bound, bound)
        else:
            vmax = float(np.percentile(vals, _GLOBAL_CLIM_PERCENTILE))
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = 1.0
            clims[key] = (0.0, vmax)

    # ref_mag is common; if missing, fall back to local behavior.
    return clims


# -----------------------------------------------------------------------------
# Metric plots across basis sets
# -----------------------------------------------------------------------------


_PLOT_FAMILY_ORDER: List[str] = [
    "aug-cc-pVnZ",
    "aug-cc-pCVnZ",
    "d-aug-cc-pVnZ",
    "d-aug-cc-pCVnZ",
]

_PLOT_CARD_ORDER: List[str] = ["DZ", "TZ", "QZ"]


def _basis_family_and_cardinal(basis: str) -> Tuple[str | None, str | None]:
    """Return (family, cardinal) for a basis name.

    Supported families:
      aug-cc-pVnZ, aug-cc-pCVnZ, d-aug-cc-pVnZ, d-aug-cc-pCVnZ
    Supported cardinals:
      DZ/TZ/QZ
    """
    b = re.sub(r"\s+", "", str(basis)).strip()

    # pV family: aug-cc-pVDZ / aug-cc-pVTZ / aug-cc-pVQZ (and d-aug-...)
    m = re.match(r"^(d-aug|aug)-cc-pV(D|T|Q)Z$", b)
    if m:
        aug_prefix = "d-aug" if m.group(1) == "d-aug" else "aug"
        card = f"{m.group(2)}Z"
        return f"{aug_prefix}-cc-pVnZ", card

    # pCV family: aug-cc-pCVDZ / aug-cc-pCVTZ / aug-cc-pCVQZ (and d-aug-...)
    m = re.match(r"^(d-aug|aug)-cc-pCV(D|T|Q)Z$", b)
    if m:
        aug_prefix = "d-aug" if m.group(1) == "d-aug" else "aug"
        card = f"{m.group(2)}Z"
        return f"{aug_prefix}-cc-pCVnZ", card

    return None, None


PLOT_METRIC_CHOICES: List[Tuple[str, str]] = [
    ("L2 relative error (rel_L2)", "rel_L2"),
    ("Mean signed rel parallel (mean_par_rel_signed)", "mean_par_rel_signed"),
    (
        "Mean + / Mean − signed rel parallel (posneg_par_rel_signed)",
        "posneg_par_rel_signed",
    ),
    ("P95 signed rel parallel (p95_par_rel_signed)", "p95_par_rel_signed"),
    ("Mean angle [deg] (mean_ang)", "mean_ang"),
    ("P95 angle [deg] (p95_ang)", "p95_ang"),
]


def _safe_weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).reshape(-1)
    w = np.asarray(w, dtype=float).reshape(-1)
    if x.size != w.size or x.size == 0:
        return float("nan")
    finite = np.isfinite(x) & np.isfinite(w)
    if not np.any(finite):
        return float("nan")
    return float(np.sum(w[finite] * x[finite]) / (4.0 * np.pi))


def _safe_weighted_l2(x: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float).reshape(-1)
    w = np.asarray(w, dtype=float).reshape(-1)
    if x.size != w.size or x.size == 0:
        return float("nan")
    finite = np.isfinite(x) & np.isfinite(w)
    if not np.any(finite):
        return float("nan")
    val = float(np.sum(w[finite] * (x[finite] ** 2)) / (4.0 * np.pi))
    return float(np.sqrt(max(0.0, val)))


def _derived_error_metrics(arrays: Dict[str, np.ndarray], w: np.ndarray) -> Dict[str, float]:
    rel_err = np.asarray(arrays.get("rel_err", np.asarray([])), dtype=float)
    par_abs = np.asarray(arrays.get("par_abs", np.asarray([])), dtype=float)
    abs_err = np.asarray(arrays.get("abs_err", np.asarray([])), dtype=float)
    ref_mag = np.asarray(arrays.get("ref_mag", np.asarray([])), dtype=float)
    ang_err = np.asarray(arrays.get("ang_err", np.asarray([])), dtype=float)
    par_rel_signed = np.asarray(arrays.get("par_rel_signed", np.asarray([])), dtype=float)
    dv_tan = np.asarray(arrays.get("dv_tan", np.asarray([])), dtype=float)

    rel_finite = rel_err[np.isfinite(rel_err)]
    ang_finite = ang_err[np.isfinite(ang_err)]
    par_rel_finite = par_rel_signed[np.isfinite(par_rel_signed)]

    rel_l2 = float("nan")
    ref_l2 = _safe_weighted_l2(ref_mag, w)
    err_l2 = _safe_weighted_l2(abs_err, w)
    if np.isfinite(ref_l2) and ref_l2 > 0 and np.isfinite(err_l2):
        rel_l2 = float(err_l2 / ref_l2)

    l2_par = _safe_weighted_l2(par_abs, w)
    if dv_tan.size:
        l2_perp = _safe_weighted_l2(dv_tan, w)
    else:
        if np.isfinite(err_l2) and np.isfinite(l2_par):
            l2_perp = float(np.sqrt(max(err_l2 * err_l2 - l2_par * l2_par, 0.0)))
        else:
            l2_perp = float("nan")

    ratio_perp_par = float("nan")
    if np.isfinite(l2_par) and l2_par > 0 and np.isfinite(l2_perp):
        ratio_perp_par = float(l2_perp / l2_par)

    if "mask" in arrays:
        mask = np.asarray(arrays["mask"], dtype=float).reshape(-1) > 0.5
        masked_rel = rel_err[mask] if rel_err.size == mask.size else np.asarray([])
        masked_w = np.asarray(w, dtype=float).reshape(-1)[mask] if rel_err.size == mask.size else np.asarray([])
        masked_rel_l2 = _safe_weighted_l2(masked_rel, masked_w)
    else:
        masked_rel_l2 = _safe_weighted_l2(rel_err, w)

    bias_ratio = float("nan")
    if par_rel_finite.size:
        pos = np.mean(par_rel_finite[par_rel_finite > 0]) if np.any(par_rel_finite > 0) else np.nan
        neg = np.mean(np.abs(par_rel_finite[par_rel_finite < 0])) if np.any(par_rel_finite < 0) else np.nan
        if np.isfinite(pos) and np.isfinite(neg) and neg > 0:
            bias_ratio = float(pos / neg)

    out = {
        "rel_L2": rel_l2,
        "L2_par": l2_par,
        "L2_perp": l2_perp,
        "ratio_perp_par": ratio_perp_par,
        "max_rel": float(np.max(rel_finite)) if rel_finite.size else float("nan"),
        "p95_rel": float(np.percentile(rel_finite, 95)) if rel_finite.size else float("nan"),
        "mean_ang": _safe_weighted_mean(ang_err, w),
        "p95_ang": float(np.percentile(ang_finite, 95)) if ang_finite.size else float("nan"),
        "bias_ratio": bias_ratio,
        "masked_rel_L2": masked_rel_l2,
        "mean_par_rel_signed": _safe_weighted_mean(par_rel_signed, w),
        "p95_par_rel_signed": float(np.percentile(par_rel_finite, 95)) if par_rel_finite.size else float("nan"),
    }
    return out


def _posneg_l2_from_signed_field(x: np.ndarray, w: np.ndarray) -> Tuple[float, float]:
    """Return (L2_positive, L2_negative) for a signed field over the sphere.

    L2 norms are computed using the Lebedev weights with the convention:
      L2(x) = sqrt( (1/(4π)) * Σ w_i * x_i^2 )

    The negative contribution is returned as a negative value so it plots below
    zero: (-L2(|min(x,0)|)).
    """

    x = np.asarray(x, dtype=float).reshape(-1)
    w = np.asarray(w, dtype=float).reshape(-1)
    if x.size != w.size:
        raise ValueError("Signed field and weights must have the same length")

    finite = np.isfinite(x) & np.isfinite(w)
    if not np.any(finite):
        return (float("nan"), float("nan"))

    xf = x[finite]
    wf = w[finite]
    pos_sq = float(np.sum(wf * (np.maximum(xf, 0.0) ** 2)) / (4.0 * np.pi))
    neg_sq = float(np.sum(wf * (np.minimum(xf, 0.0) ** 2)) / (4.0 * np.pi))
    pos = float(np.sqrt(max(0.0, pos_sq)))
    neg = -float(np.sqrt(max(0.0, neg_sq)))
    return pos, neg


def _metric_plot_value(
    metric_key: str,
    *,
    arrays: Dict[str, np.ndarray],
    metrics: Dict[str, Any],
    w: np.ndarray,
) -> float:
    key = str(metric_key)
    if key in metrics:
        return float(metrics[key])
    derived = _derived_error_metrics(arrays, w)
    if key in derived:
        return float(derived[key])
    if key == "mean_par_rel_signed":
        x = np.asarray(arrays["par_rel_signed"], dtype=float)
        return float(np.sum(w * x) / (4.0 * np.pi))
    if key == "p95_par_rel_signed":
        x = np.asarray(arrays["par_rel_signed"], dtype=float)
        x = x[np.isfinite(x)]
        return float(np.percentile(x, 95)) if x.size else float("nan")
    raise KeyError(f"Unknown plot metric '{metric_key}'")


@lru_cache(maxsize=32)
def _posneg_par_rel_signed_across_bases(
    mol: str,
    omega: float,
    ref_basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
) -> Dict[str, Tuple[float, float]]:
    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()

    beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))

    R = None
    if coord_system == "principal_axes":
        mol_obj = _load_molecule(mol)
        R = _principal_axes_rotation(mol_obj)
        beta_ref = _rotate_beta(beta_ref, R)

    n_hat = np.asarray(grid.n_hat, dtype=float)
    w = np.asarray(grid.w, dtype=float)
    v_ref = evaluate_field(beta_ref, n_hat)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )

    out: Dict[str, Tuple[float, float]] = {}
    for basis in _basis_list():
        family, card = _basis_family_and_cardinal(str(basis))
        if family is None or card is None:
            continue
        if family not in _PLOT_FAMILY_ORDER:
            continue
        if card not in _PLOT_CARD_ORDER:
            continue

        beta_bas = _tensor_from_long(df, mol, str(basis), float(omega))
        if R is not None:
            beta_bas = _rotate_beta(beta_bas, R)
        v_bas = evaluate_field(beta_bas, n_hat)
        arrays, _metrics = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)

        pos, neg = _posneg_l2_from_signed_field(arrays["par_rel_signed"], w)
        out[str(basis)] = (pos, neg)

    return out


@lru_cache(maxsize=128)
def _posneg_par_rel_signed_across_omegas_for_basis(
    mol: str,
    ref_basis: str,
    basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
) -> Dict[float, Tuple[float, float]]:
    """Return {omega: (L2_pos, L2_neg)} for a single basis across all omegas."""
    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()

    R = None
    if coord_system == "principal_axes":
        mol_obj = _load_molecule(mol)
        R = _principal_axes_rotation(mol_obj)

    n_hat = np.asarray(grid.n_hat, dtype=float)
    w = np.asarray(grid.w, dtype=float)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )

    out: Dict[float, Tuple[float, float]] = {}
    for omega in _omega_list(mol):
        try:
            beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))
            beta_bas = _tensor_from_long(df, mol, basis, float(omega))
        except KeyError:
            continue
        if R is not None:
            beta_ref = _rotate_beta(beta_ref, R)
            beta_bas = _rotate_beta(beta_bas, R)
        v_ref = evaluate_field(beta_ref, n_hat)
        v_bas = evaluate_field(beta_bas, n_hat)
        arrays, _metrics = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)
        out[float(omega)] = _posneg_l2_from_signed_field(arrays["par_rel_signed"], w)

    return out


@lru_cache(maxsize=32)
def _metric_values_across_bases(
    mol: str,
    omega: float,
    ref_basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
    metric_key: str,
) -> Dict[str, float]:
    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()

    beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))

    R = None
    if coord_system == "principal_axes":
        mol_obj = _load_molecule(mol)
        R = _principal_axes_rotation(mol_obj)
        beta_ref = _rotate_beta(beta_ref, R)

    n_hat = np.asarray(grid.n_hat, dtype=float)
    w = np.asarray(grid.w, dtype=float)
    v_ref = evaluate_field(beta_ref, n_hat)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )

    values: Dict[str, float] = {}
    for basis in _basis_list():
        family, card = _basis_family_and_cardinal(str(basis))
        if family is None or card is None:
            continue
        if family not in _PLOT_FAMILY_ORDER:
            continue
        if card not in _PLOT_CARD_ORDER:
            continue
        beta_bas = _tensor_from_long(df, mol, str(basis), float(omega))
        if R is not None:
            beta_bas = _rotate_beta(beta_bas, R)
        v_bas = evaluate_field(beta_bas, n_hat)
        arrays, metrics = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)
        values[str(basis)] = _metric_plot_value(metric_key, arrays=arrays, metrics=metrics, w=w)

    return values


@lru_cache(maxsize=128)
def _metric_values_across_omegas_for_basis(
    mol: str,
    ref_basis: str,
    basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
    metric_key: str,
) -> Dict[float, float]:
    """Return {omega: metric_value} for a single basis across all omegas."""
    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()

    R = None
    if coord_system == "principal_axes":
        mol_obj = _load_molecule(mol)
        R = _principal_axes_rotation(mol_obj)

    n_hat = np.asarray(grid.n_hat, dtype=float)
    w = np.asarray(grid.w, dtype=float)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )

    out: Dict[float, float] = {}
    for omega in _omega_list(mol):
        try:
            beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))
            beta_bas = _tensor_from_long(df, mol, basis, float(omega))
        except KeyError:
            continue
        if R is not None:
            beta_ref = _rotate_beta(beta_ref, R)
            beta_bas = _rotate_beta(beta_bas, R)
        v_ref = evaluate_field(beta_ref, n_hat)
        v_bas = evaluate_field(beta_bas, n_hat)
        arrays, metrics = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)
        out[float(omega)] = _metric_plot_value(metric_key, arrays=arrays, metrics=metrics, w=w)

    return out


@lru_cache(maxsize=32)
def _metric_plot_data_url(
    mol: str,
    omega: float,
    ref_basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
    metric_key: str,
    plot_x: str,
    omega_agg: str,
) -> str:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    plot_x = str(plot_x)
    omega_agg = str(omega_agg)
    yscale = str(state.plot_yscale)

    metric_key = str(metric_key)

    if plot_x == "omega":
        fig, ax = plt.subplots(figsize=(9.0, 5.0), dpi=160)
        omegas = [float(x) for x in _omega_list(mol)]

        # Plot up to 12 lines: each family/cardinal (DZ/TZ/QZ) across omega.
        # For the special pos/neg plot, each basis becomes two lines (+ solid, − dashed).
        bases = _basis_list()
        for fam in _PLOT_FAMILY_ORDER:
            for card in _PLOT_CARD_ORDER:
                # Find the basis string for this family/card.
                basis_name = None
                for b in bases:
                    f, c = _basis_family_and_cardinal(b)
                    if f == fam and c == card:
                        basis_name = b
                        break
                if basis_name is None:
                    continue

                if metric_key == "posneg_par_rel_signed":
                    series = _posneg_par_rel_signed_across_omegas_for_basis(
                        mol,
                        ref_basis,
                        str(basis_name),
                        int(lebedev_order),
                        coord_system,
                        rel_norm,
                    )
                    xs = [w for w in omegas if w in series]
                    ys_pos = [float(series[w][0]) for w in xs]
                    ys_neg = [float(series[w][1]) for w in xs]
                    if not xs:
                        continue
                    ax.plot(xs, ys_pos, marker="o", linewidth=1.5, label=f"{fam}:{card} (+)")
                    ax.plot(
                        xs,
                        ys_neg,
                        marker="o",
                        linewidth=1.5,
                        linestyle="--",
                        label=f"{fam}:{card} (-)",
                    )
                else:
                    series = _metric_values_across_omegas_for_basis(
                        mol,
                        ref_basis,
                        str(basis_name),
                        int(lebedev_order),
                        coord_system,
                        rel_norm,
                        metric_key,
                    )
                    xs = [w for w in omegas if w in series]
                    ys = [float(series[w]) for w in xs]
                    if not xs:
                        continue
                    ax.plot(xs, ys, marker="o", linewidth=1.5, label=f"{fam}:{card}")

        ax.set_xlabel("omega")
        ax.set_ylabel("mean(par_rel_signed)+ / mean(par_rel_signed)-" if metric_key == "posneg_par_rel_signed" else str(metric_key))
        ax.grid(True, alpha=0.25)
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=1e-3)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()

    else:
        # Default: x-axis is basis set blocks.
        # Optionally show omega-variation as error bars (mean ± std across all ω).

        # Arrange x as blocks: each family gets DZ/TZ/QZ slots.
        x_positions: Dict[Tuple[str, str], float] = {}
        tick_pos: List[float] = []
        tick_lab: List[str] = []
        family_centers: Dict[str, float] = {}
        x = 0.0
        gap = 0.8
        step = 1.0
        for fam in _PLOT_FAMILY_ORDER:
            start = x
            for card in _PLOT_CARD_ORDER:
                x_positions[(fam, card)] = x
                tick_pos.append(x)
                tick_lab.append(card)
                x += step
            end = x - step
            family_centers[fam] = 0.5 * (start + end)
            x += gap

        fig, ax = plt.subplots(figsize=(9.0, 4.6), dpi=160)

        bases = _basis_list()
        for fam in _PLOT_FAMILY_ORDER:
            xs: List[float] = []
            ys: List[float] = []
            yerr: List[float] = []

            xs2: List[float] = []
            ys2: List[float] = []
            yerr2: List[float] = []

            for card in _PLOT_CARD_ORDER:
                basis_name = None
                for b in bases:
                    f, c = _basis_family_and_cardinal(b)
                    if f == fam and c == card:
                        basis_name = b
                        break
                if basis_name is None:
                    continue

                x = float(x_positions[(fam, card)])
                if omega_agg == "mean_std":
                    if metric_key == "posneg_par_rel_signed":
                        series = _posneg_par_rel_signed_across_omegas_for_basis(
                            mol,
                            ref_basis,
                            str(basis_name),
                            int(lebedev_order),
                            coord_system,
                            rel_norm,
                        )
                        vals = np.asarray(list(series.values()), dtype=float)
                        if vals.size == 0:
                            continue
                        pos = vals[:, 0]
                        neg = vals[:, 1]
                        pos = pos[np.isfinite(pos)]
                        neg = neg[np.isfinite(neg)]
                        if pos.size == 0 or neg.size == 0:
                            continue
                        y = float(np.mean(pos))
                        e = float(np.std(pos))
                        y2 = float(np.mean(neg))
                        e2 = float(np.std(neg))
                    else:
                        series = _metric_values_across_omegas_for_basis(
                            mol,
                            ref_basis,
                            str(basis_name),
                            int(lebedev_order),
                            coord_system,
                            rel_norm,
                            metric_key,
                        )
                        vals = np.asarray(list(series.values()), dtype=float)
                        vals = vals[np.isfinite(vals)]
                        if vals.size == 0:
                            continue
                        y = float(np.mean(vals))
                        e = float(np.std(vals))
                else:
                    # current ω
                    if metric_key == "posneg_par_rel_signed":
                        values = _posneg_par_rel_signed_across_bases(
                            mol,
                            float(omega),
                            ref_basis,
                            int(lebedev_order),
                            coord_system,
                            rel_norm,
                        )
                        if str(basis_name) not in values:
                            continue
                        y, y2 = values[str(basis_name)]
                        e, e2 = 0.0, 0.0
                    else:
                        values = _metric_values_across_bases(
                            mol,
                            float(omega),
                            ref_basis,
                            int(lebedev_order),
                            coord_system,
                            rel_norm,
                            metric_key,
                        )
                        if str(basis_name) not in values:
                            continue
                        y = float(values[str(basis_name)])
                        e = 0.0

                xs.append(x)
                ys.append(y)
                yerr.append(e)

                if metric_key == "posneg_par_rel_signed":
                    xs2.append(x)
                    ys2.append(float(y2))
                    yerr2.append(float(e2))

            if not xs:
                continue
            order = np.argsort(np.asarray(xs))
            xs = [xs[i] for i in order]
            ys = [ys[i] for i in order]
            yerr = [yerr[i] for i in order]

            if metric_key == "posneg_par_rel_signed" and xs2:
                order2 = np.argsort(np.asarray(xs2))
                xs2 = [xs2[i] for i in order2]
                ys2 = [ys2[i] for i in order2]
                yerr2 = [yerr2[i] for i in order2]

            if metric_key == "posneg_par_rel_signed":
                if omega_agg == "mean_std":
                    ax.errorbar(
                        xs,
                        ys,
                        yerr=yerr,
                        fmt="o-",
                        capsize=3,
                        linewidth=1.5,
                        label=f"{fam} (+)",
                    )
                    ax.errorbar(
                        xs2,
                        ys2,
                        yerr=yerr2,
                        fmt="^-",
                        capsize=3,
                        linewidth=1.5,
                        linestyle="--",
                        label=f"{fam} (-)",
                    )
                else:
                    ax.plot(xs, ys, marker="o", linewidth=1.5, label=f"{fam} (+)")
                    ax.plot(
                        xs2,
                        ys2,
                        marker="^",
                        linewidth=1.5,
                        linestyle="--",
                        label=f"{fam} (-)",
                    )
            else:
                if omega_agg == "mean_std":
                    ax.errorbar(
                        xs,
                        ys,
                        yerr=yerr,
                        fmt="o-",
                        capsize=3,
                        linewidth=1.5,
                        label=fam,
                    )
                else:
                    ax.plot(xs, ys, marker="o", linewidth=1.5, label=fam)

        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lab)
        if tick_pos:
            ax.set_xlim(min(tick_pos) - 0.5, max(tick_pos) + 0.5)
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_ylabel("mean(par_rel_signed)+ / mean(par_rel_signed)-" if metric_key == "posneg_par_rel_signed" else str(metric_key))
        if yscale == "symlog":
            ax.set_yscale("symlog", linthresh=1e-3)

        # Family labels over each block.
        if tick_pos:
            ylim = ax.get_ylim()
            y_text = ylim[1]
            for fam, cx in family_centers.items():
                ax.text(cx, y_text, fam, ha="center", va="bottom", fontsize=8)

        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    data = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{data}"


def _triple_metric_plot_data_url(
    mol: str,
    omega: float,
    ref_basis: str,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
    *,
    posneg_yscale: str,
) -> str:
    """Return a data URL for a fixed 3-panel summary plot.

    Panels (top -> bottom):
      1) rel_L2 across basis sets (linear)
      2) mean positive / mean negative par_rel_signed across basis sets (linear|symlog)
      3) mean_ang across basis sets (linear)
    """

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    # Arrange x as blocks: each family gets DZ/TZ/QZ slots.
    x_positions: Dict[Tuple[str, str], float] = {}
    tick_pos: List[float] = []
    tick_lab: List[str] = []
    family_centers: Dict[str, float] = {}
    x = 0.0
    gap = 0.8
    step = 1.0
    for fam in _PLOT_FAMILY_ORDER:
        start = x
        for card in _PLOT_CARD_ORDER:
            x_positions[(fam, card)] = x
            tick_pos.append(x)
            tick_lab.append(card)
            x += step
        end = x - step
        family_centers[fam] = 0.5 * (start + end)
        x += gap

    fig, (ax1, ax2, ax3) = plt.subplots(
        3,
        1,
        figsize=(9.0, 8.8),
        dpi=160,
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.15, 1.0]},
    )

    bases = _basis_list()

    # 1) rel_L2 (linear)
    values_l2 = _metric_values_across_bases(
        mol,
        float(omega),
        ref_basis,
        int(lebedev_order),
        coord_system,
        rel_norm,
        "rel_L2",
    )
    for fam in _PLOT_FAMILY_ORDER:
        xs: List[float] = []
        ys: List[float] = []
        for card in _PLOT_CARD_ORDER:
            basis_name = None
            for b in bases:
                f, c = _basis_family_and_cardinal(b)
                if f == fam and c == card:
                    basis_name = b
                    break
            if basis_name is None:
                continue
            if str(basis_name) not in values_l2:
                continue
            xs.append(float(x_positions[(fam, card)]))
            ys.append(float(values_l2[str(basis_name)]))
        if xs:
            order = np.argsort(np.asarray(xs))
            xs = [xs[i] for i in order]
            ys = [ys[i] for i in order]
            ax1.plot(xs, ys, marker="o", linewidth=1.5, label=fam)
    ax1.set_ylabel("rel_L2")
    ax1.grid(True, axis="y", alpha=0.25)
    ax1.legend(loc="best", fontsize=8)

    # 2) mean positive / mean negative par_rel_signed (linear|symlog)
    values_posneg = _posneg_par_rel_signed_across_bases(
        mol,
        float(omega),
        ref_basis,
        int(lebedev_order),
        coord_system,
        rel_norm,
    )
    for fam in _PLOT_FAMILY_ORDER:
        xs: List[float] = []
        ys_pos: List[float] = []
        ys_neg: List[float] = []
        for card in _PLOT_CARD_ORDER:
            basis_name = None
            for b in bases:
                f, c = _basis_family_and_cardinal(b)
                if f == fam and c == card:
                    basis_name = b
                    break
            if basis_name is None:
                continue
            if str(basis_name) not in values_posneg:
                continue
            pos, neg = values_posneg[str(basis_name)]
            xs.append(float(x_positions[(fam, card)]))
            ys_pos.append(float(pos))
            ys_neg.append(float(neg))
        if xs:
            order = np.argsort(np.asarray(xs))
            xs = [xs[i] for i in order]
            ys_pos = [ys_pos[i] for i in order]
            ys_neg = [ys_neg[i] for i in order]
            ax2.plot(xs, ys_pos, marker="o", linewidth=1.5, label=f"{fam} (+)")
            ax2.plot(
                xs,
                ys_neg,
                marker="^",
                linewidth=1.5,
                linestyle="--",
                label=f"{fam} (-)",
            )
    ax2.set_ylabel("L2 par_rel_signed (+ / -)")
    ax2.grid(True, axis="y", alpha=0.25)
    if str(posneg_yscale) == "symlog":
        ax2.set_yscale("symlog", linthresh=1e-3)
    ax2.legend(loc="best", fontsize=8)

    # 3) mean angle (linear)
    values_ang = _metric_values_across_bases(
        mol,
        float(omega),
        ref_basis,
        int(lebedev_order),
        coord_system,
        rel_norm,
        "mean_ang",
    )
    for fam in _PLOT_FAMILY_ORDER:
        xs: List[float] = []
        ys: List[float] = []
        for card in _PLOT_CARD_ORDER:
            basis_name = None
            for b in bases:
                f, c = _basis_family_and_cardinal(b)
                if f == fam and c == card:
                    basis_name = b
                    break
            if basis_name is None:
                continue
            if str(basis_name) not in values_ang:
                continue
            xs.append(float(x_positions[(fam, card)]))
            ys.append(float(values_ang[str(basis_name)]))
        if xs:
            order = np.argsort(np.asarray(xs))
            xs = [xs[i] for i in order]
            ys = [ys[i] for i in order]
            ax3.plot(xs, ys, marker="o", linewidth=1.5, label=fam)
    ax3.set_ylabel("mean_ang [deg]")
    ax3.grid(True, axis="y", alpha=0.25)
    ax3.legend(loc="best", fontsize=8)

    ax3.set_xticks(tick_pos)
    ax3.set_xticklabels(tick_lab)
    if tick_pos:
        ax3.set_xlim(min(tick_pos) - 0.5, max(tick_pos) + 0.5)

    # Family labels over each block (top plot, for compactness).
    if tick_pos:
        ylim = ax1.get_ylim()
        y_text = ylim[1]
        for fam, cx in family_centers.items():
            ax1.text(cx, y_text, fam, ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    data = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{data}"


def _update_metric_plot() -> None:
    try:
        src = _triple_metric_plot_data_url(
            str(state.mol),
            float(state.omega),
            str(state.ref_basis),
            int(state.lebedev_order),
            str(state.coord_system),
            str(state.rel_norm),
            posneg_yscale=str(state.posneg_yscale),
        )
        state.metric_plot_src = src
        state.metric_plot_error = "" if src else "Plot unavailable (missing matplotlib?)"
    except Exception as exc:
        state.metric_plot_src = ""
        state.metric_plot_error = f"Plot error: {type(exc).__name__}: {exc}"


@lru_cache(maxsize=1)
def _available_tensor_keys() -> set[Tuple[str, str, float]]:
    return set(_tensor_lookup().keys())


def _has_tensor(mol: str, basis: str, omega: float | int) -> bool:
    return (str(mol), str(basis), float(omega)) in _available_tensor_keys()


@lru_cache(maxsize=4096)
def _pair_metric_summary(
    mol: str,
    basis: str,
    ref_basis: str,
    omega: float,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
) -> Dict[str, float]:
    if str(basis) == str(ref_basis):
        return {}
    if not _has_tensor(mol, ref_basis, omega):
        return {}
    if not _has_tensor(mol, basis, omega):
        return {}

    grid = load_lebedev_grid(int(lebedev_order))
    df = _data()
    beta_ref = _tensor_from_long(df, mol, ref_basis, float(omega))
    beta_bas = _tensor_from_long(df, mol, basis, float(omega))

    R = None
    if coord_system == "principal_axes":
        mol_obj = _load_molecule(mol)
        R = _principal_axes_rotation(mol_obj)
        beta_ref = _rotate_beta(beta_ref, R)
        beta_bas = _rotate_beta(beta_bas, R)

    n_hat = np.asarray(grid.n_hat, dtype=float)
    w = np.asarray(grid.w, dtype=float)
    v_ref = evaluate_field(beta_ref, n_hat)
    v_bas = evaluate_field(beta_bas, n_hat)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(rel_norm),
    )
    arrays, metrics = compute_error_fields(v_ref, v_bas, n_hat, w, settings=settings)
    out = dict(metrics)
    out.update(_derived_error_metrics(arrays, w))
    return out


def _build_global_rows(
    *,
    metric_key: str,
    ref_basis: str,
    omega: float,
    lebedev_order: int,
    coord_system: str,
    rel_norm: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for mol in _molecule_list():
        if not _has_tensor(str(mol), str(ref_basis), float(omega)):
            continue
        for basis in _basis_list():
            if str(basis) == str(ref_basis):
                continue
            summary = _pair_metric_summary(
                str(mol),
                str(basis),
                str(ref_basis),
                float(omega),
                int(lebedev_order),
                str(coord_system),
                str(rel_norm),
            )
            if not summary:
                continue
            value = float(summary.get(str(metric_key), float("nan")))
            if not np.isfinite(value):
                continue
            family, cardinal = _basis_family_and_cardinal(str(basis))
            rows.append(
                {
                    "molecule": str(mol),
                    "basis": str(basis),
                    "family": str(family) if family is not None else "other",
                    "cardinal": str(cardinal) if cardinal is not None else "-",
                    "omega": float(omega),
                    "metric_value": value,
                }
            )
    return rows


def _kmeans_labels(X: np.ndarray, k: int, *, max_iter: int = 50) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    n = int(X.shape[0])
    if n == 0:
        return np.asarray([], dtype=int)
    k = max(1, min(int(k), n))
    # Deterministic initialization: evenly spaced row picks.
    seed_idx = np.linspace(0, n - 1, k, dtype=int)
    centers = X[seed_idx, :].copy()

    labels = np.zeros(n, dtype=int)
    for _ in range(max_iter):
        d2 = np.sum((X[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(d2, axis=1).astype(int)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for cid in range(k):
            mask = labels == cid
            if np.any(mask):
                centers[cid, :] = np.mean(X[mask, :], axis=0)
    return labels


def _assign_molecule_clusters(rows: List[Dict[str, Any]], *, cluster_k: int) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    if not rows:
        return rows, []

    molecules = sorted(set(str(r["molecule"]) for r in rows))
    basis_order = sorted(set(str(r["basis"]) for r in rows))
    by_mol_basis: Dict[Tuple[str, str], float] = {
        (str(r["molecule"]), str(r["basis"])): float(r["metric_value"]) for r in rows
    }

    X = np.full((len(molecules), len(basis_order)), np.nan, dtype=float)
    for i, mol in enumerate(molecules):
        for j, basis in enumerate(basis_order):
            if (mol, basis) in by_mol_basis:
                X[i, j] = float(by_mol_basis[(mol, basis)])

    col_mean = np.nanmean(X, axis=0)
    col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    inds = np.where(~np.isfinite(X))
    if inds[0].size:
        X[inds] = col_mean[inds[1]]
    col_std = np.std(X, axis=0)
    col_std = np.where(col_std > 0, col_std, 1.0)
    Xz = (X - np.mean(X, axis=0)) / col_std

    labels = _kmeans_labels(Xz, int(cluster_k))
    cluster_map = {mol: int(labels[i]) + 1 for i, mol in enumerate(molecules)}

    out_rows: List[Dict[str, Any]] = []
    for row in rows:
        row2 = dict(row)
        cid = cluster_map.get(str(row["molecule"]), 1)
        row2["cluster"] = f"C{cid}"
        out_rows.append(row2)

    cluster_rows: List[Dict[str, Any]] = []
    for cid in sorted(set(cluster_map.values())):
        mols = [m for m, c in cluster_map.items() if c == cid]
        vals = [float(r["metric_value"]) for r in rows if str(r["molecule"]) in mols]
        metric_mean = float(np.mean(vals)) if vals else float("nan")
        cluster_rows.append(
            {
                "cluster": f"C{cid}",
                "molecules": len(mols),
                "mean_metric": metric_mean,
            }
        )
    return out_rows, cluster_rows


def _global_trend_plot_data_url(rows: List[Dict[str, Any]], *, metric_key: str) -> str:
    if not rows:
        return ""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    grouped: Dict[Tuple[str, str], List[float]] = {}
    for row in rows:
        fam = str(row.get("family", "other"))
        card = str(row.get("cardinal", "-"))
        grouped.setdefault((fam, card), []).append(float(row["metric_value"]))

    x_positions: Dict[Tuple[str, str], float] = {}
    tick_pos: List[float] = []
    tick_lab: List[str] = []
    family_centers: Dict[str, float] = {}
    x = 0.0
    gap = 0.8
    step = 1.0
    for fam in _PLOT_FAMILY_ORDER:
        start = x
        for card in _PLOT_CARD_ORDER:
            x_positions[(fam, card)] = x
            tick_pos.append(x)
            tick_lab.append(card)
            x += step
        end = x - step
        family_centers[fam] = 0.5 * (start + end)
        x += gap

    fig, ax = plt.subplots(figsize=(9.0, 3.9), dpi=160)
    for fam in _PLOT_FAMILY_ORDER:
        xs: List[float] = []
        ys: List[float] = []
        yerr: List[float] = []
        for card in _PLOT_CARD_ORDER:
            vals = np.asarray(grouped.get((fam, card), []), dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            xs.append(float(x_positions[(fam, card)]))
            ys.append(float(np.mean(vals)))
            yerr.append(float(np.std(vals)))
        if xs:
            ax.errorbar(xs, ys, yerr=yerr, fmt="o-", capsize=3, linewidth=1.5, label=fam)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab)
    if tick_pos:
        ax.set_xlim(min(tick_pos) - 0.5, max(tick_pos) + 0.5)
    vals_all = np.asarray([float(r["metric_value"]) for r in rows], dtype=float)
    vals_all = vals_all[np.isfinite(vals_all)]
    if vals_all.size:
        ymax = float(np.percentile(np.abs(vals_all), _GLOBAL_TREND_PERCENTILE))
        if ymax > 0:
            if np.nanmin(vals_all) < 0:
                ax.set_ylim(-1.1 * ymax, 1.1 * ymax)
            else:
                ax.set_ylim(0.0, 1.1 * ymax)
    ax.grid(True, axis="y", alpha=0.25)
    ax.set_ylabel(metric_key)
    ax.legend(loc="best", fontsize=8)
    if tick_pos:
        ylim = ax.get_ylim()
        y_text = ylim[1]
        for fam, cx in family_centers.items():
            ax.text(cx, y_text, fam, ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    data = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{data}"


state.setdefault("mol", "H2O")
state.setdefault("omega", 0.0)
state.setdefault("ref_basis", "MRA")
state.setdefault("bas_basis_a", "aug-cc-pVDZ")
state.setdefault("bas_basis_b", "d-aug-cc-pVDZ")
state.setdefault("lebedev_order", 29)
state.setdefault("field_mode", "bas")
state.setdefault("scalar_mode", "par_rel_signed")
state.setdefault("glyph_factor", 0.12)
state.setdefault("molecule_scale", 0.7)
state.setdefault("atom_radius", 0.25)
state.setdefault("coord_system", "as_is")
state.setdefault("view_layout", "compare_bases")
state.setdefault("scale_by_metric", False)
state.setdefault("rel_norm", "global_rms")
state.setdefault("show_axes", True)
state.setdefault("metrics", {})
state.setdefault("metric_rows", [])
state.setdefault("posneg_yscale", "linear")
state.setdefault("metric_plot_src", "")
state.setdefault("metric_plot_error", "")
state.setdefault(
    "metric_headers",
    [
        {"text": "Metric", "value": "name"},
        {"text": "A", "value": "a"},
        {"text": "B", "value": "b"},
    ],
)
state.setdefault("dataset_db_dir", str(_SHG_DB_DIR) if _SHG_DB_DIR is not None else "")
state.setdefault("active_data_source", _active_data_source_label())
state.setdefault("dataset_status", "ready")
state.setdefault("dataset_status_color", "blue")
state.setdefault("dashboard_assumptions", "")
state.setdefault("molecule_items", _molecule_list())
state.setdefault("basis_items", _basis_list())
state.setdefault("omega_items", _omega_list(str(state.mol)))
state.setdefault("global_metric", "rel_L2")
state.setdefault("global_cluster_k", 3)
state.setdefault("global_auto_update", False)
state.setdefault("global_stale_note", "")
state.setdefault("global_last_runtime", "")
state.setdefault("global_kpi_molecules", "0")
state.setdefault("global_kpi_comparisons", "0")
state.setdefault("global_kpi_mean", "n/a")
state.setdefault("global_kpi_worst", "n/a")
state.setdefault("global_trend_src", "")
state.setdefault("global_trend_error", "")
state.setdefault("global_breakdown_rows", [])
state.setdefault(
    "global_breakdown_headers",
    [
        {"text": "Molecule", "value": "molecule"},
        {"text": "Basis", "value": "basis"},
        {"text": "Metric", "value": "metric_value"},
        {"text": "Cluster", "value": "cluster"},
        {"text": "Family", "value": "family"},
        {"text": "Cardinal", "value": "cardinal"},
    ],
)
state.setdefault(
    "global_cluster_headers",
    [
        {"text": "Cluster", "value": "cluster"},
        {"text": "Molecules", "value": "molecules"},
        {"text": "Mean metric", "value": "mean_metric"},
    ],
)
state.setdefault("global_cluster_rows", [])


@state.change("show_axes")
def _on_show_axes(show_axes, **_):
    vis = bool(show_axes)
    _scene_left["axes_actor"].SetVisibility(1 if vis else 0)
    _scene_mid["axes_actor"].SetVisibility(1 if vis else 0)
    _scene_right["axes_actor"].SetVisibility(1 if vis else 0)
    ctrl.view_update()


def _clear_cached_data_views() -> None:
    for fn in (
        _shg_long,
        _geometry_map,
        _data,
        _tensor_lookup,
        _available_tensor_keys,
        _load_molecule,
        _tensor_from_long,
        _global_metric_clims,
        _posneg_par_rel_signed_across_bases,
        _posneg_par_rel_signed_across_omegas_for_basis,
        _metric_values_across_bases,
        _metric_values_across_omegas_for_basis,
        _metric_plot_data_url,
        _triple_metric_plot_data_url,
        _pair_metric_summary,
    ):
        try:
            fn.cache_clear()
        except Exception:
            pass


def _refresh_selector_items() -> None:
    state.molecule_items = _molecule_list()
    state.basis_items = _basis_list()
    if state.molecule_items and str(state.mol) not in set(state.molecule_items):
        state.mol = str(state.molecule_items[0])
    omegas = _omega_list(str(state.mol))
    state.omega_items = [float(x) for x in omegas]
    if omegas and float(state.omega) not in {float(x) for x in omegas}:
        state.omega = float(omegas[0])


def _metric_rank_value(metric_key: str, metric_value: float) -> float:
    if not np.isfinite(metric_value):
        return float("-inf")
    if str(metric_key) in ("mean_par_rel_signed", "p95_par_rel_signed"):
        return float(abs(metric_value))
    return float(metric_value)


def _update_global_dashboard() -> None:
    t0 = time.perf_counter()
    state.global_stale_note = ""
    metric_key = str(state.global_metric)
    rows: List[Dict[str, Any]] = []
    try:
        rows = _build_global_rows(
            metric_key=metric_key,
            ref_basis=str(state.ref_basis),
            omega=float(state.omega),
            lebedev_order=int(state.lebedev_order),
            coord_system=str(state.coord_system),
            rel_norm=str(state.rel_norm),
        )
    except Exception as exc:
        state.global_kpi_molecules = "0"
        state.global_kpi_comparisons = "0"
        state.global_kpi_mean = "n/a"
        state.global_kpi_worst = "n/a"
        state.global_trend_src = ""
        state.global_trend_error = f"Global dashboard error: {type(exc).__name__}: {exc}"
        state.global_breakdown_rows = []
        state.global_cluster_rows = []
        state.global_last_runtime = f"{(time.perf_counter() - t0):.2f}s"
        return

    if not rows:
        state.global_kpi_molecules = "0"
        state.global_kpi_comparisons = "0"
        state.global_kpi_mean = "n/a"
        state.global_kpi_worst = "n/a"
        state.global_trend_src = ""
        state.global_trend_error = "No comparisons available for current selections."
        state.global_breakdown_rows = []
        state.global_cluster_rows = []
        state.dashboard_assumptions = (
            "Assumptions: global view uses the selected omega only; clustering is heuristic k-means on "
            "per-molecule basis-metric vectors with missing entries imputed by basis-wise means."
        )
        state.global_last_runtime = f"{(time.perf_counter() - t0):.2f}s"
        return

    rows, cluster_rows = _assign_molecule_clusters(
        rows,
        cluster_k=max(1, min(int(state.global_cluster_k), _GLOBAL_CLUSTER_MAX)),
    )
    vals = np.asarray([float(r["metric_value"]) for r in rows], dtype=float)
    vals = vals[np.isfinite(vals)]
    molecules = sorted(set(str(r["molecule"]) for r in rows))

    state.global_kpi_molecules = str(len(molecules))
    state.global_kpi_comparisons = str(len(rows))
    state.global_kpi_mean = f"{float(np.mean(vals)):.4g}" if vals.size else "n/a"

    worst = max(rows, key=lambda r: _metric_rank_value(metric_key, float(r["metric_value"])))
    state.global_kpi_worst = f"{worst['molecule']} / {worst['basis']} ({float(worst['metric_value']):.4g})"

    try:
        state.global_trend_src = _global_trend_plot_data_url(rows, metric_key=metric_key)
        state.global_trend_error = "" if state.global_trend_src else "Trend plot unavailable (missing matplotlib?)"
    except Exception as exc:
        state.global_trend_src = ""
        state.global_trend_error = f"Trend plot error: {type(exc).__name__}: {exc}"

    ranked = sorted(
        rows,
        key=lambda r: _metric_rank_value(metric_key, float(r["metric_value"])),
        reverse=True,
    )
    state.global_breakdown_rows = ranked[:_GLOBAL_BREAKDOWN_MAX_ROWS]
    state.global_cluster_rows = cluster_rows
    state.dashboard_assumptions = (
        "Assumptions: global metric rows use selected omega and reference basis; clustering is a simple "
        "k-means grouping on per-molecule basis-metric vectors (easy to replace with domain-specific models)."
    )
    state.global_last_runtime = f"{(time.perf_counter() - t0):.2f}s"


def _maybe_update_global_dashboard(*, force: bool = False) -> None:
    if force or bool(state.global_auto_update):
        _update_global_dashboard()
        return
    state.global_stale_note = "Global dashboard paused for responsiveness. Click Refresh Global."


def _refresh_global_dashboard() -> None:
    _update_global_dashboard()


ctrl.refresh_global_dashboard = _refresh_global_dashboard


def _load_dataset_db() -> None:
    global _SHG_DB_DIR
    global _SHG_CSV_PATH
    global _BUNDLE_DIR
    global _GEOM_MAP_PATH

    raw = str(state.dataset_db_dir).strip()
    if not raw:
        state.dataset_status = "Enter a database directory path."
        state.dataset_status_color = "red"
        return
    db_dir = Path(raw).expanduser().resolve()
    if not db_dir.exists() or not db_dir.is_dir():
        state.dataset_status = f"Invalid directory: {db_dir}"
        state.dataset_status_color = "red"
        return

    try:
        _SHG_DB_DIR = db_dir
        _SHG_CSV_PATH = None
        _BUNDLE_DIR = None
        _GEOM_MAP_PATH = None
        _clear_cached_data_views()
        _refresh_selector_items()
        _ensure_basis_selections()
        state.active_data_source = _active_data_source_label()
        state.dataset_status = f"Loaded DB: {db_dir}"
        state.dataset_status_color = "green"
        _rebuild_both()
        _maybe_update_global_dashboard(force=False)
    except Exception as exc:
        state.dataset_status = f"Load failed: {type(exc).__name__}: {exc}"
        state.dataset_status_color = "red"


ctrl.load_dataset_db = _load_dataset_db


def _ensure_basis_selections() -> None:
    bases = _basis_list()
    if not bases:
        return

    if str(state.ref_basis) not in bases:
        state.ref_basis = str(bases[0])

    if str(state.bas_basis_a) not in bases:
        preferred = "aug-cc-pVDZ"
        state.bas_basis_a = preferred if preferred in bases else str(bases[0])

    if str(state.bas_basis_b) not in bases:
        preferred = "d-aug-cc-pVDZ"
        if preferred in bases:
            state.bas_basis_b = preferred
        else:
            fallback = "aug-cc-pVDZ"
            state.bas_basis_b = fallback if fallback in bases else str(bases[0])


def _apply_view_layout() -> None:
    layout = str(state.view_layout)
    if layout == "reference_only":
        _scene_left["renderer"].SetViewport(0.0, 0.0, 1.0, 1.0)
        _scene_mid["renderer"].SetViewport(0.0, 0.0, 0.0, 0.0)
        _scene_right["renderer"].SetViewport(0.0, 0.0, 0.0, 0.0)
        _scene_mid["scalar_bar"].VisibilityOff()
        _scene_right["scalar_bar"].VisibilityOff()
    elif layout == "basis_only":
        _scene_left["renderer"].SetViewport(0.0, 0.0, 0.0, 0.0)
        _scene_left["scalar_bar"].VisibilityOff()
        _scene_mid["renderer"].SetViewport(0.0, 0.0, 1.0, 1.0)
        _scene_right["renderer"].SetViewport(0.0, 0.0, 0.0, 0.0)
        _scene_right["scalar_bar"].VisibilityOff()
    elif layout == "compare_bases":
        _scene_left["renderer"].SetViewport(0.0, 0.0, 1.0 / 3.0, 1.0)
        _scene_mid["renderer"].SetViewport(1.0 / 3.0, 0.0, 2.0 / 3.0, 1.0)
        _scene_right["renderer"].SetViewport(2.0 / 3.0, 0.0, 1.0, 1.0)
    else:
        # side_by_side: reference + basis A
        _scene_left["renderer"].SetViewport(0.0, 0.0, 0.5, 1.0)
        _scene_mid["renderer"].SetViewport(0.5, 0.0, 1.0, 1.0)
        _scene_right["renderer"].SetViewport(0.0, 0.0, 0.0, 0.0)
        _scene_right["scalar_bar"].VisibilityOff()


@state.change("view_layout")
def _on_layout_change(**_):
    _apply_view_layout()
    _update_corner_annotations()
    ctrl.view_update()


def _clear_scene(scene: Dict[str, Any]) -> None:
    renderer = scene["renderer"]
    if scene.get("glyph_actor") is not None:
        renderer.RemoveActor(scene["glyph_actor"])
        scene["glyph_actor"] = None
    for a in scene.get("mol_actors", []):
        renderer.RemoveActor(a)
    scene["mol_actors"] = []


def _set_scalar_bar(
    scene: Dict[str, Any],
    *,
    actor: vtkActor,
    title: str,
    mode: str,
    clim: Tuple[float, float],
) -> None:
    scalar_bar = scene["scalar_bar"]
    glyph_mapper = actor.GetMapper()
    lut = glyph_mapper.GetLookupTable()
    scalar_bar.SetLookupTable(lut)

    scalar_bar.SetTitle(str(title))

    # Standard scalar bar behavior.
    scalar_bar.SetUseCustomLabels(False)
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetLabelFormat("%-#6.3g")

    scalar_bar.VisibilityOn()


def _coord_label() -> str:
    cs = str(state.coord_system)
    if cs == "principal_axes":
        return "principal_axes"
    return "as_is"


def _layout_label() -> str:
    layout = str(state.view_layout)
    if layout == "reference_only":
        return "ref_only"
    if layout == "basis_only":
        return "basis_only"
    if layout == "compare_bases":
        return "compare"
    return "side_by_side"


def _update_corner_annotations() -> None:
    # Titles: upper-left. CornerAnnotation is reliably supported in VTK.js.
    _scene_left["corner"].SetText(2, f"Reference ({state.ref_basis})")
    _scene_mid["corner"].SetText(2, f"Basis A ({state.bas_basis_a})")
    _scene_right["corner"].SetText(2, f"Basis B ({state.bas_basis_b})")

    # Context: compact upper-right.
    base = f"Coord: {_coord_label()}  |  View: {_layout_label()}  |  RelNorm: {_rel_norm_label()}"
    _scene_left["corner"].SetText(3, base)
    _scene_mid["corner"].SetText(3, base)
    _scene_right["corner"].SetText(3, base)

    # Best-effort: also set TextActor titles (may not render in all VTK.js builds).
    if _scene_left.get("title_actor") is not None:
        _scene_left["title_actor"].SetInput(f"Reference ({state.ref_basis})")
    if _scene_mid.get("title_actor") is not None:
        _scene_mid["title_actor"].SetInput(f"Basis A ({state.bas_basis_a})")
    if _scene_right.get("title_actor") is not None:
        _scene_right["title_actor"].SetInput(f"Basis B ({state.bas_basis_b})")


def _compute_bundle():
    grid = load_lebedev_grid(int(state.lebedev_order))
    df = _data()

    beta_ref = _tensor_from_long(df, state.mol, state.ref_basis, float(state.omega))
    beta_a = _tensor_from_long(df, state.mol, state.bas_basis_a, float(state.omega))
    beta_b = _tensor_from_long(df, state.mol, state.bas_basis_b, float(state.omega))

    R = None
    if state.coord_system == "principal_axes":
        mol = _load_molecule(state.mol)
        R = _principal_axes_rotation(mol)
        beta_ref = _rotate_beta(beta_ref, R)
        beta_a = _rotate_beta(beta_a, R)
        beta_b = _rotate_beta(beta_b, R)

    # Support both: LebedevGrid with attributes (.n_hat/.w) and legacy tuple (n_hat, w)
    if isinstance(grid, tuple) and len(grid) == 2:
        n_hat = np.asarray(grid[0], dtype=float)
        w = np.asarray(grid[1], dtype=float)
    else:
        n_hat = np.asarray(grid.n_hat, dtype=float)
        w = np.asarray(grid.w, dtype=float)

    v_ref = evaluate_field(beta_ref, n_hat)
    v_a = evaluate_field(beta_a, n_hat)
    v_b = evaluate_field(beta_b, n_hat)

    settings = ErrorSettings(
        enable_mask=True,
        mask_threshold_fraction=0.01,
        component_metrics=True,
        radial_tangential_metrics=True,
        rel_norm=str(state.rel_norm),
    )

    arrays_a, metrics_a = compute_error_fields(v_ref, v_a, n_hat, w, settings=settings)
    arrays_b, metrics_b = compute_error_fields(v_ref, v_b, n_hat, w, settings=settings)
    metrics_a = {**metrics_a, **_derived_error_metrics(arrays_a, w)}
    metrics_b = {**metrics_b, **_derived_error_metrics(arrays_b, w)}

    return n_hat, w, v_ref, v_a, v_b, arrays_a, arrays_b, metrics_a, metrics_b, R


def _rebuild_scene(
    scene: Dict[str, Any],
    *,
    n_hat: np.ndarray,
    vectors: np.ndarray,
    scalars: np.ndarray,
    scale: np.ndarray,
    clim: Tuple[float, float],
    scalar_title: str,
    scalar_mode: str,
    lut_kind: str,
    R: np.ndarray | None,
 ) -> None:
    renderer = scene["renderer"]
    renderWindow = scene["renderWindow"]
    _clear_scene(scene)

    vectors = np.nan_to_num(vectors, nan=0.0, posinf=0.0, neginf=0.0)
    scalars = np.nan_to_num(scalars, nan=0.0, posinf=0.0, neginf=0.0)
    scale = np.nan_to_num(scale, nan=0.0, posinf=0.0, neginf=0.0)

    glyph_actor = _build_glyph_actor(
        n_hat=n_hat,
        vectors=vectors,
        scalars=scalars,
        scale=scale,
        glyph_factor=float(state.glyph_factor),
        clim=clim,
        lut_kind=lut_kind,
        scalar_mode=str(scalar_mode),
    )
    renderer.AddActor(glyph_actor)
    scene["glyph_actor"] = glyph_actor
    _set_scalar_bar(scene, actor=glyph_actor, title=scalar_title, mode=str(scalar_mode), clim=clim)

    mol = _load_molecule(state.mol)
    if mol is None:
        return

    symbols = list(mol.symbols)
    coords = np.asarray(mol.geometry, dtype=float)
    coords = coords - np.mean(coords, axis=0, keepdims=True)
    if R is not None:
        coords = coords @ R
    max_norm = float(np.max(np.linalg.norm(coords, axis=1))) if coords.size else 1.0
    if max_norm > 0:
        coords = coords / max_norm
    coords = coords * float(state.molecule_scale)

    for sym in sorted(set(symbols)):
        idx = [i for i, s in enumerate(symbols) if s == sym]
        pts = coords[idx]
        actor = _build_atom_actor(pts, radius=float(state.atom_radius), rgb=_element_rgb(sym))
        renderer.AddActor(actor)
        scene["mol_actors"].append(actor)

    # With VTK.js client-side rendering, we do not need server-side OpenGL renders.
    # Avoid calling Render() to prevent headless/OpenGL stalls.


def _rebuild_both() -> None:
    import traceback
    global _camera_initialized

    try:
        print("[trame] rebuilding scenes ...", flush=True)
        state.status = "building scene"
        state.status_color = "orange"
        state.last_error = ""

        # Apply the current layout (important for initial load).
        _apply_view_layout()

        n_hat, w, v_ref, v_a, v_b, arrays_a, arrays_b, metrics_a, metrics_b, R = _compute_bundle()
        _update_corner_annotations()

        # LEFT: fixed to reference magnitude, colored by |v_ref|
        left_vectors = v_ref
        left_scalars = arrays_a["ref_mag"]
        left_scale = arrays_a["ref_mag"]
        left_clim = _auto_clim(left_scalars, mode="ref_mag")
        _rebuild_scene(
            _scene_left,
            n_hat=n_hat,
            vectors=left_vectors,
            scalars=left_scalars,
            scale=left_scale,
            clim=left_clim,
            scalar_title=_scalar_mode_label("ref_mag"),
            scalar_mode="ref_mag",
            lut_kind="sequential",
            R=R,
        )

        def _vectors_and_scale(field_mode: str, *, v_bas: np.ndarray, arrays: Dict[str, np.ndarray]):
            if field_mode == "ref":
                return v_ref, arrays["ref_mag"]
            if field_mode == "bas":
                return v_bas, arrays["bas_mag"]
            return arrays["dv"], arrays["abs_err"]

        mode = str(state.scalar_mode)
        if mode not in arrays_a:
            raise KeyError(
                f"Unknown scalar_mode={mode}. Available: {sorted(arrays_a.keys())}"
            )
        style = _metric_style(mode)
        lut_kind = "diverging" if style["signed"] else "sequential"

        global_clims = _global_metric_clims(
            str(state.mol),
            float(state.omega),
            str(state.ref_basis),
            int(state.lebedev_order),
            str(state.coord_system),
            str(state.rel_norm),
        )

        # BASIS A scene
        a_vectors, a_scale = _vectors_and_scale(str(state.field_mode), v_bas=v_a, arrays=arrays_a)
        a_scalars = arrays_a[mode]
        a_clim = global_clims.get(mode, _auto_clim(a_scalars, mode=mode))
        if bool(state.scale_by_metric):
            a_scale = np.abs(a_scalars) if style["signed"] else a_scalars
        _rebuild_scene(
            _scene_mid,
            n_hat=n_hat,
            vectors=a_vectors,
            scalars=a_scalars,
            scale=a_scale,
            clim=a_clim,
            scalar_title=_scalar_mode_label(mode),
            scalar_mode=mode,
            lut_kind=lut_kind,
            R=R,
        )

        # BASIS B scene
        b_vectors, b_scale = _vectors_and_scale(str(state.field_mode), v_bas=v_b, arrays=arrays_b)
        b_scalars = arrays_b[mode]
        b_clim = global_clims.get(mode, _auto_clim(b_scalars, mode=mode))
        if bool(state.scale_by_metric):
            b_scale = np.abs(b_scalars) if style["signed"] else b_scalars
        _rebuild_scene(
            _scene_right,
            n_hat=n_hat,
            vectors=b_vectors,
            scalars=b_scalars,
            scale=b_scale,
            clim=b_clim,
            scalar_title=_scalar_mode_label(mode),
            scalar_mode=mode,
            lut_kind=lut_kind,
            R=R,
        )

        # One shared camera for all renderers.
        if not _camera_initialized:
            _scene_left["renderer"].ResetCamera()
            _camera_initialized = True

        metric_keys = (
            "rel_L2",
            "L2_par",
            "L2_perp",
            "ratio_perp_par",
            "max_rel",
            "p95_rel",
            "mean_ang",
            "p95_ang",
            "bias_ratio",
            "masked_rel_L2",
        )
        state.metric_rows = [
            {"name": k, "a": metrics_a.get(k), "b": metrics_b.get(k)} for k in metric_keys
        ]
        state.metric_headers = [
            {"text": "Metric", "value": "name"},
            {"text": f"A: {state.bas_basis_a}", "value": "a"},
            {"text": f"B: {state.bas_basis_b}", "value": "b"},
        ]

        _update_metric_plot()

        print("[trame] scenes rebuilt", flush=True)
        state.status = "ready"
        state.status_color = "green"

        ctrl.view_update()

    except Exception as exc:
        tb = traceback.format_exc()
        print("[trame] rebuild failed:\n" + tb, flush=True)
        state.last_error = f"{type(exc).__name__}: {exc}"
        state.status = "error"
        state.status_color = "red"
        try:
            _scene_left["scalar_bar"].VisibilityOff()
            _scene_mid["scalar_bar"].VisibilityOff()
            _scene_right["scalar_bar"].VisibilityOff()
            ctrl.view_update()
        except Exception:
            pass


@state.change("mol")
def _on_mol_change(mol, **_):
    global _camera_initialized
    omegas = _omega_list(mol)
    state.omega_items = [float(x) for x in omegas]
    if omegas and float(state.omega) not in {float(x) for x in omegas}:
        state.omega = float(omegas[0])
    _camera_initialized = False
    _rebuild_both()


@state.change(
    "omega",
    "ref_basis",
    "bas_basis_a",
    "bas_basis_b",
    "lebedev_order",
    "field_mode",
    "scalar_mode",
    "glyph_factor",
    "molecule_scale",
    "atom_radius",
    "coord_system",
    "rel_norm",
    "scale_by_metric",
)
def _on_params_change(**_):
    _rebuild_both()


@state.change(
    "omega",
    "ref_basis",
    "lebedev_order",
    "coord_system",
    "rel_norm",
    "global_metric",
    "global_cluster_k",
)
def _on_global_dashboard_params_change(**_):
    _maybe_update_global_dashboard(force=False)


@state.change("global_auto_update")
def _on_global_auto_update(global_auto_update, **_):
    if bool(global_auto_update):
        _update_global_dashboard()


@state.change("posneg_yscale")
def _on_posneg_yscale_change(**_):
    _update_metric_plot()
    ctrl.view_update()


def standard_buttons():
    vuetify.VCheckbox(
        v_model="$vuetify.theme.dark",
        on_icon="mdi-lightbulb-off-outline",
        off_icon="mdi-lightbulb-outline",
        classes="mx-1",
        hide_details=True,
        dense=True,
    )
    with vuetify.VBtn(icon=True, click=ctrl.reset_camera):
        vuetify.VIcon("mdi-crop-free")


def _metrics_table():
    vuetify.VDataTable(
        dense=True,
        disable_pagination=True,
        hide_default_footer=True,
        items=("metric_rows", []),
        headers=("metric_headers", []),
    )


def _global_breakdown_table():
    vuetify.VDataTable(
        dense=True,
        disable_pagination=False,
        hide_default_footer=False,
        items=("global_breakdown_rows", []),
        headers=("global_breakdown_headers", []),
        items_per_page=12,
        classes="mb-2",
    )


def _global_cluster_table():
    vuetify.VDataTable(
        dense=True,
        disable_pagination=True,
        hide_default_footer=True,
        items=("global_cluster_rows", []),
        headers=("global_cluster_headers", []),
        classes="mb-2",
    )


with SinglePageWithDrawerLayout(server) as layout:
    layout.title.set_text("β unit-sphere viewer")

    with layout.toolbar:
        vuetify.VSpacer()
        vuetify.VDivider(vertical=True, classes="mx-2")
        vuetify.VChip(
            "{{ status }}",
            small=True,
            outlined=True,
            color=("status_color", "orange"),
            classes="mr-2",
        )
        standard_buttons()

    with layout.drawer as drawer:
        drawer.width = 360
        vuetify.VContainer(fluid=True, classes="pa-2")

        vuetify.VAlert(
            "{{ last_error }}",
            type="error",
            dense=True,
            outlined=True,
            prominent=False,
            v_if="last_error",
            classes="mb-2",
        )

        vuetify.VSubheader("Dataset")
        vuetify.VTextField(
            label="Database directory",
            v_model=("dataset_db_dir", state.dataset_db_dir),
            dense=True,
            outlined=True,
            hide_details=True,
            placeholder="/path/to/database/root",
            classes="mb-1",
        )
        with vuetify.VRow(no_gutters=True, classes="mb-2"):
            with vuetify.VCol(cols=7, classes="pr-1"):
                vuetify.VChip(
                    "{{ dataset_status }}",
                    small=True,
                    outlined=True,
                    color=("dataset_status_color", "blue"),
                )
            with vuetify.VCol(cols=5, classes="pl-1"):
                vuetify.VBtn(
                    "Load DB",
                    small=True,
                    block=True,
                    color="primary",
                    click=ctrl.load_dataset_db,
                )
        vuetify.VTextField(
            label="Active source",
            v_model=("active_data_source", state.active_data_source),
            dense=True,
            outlined=True,
            hide_details=True,
            readonly=True,
            classes="mb-2",
        )
        vuetify.VDivider(classes="my-2")
        vuetify.VSubheader("Sphere Controls")

        vuetify.VSelect(
            label="Molecule",
            v_model=("mol", state.mol),
            items=("molecule_items", state.molecule_items),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Coordinate system",
            v_model=("coord_system", state.coord_system),
            items=(
                "coord_items",
                [
                    {"text": "As-is (tensor frame)", "value": "as_is"},
                    {"text": "Principal axes (inertia)", "value": "principal_axes"},
                ],
            ),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VCheckbox(
            label="Show axes",
            v_model=("show_axes", state.show_axes),
            dense=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="View",
            v_model=("view_layout", state.view_layout),
            items=(
                "view_items",
                [
                    {"text": "Side-by-side", "value": "side_by_side"},
                    {"text": "Compare bases (A vs B)", "value": "compare_bases"},
                    {"text": "Reference only", "value": "reference_only"},
                    {"text": "Basis/errors only", "value": "basis_only"},
                ],
            ),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Omega",
            v_model=("omega", state.omega),
            items=("omega_items", state.omega_items),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Reference basis",
            v_model=("ref_basis", state.ref_basis),
            items=("basis_items", state.basis_items),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Basis A",
            v_model=("bas_basis_a", state.bas_basis_a),
            items=("basis_items", state.basis_items),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Basis B",
            v_model=("bas_basis_b", state.bas_basis_b),
            items=("basis_items", state.basis_items),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Lebedev order",
            v_model=("lebedev_order", state.lebedev_order),
            items=("lebedev_items", [3, 5, 7, 9, 11, 13, 15, 17,
                19, 21, 23, 25, 27, 29, 31, 35,
                41, 47, 53, 59, 65, 71, 77, 83,
                89, 95, 101, 107, 113, 119, 125, 131]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Vector field",
            v_model=("field_mode", state.field_mode),
            items=("field_items", [{"text": t, "value": v} for t, v in FIELD_CHOICES]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSelect(
            label="Color by",
            v_model=("scalar_mode", state.scalar_mode),
            items=("metric_items", [{"text": t, "value": v} for t, v in METRIC_CHOICES]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSelect(
            label="Relative normalization",
            v_model=("rel_norm", state.rel_norm),
            items=(
                "rel_norm_items",
                [
                    {"text": "Global RMS (recommended)", "value": "global_rms"},
                    {"text": "Global max", "value": "global_max"},
                    {"text": "Global mean", "value": "global_mean"},
                    {"text": "Pointwise |v_ref(n)| (old)", "value": "pointwise"},
                ],
            ),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VDivider(classes="my-2")
        vuetify.VSubheader("Global Dashboard")
        vuetify.VSelect(
            label="Global metric",
            v_model=("global_metric", state.global_metric),
            items=("global_metric_items", [{"text": t, "value": v} for t, v in GLOBAL_DASH_METRIC_CHOICES]),
            dense=True,
            outlined=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSlider(
            label="Molecule clusters (heuristic)",
            v_model=("global_cluster_k", state.global_cluster_k),
            min=1,
            max=_GLOBAL_CLUSTER_MAX,
            step=1,
            dense=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VCheckbox(
            label="Auto-update global dashboard (slower)",
            v_model=("global_auto_update", state.global_auto_update),
            dense=True,
            hide_details=True,
            classes="mb-1",
        )
        with vuetify.VRow(no_gutters=True, classes="mb-2"):
            with vuetify.VCol(cols=7, classes="pr-1"):
                vuetify.VTextField(
                    label="Global runtime",
                    v_model=("global_last_runtime", state.global_last_runtime),
                    dense=True,
                    outlined=True,
                    hide_details=True,
                    readonly=True,
                )
            with vuetify.VCol(cols=5, classes="pl-1"):
                vuetify.VBtn(
                    "Refresh",
                    small=True,
                    block=True,
                    color="primary",
                    click=ctrl.refresh_global_dashboard,
                )

        vuetify.VCheckbox(
            label="Scale arrows by color metric",
            v_model=("scale_by_metric", state.scale_by_metric),
            dense=True,
            hide_details=True,
            classes="mb-2",
        )

        vuetify.VSlider(
            label="Arrow scale",
            v_model=("glyph_factor", state.glyph_factor),
            min=0.01,
            max=0.5,
            step=0.01,
            dense=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSlider(
            label="Molecule scale",
            v_model=("molecule_scale", state.molecule_scale),
            min=0.2,
            max=1.5,
            step=0.05,
            dense=True,
            hide_details=True,
            classes="mb-2",
        )
        vuetify.VSlider(
            label="Atom radius",
            v_model=("atom_radius", state.atom_radius),
            min=0.05,
            max=0.6,
            step=0.01,
            dense=True,
            hide_details=True,
            classes="mb-4",
        )

        vuetify.VDivider(classes="my-2")
        vuetify.VSubheader("Selected Molecule Metrics")
        _metrics_table()

    with layout.content:
        with vuetify.VContainer(fluid=True, classes="pa-0 fill-height"):
            with vuetify.VRow(no_gutters=True, classes="fill-height"):
                with vuetify.VCol(cols=8, classes="pa-0 fill-height"):
                    # Use client-side VTK.js rendering to avoid server OpenGL/EGL issues.
                    # Left/right are renderers in a shared renderWindow, with a shared camera.
                    view = vtk.VtkLocalView(
                        renderWindow,
                        ref="view",
                        style="width: 100%; height: 100%;",
                    )
                    ctrl.view_update = view.update

                with vuetify.VCol(cols=4, classes="pa-2", style="height: 100%; overflow-y: auto;"):
                    vuetify.VSubheader("Global Dashboard")
                    vuetify.VAlert(
                        "{{ dashboard_assumptions }}",
                        type="info",
                        dense=True,
                        outlined=True,
                        prominent=False,
                        classes="mb-2",
                    )
                    vuetify.VAlert(
                        "{{ global_stale_note }}",
                        type="warning",
                        dense=True,
                        outlined=True,
                        prominent=False,
                        v_if="global_stale_note",
                        classes="mb-2",
                    )
                    with vuetify.VRow(no_gutters=True, classes="mb-2"):
                        with vuetify.VCol(cols=6, classes="pr-1 pb-1"):
                            with vuetify.VCard(outlined=True):
                                vuetify.VCardText(
                                    "Molecules: {{ global_kpi_molecules }}",
                                    classes="py-2",
                                )
                        with vuetify.VCol(cols=6, classes="pl-1 pb-1"):
                            with vuetify.VCard(outlined=True):
                                vuetify.VCardText(
                                    "Comparisons: {{ global_kpi_comparisons }}",
                                    classes="py-2",
                                )
                        with vuetify.VCol(cols=6, classes="pr-1 pt-1"):
                            with vuetify.VCard(outlined=True):
                                vuetify.VCardText(
                                    "Mean metric: {{ global_kpi_mean }}",
                                    classes="py-2",
                                )
                        with vuetify.VCol(cols=6, classes="pl-1 pt-1"):
                            with vuetify.VCard(outlined=True):
                                vuetify.VCardText(
                                    "Worst: {{ global_kpi_worst }}",
                                    classes="py-2",
                                )
                    vuetify.VAlert(
                        "{{ global_trend_error }}",
                        type="warning",
                        dense=True,
                        outlined=True,
                        prominent=False,
                        v_if="global_trend_error",
                        classes="mb-2",
                    )
                    vuetify.VImg(
                        src=("global_trend_src", ""),
                        contain=True,
                        v_if="global_trend_src",
                        classes="mb-2",
                        style="width: 100%;",
                    )
                    _global_breakdown_table()
                    vuetify.VSubheader("Cluster Summary")
                    _global_cluster_table()
                    vuetify.VDivider(classes="my-2")
                    vuetify.VSubheader("Selected Molecule Trend")
                    vuetify.VSelect(
                        label="Signed pos/neg y-scale",
                        v_model=("posneg_yscale", state.posneg_yscale),
                        items=(
                            "posneg_yscale_items",
                            [
                                {"text": "Linear", "value": "linear"},
                                {"text": "SymLog", "value": "symlog"},
                            ],
                        ),
                        dense=True,
                        outlined=True,
                        hide_details=True,
                        classes="mb-2",
                    )
                    vuetify.VAlert(
                        "{{ metric_plot_error }}",
                        type="warning",
                        dense=True,
                        outlined=True,
                        prominent=False,
                        v_if="metric_plot_error",
                        classes="mb-2",
                    )
                    vuetify.VImg(
                        src=("metric_plot_src", ""),
                        contain=True,
                        v_if="metric_plot_src",
                        classes="mb-2",
                        style="width: 100%;",
                    )

            # Now that we have live views, build the initial scenes.
            _rebuild_both()
            _maybe_update_global_dashboard(force=False)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host interface to bind (default: 127.0.0.1; use 0.0.0.0 for remote access)",
    )
    parser.add_argument(
        "--port",
        default=9010,
        type=int,
        help="Port to bind (default: 9010; use 0 to auto-pick a free port)",
    )
    parser.add_argument(
        "--server",
        action="store_true",
        help="Compatibility flag for headless runs; browser auto-open is disabled by default.",
    )
    parser.add_argument(
        "--mol-file",
        dest="mol_file",
        default=None,
        help="Path to a single .mol file used when geometry is missing",
    )
    parser.add_argument(
        "--mol-dir",
        dest="mol_dir",
        default=None,
        help="Directory of .mol files (matched by molecule label)",
    )
    parser.add_argument(
        "--mol-map",
        dest="mol_map",
        default=None,
        help="JSON map of molecule label -> .mol path",
    )
    parser.add_argument(
        "--shg-csv",
        dest="shg_csv",
        default=None,
        help="Path to shg_ijk.csv (default: data/csv_data/shg_ijk.csv)",
    )
    parser.add_argument(
        "--db-dir",
        dest="db_dir",
        default=None,
        help="Path to a calculation database directory (builds SHG table on the fly)",
    )
    parser.add_argument(
        "--bundle-dir",
        dest="bundle_dir",
        default=None,
        help="Directory containing shg_ijk.csv and optional geometries.json",
    )
    parser.add_argument(
        "--geom-json",
        dest="geom_json",
        default=None,
        help="Path to a geometry map JSON (label -> {symbols, geometry})",
    )
    parser.add_argument(
        "--write-bundle",
        dest="write_bundle",
        default=None,
        help="If set with --db-dir, writes a bundle directory (CSV + geometries.json)",
    )
    return parser.parse_args(argv)


def _select_free_port(host: str, preferred_port: int, *, max_tries: int = 50) -> int:
    import errno
    import socket

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
    global _MOL_RESOLVER
    global _SHG_CSV_PATH
    global _SHG_DB_DIR
    global _GEOM_MAP_PATH
    global _BUNDLE_DIR
    global _WRITE_BUNDLE_DIR

    args = _parse_args(argv)
    if args.bundle_dir:
        _BUNDLE_DIR = Path(args.bundle_dir).expanduser().resolve()
        _SHG_CSV_PATH = _BUNDLE_DIR / "shg_ijk.csv"
    elif args.db_dir:
        _SHG_DB_DIR = Path(args.db_dir).expanduser().resolve()
    elif args.shg_csv:
        _SHG_CSV_PATH = Path(args.shg_csv).expanduser().resolve()

    if args.geom_json:
        _GEOM_MAP_PATH = Path(args.geom_json).expanduser().resolve()

    if args.write_bundle:
        _WRITE_BUNDLE_DIR = Path(args.write_bundle).expanduser().resolve()
    _MOL_RESOLVER = None

    # Re-bind dataset after CLI args are parsed (module import may have populated caches).
    _clear_cached_data_views()
    _refresh_selector_items()
    _ensure_basis_selections()
    state.active_data_source = _active_data_source_label()
    if _SHG_DB_DIR is not None:
        state.dataset_db_dir = str(_SHG_DB_DIR)
    _rebuild_both()
    _maybe_update_global_dashboard(force=False)

    host = str(args.host)
    port = _select_free_port(host, int(args.port))

    url_host = host
    if url_host in ("0.0.0.0", "::"):
        url_host = "127.0.0.1"
    print(f"Trame server running at http://{url_host}:{port}/")
    try:
        server.start(host=host, port=port, open_browser=False)
    except TypeError:
        # Backward compatibility with older Trame versions.
        server.start(host=host, port=port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
