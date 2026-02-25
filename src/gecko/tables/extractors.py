from __future__ import annotations

from typing import Any

import numpy as np

from gecko.core.model import Calculation


def make_envelope(calc: Calculation) -> dict[str, Any]:
    mol_id = calc.meta.get("mol_id") or calc.meta.get("molecule_id")
    return {
        "calc_id": calc.meta.get("calc_id"),
        "geom_id": calc.meta.get("geom_id"),
        "mol_id": mol_id,
        "molecule_id": calc.meta.get("molecule_id"),
        "label": calc.meta.get("label"),
        "code": calc.code,
        "root": str(calc.root),
        "basis": getattr(calc, "basis", None) or calc.meta.get("basis"),
        "method": calc.meta.get("method"),
    }


def _require_geom(env: dict[str, Any]) -> bool:
    return env.get("geom_id") is not None


def extract_beta(calc: Calculation) -> list[dict[str, Any]]:
    beta = calc.data.get("beta") or {}
    if not beta:
        return []
    if not all(k in beta for k in ("omega", "components", "values")):
        return []

    env = make_envelope(calc)

    omega = np.asarray(beta["omega"], dtype=float)
    comps = list(beta["components"])
    vals = np.asarray(beta["values"], dtype=float)

    rows: list[dict[str, Any]] = []
    for i in range(vals.shape[0]):
        omegaA, omegaB, omegaC = map(float, omega[i])
        for j, ijk in enumerate(comps):
            rows.append(
                {
                    **env,
                    "omegaA": omegaA,
                    "omegaB": omegaB,
                    "omegaC": omegaC,
                    "ijk": str(ijk).lower(),
                    "value": float(vals[i, j]),
                }
            )
    return rows


def extract_alpha(calc: Calculation) -> list[dict[str, Any]]:
    alpha = calc.data.get("alpha")
    if not alpha:
        return []

    env = make_envelope(calc)

    omega = np.asarray(alpha.get("omega", []), dtype=float).reshape(-1)
    comps = list(alpha.get("components", []))
    vals = np.asarray(alpha.get("values", []), dtype=float)

    rows: list[dict[str, Any]] = []
    for i, om in enumerate(omega):
        for j, ij in enumerate(comps):
            rows.append({**env, "omega": float(om), "ij": str(ij).lower(), "value": float(vals[i, j])})
    return rows


def _find_task(raw_json: dict[str, Any], task_type: str) -> dict[str, Any] | None:
    tasks = raw_json.get("tasks", [])
    if isinstance(tasks, list):
        for t in tasks:
            if isinstance(t, dict) and t.get("type") == task_type:
                return t
    return None


def extract_energy(calc: Calculation) -> list[dict[str, Any]]:
    env = make_envelope(calc)

    energy = calc.meta.get("ground_state_energy")
    if energy is None:
        raw = calc.data.get("raw_json")
        if isinstance(raw, dict):
            scf = _find_task(raw, "scf")
            if isinstance(scf, dict):
                energy = scf.get("energy")

    if energy is None:
        return []

    return [{**env, "energy": float(energy)}]


def extract_dipole(calc: Calculation) -> list[dict[str, Any]]:
    env = make_envelope(calc)
    if not _require_geom(env):
        return []

    raw = calc.data.get("raw_json")
    if not isinstance(raw, dict):
        return []

    scf = _find_task(raw, "scf")
    if not isinstance(scf, dict):
        return []

    dip = scf.get("dipole", {})
    vals = dip.get("vals") if isinstance(dip, dict) else None
    if not isinstance(vals, list) or len(vals) < 3:
        return []

    rows = []
    for comp, value in zip(["x", "y", "z"], vals, strict=False):
        rows.append({**env, "i": comp, "value": float(value)})
    return rows


def _alpha_iso_map(alpha: dict[str, Any]) -> dict[float, float]:
    omega = np.asarray(alpha.get("omega", []), dtype=float).reshape(-1)
    comps = list(alpha.get("components", []))
    vals = np.asarray(alpha.get("values", []), dtype=float)
    if omega.size == 0 or vals.size == 0 or not comps:
        return {}

    try:
        idx_xx = comps.index("xx")
        idx_yy = comps.index("yy")
        idx_zz = comps.index("zz")
    except ValueError:
        return {}

    iso = (vals[:, idx_xx] + vals[:, idx_yy] + vals[:, idx_zz]) / 3.0
    return {float(omega[i]): float(iso[i]) for i in range(len(omega))}


def _match_omega(target: float, omega_vals: list[float], *, tol: float = 1e-8) -> float | None:
    if not omega_vals:
        return None
    diffs = np.abs(np.asarray(omega_vals, dtype=float) - float(target))
    idx = int(np.argmin(diffs))
    if diffs[idx] <= tol:
        return float(omega_vals[idx])
    return None


def _iso_derivatives_by_mode(
    pol_freqs: list[float],
    derivs_by_mode: Any,
) -> dict[float, dict[int, float]]:
    if derivs_by_mode is None:
        return {}

    items: list[tuple[float, np.ndarray]] = []
    if isinstance(derivs_by_mode, dict):
        for k, v in derivs_by_mode.items():
            items.append((float(k), np.asarray(v, dtype=float)))
        items.sort(key=lambda x: x[0])
    elif isinstance(derivs_by_mode, list):
        for i, arr in enumerate(derivs_by_mode):
            if i >= len(pol_freqs):
                break
            items.append((float(pol_freqs[i]), np.asarray(arr, dtype=float)))
    else:
        return {}

    iso_map: dict[float, dict[int, float]] = {}
    for freq, arr in items:
        if arr.ndim == 3 and arr.shape[0] == 3 and arr.shape[1] == 3:
            iso = (arr[0, 0, :] + arr[1, 1, :] + arr[2, 2, :]) / 3.0
        elif arr.ndim == 2 and arr.shape[0] == 9:
            iso = (arr[0, :] + arr[4, :] + arr[8, :]) / 3.0
        elif arr.ndim == 2 and arr.shape[1] == 9:
            iso = (arr[:, 0] + arr[:, 4] + arr[:, 8]) / 3.0
        else:
            continue

        iso_map[freq] = {int(i + 1): float(val) for i, val in enumerate(iso)}

    return iso_map


def extract_raman(calc: Calculation) -> list[dict[str, Any]]:
    raman = calc.data.get("raman")
    if not raman:
        return []

    env = make_envelope(calc)

    raman_by_freq = raman.get("raman_by_freq") or {}
    if not isinstance(raman_by_freq, dict) or not raman_by_freq:
        return []

    pol_freqs = list(np.asarray(raman.get("polarization_frequencies", []), dtype=float).reshape(-1))
    alpha = calc.data.get("alpha") or {}
    alpha_iso = _alpha_iso_map(alpha)
    alpha_keys = sorted(alpha_iso.keys())

    dalpha_iso = _iso_derivatives_by_mode(
        pol_freqs,
        raman.get("polarizability_derivatives_by_mode"),
    )

    rows: list[dict[str, Any]] = []
    for freq_key, items in raman_by_freq.items():
        try:
            omega_pol = float(freq_key)
        except Exception:
            continue

        omega_match = _match_omega(omega_pol, alpha_keys) if alpha_keys else None
        alpha_iso_val = alpha_iso.get(omega_match) if omega_match is not None else None

        for row in items:
            mode = int(row.get("mode"))
            rows.append(
                {
                    **env,
                    "omega_pol": omega_pol,
                    "mode": mode,
                    "freq_cm1": float(row.get("freq_cm1")),
                    "alpha2": float(row.get("alpha2")),
                    "beta2": float(row.get("beta2")),
                    "pol_int": float(row.get("pol_int")),
                    "depol_int": float(row.get("depol_int")),
                    "dep_ratio": float(row.get("dep_ratio")),
                    "alpha_iso": float(alpha_iso_val) if alpha_iso_val is not None else None,
                    "dalpha_iso": dalpha_iso.get(omega_pol, {}).get(mode),
                }
            )

    return rows


def _to_float_or_none(value: Any) -> float | None:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int_or_none(value: Any) -> int | None:
    try:
        if value is None:
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _to_bool_or_none(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        lower = value.strip().lower()
        if lower in ("true", "t", "yes", "y", "1"):
            return True
        if lower in ("false", "f", "no", "n", "0"):
            return False
    return None


def extract_timing_points(calc: Calculation) -> list[dict[str, Any]]:
    timings = calc.data.get("timings")
    if not isinstance(timings, dict):
        return []

    point_rows = timings.get("point_rows")
    if not isinstance(point_rows, list) or not point_rows:
        return []

    env = make_envelope(calc)
    rows: list[dict[str, Any]] = []
    for item in point_rows:
        if not isinstance(item, dict):
            continue
        rows.append(
            {
                **env,
                "timing_kind": str(item["timing_kind"])
                if item.get("timing_kind") is not None
                else None,
                "state_id": str(item["state_id"]) if item.get("state_id") is not None else None,
                "protocol": str(item["protocol"]) if item.get("protocol") is not None else None,
                "frequency": _to_float_or_none(item.get("frequency")),
                "cpu_seconds": _to_float_or_none(item.get("cpu_seconds")),
                "wall_seconds": _to_float_or_none(item.get("wall_seconds")),
                "converged": _to_bool_or_none(item.get("converged")),
                "saved": _to_bool_or_none(item.get("saved")),
                "restart_kind": str(item["restart_kind"])
                if item.get("restart_kind") is not None
                else None,
                "restart_loaded_from_disk": _to_bool_or_none(
                    item.get("restart_loaded_from_disk")
                ),
                "restart_promoted_from_static": _to_bool_or_none(
                    item.get("restart_promoted_from_static")
                ),
                "restart_source_frequency": _to_float_or_none(
                    item.get("restart_source_frequency")
                ),
                "restart_source_protocol": str(item["restart_source_protocol"])
                if item.get("restart_source_protocol") is not None
                else None,
                "derived_state_id": str(item["derived_state_id"])
                if item.get("derived_state_id") is not None
                else None,
                "owner_group": _to_int_or_none(item.get("owner_group")),
                "success": _to_bool_or_none(item.get("success")),
            }
        )
    return rows


def extract_timing_summary(calc: Calculation) -> list[dict[str, Any]]:
    timings = calc.data.get("timings")
    if not isinstance(timings, dict):
        return []
    summary = timings.get("summary")
    if not isinstance(summary, dict) or not summary:
        return []
    return [{**make_envelope(calc), **summary}]
