from __future__ import annotations


from pathlib import Path
from typing import Any
import logging


import numpy as np

from gecko.core.model import Calculation
from gecko.molecule.canonical import canonicalize_atom_order

logger = logging.getLogger(__name__)

_HARTREE_TO_EV = 27.211386245988


def _beta_df_to_tensor(beta_df) -> dict[str, Any]:
    """
    Convert MADNESS df_pivot (pandas DataFrame) into a tensor-first representation.

    Expects:
      - beta_df.index = MultiIndex [omegaA, omegaB, omegaC]
      - beta_df.columns = strings like "xyz", "xxy", etc.
      - beta_df.values = scalar (float or complex)

    Returns:
      {
        "omega": np.ndarray shape (n_freq, num_compents) -> example for beta [omegaA, omegaB, omegaC],
        "components": list[str],
        "values": np.ndarray shape (n_freq, n_comp),
        "shape": ("freq", "component"),
      }
    """
    import pandas as pd  # lazy import

    if beta_df is None:
        return {}

    if not isinstance(beta_df, pd.DataFrame):
        raise TypeError(f"Expected pandas DataFrame, got {type(beta_df)}")

    beta_df = beta_df.sort_index()
    beta_df = beta_df.reindex(sorted(beta_df.columns), axis=1)

    omega = np.asarray(beta_df.index.to_list(), dtype=float)  # (n_freq, 3)
    components = [str(c) for c in beta_df.columns.to_list()]
    values = beta_df.to_numpy(dtype=float)  # (n_freq, n_comp)
    # zero out tiny values for cleaner output
    #values[np.abs(values) < 1e-3] = 0.0

    return {"omega": omega, "components": components, "values": values, "shape": ("freq", "component")}


def _tensor_has_rows(tensor: Any) -> bool:
    if not isinstance(tensor, dict):
        return False
    if not all(k in tensor for k in ("omega", "components", "values")):
        return False
    values = np.asarray(tensor.get("values"))
    components = tensor.get("components") or []
    return values.ndim == 2 and values.shape[0] > 0 and values.shape[1] > 0 and len(components) > 0


def _legacy_alpha_to_tensor(raw_json: dict[str, Any]) -> dict[str, Any]:
    response = raw_json.get("response")
    if not isinstance(response, dict):
        return {}
    alpha = response.get("alpha")
    if not isinstance(alpha, dict):
        return {}

    values = alpha.get("alpha")
    components = alpha.get("ij")
    omega = alpha.get("omega")
    if not (isinstance(values, list) and isinstance(components, list) and isinstance(omega, list)):
        return {}
    if not values or len(values) != len(components) or len(values) != len(omega):
        return {}

    by_freq: dict[float, dict[str, float]] = {}
    for freq_raw, comp_raw, value_raw in zip(omega, components, values, strict=False):
        try:
            freq = float(freq_raw)
            comp = str(comp_raw).strip().lower()
            value = float(value_raw)
        except Exception:
            continue
        if not comp:
            continue
        by_freq.setdefault(freq, {})[comp] = value

    if not by_freq:
        return {}

    freqs = sorted(by_freq.keys())
    comps = sorted({comp for comp_map in by_freq.values() for comp in comp_map.keys()})
    arr = np.full((len(freqs), len(comps)), np.nan, dtype=float)
    for i, freq in enumerate(freqs):
        for j, comp in enumerate(comps):
            if comp in by_freq[freq]:
                arr[i, j] = by_freq[freq][comp]

    return {"omega": np.asarray(freqs, dtype=float), "components": comps, "values": arr, "shape": ("freq", "component")}


def _legacy_beta_to_tensor(raw_json: dict[str, Any]) -> dict[str, Any]:
    hyper = raw_json.get("hyper")
    if not isinstance(hyper, dict):
        return {}
    beta = hyper.get("beta")
    if not isinstance(beta, dict):
        return {}

    a = beta.get("A")
    b = beta.get("B")
    c = beta.get("C")
    a_freq = beta.get("Afreq")
    b_freq = beta.get("Bfreq")
    c_freq = beta.get("Cfreq")
    values = beta.get("Beta")

    if not all(isinstance(x, list) for x in (a, b, c, a_freq, b_freq, c_freq, values)):
        return {}
    n = len(values)
    if n == 0:
        return {}
    if not all(len(x) == n for x in (a, b, c, a_freq, b_freq, c_freq)):
        return {}

    by_freq: dict[tuple[float, float, float], dict[str, float]] = {}
    for a_raw, b_raw, c_raw, af_raw, bf_raw, cf_raw, val_raw in zip(
        a, b, c, a_freq, b_freq, c_freq, values, strict=False
    ):
        try:
            comp = f"{str(a_raw).strip()}{str(b_raw).strip()}{str(c_raw).strip()}".lower()
            freq = (float(af_raw), float(bf_raw), float(cf_raw))
            val = float(val_raw)
        except Exception:
            continue
        if len(comp) != 3:
            continue
        by_freq.setdefault(freq, {})[comp] = val

    if not by_freq:
        return {}

    freqs = sorted(by_freq.keys())
    comps = sorted({comp for comp_map in by_freq.values() for comp in comp_map.keys()})
    arr = np.full((len(freqs), len(comps)), np.nan, dtype=float)
    for i, freq in enumerate(freqs):
        for j, comp in enumerate(comps):
            if comp in by_freq[freq]:
                arr[i, j] = by_freq[freq][comp]

    return {
        "omega": np.asarray(freqs, dtype=float),
        "components": comps,
        "values": arr,
        "shape": ("freq", "component"),
    }


def _format_mra_threshold(prefix: str, value: float) -> str:
    try:
        import math

        if value <= 0:
            return "mra"
        exp = int(round(-math.log10(float(value))))
        if exp < 0:
            return "mra"
        return f"mra-{prefix}{exp:02d}"
    except Exception:
        return "mra"


def _infer_mra_basis_from_obj(obj: Any) -> str | None:
    if not isinstance(obj, (dict, list)):
        return None

    # Depth-first search for known keys.
    if isinstance(obj, dict):
        # Direct dconv
        for k in ("dconv", "converged_for_dconv"):
            v = obj.get(k)
            if isinstance(v, (int, float)) and float(v) > 0:
                return _format_mra_threshold("d", float(v))

        # Protocol (list)
        v = obj.get("protocol")
        if isinstance(v, list) and v:
            last = v[-1]
            if isinstance(last, (int, float)) and float(last) > 0:
                return _format_mra_threshold("p", float(last))

        # Some calc-info payloads store convergence thresholds under different keys.
        for k in ("converged_for_thresh", "thresh"):
            v = obj.get(k)
            if isinstance(v, (int, float)) and float(v) > 0:
                return _format_mra_threshold("p", float(v))

        for v in obj.values():
            out = _infer_mra_basis_from_obj(v)
            if out is not None:
                return out
        return None

    # list
    for it in obj:
        out = _infer_mra_basis_from_obj(it)
        if out is not None:
            return out
    return None


def _infer_mra_basis_from_input_in_text(text: str) -> str | None:
    import re

    # dconv 1e-6 (or dconv=1e-6)
    m = re.search(r"(?im)^\s*dconv\s*(?:=)?\s*([0-9.+-eE]+)\s*$", text)
    if m:
        try:
            return _format_mra_threshold("d", float(m.group(1)))
        except Exception:
            pass

    # protocol [0.0001,1e-06,1e-7]
    m = re.search(r"(?im)^\s*protocol\s*(?:=)?\s*\[(.*?)\]\s*$", text)
    if m:
        raw = m.group(1)
        nums = re.findall(r"[0-9.]+(?:[eE][+-]?[0-9]+)?", raw)
        if nums:
            try:
                return _format_mra_threshold("p", float(nums[-1]))
            except Exception:
                pass

    # protocol 1e-4 1e-6 1e-8 (space-separated)
    m = re.search(r"(?im)^\s*protocol\s+(.+?)\s*$", text)
    if m:
        nums = re.findall(r"[0-9.]+(?:[eE][+-]?[0-9]+)?", m.group(1))
        if nums:
            try:
                return _format_mra_threshold("p", float(nums[-1]))
            except Exception:
                pass

    return None


def _infer_method_from_input_in_text(text: str) -> str | None:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    lower_lines = [ln.lower() for ln in lines]

    for ln in lower_lines:
        if ln.startswith("xc"):
            parts = ln.split()
            if len(parts) >= 2:
                return parts[1].upper()
            return None

    if any(ln.startswith("mp2") for ln in lower_lines):
        return "MP2"
    if any(ln.startswith("hf") for ln in lower_lines):
        return "HF"
    if any(ln.startswith("dft") for ln in lower_lines):
        return "HF"
    return None


def _as_float(value: Any) -> float | None:
    try:
        if value is None:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value: Any) -> int | None:
    try:
        if value is None:
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _as_bool(value: Any) -> bool | None:
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


def _extract_excited_states(raw_json: dict[str, Any]) -> list[dict[str, Any]]:
    tasks = raw_json.get("tasks")
    if not isinstance(tasks, list):
        return []

    rows: list[dict[str, Any]] = []
    for task_index, task in enumerate(tasks):
        if not isinstance(task, dict):
            continue

        excitations = task.get("excitations")
        if not isinstance(excitations, list):
            continue

        task_model = str(task["model"]) if task.get("model") is not None else None
        task_type = str(task["type"]) if task.get("type") is not None else None
        task_nfreeze = _as_int(task.get("nfreeze"))

        for excitation_index, excitation in enumerate(excitations):
            if not isinstance(excitation, dict):
                continue

            row: dict[str, Any] = {
                "task_index": task_index,
                "task_model": task_model,
                "task_type": task_type,
                "task_nfreeze": task_nfreeze,
                "excitation_index": excitation_index,
            }
            for key, value in excitation.items():
                if isinstance(value, (bool, int, float, str)) or value is None:
                    row[str(key)] = value

            omega_au = _as_float(excitation.get("omega"))
            row["omega_au"] = omega_au
            row["omega_ev"] = omega_au * _HARTREE_TO_EV if omega_au is not None else None

            irrep = row.get("irrep")
            if irrep is not None:
                row["irrep"] = str(irrep)

            rows.append(row)

    return rows


def _lookup_frequency_value(mapping: Any, freq_key: Any) -> Any:
    if not isinstance(mapping, dict):
        return None
    if freq_key in mapping:
        return mapping[freq_key]

    key_str = str(freq_key)
    if key_str in mapping:
        return mapping[key_str]

    target = _as_float(freq_key)
    if target is None:
        return None
    for k, v in mapping.items():
        kval = _as_float(k)
        if kval is None:
            continue
        if abs(kval - target) <= 1e-10:
            return v
    return None


def _sum_metric(rows: list[dict[str, Any]], key: str) -> float:
    total = 0.0
    for row in rows:
        value = row.get(key)
        if isinstance(value, (int, float)):
            total += float(value)
    return float(total)


def _extract_response_metadata(raw_json: dict[str, Any]) -> dict[str, Any] | None:
    tasks = raw_json.get("tasks")
    if isinstance(tasks, list):
        for task in tasks:
            if not isinstance(task, dict):
                continue
            if task.get("type") != "response":
                continue
            metadata = task.get("metadata")
            if isinstance(metadata, dict):
                return metadata

    top_metadata = raw_json.get("metadata")
    if isinstance(top_metadata, dict) and isinstance(top_metadata.get("states"), dict):
        return top_metadata

    if isinstance(raw_json.get("states"), dict):
        return raw_json

    return None


def _extract_state_point_timing_rows(metadata: dict[str, Any]) -> list[dict[str, Any]]:
    states = metadata.get("states")
    if not isinstance(states, dict):
        return []

    rows: list[dict[str, Any]] = []
    for state_id, state_data in states.items():
        if not isinstance(state_data, dict):
            continue
        protocols = state_data.get("protocols")
        if not isinstance(protocols, dict):
            continue

        for protocol_id, protocol_data in protocols.items():
            if not isinstance(protocol_data, dict):
                continue

            timings = protocol_data.get("timings")
            if not isinstance(timings, dict):
                continue

            converged = protocol_data.get("converged")
            saved = protocol_data.get("saved")
            restart_provenance = protocol_data.get("restart_provenance")

            for freq_key, timing in timings.items():
                if not isinstance(timing, dict):
                    continue
                cpu_seconds = _as_float(timing.get("cpu_seconds"))
                wall_seconds = _as_float(timing.get("wall_seconds"))
                if cpu_seconds is None and wall_seconds is None:
                    continue

                restart_data = _lookup_frequency_value(restart_provenance, freq_key)
                if not isinstance(restart_data, dict):
                    restart_data = {}

                rows.append(
                    {
                        "timing_kind": "state_point",
                        "state_id": str(state_id),
                        "protocol": str(protocol_id),
                        "frequency": _as_float(freq_key),
                        "cpu_seconds": cpu_seconds,
                        "wall_seconds": wall_seconds,
                        "converged": _as_bool(_lookup_frequency_value(converged, freq_key)),
                        "saved": _as_bool(_lookup_frequency_value(saved, freq_key)),
                        "restart_kind": str(restart_data["kind"])
                        if restart_data.get("kind") is not None
                        else None,
                        "restart_loaded_from_disk": _as_bool(
                            restart_data.get("loaded_from_disk")
                        ),
                        "restart_promoted_from_static": _as_bool(
                            restart_data.get("promoted_from_static")
                        ),
                        "restart_source_frequency": _as_float(
                            restart_data.get("source_frequency")
                        ),
                        "restart_source_protocol": str(restart_data["source_protocol"])
                        if restart_data.get("source_protocol") is not None
                        else None,
                        "derived_state_id": None,
                        "owner_group": None,
                        "success": None,
                    }
                )
    return rows


def _extract_derived_request_timing_rows(metadata: dict[str, Any]) -> list[dict[str, Any]]:
    derived = metadata.get("derived_state_planner")
    if not isinstance(derived, dict):
        return []
    execution = derived.get("execution")
    if not isinstance(execution, dict):
        return []

    request_timings = execution.get("request_timings")
    if not isinstance(request_timings, list):
        return []

    rows: list[dict[str, Any]] = []
    for entry in request_timings:
        if not isinstance(entry, dict):
            continue
        cpu_seconds = _as_float(entry.get("cpu_seconds"))
        wall_seconds = _as_float(entry.get("wall_seconds"))
        if cpu_seconds is None and wall_seconds is None:
            continue

        rows.append(
            {
                "timing_kind": "derived_request",
                "state_id": None,
                "protocol": None,
                "frequency": None,
                "cpu_seconds": cpu_seconds,
                "wall_seconds": wall_seconds,
                "converged": None,
                "saved": None,
                "restart_kind": None,
                "restart_loaded_from_disk": None,
                "restart_promoted_from_static": None,
                "restart_source_frequency": None,
                "restart_source_protocol": None,
                "derived_state_id": str(entry["derived_state_id"])
                if entry.get("derived_state_id") is not None
                else None,
                "owner_group": _as_int(entry.get("owner_group")),
                "success": _as_bool(entry.get("success")),
            }
        )
    return rows


def _build_timing_summary(
    metadata: dict[str, Any],
    state_rows: list[dict[str, Any]],
    derived_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    state_runtime = metadata.get("state_parallel_runtime")
    if not isinstance(state_runtime, dict):
        state_runtime = {}

    state_planner = metadata.get("state_parallel_planner")
    if not isinstance(state_planner, dict):
        state_planner = {}

    derived_planner = metadata.get("derived_state_planner")
    if not isinstance(derived_planner, dict):
        derived_planner = {}
    derived_execution = derived_planner.get("execution")
    if not isinstance(derived_execution, dict):
        derived_execution = {}

    effective_groups = state_runtime.get("effective_point_groups")
    if effective_groups is None:
        effective_groups = state_planner.get("effective_point_groups")

    start_protocol_index = state_runtime.get("effective_point_parallel_start_protocol_index")
    if start_protocol_index is None:
        start_protocol_index = state_planner.get("point_parallel_start_protocol_index")

    effective_mode = state_planner.get("effective_mode")
    if effective_mode is None:
        effective_mode = state_planner.get("requested_mode")

    return {
        "state_point_count": len(state_rows),
        "state_point_cpu_seconds": _sum_metric(state_rows, "cpu_seconds"),
        "state_point_wall_seconds": _sum_metric(state_rows, "wall_seconds"),
        "derived_request_count": len(derived_rows),
        "derived_request_cpu_seconds": _sum_metric(derived_rows, "cpu_seconds"),
        "derived_request_wall_seconds": _sum_metric(derived_rows, "wall_seconds"),
        "derived_request_success_count": sum(1 for row in derived_rows if row.get("success") is True),
        "derived_request_failed_count": sum(1 for row in derived_rows if row.get("success") is False),
        "derived_execution_mode": str(derived_execution["mode"])
        if derived_execution.get("mode") is not None
        else None,
        "derived_execution_attempted": _as_bool(derived_execution.get("attempted")),
        "derived_execution_groups": _as_int(derived_execution.get("execution_groups")),
        "derived_execution_total_cpu_seconds": _as_float(
            derived_execution.get("total_cpu_seconds")
        ),
        "derived_execution_total_wall_seconds": _as_float(
            derived_execution.get("total_wall_seconds")
        ),
        "derived_execution_completed_requests": _as_int(
            derived_execution.get("completed_requests")
        ),
        "derived_execution_failed_requests": _as_int(derived_execution.get("failed_requests")),
        "derived_execution_blocked_requests": _as_int(
            derived_execution.get("blocked_requests")
        ),
        "derived_execution_ready_requests": _as_int(derived_execution.get("ready_requests")),
        "state_parallel_effective_groups": _as_int(effective_groups),
        "state_parallel_start_protocol_index": _as_int(start_protocol_index),
        "state_parallel_effective_mode": str(effective_mode)
        if effective_mode is not None
        else None,
        "state_parallel_frequency_partition_policy": str(
            state_planner["frequency_partition_policy"]
        )
        if state_planner.get("frequency_partition_policy") is not None
        else None,
        "state_parallel_requested_groups": _as_int(state_planner.get("requested_groups")),
        "state_parallel_world_size": _as_int(state_planner.get("world_size")),
        "state_parallel_execution_enabled": _as_bool(state_planner.get("execution_enabled")),
        "state_parallel_subgroup_parallel_enabled": _as_bool(
            state_planner.get("subgroup_parallel_enabled")
        ),
        "state_parallel_restart_point_parallel_promoted": _as_bool(
            state_runtime.get("restart_point_parallel_promoted")
        ),
        "state_parallel_restart_protocol0_saved_complete": _as_bool(
            state_runtime.get("restart_protocol0_saved_complete")
        ),
    }


def _parse_madness_timings(raw_json: dict[str, Any]) -> dict[str, Any]:
    metadata = _extract_response_metadata(raw_json)
    if not isinstance(metadata, dict):
        return {
            "schema_version": 1,
            "point_rows": [],
            "state_point_rows": [],
            "derived_request_rows": [],
            "summary": {},
        }

    state_rows = _extract_state_point_timing_rows(metadata)
    derived_rows = _extract_derived_request_timing_rows(metadata)
    point_rows = [*state_rows, *derived_rows]
    summary = _build_timing_summary(metadata, state_rows, derived_rows)

    return {
        "schema_version": 1,
        "point_rows": point_rows,
        "state_point_rows": state_rows,
        "derived_request_rows": derived_rows,
        "summary": summary,
    }


def _merge_timing_payloads(
    primary: dict[str, Any],
    secondary: dict[str, Any],
) -> dict[str, Any]:
    def _rows(payload: dict[str, Any], key: str) -> list[dict[str, Any]]:
        rows = payload.get(key)
        if isinstance(rows, list):
            return [row for row in rows if isinstance(row, dict)]
        return []

    def _summary(payload: dict[str, Any]) -> dict[str, Any]:
        summary = payload.get("summary")
        if isinstance(summary, dict):
            return summary
        return {}

    primary_state = _rows(primary, "state_point_rows")
    secondary_state = _rows(secondary, "state_point_rows")
    primary_derived = _rows(primary, "derived_request_rows")
    secondary_derived = _rows(secondary, "derived_request_rows")

    state_rows = primary_state if primary_state else secondary_state
    derived_rows = primary_derived if primary_derived else secondary_derived
    point_rows = [*state_rows, *derived_rows]

    primary_summary = _summary(primary)
    secondary_summary = _summary(secondary)
    summary = dict(primary_summary)
    for k, v in secondary_summary.items():
        if k not in summary or summary[k] is None:
            summary[k] = v

    return {
        "schema_version": 1,
        "point_rows": point_rows,
        "state_point_rows": state_rows,
        "derived_request_rows": derived_rows,
        "summary": summary,
    }


def parse_run(calc: Calculation) -> None:
    """
    Populate calc.data for MADNESS runs.

    Supports:
      - MADQC style: *.calc_info.json
      - Legacy molresponse style: output.json

    Both are parsed by the same legacy madqc_parser (schema is nearly the same),
    with component label normalization handled inside the legacy parser.
    """
    json_path, style = _select_input_json(calc)
    if json_path is None:
        meta_path = calc.artifacts.get("responses_metadata_json")
        if isinstance(meta_path, Path) and meta_path.exists():
            calc.meta["style"] = "responses_metadata"
            calc.data["timings"] = _parse_madness_timings(_read_json(meta_path))
            if calc.basis is None:
                calc.basis = "mra"
            calc.meta.setdefault("basis", calc.basis)
        return

    calc.meta["style"] = style

    from gecko.plugins.madness.legacy.madness_data import madqc_parser

    obj = madqc_parser(json_path)

    # Keep raw json as source of truth during migration
    calc.data["raw_json"] = _read_json(json_path)
    calc.data["timings"] = _parse_madness_timings(calc.data["raw_json"])
    calc.data["excited_states"] = _extract_excited_states(calc.data["raw_json"])
    meta_path = calc.artifacts.get("responses_metadata_json")
    if isinstance(meta_path, Path) and meta_path.exists():
        calc.data["timings"] = _merge_timing_payloads(
            calc.data["timings"],
            _parse_madness_timings(_read_json(meta_path)),
        )

    # Basis label (MRA) inference.
    # Priority: paired input .in (MADQC) -> input.json -> parsed payload.
    basis = None
    input_in = calc.artifacts.get("input_in")
    if isinstance(input_in, Path) and input_in.exists():
        try:
            basis = _infer_mra_basis_from_input_in_text(
                input_in.read_text(encoding="utf-8", errors="ignore")
            )
        except Exception:
            basis = None
        if calc.meta.get("method") is None:
            try:
                calc.meta["method"] = _infer_method_from_input_in_text(
                    input_in.read_text(encoding="utf-8", errors="ignore")
                )
            except Exception:
                pass

    input_json = calc.artifacts.get("input_json")
    if isinstance(input_json, Path) and input_json.exists():
        try:
            basis = _infer_mra_basis_from_obj(_read_json(input_json))
        except Exception:
            basis = None
    if basis is None:
        basis = _infer_mra_basis_from_obj(calc.data.get("raw_json"))
    calc.basis = basis or "mra"
    calc.meta["basis"] = calc.basis

    # Tensor-first hyperpolarizability
    calc.data["beta"] = _beta_df_to_tensor(obj.beta_pivot)
    if not _tensor_has_rows(calc.data["beta"]):
        calc.data["beta"] = _legacy_beta_to_tensor(calc.data["raw_json"])

    calc.data["alpha"] = _beta_df_to_tensor(obj.alpha_pivot)
    if not _tensor_has_rows(calc.data["alpha"]):
        calc.data["alpha"] = _legacy_alpha_to_tensor(calc.data["raw_json"])

    calc.data["raman"] = {"polarization_frequencies": obj.polarization_frequencies,
                          "vibrational_frequencies": obj.vibrational_frequencies,
                          "polarizability_derivatives": obj.polarizability_derivatives,
                          "polarizability_derivatives_by_mode": obj.polarizability_derivatives,
                            "raman_by_freq": obj.raman_by_freq}


    # Keep other useful arrays (best-effort; may be None for some runs)
    calc.data["orbital_energies"] = obj.orbital_energies
    calc.data["hessian"] = obj.hessian
    calc.data["normal_modes"] = obj.normal_modes
    calc.data["vibrational_frequencies"] = obj.vibrational_frequencies
    calc.data["polarization_frequencies"] = getattr(obj, "polarization_frequencies", None)

    calc.data["molecule"] = obj.molecule
    calc.molecule = obj.molecule
    calc.meta["ground_state_energy"] = obj.ground_state_energy

    if obj.raman_by_freq is not None:
        for _, rows in obj.raman_by_freq.items():
            if isinstance(rows, list):
                rows.sort(key=lambda r: (float(r.get("freq_cm1", 0.0)), int(r.get("mode", 0))))
        calc.data["raman"] = {
            "polarization_frequencies": np.asarray(
                obj.polarization_frequencies, dtype=float
            )
            if obj.polarization_frequencies is not None
            else np.asarray([], dtype=float),
            "vibrational_frequencies": np.asarray(
                obj.vibrational_frequencies, dtype=float
            )
            if obj.vibrational_frequencies is not None
            else np.asarray([], dtype=float),
            "polarizability_derivatives": obj.polarizability_derivatives,
            "polarizability_derivatives_by_mode": obj.polarizability_derivatives_normal_modes,
            "raman_by_freq": obj.raman_by_freq,
        }

    if calc.molecule is None:
        calc.molecule = _load_molecule_from_input(calc.artifacts.get("input_json"))
        if calc.molecule is not None:
            calc.meta.setdefault("molecule_source", "input.json")
        else:
            calc.meta.setdefault("warnings", []).append(
                "MADNESS output missing molecule; input.json not found or invalid."
            )



def _select_input_json(calc: Calculation) -> tuple[Path | None, str | None]:
    """
    Decide which JSON file to parse for MADNESS.

    Priority:
      1) MADQC: calc_info_json
      2) legacy: output_json
    """
    ci_path = calc.artifacts.get("calc_info_json")
    if ci_path and ci_path.exists():
        return ci_path, "madqc"

    out_path = calc.artifacts.get("output_json")
    if out_path and out_path.exists():
        return out_path, "molresponse"

    return None, None


def _read_json(path: Path) -> dict[str, Any]:
    import json
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _load_molecule_from_input(path: Path | None):
    if path is None or not path.exists():
        return None
    data = _read_json(path)
    mol = data.get("molecule")
    if not isinstance(mol, dict):
        return None
    symbols = mol.get("symbols")
    geometry = mol.get("geometry")
    if symbols is None or geometry is None:
        return None

    import numpy as np
    import qcelemental as qcel

    units = mol.get("units") or mol.get("parameters", {}).get("units")
    coords = np.asarray(geometry, dtype=float)
    if isinstance(units, str) and units.lower() in ("bohr", "atomic", "au"):
        coords = coords * qcel.constants.bohr2angstroms
    symbols_sorted, coords_sorted = canonicalize_atom_order(list(symbols), coords, decimals=10)
    kwargs = {
        "symbols": symbols_sorted,
        "geometry": coords_sorted,
    }
    mol_charge = mol.get("charge")
    mol_mult = mol.get("multiplicity")
    if mol_charge is not None:
        kwargs["molecular_charge"] = mol_charge
    if mol_mult is not None:
        kwargs["molecular_multiplicity"] = mol_mult
    return qcel.models.Molecule(**kwargs)
