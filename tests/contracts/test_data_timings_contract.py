from __future__ import annotations

import json
from pathlib import Path

import pytest

from gecko.core.load import load_calc


def _write_timing_fixture(tmp_path: Path) -> Path:
    root = tmp_path / "06_mra-timing_h2o"
    root.mkdir(parents=True, exist_ok=True)

    stem = "mad.timing_h2o"
    (root / f"{stem}.in").write_text("dconv 1e-6\nhf\n", encoding="utf-8")

    payload = {
        "tasks": [
            {
                "type": "response",
                "metadata": {
                    "states": {
                        "Dipole_x": {
                            "protocols": {
                                "1e-04": {
                                    "converged": {"0.000": True},
                                    "saved": {"0.000": True},
                                    "restart_provenance": {
                                        "0.000": {
                                            "kind": "initial_guess",
                                            "loaded_from_disk": False,
                                            "promoted_from_static": False,
                                            "source_frequency": None,
                                            "source_protocol": None,
                                        }
                                    },
                                    "timings": {
                                        "0.000": {
                                            "cpu_seconds": 10.5,
                                            "wall_seconds": 10.0,
                                        }
                                    },
                                }
                            }
                        },
                        "Dipole_y": {
                            "protocols": {
                                "1e-04": {
                                    "converged": {"0.000": False},
                                    "saved": {"0.000": False},
                                    "timings": {
                                        "0.000": {
                                            "cpu_seconds": 11.5,
                                            "wall_seconds": 11.0,
                                        }
                                    },
                                }
                            }
                        },
                    },
                    "state_parallel_runtime": {
                        "effective_point_groups": 2,
                        "effective_point_parallel_start_protocol_index": 1,
                        "restart_point_parallel_promoted": True,
                        "restart_protocol0_saved_complete": False,
                    },
                    "state_parallel_planner": {
                        "effective_mode": "point_parallel",
                        "frequency_partition_policy": "block",
                        "requested_groups": 2,
                        "world_size": 8,
                        "execution_enabled": True,
                        "subgroup_parallel_enabled": False,
                    },
                    "derived_state_planner": {
                        "execution": {
                            "mode": "owner_group_subworld",
                            "attempted": True,
                            "completed_requests": 1,
                            "failed_requests": 0,
                            "blocked_requests": 0,
                            "ready_requests": 1,
                            "execution_groups": 2,
                            "total_cpu_seconds": 1.2,
                            "total_wall_seconds": 1.1,
                            "request_timings": [
                                {
                                    "derived_state_id": "Derived_Dipole_x__Dipole_x_Dipole_x_f0.000_0.000",
                                    "owner_group": 1,
                                    "success": True,
                                    "cpu_seconds": 1.2,
                                    "wall_seconds": 1.1,
                                }
                            ],
                        }
                    },
                },
                "properties": {"response_properties": []},
            }
        ]
    }
    (root / f"{stem}.calc_info.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    return root


def _write_responses_metadata_fixture(tmp_path: Path) -> Path:
    root = tmp_path / "07_mra-timing_metadata_only"
    (root / "responses").mkdir(parents=True, exist_ok=True)

    payload = {
        "states": {
            "Dipole_z": {
                "protocols": {
                    "1e-05": {
                        "converged": {"0.000": True},
                        "saved": {"0.000": True},
                        "timings": {"0.000": {"cpu_seconds": 5.0, "wall_seconds": 4.5}},
                    }
                }
            }
        }
    }
    (root / "responses" / "metadata.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    return root


def test_timings_data_contract(tmp_path: Path) -> None:
    calc = load_calc(_write_timing_fixture(tmp_path))
    assert calc.code == "madness"

    timings = calc.data.get("timings")
    assert isinstance(timings, dict)
    assert set(timings.keys()) >= {
        "schema_version",
        "point_rows",
        "state_point_rows",
        "derived_request_rows",
        "summary",
    }
    assert timings["schema_version"] == 1

    point_rows = timings["point_rows"]
    state_rows = timings["state_point_rows"]
    derived_rows = timings["derived_request_rows"]
    assert isinstance(point_rows, list) and len(point_rows) == 3
    assert isinstance(state_rows, list) and len(state_rows) == 2
    assert isinstance(derived_rows, list) and len(derived_rows) == 1

    required_point_keys = {
        "timing_kind",
        "state_id",
        "protocol",
        "frequency",
        "cpu_seconds",
        "wall_seconds",
        "converged",
        "saved",
        "restart_kind",
        "restart_loaded_from_disk",
        "restart_promoted_from_static",
        "restart_source_frequency",
        "restart_source_protocol",
        "derived_state_id",
        "owner_group",
        "success",
    }
    for row in point_rows:
        assert isinstance(row, dict)
        assert set(row.keys()) >= required_point_keys

    summary = timings["summary"]
    assert isinstance(summary, dict)
    assert summary["state_point_count"] == 2
    assert summary["derived_request_count"] == 1
    assert summary["state_point_cpu_seconds"] == pytest.approx(22.0)
    assert summary["state_point_wall_seconds"] == pytest.approx(21.0)
    assert summary["derived_request_cpu_seconds"] == pytest.approx(1.2)
    assert summary["derived_request_wall_seconds"] == pytest.approx(1.1)
    assert summary["derived_execution_mode"] == "owner_group_subworld"
    assert summary["state_parallel_effective_mode"] == "point_parallel"
    assert summary["state_parallel_effective_groups"] == 2


def test_timings_data_from_responses_metadata_only(tmp_path: Path) -> None:
    calc = load_calc(_write_responses_metadata_fixture(tmp_path))
    assert calc.code == "madness"
    assert calc.meta.get("style") == "responses_metadata"

    timings = calc.data.get("timings")
    assert isinstance(timings, dict)

    state_rows = timings.get("state_point_rows")
    assert isinstance(state_rows, list)
    assert len(state_rows) == 1
    row = state_rows[0]
    assert row.get("timing_kind") == "state_point"
    assert row.get("state_id") == "Dipole_z"
    assert row.get("protocol") == "1e-05"
    assert row.get("cpu_seconds") == pytest.approx(5.0)
    assert row.get("wall_seconds") == pytest.approx(4.5)
