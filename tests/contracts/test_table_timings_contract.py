from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from gecko.core.load import load_calc
from gecko.tables.builder import TableBuilder


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


def test_table_timing_points_contract(tmp_path: Path) -> None:
    calc = load_calc(_write_timing_fixture(tmp_path))
    tb = TableBuilder([calc])

    df = tb.build_timing_points()
    assert not df.empty
    assert len(df) == 3

    required_cols = {
        "calc_id",
        "geom_id",
        "mol_id",
        "molecule_id",
        "label",
        "code",
        "root",
        "basis",
        "method",
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
    assert set(df.columns) >= required_cols

    assert set(df["timing_kind"].dropna().unique()) == {"state_point", "derived_request"}

    state_rows = df[df["timing_kind"] == "state_point"]
    derived_rows = df[df["timing_kind"] == "derived_request"]
    assert len(state_rows) == 2
    assert len(derived_rows) == 1

    assert np.isfinite(state_rows["cpu_seconds"].astype(float)).all()
    assert np.isfinite(state_rows["wall_seconds"].astype(float)).all()
    assert np.isfinite(derived_rows["cpu_seconds"].astype(float)).all()
    assert np.isfinite(derived_rows["wall_seconds"].astype(float)).all()


def test_table_timing_summary_contract(tmp_path: Path) -> None:
    calc = load_calc(_write_timing_fixture(tmp_path))
    tb = TableBuilder([calc])

    df = tb.build_timing_summary()
    assert not df.empty
    assert len(df) == 1

    required_cols = {
        "calc_id",
        "geom_id",
        "mol_id",
        "molecule_id",
        "label",
        "code",
        "root",
        "basis",
        "method",
        "state_point_count",
        "state_point_cpu_seconds",
        "state_point_wall_seconds",
        "derived_request_count",
        "derived_request_cpu_seconds",
        "derived_request_wall_seconds",
        "derived_execution_mode",
        "state_parallel_effective_mode",
        "state_parallel_effective_groups",
    }
    assert set(df.columns) >= required_cols

    row = df.iloc[0]
    assert int(row["state_point_count"]) == 2
    assert int(row["derived_request_count"]) == 1
    assert float(row["state_point_cpu_seconds"]) == pytest.approx(22.0)
    assert float(row["state_point_wall_seconds"]) == pytest.approx(21.0)
    assert float(row["derived_request_cpu_seconds"]) == pytest.approx(1.2)
    assert float(row["derived_request_wall_seconds"]) == pytest.approx(1.1)
    assert str(row["derived_execution_mode"]) == "owner_group_subworld"
