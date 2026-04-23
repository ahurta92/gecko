"""Tests for gecko.madnessproject — input generation and output loading."""

from __future__ import annotations

from pathlib import Path

import pytest

from gecko.madnessproject import (
    CalculationParameters,
    Molecule,
    ResponseParameters,
    load_output,
    madness_input,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "load_calc"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@pytest.fixture()
def water_atomic() -> Molecule:
    """H2O in atomic units matching the raman fixture geometry."""
    return Molecule(
        atoms="O 0.0 0.0 0.21293780; H 0.0 1.42131521 -0.85175395; H 0.0 -1.42131521 -0.85175395",
        units="atomic",
        eprec=1e-6,
        no_orient=True,
    )


# ---------------------------------------------------------------------------
# Input generation — block order and content
# ---------------------------------------------------------------------------

class TestMadnessInput:
    def test_block_order_dft_molecule_response(self, water_atomic):
        """dft block must come before molecule, molecule before response."""
        calc = CalculationParameters(dipole=True, econv=1e-7, maxiter=20, gopt=True,
                                     protocol=[1e-4, 1e-6, 1e-7])
        resp = ResponseParameters(
            dipole=True,
            frequencies=[0.0, 0.02, 0.04, 0.06, 0.08, 0.1],
            quadratic=True,
            nuclear=True,
            nuclear_directions="xyz",
            nuclear_frequencies=0.0,
            nuclear_atom_indices=[0, 1, 2],
        )
        content = madness_input(calc, water_atomic, resp)
        dft_pos = content.index("dft\n")
        mol_pos = content.index("molecule\n")
        resp_pos = content.index("response\n")
        assert dft_pos < mol_pos < resp_pos, (
            f"Expected dft < molecule < response, got positions {dft_pos}, {mol_pos}, {resp_pos}"
        )

    def test_dft_block_contains_expected_keys(self, water_atomic):
        calc = CalculationParameters(dipole=True, econv=1e-7, maxiter=20, gopt=True)
        content = madness_input(calc, water_atomic)
        assert "dipole true" in content
        assert "econv 1e-07" in content
        assert "maxiter 20" in content
        assert "gopt true" in content

    def test_molecule_block_contains_atoms(self, water_atomic):
        calc = CalculationParameters()
        content = madness_input(calc, water_atomic)
        assert "O " in content
        assert "H " in content
        assert "units atomic" in content

    def test_response_block_contains_frequencies(self, water_atomic):
        calc = CalculationParameters()
        resp = ResponseParameters(frequencies=[0.0, 0.02, 0.04])
        content = madness_input(calc, water_atomic, resp)
        assert "dipole.frequencies" in content
        assert "0.0,0.02,0.04" in content

    def test_no_response_omits_response_block(self, water_atomic):
        calc = CalculationParameters()
        content = madness_input(calc, water_atomic)
        assert "response" not in content

    def test_emit_property_flag(self, water_atomic):
        calc = CalculationParameters()
        resp = ResponseParameters(emit_property=True)
        content = madness_input(calc, water_atomic, resp)
        assert "property true" in content

    def test_suppress_property_flag(self, water_atomic):
        calc = CalculationParameters()
        resp = ResponseParameters(emit_property=False)
        content = madness_input(calc, water_atomic, resp)
        # "property true" should not appear; "property" key absent
        assert "property true" not in content


# ---------------------------------------------------------------------------
# Output loading — round-trip against real fixtures
# ---------------------------------------------------------------------------

class TestLoadOutput:
    def test_n2_energy_not_none(self):
        """load_output on n2 output.json must return a non-None energy."""
        path = FIXTURES / "01_mra-d04_n2" / "output.json"
        result = load_output(path)
        assert result.energy is not None, "energy should not be None for n2 fixture"

    def test_n2_energy_value(self):
        """n2 SCF energy matches fixture value."""
        path = FIXTURES / "01_mra-d04_n2" / "output.json"
        result = load_output(path)
        assert result.energy == pytest.approx(-413.6193973843323, rel=1e-6)

    def test_bh2cl_beta_not_none(self):
        """load_output on bh2cl outputs.json must return non-None beta."""
        path = FIXTURES / "00_mra-d06_bh2cl" / "outputs.json"
        result = load_output(path)
        assert result.beta is not None, "beta should not be None for bh2cl fixture"

    def test_calc_info_loads_without_error(self):
        """load_output on calc_info.json format must not raise."""
        path = FIXTURES / "03_mra-raman_h2o" / "mad.h2o_gopt.calc_info.json"
        result = load_output(path)
        assert result is not None
