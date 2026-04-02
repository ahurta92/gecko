"""Parameter dataclasses for MADNESS input file sections.

Each class mirrors a C++ parameter struct in the MADNESS source:

  DFTParams      → src/madness/chem/CalculationParameters.h
  MoleculeParams → src/madness/chem/molecule.h  (GeometryParameters)
  ResponseParams → src/apps/molresponse/response_parameters.h

All fields default to None.  None means "omit from the .in file and let
MADNESS use its own internal default."  Set a field to override it.
"""

from __future__ import annotations

import dataclasses
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Value renderer
# ---------------------------------------------------------------------------


def _render_value(v: object) -> str:
    """Convert a Python value to MADNESS .in file string representation."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list):
        return "[" + ",".join(_render_value(x) for x in v) + "]"
    if isinstance(v, float) and v != 0.0 and abs(v) < 1e-2:
        return f"{v:.1e}"
    return str(v)


# ---------------------------------------------------------------------------
# DFTParams  (dft … end section)
# ---------------------------------------------------------------------------


@dataclass
class DFTParams:
    """Parameters for the ``dft`` section of a MADNESS .in file.

    Notable MADNESS defaults (as of current source):
        xc="hf", k=-1 (auto from thresh), econv=1e-5, dconv=1e-4,
        maxiter=25, protocol=[1e-4,1e-6], localize="new", save=True,
        print_level=3
    """

    xc: Optional[str] = None
    """Exchange-correlation functional: "hf", "lda", "b3lyp", "pbe0", …"""

    k: Optional[int] = None
    """Wavelet order (6–10 are typical; -1 = auto from thresh)."""

    l: Optional[float] = None
    """User coordinate box size in atomic units (default 20)."""

    charge: Optional[float] = None
    """Total molecular charge."""

    econv: Optional[float] = None
    """Energy convergence threshold."""

    dconv: Optional[float] = None
    """Density convergence threshold."""

    maxiter: Optional[int] = None
    """Maximum SCF iterations."""

    protocol: Optional[list[float]] = None
    """Multi-level convergence protocol thresholds, e.g. [1e-4, 1e-6]."""

    restart: Optional[bool] = None
    """Restart from orbitals on disk."""

    save: Optional[bool] = None
    """Save orbitals to disk (default True)."""

    localize: Optional[str] = None
    """Localization method: "pm", "boys", "new", "canon"."""

    pointgroup: Optional[str] = None
    """Point group: "c1", "c2v", "d2h", …"""

    spin_restricted: Optional[bool] = None
    """Spin-restricted calculation (default True)."""

    nopen: Optional[int] = None
    """Number of unpaired electrons (nalpha - nbeta)."""

    dipole: Optional[bool] = None
    """Calculate dipole moment at each SCF step."""

    derivatives: Optional[bool] = None
    """Calculate nuclear derivatives."""

    gopt: Optional[bool] = None
    """Geometry optimization."""

    gtol: Optional[float] = None
    """Geometry optimization gradient tolerance."""

    gmaxiter: Optional[int] = None
    """Geometry optimization max iterations."""

    algopt: Optional[str] = None
    """Optimization algorithm: "bfgs", "cg"."""

    print_level: Optional[int] = None
    """Verbosity: 0=none, 1=final, 2=iter, 3=timings, 10=debug."""

    maxsub: Optional[int] = None
    """DIIS subspace size (0/1 = disable)."""

    aobasis: Optional[str] = None
    """AO basis for initial guess: "6-31g", "3-21g", "sto-6g", …"""

    lo: Optional[float] = None
    """Smallest length scale to resolve."""

    maxrotn: Optional[float] = None
    """Max orbital rotation per SCF step."""


# ---------------------------------------------------------------------------
# MoleculeParams  (molecule … end section)
# ---------------------------------------------------------------------------


@dataclass
class MoleculeParams:
    """Parameters for the ``molecule`` section of a MADNESS .in file.

    Notable MADNESS defaults:
        eprec=1e-4, units="atomic", no_orient=False
    """

    eprec: Optional[float] = None
    """Smoothing parameter for the nuclear potential."""

    units: Optional[str] = None
    """Coordinate units: "atomic" (Bohr) or "angstrom"."""

    no_orient: Optional[bool] = None
    """If True, do not reorient/symmetrize the molecule."""

    field: Optional[list[float]] = None
    """External electric field [fx, fy, fz] in atomic units."""

    symtol: Optional[float] = None
    """Distance threshold for symmetry detection."""

    core_type: Optional[str] = None
    """Core potential type: "none" or "mcp"."""


# ---------------------------------------------------------------------------
# ResponseParams  (response … end section)
# ---------------------------------------------------------------------------

# Mapping from Python attribute names to MADNESS dot-notation keys
_RESPONSE_FIELD_MAP: dict[str, str] = {
    "dipole_frequencies": "dipole.frequencies",
    "dipole_directions": "dipole.directions",
    "nuclear_directions": "nuclear.directions",
    "nuclear_frequencies": "nuclear.frequencies",
    "nuclear_atom_indices": "nuclear.atom_indices",
}


@dataclass
class ResponseParams:
    """Parameters for the ``response`` section of a MADNESS .in file.

    Notable MADNESS defaults:
        maxiter=5, dconv=1e-4, protocol=[1e-4,1e-6], kain=False,
        maxsub=5, save=False, restart=False, localize="canon",
        print_level=3
    """

    prefix: Optional[str] = None
    """Prefix for output files."""

    archive: Optional[str] = None
    """Path to ground-state restart data (default "../moldft.restartdata")."""

    maxiter: Optional[int] = None
    """Maximum response iterations."""

    dconv: Optional[float] = None
    """Density convergence threshold."""

    protocol: Optional[list[float]] = None
    """Multi-level convergence protocol thresholds."""

    kain: Optional[bool] = None
    """Krylov Accelerated Inexact Newton solver."""

    maxsub: Optional[int] = None
    """KAIN subspace size (0/1 = disable)."""

    maxrotn: Optional[float] = None
    """Max orbital rotation per response step."""

    save: Optional[bool] = None
    """Save response orbitals to disk."""

    restart: Optional[bool] = None
    """Restart from saved response orbitals."""

    localize: Optional[str] = None
    """Localization: "pm", "boys", "new", "canon"."""

    print_level: Optional[int] = None
    """Verbosity: 0=none, 1=final, 2=iter, 3=timings, 10=debug."""

    # --- perturbation / property fields ---

    dipole_frequencies: Optional[list[float]] = None
    """Frequencies for dipole response (written as dipole.frequencies)."""

    dipole_directions: Optional[str] = None
    """Directions for dipole response: "xyz", "x", "y", "z"."""

    quadratic: Optional[bool] = None
    """Compute quadratic (hyperpolarizability) response."""

    nuclear: Optional[bool] = None
    """Compute nuclear (Raman) response."""

    nuclear_directions: Optional[str] = None
    """Directions for nuclear response (written as nuclear.directions)."""

    nuclear_frequencies: Optional[float] = None
    """Frequency for nuclear response (written as nuclear.frequencies)."""

    nuclear_atom_indices: Optional[list[int]] = None
    """Atom indices for nuclear response (written as nuclear.atom_indices)."""

    requested_properties: Optional[list[str]] = None
    """Properties to output, e.g. ["polarizability", "hyperpolarizability"]."""

    property: Optional[bool] = None
    """Boolean flag to enable response property output."""


# ---------------------------------------------------------------------------
# Numerical tier loader
# ---------------------------------------------------------------------------

_FIXTURES_DIR = Path(
    os.environ.get(
        "GECKO_FIXTURES_DIR",
        "/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next"
        "/src/apps/molresponse_v3/tests/fixtures",
    )
)


def _load_tier(tier: str) -> tuple[DFTParams | None, MoleculeParams | None]:
    """Load DFTParams and MoleculeParams for the given accuracy tier.

    Reads numerical_settings.json from GECKO_FIXTURES_DIR.
    Returns (DFTParams, MoleculeParams) or (None, None) if tier is "none"
    or the settings file cannot be found.
    """
    if tier == "none":
        return None, None

    settings_path = _FIXTURES_DIR / "numerical_settings.json"
    if not settings_path.exists():
        raise FileNotFoundError(
            f"numerical_settings.json not found at {settings_path}. "
            f"Set GECKO_FIXTURES_DIR to the fixtures directory."
        )

    with open(settings_path) as f:
        data = json.load(f)

    tiers = data.get("tiers", {})
    if tier not in tiers:
        valid = ", ".join(tiers.keys())
        raise ValueError(f"Unknown tier {tier!r}. Valid tiers: {valid}")

    tier_data = tiers[tier]
    dft_raw = tier_data.get("dft", {})
    mol_raw = tier_data.get("molecule", {})

    dft_params = DFTParams(
        dconv=dft_raw.get("dconv"),
        econv=dft_raw.get("econv"),
        protocol=dft_raw.get("protocol"),
        l=dft_raw.get("l"),
        maxiter=dft_raw.get("maxiter"),
        localize=dft_raw.get("localize"),
        dipole=dft_raw.get("dipole"),
    )
    mol_params = MoleculeParams(
        eprec=mol_raw.get("eprec"),
    )
    return dft_params, mol_params
