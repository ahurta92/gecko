from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Literal


JobType = Literal["optimize", "raman"]
Program = Literal["madness", "dalton"]


@dataclass(frozen=True)
class RamanBenchDB:
    root: Path

    def __post_init__(self):
        object.__setattr__(self, "root", self.root.expanduser().resolve())

    # ---------------- discovery ----------------
    def molecules(self) -> list[str]:
        if not self.root.exists():
            return []
        return sorted([p.name for p in self.root.iterdir() if p.is_dir()])

    # ---------------- dirs ----------------
    def molecule_dir(self, mol: str) -> Path:
        return self.root / mol

    def madness_dir(self, mol: str) -> Path:
        return self.molecule_dir(mol) / "madness"

    def dalton_dir(self, mol: str) -> Path:
        return self.molecule_dir(mol) / "dalton"

    # ---------------- inputs ----------------
    def madness_input(self, mol: str) -> Path:
        return self.madness_dir(mol) / "mad.raman.in"

    def dalton_opt_input(self, mol: str) -> Path:
        return self.dalton_dir(mol) / "optimize.dal"

    def dalton_raman_input(self, mol: str) -> Path:
        return self.dalton_dir(mol) / "raman.dal"

    def dalton_mol(self, mol: str, basis: str) -> Path:
        return self.dalton_dir(mol) / f"{mol}_{basis}.mol"

    def dalton_opt_mol(self, mol: str, basis: str) -> Path:
        return self.dalton_dir(mol) / f"{mol}_{basis}_opt.mol"

    # ---------------- helpers ----------------
    @staticmethod
    def _latest_by_mtime(paths: Iterable[Path]) -> Optional[Path]:
        paths = [p for p in paths if p.exists()]
        return max(paths, key=lambda p: p.stat().st_mtime) if paths else None

    def find_latest(self, mol: str, program: Program, pattern: str) -> Optional[Path]:
        base = self.molecule_dir(mol) / program
        if not base.is_dir():
            return None
        return self._latest_by_mtime(base.glob(pattern))

    # ---------------- dalton outputs (basis-aware) ----------------
    def dalton_output(self, mol: str, job: JobType, basis: str) -> Optional[Path]:
        """
        Preferred Dalton output name convention:
          optimize_{basis}.out
          raman_{basis}.out

        Falls back to globs if the exact name isn't present.
        """
        ddir = self.dalton_dir(mol)
        if not ddir.is_dir():
            return None

        exact_candidates = [
            ddir / f"{job}_{basis}.out",
        ]
        for p in exact_candidates:
            if p.exists():
                return p

        # Fall back to common alternates
        # (useful if you ever change naming or have legacy files)
        candidates = [
            *ddir.glob(f"{job}_{basis}*.out"),
            *ddir.glob(f"{job}_opt_{basis}*.out"),
            *ddir.glob(f"{job}*{basis}*.out")
        ]
        # also consider jobs where basis isn't in filename
        candidates.extend(ddir.glob(f"{job}.out"))
        return self._latest_by_mtime(candidates)

    def dalton_error(self, mol: str, job: JobType, basis: str) -> Optional[Path]:
        """Same logic as dalton_output, but for .err files."""
        ddir = self.dalton_dir(mol)
        if not ddir.is_dir():
            return None

        exact_candidates = [
            ddir / f"{job}_{basis}.err",
            ddir / f"{job}_opt_{basis}.err",
            ddir / f"{job}_{mol}_{basis}.err",
            ddir / f"{job}_{mol}_opt_{basis}.err",
            ddir / f"{job}.err",
        ]
        for p in exact_candidates:
            if p.exists():
                return p

        candidates = [
            *ddir.glob(f"{job}_{basis}*.err"),
            *ddir.glob(f"{job}_opt_{basis}*.err"),
            *ddir.glob(f"{job}*{basis}*.err"),
            *ddir.glob(f"{job}*.err"),
        ]
        return self._latest_by_mtime(candidates)

    # ---------------- madness outputs ----------------
    def madness_stdout(self, mol: str) -> Optional[Path]:
        # if you later adopt a naming convention, we can tighten this
        out = self.find_latest(mol, "madness", "*.out")
        return out

    def madness_calc_info_json(self, mol: str) -> Optional[Path]:
        p = self.madness_dir(mol) / "mad.raman.calc_info.json"
        print(p)
        return p if p.exists() else None
