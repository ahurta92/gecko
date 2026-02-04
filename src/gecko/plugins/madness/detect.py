from __future__ import annotations

from pathlib import Path


def _madqc_stem(calc_info_path: Path) -> str:
    name = calc_info_path.name
    suffix = ".calc_info.json"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    # Fallback: best-effort (unlikely to be used)
    return calc_info_path.stem.replace(".calc_info", "")


def can_load(path: Path) -> bool:
    """
    Detect MADNESS run directories.

    Supports:
      - MADQC style: *.calc_info.json
      - Legacy molresponse style: output.json / outputs.json
      - Optional: responses/metadata.json
    """
    if not path.exists() or not path.is_dir():
        return False

    # MADQC marker
    for p in path.glob("*.calc_info.json"):
        stem = _madqc_stem(p)
        if (path / f"{stem}.in").exists():
            return True

    # Legacy molresponse marker
    if (path / "output.json").exists() and (path / "input.json").exists():
        return True
    if (path / "outputs.json").exists() and (path / "input.json").exists():
        return True

    # Optional common marker
    if (path / "responses" / "metadata.json").exists():
        return True

    return False
