"""Common validation and output folders for BROADNESS figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


def _valid_name(value: object) -> str:
    name = re.sub(r"\W", "_", str(value).strip())
    if not name:
        return "Figure"
    return f"x{name}" if name[0].isdigit() else name


@dataclass(frozen=True, slots=True)
class BROADNESSFigureSettings:
    Mode: str
    Layout: str
    Formats: tuple[str, ...]
    Prefix: str
    Show: bool
    Save: bool
    MakeIndividual: bool
    MakeSummary: bool
    OutputFolder: Path | None


def BROADNESS_FigureSettings(
    figureMode: str = "off",
    figureLayout: str = "individual",
    outputPath: str | Path | None = None,
    figureFormats: str | tuple[str, ...] | list[str] = ("png",),
    figurePrefix: str = "",
    functionFolder: str = "Figures",
) -> BROADNESSFigureSettings:
    """Prepare the same mode/layout/export choices used by the MATLAB toolbox.

    The ``fig`` format stores a Python Matplotlib figure with pickle. It is
    intended for trusted local files and is not a MATLAB ``.fig`` file.
    """
    mode = str(figureMode).lower()
    layout = str(figureLayout).lower()
    if mode not in {"off", "show", "save", "both"}:
        raise ValueError("figureMode must be 'off', 'show', 'save', or 'both'")
    if layout not in {"individual", "summary", "both"}:
        raise ValueError("figureLayout must be 'individual', 'summary', or 'both'")
    raw_formats = (figureFormats,) if isinstance(figureFormats, str) else tuple(figureFormats)
    formats = tuple(dict.fromkeys(str(item).lower() for item in raw_formats))
    if not formats or any(item not in {"png", "pdf", "fig"} for item in formats):
        raise ValueError("figureFormats may contain only 'png', 'pdf', and 'fig'")
    prefix = f"{_valid_name(figurePrefix)}_" if figurePrefix else ""
    show = mode in {"show", "both"}
    save = mode in {"save", "both"}
    folder = None
    if save:
        if outputPath is None or str(outputPath) == "":
            raise ValueError("outputPath is required when figureMode is 'save' or 'both'")
        folder = Path(outputPath).expanduser() / "BROADNESS_Figures" / _valid_name(functionFolder)
        folder.mkdir(parents=True, exist_ok=True)
    return BROADNESSFigureSettings(
        Mode=mode, Layout=layout, Formats=formats, Prefix=prefix,
        Show=show, Save=save, MakeIndividual=layout in {"individual", "both"},
        MakeSummary=layout in {"summary", "both"}, OutputFolder=folder,
    )


__all__ = ["BROADNESSFigureSettings", "BROADNESS_FigureSettings"]
