"""Save and display Matplotlib figures with BROADNESS conventions."""

from __future__ import annotations

from io import BytesIO
import pickle
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.figure import Figure

from .BROADNESS_FigureSettings import BROADNESSFigureSettings, _valid_name


def BROADNESS_FinalizeFigure(
    fig: Figure | None, Settings: BROADNESSFigureSettings, fileName: str
) -> list[Path]:
    """Export PNG at 300 dpi, vector PDF, or a Python-native editable figure.

    ``.fig`` is a pickled Matplotlib figure, loadable only from trusted files
    with Python. It is deliberately not presented as MATLAB-compatible.
    Save-only mode closes the figure; show/both mode retains its handle.
    """
    if fig is None:
        return []
    from matplotlib import pyplot as plt
    from matplotlib.figure import Figure

    if not isinstance(fig, Figure):
        raise TypeError("fig must be a Matplotlib Figure")
    fig.patch.set_facecolor("white")
    paths: list[Path] = []
    if Settings.Save:
        assert Settings.OutputFolder is not None
        for extension in Settings.Formats:
            path = Settings.OutputFolder / f"{Settings.Prefix}{_valid_name(fileName)}.{extension}"
            if extension == "fig":
                with path.open("wb") as stream:
                    pickle.dump(fig, stream, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                fig.savefig(path, dpi=300 if extension == "png" else None,
                            bbox_inches="tight", facecolor="white")
            paths.append(path)
    if Settings.Show:
        try:
            from IPython import get_ipython
        except ImportError:
            shell = None
        else:
            shell = get_ipython()
        if shell is not None and shell.__class__.__name__ == "ZMQInteractiveShell":
            from IPython.display import Image, display

            preview = BytesIO()
            fig.savefig(preview, format="png", dpi=110, bbox_inches="tight")
            display(Image(data=preview.getvalue(), width=min(760, int(fig.get_figwidth() * 110))))
            plt.close(fig)
        elif plt.get_backend().lower() != "agg":
            fig.show()
    else:
        plt.close(fig)
    return paths


__all__ = ["BROADNESS_FinalizeFigure"]
