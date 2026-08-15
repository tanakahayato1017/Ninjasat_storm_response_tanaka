"""Shared Okabe-Ito colors and plotting defaults for repository figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DATA = ROOT / "data" / "figure_data"
FIGURES = ROOT / "figures"

C_BLUE = "#0072B2"
C_VERMILLION = "#D55E00"
C_GREEN = "#009E73"
C_SKY = "#56B4E9"
C_PURPLE = "#CC79A7"
C_ORANGE = "#E69F00"
C_GRAY = "#7F7F7F"

EVENT_COLORS = {"May": C_BLUE, "Oct": C_VERMILLION}

plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 150,
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 1.4,
    "legend.frameon": False,
})

def save(fig, stem: str) -> None:

    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        path = FIGURES / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight")
        print(f"wrote {path}")
