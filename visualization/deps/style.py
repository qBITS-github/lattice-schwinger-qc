"""Shared figure style and output paths for the paper plots.

Each paper figure has a directory under ``figures/``. Running the
visualization scripts writes two files there: ``reproduction.png`` is the
paper figure drawn from this repository, and ``overlay.png`` places that
result on the published figure.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

# Okabe-Ito, one color per quench field in the order stored in the npz files.
COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
MARKERS = ("o", "s", "^", "D", "v", "P")

FIGURE_DIRS = {
    1: "fig1_condensate_and_spectrum",
    2: "fig2_critical_field",
    3: "fig3_total_charge_and_energy",
    4: "fig4_spatial_charge",
    5: "fig5_electric_field_energy",
    6: "fig6_vacuum_fidelity",
    7: "fig7_decay_rate",
    8: "fig8_site_charge",
}


def apply_style():
    """Matplotlib defaults for the comparison overlays."""
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 160,
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 11,
            "legend.fontsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.35,
            "lines.linewidth": 1.8,
        }
    )


def paper_style():
    """Matplotlib defaults for a reproduction of a paper figure."""
    apply_style()
    plt.rcParams.update(
        {
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.frameon": False,
            "lines.linewidth": 1.35,
        }
    )


# Matplotlib tab10 colors and marker shapes used by the paper, keyed by ε.
PAPER_EPS_COLOR = {
    0.0: "#7f7f7f",
    0.5: "#1f77b4",
    1.0: "#ff7f0e",
    1.5: "#2ca02c",
    2.0: "#d62728",
    2.5: "#9467bd",
    3.0: "#8c564b",
}
PAPER_EPS_MARKER = {
    0.0: "o",
    0.5: "o",
    1.0: "^",
    1.5: "s",
    2.0: "D",
    2.5: "v",
    3.0: "P",
}
# Dashed stroke of the paper's exact curves.
PAPER_DASH = (0, (4, 2))


def save_figure(fig, out_path):
    """Write ``fig`` to ``out_path`` and close it."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


def figure_png(number, name):
    """Path of one png for a paper figure."""
    return ROOT / "figures" / FIGURE_DIRS[number] / f"{name}.png"


def reproduction_png(number):
    """Path of the paper-figure reproduction."""
    return figure_png(number, "reproduction")


def drop_legacy_pngs(number):
    """Remove algorithm-named plots that these scripts no longer write."""
    for name in ("ED", "VQE", "VQD", "Trotter"):
        path = figure_png(number, name)
        if path.is_file():
            path.unlink()


def param_line(N, a, m, g):
    """Single-line parameter caption."""
    return rf"$N={N}$, $a={a:g}$, $m={m:g}$, $g={g:g}$"


def skip_missing(path, what):
    """Print a skip line and return True when ``path`` is not a file."""
    path = Path(path)
    if path.is_file():
        return False
    print(f"skip {what}: {path} is missing")
    return True
