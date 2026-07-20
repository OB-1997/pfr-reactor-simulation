"""Shared matplotlib style for all bucket plot scripts in this project.

Import as:
    from _common import viridis_palette  # rcParams applied on import
"""

import matplotlib.pyplot as plt

# Thesis-wide unified plot style (matches Fig 2.1, all post-2026-05-10 figures).
plt.rcParams.update({
    "font.size":        10,
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "axes.linewidth":   0.8,
    "lines.linewidth":  1.4,
    "lines.markersize": 5,
    "grid.linestyle":   ":",
    "grid.linewidth":   0.5,
    "grid.alpha":       0.4,
    "savefig.bbox":     "tight",
})


def viridis_palette(n: int) -> list:
    """n viridis colours spaced through the perceptually-uniform middle of the map."""
    cmap = plt.get_cmap("viridis")
    return [cmap(0.15 + 0.70 * i / max(1, n - 1)) for i in range(n)]
