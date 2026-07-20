"""Wall-time-vs-N_cells figure (thesis §4.3) — three solver curves on log-log
axes plus a slope-1 reference line.

Reads  results/[wall_time_data]/walltime.csv
Writes results/[wall_time_data]/walltime.png
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import viridis_palette  # noqa: E402  (applies shared rcParams on import)


ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[wall_time_data]"
CSV_PATH = ROOT / "results" / BUCKET / "walltime.csv"
PNG_PATH = ROOT / "results" / BUCKET / "walltime.png"

SOLVERS = ["explicit_euler", "implicit_euler", "odepack"]
LABELS = {
    "explicit_euler": r"Forward Euler ($\Delta t = 1$ ms)",
    "implicit_euler": r"Backward Euler + banded Newton ($\Delta t = 10$ ms)",
    "odepack":        r"DLSODE BDF (mf=25, adaptive)",
}
MARKERS = {"explicit_euler": "o", "implicit_euler": "s", "odepack": "^"}
PALETTE = viridis_palette(3)
COLOURS = {
    "explicit_euler": PALETTE[0],
    "implicit_euler": PALETTE[1],
    "odepack":        PALETTE[2],
}


def read_csv() -> dict:
    """Returns {solver: (N_array, wall_array)} sorted by N_cells."""
    by_solver: dict[str, list[tuple[int, float]]] = {s: [] for s in SOLVERS}
    with open(CSV_PATH, newline="") as f:
        for row in csv.DictReader(f):
            s = row["solver"]
            if s in by_solver:
                by_solver[s].append((int(row["N_cells"]), float(row["wall_time"])))
    return {
        s: (np.array([p[0] for p in sorted(pts)]),
            np.array([p[1] for p in sorted(pts)]))
        for s, pts in by_solver.items()
    }


def main() -> None:
    data = read_csv()

    fig, ax = plt.subplots(figsize=(7.0, 4.0))

    for s in SOLVERS:
        N, t = data[s]
        ax.loglog(N, t, marker=MARKERS[s], color=COLOURS[s],
                  label=LABELS[s], linewidth=1.4, markersize=5,
                  markerfacecolor="white", markeredgewidth=1.3)

    # Slope-1 reference (theoretical for fixed-dt methods on banded LU).
    # Anchor at the implicit_euler N=50 point so it sits visibly on the curves.
    N_ref = np.array([50, 1000])
    t_anchor = data["implicit_euler"][1][0]
    t_ref = t_anchor * (N_ref / 50.0)
    ax.loglog(N_ref, t_ref, linestyle="--", color="0.3", linewidth=0.9,
              label=r"slope $= 1$ (theoretical, fixed-$\Delta t$)", alpha=0.75)

    ax.set_xlabel(r"$N_{\mathrm{cells}}$")
    ax.set_ylabel("wall-clock time / s")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.legend(loc="upper left", frameon=False)
    ax.set_xticks([50, 100, 200, 400, 800, 1000])
    ax.set_xticklabels([50, 100, 200, 400, 800, 1000])

    fig.tight_layout()
    fig.savefig(PNG_PATH, dpi=200)
    print(f"Wrote {PNG_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
