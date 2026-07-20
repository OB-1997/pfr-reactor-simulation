"""Explicit-vs-implicit Euler CPU comparison figure (thesis §4.5).

Reads  results/[euler_comparison]/cpu_vs_n.csv
Writes results/[euler_comparison]/cpu_vs_n.png

One panel, log-log:
  - 5 implicit-Euler curves (Δt ∈ {5, 10, 25, 50, 100} ms, viridis-coloured)
  - 1 explicit-Euler reference (Δt = 1 ms, black dashed) for comparison.
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import viridis_palette  # noqa: E402 (applies shared rcParams)


ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[euler_comparison]"
CSV_PATH = ROOT / "results" / BUCKET / "cpu_vs_n.csv"
PNG_PATH = ROOT / "results" / BUCKET / "cpu_vs_n.png"


N_VALUES = [100, 200, 400, 800]


def read_csv() -> tuple[dict, list[tuple[str, int, int]]]:
    """Returns ({(solver, dt_ms): (N_array, wall_array_with_nan_gaps)}, [failures]).

    Failed runs (NST = 0, Newton diverged on step 1) are kept as NaN in the
    wall_time array at their N position. This makes matplotlib break the line
    at the gap instead of drawing a misleading interpolation across it."""
    # solver+dt -> {N -> wall_time or None}
    bins: dict[tuple[str, int], dict[int, float | None]] = {}
    failures: list[tuple[str, int, int]] = []
    with open(CSV_PATH, newline="") as f:
        for row in csv.DictReader(f):
            solver = row["solver"]
            n = int(row["N_cells"])
            dt_ms = int(row["dt_ms"])
            nst = int(row["nst"])
            key = (solver, dt_ms)
            bins.setdefault(key, {})
            if nst == 0:
                failures.append((solver, n, dt_ms))
                bins[key][n] = None
            else:
                bins[key][n] = float(row["wall_time"])

    data = {}
    for key, by_n in bins.items():
        ns = sorted(N_VALUES)
        walls = [by_n.get(n) for n in ns]
        walls_np = np.array([np.nan if w is None else w for w in walls])
        data[key] = (np.array(ns), walls_np)
    return data, failures


def main() -> None:
    data, failures = read_csv()

    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    # Implicit Euler — viridis ordered by Δt (small Δt → cool, large Δt → warm)
    implicit_dts = sorted({k[1] for k in data if k[0] == "implicit_euler"})
    palette = viridis_palette(len(implicit_dts))
    for colour, dt_ms in zip(palette, implicit_dts):
        N, t = data[("implicit_euler", dt_ms)]
        ax.loglog(
            N, t, marker="o", color=colour,
            label=fr"BE $\Delta t = {dt_ms}$ ms",
            linewidth=1.4, markersize=5,
            markerfacecolor="white", markeredgewidth=1.3,
        )

    # Explicit Euler reference — distinct style so it reads as the baseline
    if ("explicit_euler", 1) in data:
        N, t = data[("explicit_euler", 1)]
        ax.loglog(
            N, t, marker="s", color="black", linestyle="--",
            label=r"FE $\Delta t = 1$ ms (ref)",
            linewidth=1.4, markersize=5,
            markerfacecolor="white", markeredgewidth=1.3,
        )

    # Newton-divergence failures show up as line breaks in the curves above
    # (NaN-gapped). They are not marked separately in the figure; the LaTeX
    # caption documents which (N, Δt) combinations diverged.

    ax.set_xlabel(r"$N_{\mathrm{cells}}$")
    ax.set_ylabel("wall-clock time / s")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_xticks([100, 200, 400, 800])
    ax.set_xticklabels([100, 200, 400, 800])
    ax.legend(loc="upper left", frameon=False, fontsize=7,
              ncol=2, handlelength=1.6, handletextpad=0.5,
              columnspacing=1.0, borderaxespad=0.5)

    fig.tight_layout()
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, dpi=200)
    print(f"Wrote {PNG_PATH.relative_to(ROOT)}")
    if failures:
        print(f"  ⚠ {len(failures)} run(s) failed (NST=0): "
              + ", ".join(f"{s} N={n} dt={dt}ms" for s, n, dt in failures))


if __name__ == "__main__":
    main()
