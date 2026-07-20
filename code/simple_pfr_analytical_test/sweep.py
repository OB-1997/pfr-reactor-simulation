"""sweep.py - Ch 4 Test 2: relative error vs N_cells for each solver.

For each (solver, N_cells) in the sweep grid, run pfr_simple to steady state,
compute the max-norm relative error of C_A vs the analytical reference, and
plot rel-err vs N on a log-log axis with one curve per solver.

Usage:
    python3 sweep.py                  # run sims + write CSV + plot
    python3 sweep.py --plot-only      # re-plot from existing CSV (no sims)

Output:
    plots/error_vs_ncells.png  - the figure
    plots/error_vs_ncells.csv  - the underlying data
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Reuse analytical formulas + run-dir reader from validate.py.
sys.path.insert(0, str(Path(__file__).parent))
from validate import (
    analytical_dirichlet_neumann,
    analytical_wehner_wilhelm,
    _read_run,
)

# === thesis-wide unified plot style (matches Fig 2.1, rtd_cascade) ===
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
    cmap = plt.get_cmap("viridis")
    return [cmap(0.15 + 0.70 * i / max(1, n - 1)) for i in range(n)]
# === end unified style ===


PROJECT = Path(__file__).resolve().parent
BIN = PROJECT / "pfr_simple"
RUNS = PROJECT / "runs"
OUTPUT = PROJECT / "output"
PLOTS = PROJECT / "plots"

SOLVERS = ["odepack", "explicit_euler", "implicit_euler"]
N_LIST = [50, 100, 200, 400, 800, 1000]

# Physical parameters (must match the namelist + params.f90 defaults)
L_REACTOR = 20.0
U_VEL = 0.516
D_AX = 1.0e-2
K_RXN = 0.026
C_A_IN = 10.4
T_END = 200.0    # ~5 tau, well past steady state


def write_namelist(n_cells: int) -> Path:
    path = RUNS / f"sweep_N{n_cells}.nml"
    path.write_text(
        f"""&reactor
    L_reactor = {L_REACTOR}
    d_reactor = 0.60
    u_vel     = {U_VEL}
    T_react   = 1173.15
    P_react   = 101325.0
    D_ax      = {D_AX:.3e}
    N_cells   = {n_cells}
/

&time_window
    t_start = 0.0
    t_end   = {T_END}
    dt_out  = 5.0
/
"""
    )
    return path


def run_one(solver: str, n_cells: int) -> Path:
    nml = write_namelist(n_cells)
    proc = subprocess.run(
        [str(BIN), f"--solver={solver}", f"--config={nml}"],
        cwd=PROJECT, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + "\n" + proc.stderr + "\n")
        raise RuntimeError(f"pfr_simple failed for solver={solver} N={n_cells}")
    dirs = sorted(
        OUTPUT.glob(f"{solver}_*"),
        key=lambda p: int(p.name.rsplit("_", 1)[-1]),
    )
    return dirs[-1]


def rel_error_vs_analytical(run_dir: Path):
    t, C, n_cells = _read_run(run_dir)
    pe = U_VEL * L_REACTOR / D_AX
    da = K_RXN * L_REACTOR / U_VEL

    j = np.arange(1, n_cells + 1)
    x_star = (j - 0.5) / n_cells

    c_a_sim = C[-1, :, 0]
    c_a_dn = C_A_IN * analytical_dirichlet_neumann(x_star, pe, da)
    c_a_ww = C_A_IN * analytical_wehner_wilhelm(x_star, pe, da)

    return {
        "rel_dn_max": float(np.max(np.abs(c_a_sim - c_a_dn) / c_a_dn)),
        "rel_ww_max": float(np.max(np.abs(c_a_sim - c_a_ww) / c_a_ww)),
    }


def reset_outputs():
    """Wipe output/<solver>_* so the auto-incremented dirs start at _1."""
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir()


def write_csv(results: dict, csv_path: Path) -> None:
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solver", "N_cells", "rel_err_dirichlet_neumann", "rel_err_wehner_wilhelm"])
        for solver in SOLVERS:
            for n in N_LIST:
                r = results[solver][n]
                w.writerow([solver, n, r["rel_dn_max"], r["rel_ww_max"]])
    print(f"  wrote {csv_path}")


def load_csv(csv_path: Path) -> dict:
    """Read existing CSV into the {solver: {N: {rel_dn_max, rel_ww_max}}} structure."""
    results = {s: {} for s in SOLVERS}
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            s = row["solver"]
            n = int(row["N_cells"])
            results[s][n] = {
                "rel_dn_max": float(row["rel_err_dirichlet_neumann"]),
                "rel_ww_max": float(row["rel_err_wehner_wilhelm"]),
            }
    return results


def plot_results(results: dict, png_path: Path) -> None:
    """Generate the main rel-err vs N_cells figure (thesis Fig 4.2)."""
    palette = viridis_palette(3)

    # Concentric marker sizes so all three solvers remain individually visible
    # despite the curves overlapping to four significant figures (the honest
    # finding: time-integration error sits below the spatial-discretisation
    # floor for every solver at every grid).
    style = {
        "odepack":        dict(marker="o", color=palette[0], ms=10, mfc="none", mew=1.6,
                               label="DLSODE (BDF, banded FD Jac)"),
        "explicit_euler": dict(marker="s", color=palette[1], ms=6,  mfc="none", mew=1.6,
                               label=r"Forward Euler ($\Delta t = 1$ ms)"),
        "implicit_euler": dict(marker="^", color=palette[2], ms=4,  mfc=palette[2], mew=1.0,
                               label=r"Backward Euler ($\Delta t = 10$ ms)"),
    }

    fig, ax = plt.subplots(figsize=(7.0, 4.0))

    ns = np.array(N_LIST, dtype=float)
    for solver in SOLVERS:
        err = np.array([results[solver][n]["rel_dn_max"] for n in N_LIST])
        ax.loglog(ns, err, ls="-", lw=1.4, **style[solver])

    # First-order reference slope anchored at the smallest N.
    n_ref = np.array([N_LIST[0], N_LIST[-1]], dtype=float)
    err_anchor = results["odepack"][N_LIST[0]]["rel_dn_max"]
    err_first = err_anchor * (N_LIST[0] / n_ref)
    ax.loglog(n_ref, err_first, color="0.3", linestyle="--", linewidth=0.9, alpha=0.75,
              label=r"first-order reference: err $\propto N^{-1}$")

    ax.set_xlabel(r"$N_\mathrm{cells}$")
    ax.set_ylabel(r"max-norm relative error")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.legend(loc="lower left", frameon=False)

    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    print(f"  wrote {png_path}")


def main():
    PLOTS.mkdir(exist_ok=True)
    csv_path = PLOTS / "error_vs_ncells.csv"
    png_path = PLOTS / "error_vs_ncells.png"

    if "--plot-only" in sys.argv:
        if not csv_path.exists():
            sys.exit(f"--plot-only set but {csv_path} does not exist")
        results = load_csv(csv_path)
        plot_results(results, png_path)
        return

    reset_outputs()
    results = {s: {} for s in SOLVERS}
    for solver in SOLVERS:
        for n in N_LIST:
            print(f"  running solver={solver:<16} N={n:>5} ...", end=" ", flush=True)
            run_dir = run_one(solver, n)
            err = rel_error_vs_analytical(run_dir)
            results[solver][n] = err
            print(f"rel_err(D/N) = {err['rel_dn_max']:.3e}   "
                  f"rel_err(W-W) = {err['rel_ww_max']:.3e}   ({run_dir.name})")

    write_csv(results, csv_path)
    plot_results(results, png_path)


if __name__ == "__main__":
    main()
