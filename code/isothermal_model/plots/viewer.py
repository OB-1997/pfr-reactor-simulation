#!/usr/bin/env python3
"""
viewer.py — generic per-run inspector for the isothermal reactor model.

Writes three PNGs per run dir into plots/[viewer_output]/<bucket>/<run_name>/:
  - profile_final.png : spatial profiles at the final time step
  - time_midpoint.png : time series at the reactor midpoint
  - heatmaps.png      : time-space heatmaps for the key species
  - soot.png          : soot (species C) close-up

By default, scans every run dir under output/<bucket>/<run_name>/ across all
buckets. Pass a single run dir to scope to one run:
    python3 viewer.py output/[wall_time_data]/odepack_3
"""

from pathlib import Path
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

N_SPECIES = 10

SPECIES = ["CH4", "C6H6", "C9H10", "O2", "CO", "H2", "H2O", "CO2", "C", "N2"]
DISPLAY = {
    "CH4":   r"CH$_4$",   "C6H6":  r"C$_6$H$_6$",  "C9H10": r"C$_9$H$_{10}$",
    "O2":    r"O$_2$",    "CO":    "CO",           "H2":    r"H$_2$",
    "H2O":   r"H$_2$O",   "CO2":   r"CO$_2$",      "C":     "Soot (C)",
    "N2":    r"N$_2$",
}
KEY_SPECIES = ["CH4", "C6H6", "O2", "H2", "CO", "H2O"]

SOLVER_LABEL = {
    "explicit_euler": "Forward Euler",
    "implicit_euler": "Backward Euler",
    "odepack":        "DLSODE (BDF)",
}

ROOT  = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = ROOT / "output"
VIEWER_OUT  = ROOT / "plots" / "[viewer_output]"


def parse_params(params_path):
    """Pull the few scalars we need (N_cells, L, T_react) from parameters.txt."""
    text = params_path.read_text()

    def grab(pattern, cast):
        m = re.search(pattern, text)
        if not m:
            raise ValueError(f"could not find {pattern!r} in {params_path}")
        return cast(m.group(1))

    return {
        "N_cells": grab(r"N_cells\s*=\s*(\d+)", int),
        "L":       grab(r"L_reactor\s*=\s*([0-9.+\-eE]+)", float),
        "T_react": grab(r"T_react\s*=\s*([0-9.+\-eE]+)", float),
    }


def split_run_name(name):
    """`explicit_euler_3` -> (`explicit_euler`, 3); `odepack_5` -> (`odepack`, 5)."""
    m = re.match(r"^(.*)_(\d+)$", name)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def resolve_plot_dir(run_dir):
    """Map output/<bucket>/<run_name>/ → plots/[viewer_output]/<bucket>/<run_name>/."""
    try:
        rel = run_dir.relative_to(OUTPUT_ROOT)
    except ValueError:
        return run_dir  # unknown layout — write next to output.txt
    return VIEWER_OUT / rel


def list_runs():
    """All run dirs across all buckets in output/ that have output.txt."""
    return sorted(
        p for p in OUTPUT_ROOT.glob("*/*")
        if p.is_dir() and (p / "output.txt").exists()
    )


def load(run_dir, params):
    output = run_dir / "output.txt"
    if not output.exists():
        raise FileNotFoundError(f"No output.txt at {output}.")
    data = np.loadtxt(output)
    t = data[:, 0]
    n_cells = params["N_cells"]
    expected = 1 + n_cells * N_SPECIES
    if data.shape[1] != expected:
        raise ValueError(
            f"{output}: expected {expected} columns (1 + {n_cells}*{N_SPECIES}), "
            f"got {data.shape[1]}"
        )
    C = data[:, 1:].reshape(len(t), n_cells, N_SPECIES)
    L = params["L"]
    x = np.linspace(0.5 * L / n_cells, L - 0.5 * L / n_cells, n_cells)
    return t, x, C


def plot_final_profile(t, x, C, params, solver, plots_dir):
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    fig.suptitle(
        f"Spatial profiles at t = {t[-1]:g} s   "
        f"(T = {params['T_react']:g} K, N = {params['N_cells']}, {SOLVER_LABEL.get(solver, solver)})",
        fontsize=13, fontweight="bold",
    )
    for k, sp in enumerate(KEY_SPECIES):
        i  = SPECIES.index(sp)
        ax = axes[k // 2, k % 2]
        ax.plot(x, C[-1, :, i], color=f"C{i}", linewidth=1.5)
        ax.set_title(DISPLAY[sp])
        ax.set_ylabel(r"[mol/m$^3$]")
        ax.set_xlim(0.0, params["L"])
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    axes[-1, 0].set_xlabel("x [m]")
    axes[-1, 1].set_xlabel("x [m]")
    plt.tight_layout()
    fig.savefig(plots_dir / "profile_final.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_midpoint_time(t, x, C, params, solver, plots_dir):
    mid = params["N_cells"] // 2
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(
        f"Time evolution at reactor midpoint (x = {x[mid]:g} m, "
        f"T = {params['T_react']:g} K, N = {params['N_cells']}, {SOLVER_LABEL.get(solver, solver)})",
        fontsize=12, fontweight="bold",
    )
    for sp in KEY_SPECIES:
        i = SPECIES.index(sp)
        ax.plot(t, C[:, mid, i], color=f"C{i}", linewidth=1.5, label=DISPLAY[sp])
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Concentration [mol/m$^3$]")
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, alpha=0.3)
    ax.legend()
    plt.tight_layout()
    fig.savefig(plots_dir / "time_midpoint.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_heatmaps(t, x, C, params, solver, plots_dir):
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    fig.suptitle(
        f"Time-space heatmaps (T = {params['T_react']:g} K, "
        f"N = {params['N_cells']}, {SOLVER_LABEL.get(solver, solver)})",
        fontsize=13, fontweight="bold",
    )
    for k, sp in enumerate(KEY_SPECIES):
        i  = SPECIES.index(sp)
        ax = axes[k // 3, k % 3]
        im = ax.pcolormesh(x, t, C[:, :, i], cmap="inferno", shading="auto")
        ax.set_title(DISPLAY[sp])
        ax.set_xlabel("x [m]")
        if k % 3 == 0:
            ax.set_ylabel("t [s]")
        plt.colorbar(im, ax=ax, pad=0.02).set_label(r"[mol/m$^3$]")
    plt.tight_layout()
    fig.savefig(plots_dir / "heatmaps.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_soot(t, x, C, params, solver, plots_dir):
    i = SPECIES.index("C")
    mid = params["N_cells"] // 2
    color = "0.15"

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    fig.suptitle(
        f"Soot (C) — T = {params['T_react']:g} K, N = {params['N_cells']}, "
        f"{SOLVER_LABEL.get(solver, solver)}",
        fontsize=13, fontweight="bold",
    )

    ax = axes[0]
    ax.plot(x, C[-1, :, i], color=color, linewidth=1.8)
    ax.fill_between(x, 0.0, C[-1, :, i], color=color, alpha=0.15)
    ax.set_title(f"Spatial profile at t = {t[-1]:g} s")
    ax.set_xlabel("x [m]")
    ax.set_ylabel(r"[mol/m$^3$]")
    ax.set_xlim(0.0, params["L"])
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))

    ax = axes[1]
    ax.plot(t, C[:, mid, i], color=color, linewidth=1.8)
    ax.set_title(f"Time evolution at midpoint (x = {x[mid]:g} m)")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"[mol/m$^3$]")
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))

    ax = axes[2]
    im = ax.pcolormesh(x, t, C[:, :, i], cmap="inferno", shading="auto")
    ax.set_title("Time-space heatmap")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("t [s]")
    plt.colorbar(im, ax=ax, pad=0.02).set_label(r"[mol/m$^3$]")

    plt.tight_layout()
    fig.savefig(plots_dir / "soot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_one(run_dir):
    solver, idx = split_run_name(run_dir.name)
    params = parse_params(run_dir / "parameters.txt")
    t, x, C = load(run_dir, params)
    plots_dir = resolve_plot_dir(run_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    print(f"  {run_dir.name}: {len(t)} snapshots, N={params['N_cells']} -> {plots_dir.relative_to(ROOT)}")
    plot_final_profile(t, x, C, params, solver, plots_dir)
    plot_midpoint_time(t, x, C, params, solver, plots_dir)
    plot_heatmaps(t, x, C, params, solver, plots_dir)
    plot_soot(t, x, C, params, solver, plots_dir)


def main():
    if len(sys.argv) >= 2:
        candidate = Path(sys.argv[1])
        if not candidate.is_absolute():
            candidate = ROOT / candidate
        plot_one(candidate)
        return

    runs = list_runs()
    if not runs:
        raise SystemExit("No run dirs found under output/.")
    print(f"  Plotting {len(runs)} runs from {ROOT / 'output'}")
    for run_dir in runs:
        plot_one(run_dir)
    print(f"  All plots written under {PLOTS_ROOT}")


if __name__ == "__main__":
    main()
