"""Engineering-result figure (thesis §4.2.1) — concentration profiles of all
ten species along the reactor at steady state, base case (T = 1173.15 K,
u = 0.516 m/s, N_cells = 200, ODEPACK).

Reads  output/[base_case]/canonical/output.txt
Writes results/[base_case]/profiles.png
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common  # noqa: F401, E402  (applies shared rcParams on import)


# Categorical palette — each species gets a distinct tab10 colour. Used for
# species-vs-x figures where the lines are categorical labels, not ordered
# levels. Viridis stays the convention for the T- and u-sweep figures where
# the lines ARE ordered.
_TAB = plt.get_cmap("tab10")
SPECIES_COLOURS = {
    "CH4":   _TAB(0),  # blue
    "C6H6":  _TAB(1),  # orange
    "C9H10": _TAB(2),  # green
    "O2":    _TAB(3),  # red
    "CO":    _TAB(4),  # purple
    "H2":    _TAB(5),  # brown
    "H2O":   _TAB(9),  # cyan (skip pink for water-conventional blue-cyan)
    "CO2":   _TAB(7),  # grey
    "C":     "black",  # soot
    "N2":    _TAB(8),  # olive
}


ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[base_case]"
RUN_DIR = ROOT / "output" / BUCKET / "canonical"
PNG_PATH = ROOT / "results" / BUCKET / "profiles.png"

N_SPECIES = 10
SPECIES = ["CH4", "C6H6", "C9H10", "O2", "CO", "H2", "H2O", "CO2", "C", "N2"]
DISPLAY = {
    "CH4":   r"CH$_4$",   "C6H6":  r"C$_6$H$_6$",  "C9H10": r"C$_9$H$_{10}$",
    "O2":    r"O$_2$",    "CO":    "CO",           "H2":    r"H$_2$",
    "H2O":   r"H$_2$O",   "CO2":   r"CO$_2$",      "C":     "C (soot)",
    "N2":    r"N$_2$",
}
# Split by concentration magnitude so the low-magnitude species (fuels,
# injected reactants, soot, inert) get their own y-scale and remain visible.
LEFT_PANEL  = ["H2", "CO", "H2O"]                                  # high magnitude (~10–55 mol/m³)
RIGHT_PANEL = ["CH4", "C6H6", "C9H10", "O2", "CO2", "C", "N2"]     # low magnitude (~0–5 mol/m³)


def parse_params(path: Path) -> dict:
    text = path.read_text()

    def grab(pattern, cast):
        m = re.search(pattern, text)
        if not m:
            raise ValueError(f"missing {pattern!r} in {path}")
        return cast(m.group(1))

    return {
        "N_cells": grab(r"N_cells\s*=\s*(\d+)", int),
        "L":       grab(r"L_reactor\s*=\s*([0-9.+\-eE]+)", float),
        "T_react": grab(r"T_react\s*=\s*([0-9.+\-eE]+)", float),
        "u_vel":   grab(r"u_vel\s*=\s*([0-9.+\-eE]+)", float),
    }


def load_final_profile(run_dir: Path):
    params = parse_params(run_dir / "parameters.txt")
    data = np.loadtxt(run_dir / "output.txt")
    n = params["N_cells"]
    if data.shape[1] != 1 + n * N_SPECIES:
        raise ValueError(
            f"unexpected output.txt shape {data.shape}, "
            f"expected (T, 1 + {n}*{N_SPECIES})"
        )
    t_end = data[-1, 0]
    C_end = data[-1, 1:].reshape(n, N_SPECIES)        # (cells, species)
    L = params["L"]
    x = np.linspace(0.5 * L / n, L - 0.5 * L / n, n)  # cell centres
    return params, t_end, x, C_end


def plot_panel(ax, x, C_end, species_list):
    for sp in species_list:
        i = SPECIES.index(sp)
        ax.plot(x, C_end[:, i], color=SPECIES_COLOURS[sp], label=DISPLAY[sp],
                linewidth=1.5)
    ax.set_xlabel("axial position / m")
    ax.set_ylabel(r"concentration / (mol m$^{-3}$)")
    ax.set_xlim(0, x[-1] + (x[1] - x[0]) / 2)
    ax.grid(True, which="both")


def main() -> None:
    params, t_end, x, C_end = load_final_profile(RUN_DIR)

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 4.2), sharey=False)
    plot_panel(axes[0], x, C_end, LEFT_PANEL)
    plot_panel(axes[1], x, C_end, RIGHT_PANEL)

    # Single combined legend below both panels, ordered LEFT_PANEL first then
    # RIGHT_PANEL so the reader scans the high-magnitude species before the
    # low-magnitude ones.
    handles_l, labels_l = axes[0].get_legend_handles_labels()
    handles_r, labels_r = axes[1].get_legend_handles_labels()
    fig.legend(handles_l + handles_r, labels_l + labels_r,
               loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=5, frameon=False)

    fig.tight_layout(rect=[0, 0.10, 1, 1])
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, dpi=200)
    print(f"Wrote {PNG_PATH.relative_to(ROOT)}")
    print(f"  T = {params['T_react']:g} K, u = {params['u_vel']:g} m/s, "
          f"N_cells = {params['N_cells']}, t_end = {t_end:g} s")


if __name__ == "__main__":
    main()
