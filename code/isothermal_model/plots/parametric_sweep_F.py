"""Side-injection flow-rate sweep figure (thesis §4.2.4).

Reads  output/[parametric_sweep_F]/F<scale>/{output.txt,parameters.txt}
Writes results/[parametric_sweep_F]/profiles_F_sweep.png
        results/[parametric_sweep_F]/yields.csv
"""

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import viridis_palette  # noqa: E402 (applies shared rcParams)


ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[parametric_sweep_F]"
RUN_ROOT = ROOT / "output" / BUCKET
OUT_DIR  = ROOT / "results" / BUCKET

SPECIES = ["CH4", "C6H6", "C9H10", "O2", "CO", "H2", "H2O", "CO2", "C", "N2"]
DISPLAY = {
    "CH4":   r"CH$_4$",   "C6H6":  r"C$_6$H$_6$",  "C9H10": r"C$_9$H$_{10}$",
    "O2":    r"O$_2$",    "CO":    "CO",           "H2":    r"H$_2$",
    "H2O":   r"H$_2$O",   "CO2":   r"CO$_2$",      "C":     "C (soot)",
    "N2":    r"N$_2$",
}
KEY_PANEL_SPECIES = ["CO2", "C6H6", "H2", "CO"]   # matches T and u sweep figures

F_VALUES = [0.5, 1.0, 1.5]
T_BASE = 1173.15
U_BASE = 0.516

Y_IN = {"CH4": 0.51382, "C9H10": 0.39711, "H2": 0.06714, "CO": 0.01141, "CO2": 0.01053}
P_REACT = 101325.0
R_GAS   = 8.314


def bucket_name(F: float) -> str:
    return f"F{int(round(F * 100)):03d}"


def parse_params(path: Path) -> dict:
    text = path.read_text()
    def grab(pat, cast):
        m = re.search(pat, text)
        if not m:
            raise ValueError(f"missing {pat!r} in {path}")
        return cast(m.group(1))
    return {
        "N_cells": grab(r"N_cells\s*=\s*(\d+)", int),
        "L":       grab(r"L_reactor\s*=\s*([0-9.+\-eE]+)", float),
        "T_react": grab(r"T_react\s*=\s*([0-9.+\-eE]+)", float),
        "u_vel":   grab(r"u_vel\s*=\s*([0-9.+\-eE]+)", float),
    }


def load_final(run_dir: Path):
    params = parse_params(run_dir / "parameters.txt")
    data = np.loadtxt(run_dir / "output.txt")
    n = params["N_cells"]
    C_end = data[-1, 1:].reshape(n, len(SPECIES))
    L = params["L"]
    x = np.linspace(0.5 * L / n, L - 0.5 * L / n, n)
    return params, x, C_end


def plot_F_sweep(out_png: Path) -> None:
    runs = [load_final(RUN_ROOT / bucket_name(F)) for F in F_VALUES]
    palette = viridis_palette(len(F_VALUES))

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.2), sharex=True)
    for k, sp in enumerate(KEY_PANEL_SPECIES):
        ax = axes[k // 2, k % 2]
        i = SPECIES.index(sp)
        for (params, x, C_end), colour, F in zip(runs, palette, F_VALUES):
            ax.plot(x, C_end[:, i], color=colour,
                    label=rf"$F_\mathrm{{scale}}$ = {F:.1f}", linewidth=1.5)
        ax.set_ylabel(rf"{DISPLAY[sp]} / (mol m$^{{-3}}$)")
        ax.grid(True, which="both")
    axes[1, 0].set_xlabel("axial position / m")
    axes[1, 1].set_xlabel("axial position / m")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=len(F_VALUES), frameon=False)

    fig.suptitle("", y=1.0)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_png.relative_to(ROOT)}  "
          f"(T = {T_BASE} K, u = {U_BASE} m/s)")


def compute_yields() -> list[dict]:
    rows = []
    c_tot_base = P_REACT / (R_GAS * T_BASE)
    c_ch4_in = Y_IN["CH4"] * c_tot_base
    for F in F_VALUES:
        params, x, C_end = load_final(RUN_ROOT / bucket_name(F))
        out = C_end[-1, :]
        c_ch4_out  = out[SPECIES.index("CH4")]
        c_c6h6_out = out[SPECIES.index("C6H6")]
        c_h2_out   = out[SPECIES.index("H2")]
        c_co_out   = out[SPECIES.index("CO")]
        c_co2_out  = out[SPECIES.index("CO2")]
        c_soot_peak = C_end[:, SPECIES.index("C")].max()
        rows.append({
            "F_scale":      F,
            "CH4_conv_pct": 100.0 * (1.0 - c_ch4_out / c_ch4_in),
            "C6H6_out":     c_c6h6_out,
            "H2_out":       c_h2_out,
            "CO_out":       c_co_out,
            "H2_to_CO":     (c_h2_out / c_co_out) if c_co_out > 0 else float("inf"),
            "CO2_out":      c_co2_out,
            "C_soot_peak":  c_soot_peak,
        })
    return rows


def write_yields_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Wrote {path.relative_to(ROOT)}")


def main() -> None:
    plot_F_sweep(OUT_DIR / "profiles_F_sweep.png")
    rows = compute_yields()
    write_yields_csv(rows, OUT_DIR / "yields.csv")

    print()
    print(f"{'F_scale':>8} {'CH4 conv [%]':>13} {'C6H6 out':>10} {'H2 out':>8} "
          f"{'CO out':>8} {'H2/CO':>7} {'CO2 out':>9} {'C peak':>8}")
    print("-" * 80)
    for r in rows:
        print(f"{r['F_scale']:>8.2f} {r['CH4_conv_pct']:>13.3f} {r['C6H6_out']:>10.4f} "
              f"{r['H2_out']:>8.3f} {r['CO_out']:>8.3f} {r['H2_to_CO']:>7.3f} "
              f"{r['CO2_out']:>9.4f} {r['C_soot_peak']:>8.4f}")


if __name__ == "__main__":
    main()
