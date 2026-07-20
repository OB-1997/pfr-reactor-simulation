"""Parametric-sweep figures (thesis §4.2.2 / §4.2.3) for the engineering
results — 3 T × 3 u industrial operating envelope.

Reads  output/[parametric_sweep]/T<K>_u<mm-s>/{output.txt,parameters.txt}
Writes results/[parametric_sweep]/profiles_T_sweep.png
        results/[parametric_sweep]/profiles_u_sweep.png
        results/[parametric_sweep]/yields.csv
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
BUCKET = "[parametric_sweep]"
RUN_ROOT = ROOT / "output" / BUCKET
OUT_DIR  = ROOT / "results" / BUCKET

SPECIES = ["CH4", "C6H6", "C9H10", "O2", "CO", "H2", "H2O", "CO2", "C", "N2"]
DISPLAY = {
    "CH4":   r"CH$_4$",   "C6H6":  r"C$_6$H$_6$",  "C9H10": r"C$_9$H$_{10}$",
    "O2":    r"O$_2$",    "CO":    "CO",           "H2":    r"H$_2$",
    "H2O":   r"H$_2$O",   "CO2":   r"CO$_2$",      "C":     "C (soot)",
    "N2":    r"N$_2$",
}
KEY_PANEL_SPECIES = ["CO2", "C6H6", "H2", "CO"]   # 2x2 layout for the figures

# Sweep design — 3 T at u_base, 3 u at T_base (base case appears in both)
T_VALUES = [1073, 1173, 1273]   # K
U_VALUES = [0.400, 0.516, 0.700]  # m/s
T_BASE = 1173
U_BASE = 0.516

# Inlet feed mole fractions (Table 3.3 of the thesis)
Y_IN = {"CH4": 0.51382, "C9H10": 0.39711, "H2": 0.06714, "CO": 0.01141, "CO2": 0.01053}
P_REACT = 101325.0  # Pa
R_GAS   = 8.314     # J/mol/K


def bucket_name(T: int, u: float) -> str:
    return f"T{T}_u{int(round(u * 1000)):03d}"


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


def plot_sweep(values, fixed_label, fixed_value, bucket_fn, value_label,
               value_unit, palette, out_png):
    """Generic 2x2 sweep plot: one panel per KEY_PANEL_SPECIES, len(values)
    viridis lines per panel for the swept parameter."""
    runs = [load_final(RUN_ROOT / bucket_fn(v)) for v in values]

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.2), sharex=True)
    for k, sp in enumerate(KEY_PANEL_SPECIES):
        ax = axes[k // 2, k % 2]
        i = SPECIES.index(sp)
        for (params, x, C_end), colour, v in zip(runs, palette, values):
            label = f"{value_label} = {v}{value_unit}"
            ax.plot(x, C_end[:, i], color=colour, label=label, linewidth=1.5)
        ax.set_ylabel(rf"{DISPLAY[sp]} / (mol m$^{{-3}}$)")
        ax.grid(True, which="both")
    axes[1, 0].set_xlabel("axial position / m")
    axes[1, 1].set_xlabel("axial position / m")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=len(values), frameon=False)

    fig.suptitle("", y=1.0)  # no title, per convention
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_png.relative_to(ROOT)}  (fixed {fixed_label} = {fixed_value})")


def compute_yields() -> list[dict]:
    rows = []
    for T in T_VALUES:
        c_tot = P_REACT / (R_GAS * T)
        c_ch4_in = Y_IN["CH4"] * c_tot
        for u in U_VALUES:
            params, x, C_end = load_final(RUN_ROOT / bucket_name(T, u))
            out = C_end[-1, :]
            c_ch4_out = out[SPECIES.index("CH4")]
            c_c6h6_out = out[SPECIES.index("C6H6")]
            c_h2_out  = out[SPECIES.index("H2")]
            c_co_out  = out[SPECIES.index("CO")]
            c_co2_out = out[SPECIES.index("CO2")]
            c_soot_peak = C_end[:, SPECIES.index("C")].max()
            rows.append({
                "T_K":            T,
                "u_m_per_s":      u,
                "CH4_conv_pct":   100.0 * (1.0 - c_ch4_out / c_ch4_in),
                "C6H6_out":       c_c6h6_out,
                "H2_out":         c_h2_out,
                "CO_out":         c_co_out,
                "H2_to_CO":       (c_h2_out / c_co_out) if c_co_out > 0 else float("inf"),
                "CO2_out":        c_co2_out,
                "C_soot_peak":    c_soot_peak,
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
    # Figure 1 — T sweep at u_base
    plot_sweep(
        values=T_VALUES,
        fixed_label="u",
        fixed_value=f"{U_BASE} m/s",
        bucket_fn=lambda v: bucket_name(v, U_BASE),
        value_label="T",
        value_unit=" K",
        palette=viridis_palette(len(T_VALUES)),
        out_png=OUT_DIR / "profiles_T_sweep.png",
    )

    # Figure 2 — u sweep at T_base
    plot_sweep(
        values=U_VALUES,
        fixed_label="T",
        fixed_value=f"{T_BASE} K",
        bucket_fn=lambda v: bucket_name(T_BASE, v),
        value_label="u",
        value_unit=" m/s",
        palette=viridis_palette(len(U_VALUES)),
        out_png=OUT_DIR / "profiles_u_sweep.png",
    )

    # Yields table — all 9 (T, u) combinations
    rows = compute_yields()
    write_yields_csv(rows, OUT_DIR / "yields.csv")

    # Pretty print for the terminal
    print()
    print(f"{'T [K]':>6} {'u [m/s]':>8} {'CH4 conv [%]':>13} "
          f"{'C6H6 out':>10} {'H2 out':>8} {'CO out':>8} "
          f"{'H2/CO':>7} {'CO2 out':>9} {'C peak':>8}")
    print("-" * 88)
    for r in rows:
        print(f"{r['T_K']:>6} {r['u_m_per_s']:>8.3f} {r['CH4_conv_pct']:>13.3f} "
              f"{r['C6H6_out']:>10.4f} {r['H2_out']:>8.3f} {r['CO_out']:>8.3f} "
              f"{r['H2_to_CO']:>7.3f} {r['CO2_out']:>9.4f} {r['C_soot_peak']:>8.4f}")


if __name__ == "__main__":
    main()
