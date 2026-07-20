"""Collect wall-time data from output/[wall_time_data]/<solver>_<N>/performance.txt
and parameters.txt into a single CSV for plotting.

Outputs results/[wall_time_data]/walltime.csv with columns:
    solver, N_cells, wall_time, cpu_time, nst, nfe, nje, bdf_order_last, run_dir

The 18 runs themselves were generated 2026-05-10 from the namelists in
runs/[wall_time_data]/sweep_N*.nml (one namelist per N, three solvers each).
"""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[wall_time_data]"
OUTPUT_DIR = ROOT / "output" / BUCKET
CSV_PATH = ROOT / "results" / BUCKET / "walltime.csv"

PARAM_RE = re.compile(r"\s*(\w+)\s*=\s*(\S+)")
SOLVERS = ["odepack", "explicit_euler", "implicit_euler"]


def parse_kv(path: Path, key: str) -> str | None:
    """Return the first value matching `key = ...` in `path`, or None."""
    pattern = re.compile(rf"\b{re.escape(key)}\s*=\s*(\S+)")
    with open(path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                return m.group(1)
    return None


def collect(run_dir: Path) -> dict:
    """Pull the fields we need from one run folder."""
    params = run_dir / "parameters.txt"
    perf = run_dir / "performance.txt"

    n_cells = int(parse_kv(params, "N_cells"))

    wall_str = parse_kv(perf, "Wall-clock time")
    cpu_str = parse_kv(perf, "CPU time")
    nst = int(parse_kv(perf, "Internal steps (NST)"))
    nfe = int(parse_kv(perf, "RHS evaluations (NFE)"))
    nje = int(parse_kv(perf, "Jacobian evaluations (NJE)"))
    bdf = int(parse_kv(perf, "Last BDF order"))

    # Strip "s" suffix and convert
    def to_seconds(s: str) -> float:
        return float(s.rstrip("s").strip())

    return {
        "N_cells": n_cells,
        "wall_time": to_seconds(wall_str),
        "cpu_time": to_seconds(cpu_str),
        "nst": nst,
        "nfe": nfe,
        "nje": nje,
        "bdf_order_last": bdf,
    }


def main() -> None:
    rows = []
    for solver in SOLVERS:
        # Sort by run-dir trailing integer so we collect in creation order
        dirs = sorted(
            OUTPUT_DIR.glob(f"{solver}_*"),
            key=lambda p: int(p.name.split("_")[-1]),
        )
        for d in dirs:
            data = collect(d)
            data["solver"] = solver
            data["run_dir"] = f"{BUCKET}/{d.name}"
            rows.append(data)

    # Sort: solver alphabetically, then N_cells ascending
    rows.sort(key=lambda r: (r["solver"], r["N_cells"]))

    fieldnames = [
        "solver",
        "N_cells",
        "wall_time",
        "cpu_time",
        "nst",
        "nfe",
        "nje",
        "bdf_order_last",
        "run_dir",
    ]
    with open(CSV_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {len(rows)} rows to {CSV_PATH.relative_to(ROOT)}")
    print()
    print(f"{'solver':<16} {'N':>5} {'wall_t [s]':>12} {'NST':>8} {'NFE':>10} {'NJE':>8} {'BDF':>5}")
    print("-" * 70)
    for r in rows:
        print(
            f"{r['solver']:<16} {r['N_cells']:>5} {r['wall_time']:>12.4f} "
            f"{r['nst']:>8} {r['nfe']:>10} {r['nje']:>8} {r['bdf_order_last']:>5}"
        )


if __name__ == "__main__":
    main()