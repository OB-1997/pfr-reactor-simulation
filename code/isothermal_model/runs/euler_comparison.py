"""Explicit-vs-implicit Euler CPU-time comparison (thesis §4.5).

Reference baseline: explicit Euler at its hardcoded Δt = 1 ms (CFL bound).
Sweep: implicit Euler at Δt ∈ {5, 10, 25, 50, 100} ms.
Grids: N_cells ∈ {100, 200, 400, 800}.

Per (solver, N, dt): one run; output moved to
output/[euler_comparison]/<bucket>/. Also writes
results/[euler_comparison]/cpu_vs_n.csv from the performance.txt of each run.

Usage:
    python3 runs/euler_comparison.py
"""

import csv
import re
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[euler_comparison]"
NML_DIR = ROOT / "runs" / BUCKET
DEST_DIR = ROOT / "output" / BUCKET
RESULTS_DIR = ROOT / "results" / BUCKET
PROTOTYPE = ROOT / "prototype"

N_VALUES = [100, 200, 400, 800]
IMPLICIT_DT_MS = [5, 10, 25, 50, 100]
EXPLICIT_DT_MS = 1  # hardcoded in solvers.f90; informational here


def ensure_built() -> None:
    src = list((ROOT / "src").glob("*.f90"))
    if not PROTOTYPE.exists() or any(
        s.stat().st_mtime > PROTOTYPE.stat().st_mtime for s in src
    ):
        print("→ make")
        subprocess.run(["make"], cwd=ROOT, check=True)


def run(solver: str, n: int, dt_ms: int | None, bucket_name: str) -> float:
    cfg = NML_DIR / f"N{n}.nml"
    dest = DEST_DIR / bucket_name
    if dest.exists():
        shutil.rmtree(dest)

    args = [str(PROTOTYPE), f"--solver={solver}", f"--config={cfg}"]
    if dt_ms is not None:
        args.append(f"--dt-step={dt_ms / 1000.0}")

    t0 = time.time()
    subprocess.run(args, cwd=ROOT, check=True, capture_output=True)
    dt = time.time() - t0

    # The Fortran auto-incrementer writes output/<solver>_N/ at the top level.
    candidates = sorted(
        (ROOT / "output").glob(f"{solver}_*"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise SystemExit(f"No output/{solver}_*/ after {bucket_name}.")
    candidates[-1].rename(dest)
    print(f"  {bucket_name:<26} → output/{BUCKET}/   ({dt:.1f} s wall)")
    return dt


def parse_perf(run_dir: Path) -> dict:
    """Extract wall-clock time and step counters from performance.txt."""
    text = (run_dir / "performance.txt").read_text()

    def grab(pat, cast):
        m = re.search(pat, text)
        if not m:
            return None
        return cast(m.group(1))

    return {
        "wall_time": grab(r"Wall-clock time\s*=\s*([0-9.+\-eE]+)", float),
        "cpu_time":  grab(r"CPU time\s*=\s*([0-9.+\-eE]+)", float),
        "nst":       grab(r"Internal steps \(NST\)\s*=\s*(\d+)", int),
    }


def collect_csv() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for d in sorted(DEST_DIR.iterdir()):
        if not (d / "performance.txt").exists():
            continue
        # Bucket name: explicit_N100  or  implicit_N400_dt25ms
        name = d.name
        if name.startswith("explicit_"):
            solver = "explicit_euler"
            n = int(re.search(r"N(\d+)", name).group(1))
            dt_ms = EXPLICIT_DT_MS
        elif name.startswith("implicit_"):
            solver = "implicit_euler"
            n = int(re.search(r"N(\d+)", name).group(1))
            dt_ms = int(re.search(r"dt(\d+)ms", name).group(1))
        else:
            continue
        perf = parse_perf(d)
        rows.append({
            "solver":    solver,
            "N_cells":   n,
            "dt_ms":     dt_ms,
            "wall_time": perf["wall_time"],
            "cpu_time":  perf["cpu_time"],
            "nst":       perf["nst"],
            "run_dir":   f"{BUCKET}/{name}",
        })
    rows.sort(key=lambda r: (r["solver"], r["dt_ms"], r["N_cells"]))
    csv_path = RESULTS_DIR / "cpu_vs_n.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {csv_path.relative_to(ROOT)} ({len(rows)} rows)")


def main() -> None:
    ensure_built()
    DEST_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Running 4 explicit + {4 * len(IMPLICIT_DT_MS)} implicit "
          f"= {4 + 4 * len(IMPLICIT_DT_MS)} configurations:")
    t0 = time.time()

    # Explicit Euler reference (dt=1 ms hardcoded → no override needed)
    for n in N_VALUES:
        run("explicit_euler", n, None, f"explicit_N{n}")

    # Implicit Euler sweep
    for dt_ms in IMPLICIT_DT_MS:
        for n in N_VALUES:
            run("implicit_euler", n, dt_ms, f"implicit_N{n}_dt{dt_ms}ms")

    print(f"\n{4 + 4 * len(IMPLICIT_DT_MS)} runs in {time.time() - t0:.1f} s total.")

    collect_csv()


if __name__ == "__main__":
    main()
