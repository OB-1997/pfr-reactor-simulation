"""Run the 3 T × 3 u parametric sweep on the isothermal reactor model
(thesis §4.2.2 / §4.2.3 — industrial operating envelope).

For each of the 9 namelists in runs/[parametric_sweep]/, runs ODEPACK and
moves the auto-incremented output dir into output/[parametric_sweep]/<name>/.

Usage:
    python3 runs/parametric_sweep.py
"""

import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[parametric_sweep]"
NML_DIR = ROOT / "runs" / BUCKET
DEST_DIR = ROOT / "output" / BUCKET
PROTOTYPE = ROOT / "prototype"

CONFIGS = [
    "T1073_u400", "T1073_u516", "T1073_u700",
    "T1173_u400", "T1173_u516", "T1173_u700",
    "T1273_u400", "T1273_u516", "T1273_u700",
]


def ensure_built() -> None:
    src = list((ROOT / "src").glob("*.f90"))
    if not PROTOTYPE.exists() or any(
        s.stat().st_mtime > PROTOTYPE.stat().st_mtime for s in src
    ):
        print("→ make")
        subprocess.run(["make"], cwd=ROOT, check=True)


def run_one(name: str) -> None:
    cfg = NML_DIR / f"{name}.nml"
    dest = DEST_DIR / name
    if dest.exists():
        shutil.rmtree(dest)

    t0 = time.time()
    subprocess.run(
        [str(PROTOTYPE), "--solver=odepack", f"--config={cfg}"],
        cwd=ROOT, check=True, capture_output=True,
    )
    dt = time.time() - t0

    candidates = sorted(
        (ROOT / "output").glob("odepack_*"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise SystemExit(f"No new output/odepack_*/ dir after run {name}.")
    fresh = candidates[-1]
    fresh.rename(dest)
    print(f"  {name:<14} → output/{BUCKET}/{name}/   ({dt:.1f} s)")


def main() -> None:
    ensure_built()
    DEST_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Running {len(CONFIGS)} configurations:")
    t0 = time.time()
    for name in CONFIGS:
        run_one(name)
    print(f"Done — {len(CONFIGS)} runs in {time.time() - t0:.1f} s total.")


if __name__ == "__main__":
    main()
