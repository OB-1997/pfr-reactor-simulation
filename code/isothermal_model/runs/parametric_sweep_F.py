"""Side-injection flow-rate parametric sweep (thesis §4.2.4).

Three runs at base T=1173.15 K and base u=0.516 m/s, varying the global
F_scale knob added to params.f90 (Saša 2026-05-21 ask). Holds the H2O/O2
split ratio constant by multiplying every active injector uniformly.

Usage:
    python3 runs/parametric_sweep_F.py
"""

import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[parametric_sweep_F]"
NML_DIR = ROOT / "runs" / BUCKET
DEST_DIR = ROOT / "output" / BUCKET
PROTOTYPE = ROOT / "prototype"

CONFIGS = ["F050", "F100", "F150"]


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
    candidates[-1].rename(dest)
    print(f"  {name:<6} → output/{BUCKET}/{name}/   ({dt:.1f} s)")


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
