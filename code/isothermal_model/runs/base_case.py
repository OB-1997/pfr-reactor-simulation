"""Run the canonical base-case for the engineering result (thesis §4.2.1).

Builds the model if needed, runs ODEPACK at canonical conditions
(T_react = 1173.15 K, u_vel = 0.516 m/s, N_cells = 200, t_end = 50 s),
then moves the auto-incremented output into output/[base_case]/canonical/.

Usage:
    python3 runs/base_case.py
"""

import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "[base_case]"
CONFIG = ROOT / "runs" / BUCKET / "canonical.nml"
DEST = ROOT / "output" / BUCKET / "canonical"
PROTOTYPE = ROOT / "prototype"


def main() -> None:
    # Build if the binary is missing or older than any src/*.f90.
    src = list((ROOT / "src").glob("*.f90"))
    if not PROTOTYPE.exists() or any(
        s.stat().st_mtime > PROTOTYPE.stat().st_mtime for s in src
    ):
        print("→ make")
        subprocess.run(["make"], cwd=ROOT, check=True)

    # Wipe any prior canonical run.
    if DEST.exists():
        shutil.rmtree(DEST)

    print(f"→ ./prototype --solver=odepack --config={CONFIG.relative_to(ROOT)}")
    subprocess.run(
        [str(PROTOTYPE), "--solver=odepack", f"--config={CONFIG}"],
        cwd=ROOT, check=True,
    )

    # Find the most-recent odepack_N dir created at the top level.
    candidates = sorted(
        (ROOT / "output").glob("odepack_*"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise SystemExit("No new output/odepack_N/ dir found after run.")
    fresh = candidates[-1]

    print(f"→ mv {fresh.relative_to(ROOT)} {DEST.relative_to(ROOT)}")
    fresh.rename(DEST)
    print("Done.")


if __name__ == "__main__":
    main()
