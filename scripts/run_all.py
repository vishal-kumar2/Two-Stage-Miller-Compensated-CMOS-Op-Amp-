"""Reproduce every result in the repo:  python scripts/run_all.py   (needs ngspice on PATH)"""
import subprocess, sys
from pathlib import Path

here = Path(__file__).resolve().parent
for script in ["run_typical.py", "run_corners.py", "run_montecarlo.py", "draw_schematic.py"]:
    print(f"\n=== {script} ===", flush=True)
    r = subprocess.run([sys.executable, script], cwd=here)
    if r.returncode:
        sys.exit(f"{script} failed")
print("\nAll done. See results/ and docs/.")
