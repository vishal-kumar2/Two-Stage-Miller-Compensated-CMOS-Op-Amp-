"""Tiny helper layer: run an ngspice testbench and parse its outputs."""
import re, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
NETLISTS = ROOT / "netlists"
RESULTS = ROOT / "results"

def run(tb, cwd=NETLISTS):
    """Run a testbench in batch mode, return stdout (raises on ngspice errors)."""
    p = subprocess.run(["ngspice", "-b", str(tb)], cwd=cwd,
                       capture_output=True, text=True, timeout=300)
    out = p.stdout + p.stderr
    if re.search(r"fatal|no such vector|singular matrix|convergence", out, re.I):
        raise RuntimeError(f"ngspice problem in {tb}:\n{out[-800:]}")
    return out

def scalar(out, name):
    """Parse a 'name = value' line printed by ngspice."""
    m = re.search(rf"^{re.escape(name)}\s*=\s*([-+0-9.eE]+)", out, re.M)
    if not m:
        raise KeyError(name)
    return float(m.group(1))

def load_ac(path):
    """wrdata of a complex vector -> freq, complex response."""
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1] + 1j * d[:, 2]

def load_real(path, cols=1):
    d = np.loadtxt(path)
    return d[:, 0], [d[:, 1 + 2 * i] for i in range(cols)] if cols > 1 else d[:, 1]

def db(x):
    return 20 * np.log10(np.abs(x))

def crossing(x, y, level, rising=False):
    """Linear-interpolated x where y first crosses `level`."""
    s = np.sign(y - level)
    idx = np.where(np.diff(s) > 0)[0] if rising else np.where(np.diff(s) < 0)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    return x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])

def ac_metrics(f, h):
    """DC gain, UGF, phase margin, gain margin, -3 dB bandwidth from a loop response."""
    mag = db(h)
    ph = np.unwrap(np.angle(h)) * 180 / np.pi
    ph -= 360 * np.round(ph[0] / 360)          # start phase at ~0 deg
    lf = np.log10(f)
    a0 = mag[0]
    ugf = 10 ** crossing(lf, mag, 0.0)
    pm = 180 + np.interp(np.log10(ugf), lf, ph)
    f180 = 10 ** crossing(lf, ph, -180.0)
    gm = -np.interp(np.log10(f180), lf, mag) if not np.isnan(f180) else np.nan
    f3db = 10 ** crossing(lf, mag, a0 - 3.0)
    return dict(a0_db=a0, ugf_hz=ugf, pm_deg=pm, gm_db=gm, f3db_hz=f3db), mag, ph
