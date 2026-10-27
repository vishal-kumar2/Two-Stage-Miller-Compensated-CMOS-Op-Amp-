"""PVT corner sweep: 5 process x 3 supply x 3 temperature = 45 corners (open-loop AC + op)."""
import itertools, json, csv, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from spice import *

PROCESS = {  # name: (vtn, vtp, kpn, kpp)
    "TT": (0.50, -0.50, "200u", "70u"),
    "FF": (0.45, -0.45, "230u", "80u"),
    "SS": (0.55, -0.55, "170u", "60u"),
    "FS": (0.45, -0.55, "230u", "60u"),   # fast NMOS, slow PMOS
    "SF": (0.55, -0.45, "170u", "80u"),   # slow NMOS, fast PMOS
}
VDDS = [1.62, 1.8, 1.98]
TEMPS = [-40, 27, 125]
corner_file = NETLISTS / "corner.inc"
backup = corner_file.read_text()
rows = []
try:
    for (pn, (vtn, vtp, kpn, kpp)), vdd, temp in itertools.product(PROCESS.items(), VDDS, TEMPS):
        corner_file.write_text(f".param vtn={vtn} vtp={vtp} kpn={kpn} kpp={kpp}\n.param vdd_val={vdd}\n.temp {temp}\n")
        out = run("tb_ac_open_loop.cir")
        f, h = load_ac(RESULTS / "ac_open_loop.dat")
        m, _, _ = ac_metrics(f, h)
        rows.append(dict(process=pn, vdd=vdd, temp=temp, a0_db=m["a0_db"], ugf_mhz=m["ugf_hz"] / 1e6,
                         pm_deg=m["pm_deg"], idd_uA=scalar(out, "idd") * 1e6, vos_uV=scalar(out, "vos") * 1e6))
finally:
    corner_file.write_text(backup)          # always restore typical corner

with open(RESULTS / "pvt_corners.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

def col(k): return np.array([r[k] for r in rows])
summary = {}
for k, lab in [("a0_db", "A0 (dB)"), ("ugf_mhz", "UGF (MHz)"), ("pm_deg", "PM (deg)"), ("idd_uA", "Idd (uA)")]:
    v = col(k); i0, i1 = v.argmin(), v.argmax()
    summary[k] = dict(min=round(v.min(), 2), max=round(v.max(), 2),
                      min_at=f"{rows[i0]['process']}/{rows[i0]['vdd']}V/{rows[i0]['temp']}C",
                      max_at=f"{rows[i1]['process']}/{rows[i1]['vdd']}V/{rows[i1]['temp']}C")
json.dump(summary, open(RESULTS / "pvt_summary.json", "w"), indent=2)
for k, v in summary.items(): print(k, v)

# plot: PM and A0 across all 45 corners
fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
colors = {-40: "tab:blue", 27: "tab:green", 125: "tab:red"}
for ai, (k, lab) in enumerate([("pm_deg", "Phase margin (deg)"), ("ugf_mhz", "UGF (MHz)"), ("a0_db", "DC gain (dB)")]):
    for t in TEMPS:
        for vi, vdd in enumerate(VDDS):
            xs, ys = [], []
            for pi, pn in enumerate(PROCESS):
                r = [r for r in rows if r["process"] == pn and r["vdd"] == vdd and r["temp"] == t][0]
                xs.append(pi + (vi - 1) * 0.18); ys.append(r[k])
            ax[ai].scatter(xs, ys, c=colors[t], marker="o^s"[vi], s=22, label=f"{t}°C, {vdd} V" if ai == 0 else None)
    ax[ai].set_xticks(range(5)); ax[ai].set_xticklabels(PROCESS.keys()); ax[ai].set_ylabel(lab)
ax[0].axhline(60, c="k", ls="--", lw=.7); ax[0].text(0, 60.5, "PM = 60° target", fontsize=7)
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=9, fontsize=7, frameon=False)
fig.suptitle("PVT corner sweep (45 corners)", y=1.0); fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(RESULTS / "pvt_corners.png"); plt.close(fig)
