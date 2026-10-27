"""Monte Carlo mismatch: per-device threshold (Pelgrom A_VT) and beta (A_beta) variation -> input offset."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from spice import *

N, SEED = 300, 2027
AVT = 5.0      # mV*um   (representative 0.18um-class)
ABETA = 1.0    # %*um
VDD = 1.8
rng = np.random.default_rng(SEED)

# (name, drain, gate, source, bulk, type, W_um, L_um)   -- identical to netlists/opamp.cir
DEV = [("M8", "ibias", "ibias", "0", "0", "n", 5, .5), ("M5", "tail", "ibias", "0", "0", "n", 10, .5),
       ("M1", "d1", "inn", "tail", "0", "n", 5, .5), ("M2", "d2", "inp", "tail", "0", "n", 5, .5),
       ("M3", "d1", "d1", "vdd", "vdd", "p", 7, .5), ("M4", "d2", "d1", "vdd", "vdd", "p", 7, .5),
       ("M6", "out", "d2", "vdd", "vdd", "p", 70, .5), ("M7", "out", "ibias", "0", "0", "n", 50, .5)]
BASE = dict(n=dict(vto=0.50, kp=200e-6, lam=0.08, kf=4e-31), p=dict(vto=-0.50, kp=70e-6, lam=0.10, kf=4e-32))

def netlist(sample):
    L = ["* MC sample", f"Vdd vdd 0 {VDD}", "Iref vdd ibias 20u", "Vcm cm 0 0.9", "Vin inp cm dc 0", "Lfb out inn 1G", "CL out 0 2p",
         "Cc d2 cz 1p", "Rz cz out 600"]
    for (nm, d, g, s, b, t, w, l), (dvt, dbe) in zip(DEV, sample):
        p = BASE[t]; typ = "nmos" if t == "n" else "pmos"
        # |VT| shifts by dvt; beta scales by (1+dbe)
        vto = p["vto"] + (dvt if t == "n" else -dvt)
        L.append(f".model mod_{nm} {typ} level=1 vto={vto:.6f} kp={p['kp']*(1+dbe):.6e} gamma=0.4 phi=0.8 lambda={p['lam']} tox=4e-9 "
                 f"cgso=3e-10 cgdo=3e-10 cj=1e-3 cjsw=2e-10 mj=0.5 mjsw=0.33 pb=0.9 kf={p['kf']} af=1")
        L.append(f"{nm} {d} {g} {s} {b} mod_{nm} w={w}u l={l}u")
    L += [".control", "op", "let vos = v(inn) - v(inp)", "print vos", ".endc", ".end"]
    return "\n".join(L)

vos = []
for _ in range(N):
    sample = []
    for (nm, *_r, t, w, l) in DEV:
        area = w * l
        sample.append((rng.normal(0, AVT * 1e-3 / np.sqrt(2 * area)), rng.normal(0, ABETA / 100 / np.sqrt(2 * area))))
    (NETLISTS / "_mc.cir").write_text(netlist(sample))
    vos.append(scalar(run("_mc.cir"), "vos") * 1e3)       # mV
(NETLISTS / "_mc.cir").unlink()
vos = np.array(vos); mu, sd = vos.mean(), vos.std(ddof=1)

# analytic hand estimate: pair VT mismatch + load mirror VT mismatch referred to input via gm3/gm1
out = run("tb_ac_open_loop.cir")
gm1, gm3 = scalar(out, "@m.xu1.m1[gm]"), scalar(out, "@m.xu1.m3[gm]")
s_pair = AVT / np.sqrt(5 * .5)                 # sigma of delta-VT for the M1/M2 pair (mV)
s_load = AVT / np.sqrt(7 * .5)                 # M3/M4
hand = np.sqrt(s_pair ** 2 + (gm3 / gm1 * s_load) ** 2)

res = dict(samples=N, seed=SEED, mean_mV=round(mu, 3), sigma_mV=round(sd, 3), three_sigma_mV=round(3 * sd, 3),
           min_mV=round(vos.min(), 2), max_mV=round(vos.max(), 2), hand_calc_sigma_mV=round(hand, 3),
           AVT_mV_um=AVT, Abeta_pct_um=ABETA)
np.savetxt(RESULTS / "montecarlo_offsets_mV.csv", vos, header="offset_mV", fmt="%.4f")
json.dump(res, open(RESULTS / "montecarlo_summary.json", "w"), indent=2); print(res)

fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.hist(vos, bins=30, density=True, alpha=.7, label=f"{N} MC runs")
xs = np.linspace(vos.min(), vos.max(), 200)
ax.plot(xs, np.exp(-(xs - mu) ** 2 / (2 * sd ** 2)) / (sd * np.sqrt(2 * np.pi)), "r", label=f"Gaussian fit σ = {sd:.2f} mV")
ax.set_xlabel("Input-referred offset (mV)"); ax.set_ylabel("Probability density")
ax.set_title(f"Monte Carlo input offset: μ = {mu:.2f} mV, σ = {sd:.2f} mV\n(hand-calc σ ≈ {hand:.2f} mV)", fontsize=10); ax.legend()
fig.tight_layout(); fig.savefig(RESULTS / "montecarlo_offset.png"); plt.close(fig)
