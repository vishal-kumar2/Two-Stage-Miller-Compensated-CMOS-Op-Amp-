"""Run every testbench at the typical corner, extract specs, make plots."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from spice import *

plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": .3,
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
CL, CC = 2e-12, 1e-12
spec = {}

# ---------------- AC open loop + operating point -------------------------------
for tb in ["tb_slew", "tb_settle", "tb_transfer", "tb_cmrr", "tb_psrr", "tb_noise"]:
    run(tb + ".cir")
out = run("tb_ac_open_loop.cir")
idd = scalar(out, "idd"); vdd = 1.8
spec["supply_current_uA"] = idd * 1e6
spec["power_uW"] = idd * vdd * 1e6
spec["offset_uV_systematic"] = scalar(out, "vos") * 1e6
gm1 = scalar(out, "@m.xu1.m1[gm]"); gm6 = scalar(out, "@m.xu1.m6[gm]")
spec["gm1_uS"] = gm1 * 1e6; spec["gm6_uS"] = gm6 * 1e6
f, h = load_ac(RESULTS / "ac_open_loop.dat")
m, mag, ph = ac_metrics(f, h)
spec.update({k: v for k, v in m.items()})

# ---------------- CMRR / PSRR ---------------------------------------------------
fc, hc = load_ac(RESULTS / "ac_cm.dat")
fp, hp = load_ac(RESULTS / "ac_psrr.dat")
cmrr = mag - db(hc); psrr = mag - db(hp)
i1k = np.argmin(abs(f - 1e3))
spec["cmrr_dc_db"] = cmrr[0]; spec["cmrr_1k_db"] = cmrr[i1k]
spec["psrr_dc_db"] = psrr[0]; spec["psrr_1k_db"] = psrr[i1k]; spec["psrr_1M_db"] = psrr[np.argmin(abs(f-1e6))]

# ---------------- Slew rate -----------------------------------------------------
t, vo = (lambda d: (d[:, 0], d[:, 3]))(np.loadtxt(RESULTS / "tran_slew.dat"))
V0, V1 = 0.7, 1.4                # step levels, both inside the ICMR
lo, hi = V0 + 0.2 * (V1 - V0), V0 + 0.8 * (V1 - V0)   # 20%..80% points
mk = (t > 1e-6) & (t < 3e-6)
tr = crossing(t[mk], vo[mk], hi, rising=True) - crossing(t[mk], vo[mk], lo, rising=True)
spec["slew_rise_V_per_us"] = (hi - lo) / tr / 1e6
mk = t > 3e-6
tf = crossing(t[mk], vo[mk], lo) - crossing(t[mk], vo[mk], hi)
spec["slew_fall_V_per_us"] = (hi - lo) / tf / 1e6
spec["slew_theory_V_per_us"] = 40e-6 / CC / 1e6

# ---------------- Settling / overshoot -----------------------------------------
d = np.loadtxt(RESULTS / "tran_settle.dat"); ts, vs = d[:, 0], d[:, 3]
final = vs[-1]; step = 0.1
over = max(0.0, (vs.max() - final) / step * 100)
band = 0.01 * step
outside = np.where(abs(vs - final) > band)[0]
spec["overshoot_pct"] = over
spec["settling_1pct_ns"] = (ts[outside[-1]] - 100e-9) * 1e9

# ---------------- Input / output range -----------------------------------------
run("tb_icmr.cir")
vcm_pts = np.array([0.30 + 0.05 * i for i in range(30)])
raw = np.loadtxt(RESULTS / "icmr_gain.dat")
blocks = raw.reshape(len(vcm_pts), -1, raw.shape[1])           # (vcm, freq, cols)
gdb, ugf_v = [], []
for blk in blocks:
    hh = blk[:, 1] + 1j * blk[:, 2]; ff = blk[:, 0]
    gdb.append(db(hh[0]))
    mg = db(hh)
    ugf_v.append(10 ** crossing(np.log10(ff), mg, 0.0) if (mg < 0).any() and mg[0] > 0 else 0.0)
gdb, ugf_v = np.array(gdb), np.array(ugf_v)
mid = int(np.argmin(abs(vcm_pts - 0.9)))
ok = ugf_v > 0.9 * ugf_v[mid]                                  # gm1 (and hence UGF) within 10% => devices saturated
lo_i = mid
while lo_i > 0 and ok[lo_i - 1]: lo_i -= 1
hi_i = mid
while hi_i < len(ok) - 1 and ok[hi_i + 1]: hi_i += 1
spec["icmr_low_V"] = vcm_pts[lo_i]; spec["icmr_high_V"] = vcm_pts[hi_i]
d = np.loadtxt(RESULTS / "dc_transfer.dat"); vd, vt = d[:, 0], d[:, 1]
gain = np.gradient(vt, vd)
hi_gain = gain > 0.707 * 10 ** (spec["a0_db"] / 20)
spec["out_swing_low_V"] = vt[hi_gain].min(); spec["out_swing_high_V"] = vt[hi_gain].max()

# ---------------- Noise ---------------------------------------------------------
d = np.loadtxt(RESULTS / "noise_in.dat"); fn, nv = d[:, 0], d[:, 1]   # V/sqrt(Hz) (inoise_spectrum)
nv_nV = nv * 1e9
spec["noise_1k_nV_rtHz"] = float(np.interp(1e3, fn, nv_nV))
spec["noise_100k_nV_rtHz"] = float(np.interp(1e5, fn, nv_nV))
white = nv_nV[fn > 1e5][:5].mean()
spec["flicker_corner_Hz"] = float(fn[np.where(nv_nV < 1.414 * white)[0][0]])
gm3 = scalar(out, "@m.xu1.m3[gm]")
k = 1.380649e-23; T = 300.15
spec["noise_thermal_hand_nV_rtHz"] = float(np.sqrt(2 * 4 * k * T * (2/3) / gm1 * (1 + gm3 / gm1)) * 1e9)
spec["noise_10M_sim_nV_rtHz"] = float(np.interp(1e7, fn, nv_nV))
band_mk = (fn >= 10) & (fn <= 1e5)
spec["noise_10Hz_100kHz_uVrms"] = float(np.sqrt(np.trapezoid(nv[band_mk] ** 2, fn[band_mk])) * 1e6)

spec = {k: (round(float(v), 4)) for k, v in spec.items()}
json.dump(spec, open(RESULTS / "typical_summary.json", "w"), indent=2)
for k, v in spec.items():
    print(f"{k:28s} {v}")

# ============================== PLOTS ============================================
# Bode
fig, ax = plt.subplots(2, 1, figsize=(6.4, 5), sharex=True)
ax[0].semilogx(f, mag, lw=1.6); ax[0].axhline(0, c="k", lw=.7)
ax[0].plot(spec["ugf_hz"], 0, "ro"); ax[0].annotate(f"UGF = {spec['ugf_hz']/1e6:.1f} MHz", (spec["ugf_hz"], 0), (spec["ugf_hz"]*.05, 25), arrowprops=dict(arrowstyle="->"))
ax[0].set_ylabel("Gain (dB)"); ax[0].set_title(f"Open-loop response  |  A0 = {spec['a0_db']:.1f} dB, PM = {spec['pm_deg']:.1f}°")
ax[1].semilogx(f, ph, lw=1.6, c="tab:orange"); ax[1].axhline(-180, c="k", lw=.7, ls="--")
ax[1].plot(spec["ugf_hz"], -180 + spec["pm_deg"], "ro")
ax[1].set_ylabel("Phase (deg)"); ax[1].set_xlabel("Frequency (Hz)"); ax[1].set_xlim(1, 1e9)
fig.tight_layout(); fig.savefig(RESULTS / "bode_open_loop.png"); plt.close(fig)

# CMRR / PSRR
fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.semilogx(f, cmrr, label="CMRR"); ax.semilogx(f, psrr, label="PSRR (VDD)")
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("dB"); ax.legend(); ax.set_title("CMRR and PSRR vs frequency")
fig.tight_layout(); fig.savefig(RESULTS / "cmrr_psrr.png"); plt.close(fig)

# Transient
fig, ax = plt.subplots(1, 2, figsize=(8, 3.3))
d = np.loadtxt(RESULTS / "tran_slew.dat")
ax[0].plot(d[:, 0] * 1e6, d[:, 1], "--", label="Vin"); ax[0].plot(d[:, 0] * 1e6, d[:, 3], label="Vout")
ax[0].set_title(f"Slew: +{spec['slew_rise_V_per_us']:.0f} / -{spec['slew_fall_V_per_us']:.0f} V/µs"); ax[0].set_xlabel("t (µs)"); ax[0].set_ylabel("V"); ax[0].legend()
d = np.loadtxt(RESULTS / "tran_settle.dat")
ax[1].plot((d[:, 0] - 100e-9) * 1e9, d[:, 3]); ax[1].set_xlim(-20, 300)
ax[1].set_title(f"100 mV step: ts(1%) = {spec['settling_1pct_ns']:.0f} ns, OS = {spec['overshoot_pct']:.1f}%"); ax[1].set_xlabel("t (ns)")
fig.tight_layout(); fig.savefig(RESULTS / "transient.png"); plt.close(fig)

# DC
fig, ax = plt.subplots(1, 2, figsize=(8, 3.3))
ax[0].plot(vcm_pts, ugf_v / 1e6, "o-", ms=3); ax[0].axvspan(spec["icmr_low_V"], spec["icmr_high_V"], color="g", alpha=.12)
ax[0].axhline(0.9 * ugf_v[mid] / 1e6, c="r", lw=.7, ls="--")
ax[0].set_title(f"UGF vs Vcm: ICMR {spec['icmr_low_V']:.2f}–{spec['icmr_high_V']:.2f} V (UGF within 10%)", fontsize=9)
ax[0].set_xlabel("Vcm (V)"); ax[0].set_ylabel("Unity-gain freq (MHz)")
d = np.loadtxt(RESULTS / "dc_transfer.dat"); ax[1].plot((d[:, 0] - 0.9) * 1e3, d[:, 1])
ax[1].set_title(f"Open-loop transfer: swing {spec['out_swing_low_V']:.2f}–{spec['out_swing_high_V']:.2f} V", fontsize=9); ax[1].set_xlabel("Vid (mV)"); ax[1].set_ylabel("Vout (V)")
fig.tight_layout(); fig.savefig(RESULTS / "dc_characteristics.png"); plt.close(fig)

# Noise
fig, ax = plt.subplots(figsize=(6.4, 3.4))
vm = fn <= 2e7   # beyond ~UGF the input-referral divides by a rolling-off gain (not meaningful)
ax.loglog(fn[vm], nv_nV[vm]); ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("nV/√Hz")
ax.set_title(f"Input-referred noise: {spec['noise_1k_nV_rtHz']:.0f} nV/√Hz @1 kHz, corner ≈ {spec['flicker_corner_Hz']:.0f} Hz")
fig.tight_layout(); fig.savefig(RESULTS / "noise.png"); plt.close(fig)
