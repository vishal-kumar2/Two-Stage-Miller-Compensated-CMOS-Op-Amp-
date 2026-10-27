# Design notes — two-stage Miller-compensated CMOS op-amp

Everything here can be derived by hand and is cross-checked against the ngspice results.
Process: generic Level-1 models (`netlists/models.lib`): µnCox = 200 µA/V², µpCox = 70 µA/V², |VT| = 0.5 V, λn = 0.08, λp = 0.10 V⁻¹, VDD = 1.8 V.

## 1. Specifications and architecture

| Target | Value |
|---|---|
| Supply | 1.8 V |
| Load | 2 pF |
| DC gain | > 70 dB |
| Phase margin | > 60° (all PVT corners) |
| Unity-gain frequency | ~40 MHz |
| Power | < 0.5 mW |

**Topology:** NMOS differential pair (M1, M2) with PMOS current-mirror load (M3, M4) and NMOS tail source (M5), followed by a PMOS common-source second stage (M6) with an NMOS current-source load (M7). Miller capacitor Cc with nulling resistor Rz. Bias: M8 mirrors a 20 µA reference into M5 (×2 → 40 µA) and M7 (×10 → 200 µA).

Why two stages: the first stage alone gives only A1 = gm1·(ro2‖ro4) ≈ 78 V/V ≈ 38 dB with these λ values. Cascoding would raise that but costs headroom at 1.8 V; a second common-source stage adds ≈ 35 dB and a near-rail-to-rail output swing, at the cost of needing frequency compensation.

## 2. Sizing procedure

| Device | W/L (µm) | ID (µA) | Role |
|---|---|---|---|
| M8 | 5 / 0.5 | 20 | bias diode |
| M5 | 10 / 0.5 | 40 | tail (2× M8) |
| M1, M2 | 5 / 0.5 | 20 | input pair |
| M3, M4 | 7 / 0.5 | 20 | mirror load |
| M6 | 70 / 0.5 | 200 | 2nd stage |
| M7 | 50 / 0.5 | 200 | 2nd-stage load (10× M8) |

1. **Pick the first-stage gm from the target bandwidth.** GBW = gm1 / (2π·Cc). With Cc = 1 pF and gm1 ≈ 280 µS → GBW ≈ 45 MHz.
2. **Pick the first-stage current from slew rate.** SR = I_tail / Cc = 40 µA / 1 pF = 40 V/µs.
3. **Size M1/M2 for gm1:** gm = √(2·µnCox·(W/L)·ID) → W/L = gm² / (2·µnCox·ID) = (283 µ)² / (2·200 µ·20 µ) = 10.
4. **Place the second pole beyond GBW.** p2 ≈ gm6 / CL. Rule of thumb for ~60° PM: gm6 ≳ 2.2·gm1·(CL/Cc) → gm6 ≳ 1.2 mS. Choosing ID6 = 200 µA gives gm6 ≈ 2 mS (p2 ≈ 3.5× GBW).
5. **Systematic-offset condition.** For zero systematic offset, (W/L)6 / (W/L)4 = 2·(I6 / I5): (70/0.5) / (7/0.5) = 10 = 2·(200 µA / 40 µA), so M6 carries the same current density as M4 and the first-stage output sits at the voltage that biases M6 correctly. Result: simulated systematic offset only **−36 µV**.
6. **Nulling resistor.** The Miller path creates a right-half-plane zero at z = gm6/Cc. Adding Rz in series moves it to z = 1/(Cc·(1/gm6 − Rz)); choosing Rz = 600 Ω > 1/gm6 ≈ 500 Ω turns it into a **left-half-plane zero** (~1.7 GHz) that adds a few degrees of phase lead.

## 3. Hand calculation vs simulation

| Quantity | Hand calc | ngspice (TT, 27 °C) | Comment |
|---|---|---|---|
| gm1 | 283 µS | 288 µS | sim includes (1+λVds) |
| gm6 | 1.98 mS | 2.09 mS | same |
| DC gain A0 | 72.7 dB | 74.7 dB | A1 = gm1(ro2‖ro4) ≈ 78.6, A2 = gm6(ro6‖ro7) ≈ 55.0; ro = 1/(λ·ID) ignores Vds dependence |
| Dominant pole / −3 dB freq | 10.4 kHz | 8.0 kHz | f = GBW/A0 |
| Unity-gain frequency | 45 MHz | 41.3 MHz | extra poles reduce it |
| Second pole | gm6/(2π·CL) ≈ 157 MHz | — | |
| Phase margin | ≈ 75° | 70.3° | parasitic poles (Cgs6, Rz·Cgd) not in hand calc |
| Slew rate | 40 V/µs | +38.9 / −35.5 V/µs | |
| Power | 468 µW | 473 µW | 1.8 V × (20 + 40 + 200) µA |
| Thermal noise (above 1/f corner) | 11.4 nV/√Hz | 11.7 nV/√Hz | √(2·4kTγ/gm1·(1 + gm3/gm1)), γ = 2/3 |
| σ(offset), mismatch | 3.67 mV | 3.65 mV (300 MC runs) | σ² = σ²(ΔVT12) + (gm3/gm1)²·σ²(ΔVT34), AVT = 5 mV·µm |

The agreement within ~2 dB / ~10 % everywhere is expected: hand analysis uses square-law, constant-λ small-signal models.

## 4. Verification methodology (and traps I hit)

* **Open-loop AC with DC feedback.** A 1 GH inductor from output to the inverting input closes the loop at DC (so the amplifier sits at its true operating point and the offset is nulled) while leaving the AC loop open. A 1 GF capacitor grounds the inverting input at AC.
* **CMRR:** both inputs driven by the same AC source (the 1 GF capacitor couples the AC to the inverting input), Acm = vout. CMRR = Adm / Acm.
* **PSRR:** AC on VDD, inputs at AC ground, PSRR = Adm / Avdd. Uses an ideal bias current source, so this is the amplifier's intrinsic PSRR, not including a supply-dependent bias generator.
* **Slew/settling:** unity-gain follower. **Step levels must stay inside the ICMR.** My first slew testbench stepped between 0.4 V and 1.4 V; 0.4 V is below the lower input limit (~0.75 V), so the input pair shut off, the output stuck near ground, and the result was meaningless. Caught by looking at the waveform, fixed by using 0.7 V → 1.4 V.
* **ICMR definition.** Measured as the range of input common-mode where the unity-gain frequency stays within 10 % of its mid-supply value (i.e. gm1 stays constant, all devices saturated). A "follower error < 1 mV" or "gain within 3 dB" criterion is misleading: loop gain stays high even as gm1 collapses, and gm1/gds actually *rises* when the tail current falls. Result 0.75–1.6 V, matching the analytic limits: Vin,min ≈ VGS1 + VDSAT5 ≈ 0.77 V and Vin,max ≈ VDD − |VSG3| + VTn ≈ 1.6 V.
* **Output swing** is the range where the small-signal gain stays above 0.707 × A0: 0.15–1.59 V (≈ VDSAT7 to VDD − VDSAT6).

## 5. Known limitations

* **Models are generic Level-1, not a foundry PDK.** Absolute numbers (gain, speed) are not predictions of silicon; trends, trade-offs and hand-calc agreement are what this project demonstrates. See "Next steps" in the README.
* **1/f noise coefficient (`KF`) is chosen empirically** to give a ~8 kHz flicker corner. The thermal floor is model-independent and verified; the 1/f level is representative only.
* **No layout parasitics, no wire RC, ideal bias current source.** The Iref reference is not generated on-chip.
* NMOS-input pair → input range does not reach ground. A rail-to-rail input stage (complementary pairs) would be the next step if that were required.
* Corner models vary VTO and KP only; no junction-cap or resistor-tolerance corners.
