"""Draw the op-amp schematic with matplotlib -> docs/schematic.png (no external schematic lib)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT = Path(__file__).resolve().parent.parent / "docs" / "schematic.png"
fig, ax = plt.subplots(figsize=(11.5, 6.6))
ax.set_aspect("equal"); ax.axis("off")
LW, C = 1.6, "#1e293b"

def wire(*pts, **kw):
    p = np.array(pts); ax.plot(p[:, 0], p[:, 1], c=kw.get("c", C), lw=LW, solid_capstyle="round", zorder=1)
def dot(x, y): ax.add_patch(Circle((x, y), .075, color=C, zorder=3))

def mos(x, y, kind, name, flip=False, side="auto"):
    """MOSFET centred at (x,y). Gate on the left (flip=False) or right. Returns terminals (g, top, bottom)."""
    s = -1 if flip else 1
    X = lambda dx: x + s * dx
    wire((X(-.3), y - .7), (X(-.3), y + .7))                       # channel
    wire((X(-.3), y + .6), (X(.3), y + .6), (X(.3), y + 1.2))       # top terminal
    wire((X(-.3), y - .6), (X(.3), y - .6), (X(.3), y - 1.2))       # bottom terminal
    wire((X(-.45), y - .5), (X(-.45), y + .5))                      # gate plate
    gx = X(-1.2)
    if kind == "p":
        ax.add_patch(Circle((X(-.62), y), .13, fc="white", ec=C, lw=LW, zorder=3))
        wire((gx, y), (X(-.75), y))
    else:
        wire((gx, y), (X(-.45), y))
    ax.text(X(.55) + (.15 if not flip else -.15), y, name, ha="left" if not flip else "right", va="center", fontsize=10, fontweight="bold", color="#1a56db")
    return (gx, y), (X(.3), y + 1.2), (X(.3), y - 1.2)

def vdd(x, y, label=False):
    wire((x, y), (x, y + .35)); wire((x - .3, y + .35), (x + .3, y + .35))
    if label: ax.text(x + .4, y + .35, "VDD", fontsize=9, va="center")
def gnd(x, y):
    wire((x, y), (x, y - .15))
    for i, w in enumerate((.32, .2, .08)): wire((x - w, y - .15 - .12 * i), (x + w, y - .15 - .12 * i))
def cap_h(x, y, label, lp="above"):
    wire((x - .9, y), (x - .1, y)); wire((x + .1, y), (x + .9, y))
    wire((x - .1, y - .3), (x - .1, y + .3)); wire((x + .1, y - .3), (x + .1, y + .3))
    ax.text(x, y + .45, label, ha="center", fontsize=9)
def res_h(x, y, label):
    xs = np.linspace(x - .9, x + .9, 9); ys = [y, y + .22, y - .22, y + .22, y - .22, y + .22, y - .22, y + .22, y][:9]
    wire(*zip(xs, ys)); ax.text(x, y + .45, label, ha="center", fontsize=9)
def cap_v(x, y, label):
    wire((x, y + .9), (x, y + .1)); wire((x, y - .1), (x, y - .9))
    wire((x - .3, y + .1), (x + .3, y + .1)); wire((x - .3, y - .1), (x + .3, y - .1))
    ax.text(x + .45, y, label, fontsize=9, va="center")

# ============================== FIRST STAGE ===================================
g3, s3, d3 = mos(3.0, 7.5, "p", "M3", flip=True)       # diode-connected PMOS
g4, s4, d4 = mos(6.0, 7.5, "p", "M4")
vdd(*s3, True); vdd(*s4)
wire(g3, g4); dot(4.5, 7.5); wire((4.5, 7.5), (4.5, 6.3), d3)           # mirror gate bus + diode connection
g1, t1, b1 = mos(2.4, 4.2, "n", "M1")                                    # inverting input
g2, t2, b2 = mos(6.6, 4.2, "n", "M2", flip=True)                         # non-inverting input
wire(d3, t1); wire(d4, t2)
wire(b1, b2); dot(4.5, 3.0)
g5, t5, b5 = mos(4.2, 1.2, "n", "M5"); wire((4.5, 3.0), t5); gnd(*b5)
wire(g1, (.5, 4.2)); ax.text(.4, 4.2, "Vin−", ha="right", va="center", fontsize=10)
wire(g2, (8.2, 4.2)); ax.text(8.3, 4.2, "Vin+", ha="left", va="center", fontsize=10)

# ============================== BIAS ==========================================
g8, t8, b8 = mos(-1.6, 1.2, "n", "M8", flip=True)
wire(t8, (t8[0], 2.5), (-.9, 2.5), (-.9, 1.2)); dot(-.9, 1.2); wire(g8, g5); gnd(*b8)
# current source Iref
ax.add_patch(Circle((t8[0], 4.6), .5, fc="white", ec=C, lw=LW, zorder=3))
wire((t8[0], 4.1), (t8[0], 4.95)); ax.annotate("", xy=(t8[0], 4.05), xytext=(t8[0], 4.6), arrowprops=dict(arrowstyle="-|>", color=C, lw=1.3))
wire(t8, (t8[0], 4.1)); wire((t8[0], 5.1), (t8[0], 6.2)); vdd(t8[0], 6.2, True)
ax.text(t8[0] - .7, 4.6, "Iref\n20 µA", ha="right", va="center", fontsize=9)

# ============================== SECOND STAGE ==================================
g6, s6, d6 = mos(10.5, 7.5, "p", "M6"); vdd(*s6)
g7, t7, b7 = mos(10.5, 3.6, "n", "M7"); gnd(*b7)
wire(d6, t7)
out = (d6[0], (d6[1] + t7[1]) / 2); dot(*out)
# d2 -> M6 gate
d2 = (d4[0], 5.85); dot(*d2)
wire(d2, (8.7, 5.85), (8.7, 7.5), g6); dot(8.7, 7.5)
# bias to M7 gate
wire((-.9, 1.2), (-.9, -1.0), (9.3, -1.0), g7[::1] and (9.3, 3.6)); wire((9.3, 3.6), g7)
# output + load
wire(out, (13.3, out[1])); dot(13.3, out[1]); wire((13.3, out[1]), (14.3, out[1]))
ax.add_patch(Circle((14.3, out[1]), .1, fc="white", ec=C, lw=LW, zorder=3)); ax.text(14.55, out[1], "Vout", va="center", fontsize=10)
cap_v(13.3, 4.0, "CL\n2 pF"); wire((13.3, out[1]), (13.3, 4.9)); gnd(13.3, 3.1)
# Miller compensation path: M6 gate -> Cc -> Rz -> output (Rz right end lands on the Vout node at x=13.3)
wire((8.7, 7.5), (8.7, 9.9), (9.4, 9.9)); cap_h(10.3, 9.9, "Cc  1 pF")
wire((11.2, 9.9), (11.5, 9.9)); res_h(12.4, 9.9, "Rz  600 Ω")
wire((13.3, 9.9), (13.3, out[1]))

# annotations
ax.text(4.5, 9.35, "Stage 1: NMOS diff pair + PMOS mirror load", ha="center", fontsize=9, style="italic", color="#64748b")
ax.text(11.2, 11.0, "Miller compensation (Cc, Rz)", ha="center", fontsize=9, style="italic", color="#64748b")
ax.text(10.5, -1.6, "Stage 2: PMOS common-source + NMOS current-source load", ha="center", fontsize=9, style="italic", color="#64748b")
ax.text(-1.6, -1.6, "Bias mirror", ha="center", fontsize=9, style="italic", color="#64748b")
ax.set_xlim(-4.2, 15.8); ax.set_ylim(-2.2, 11.6)
fig.tight_layout(); fig.savefig(OUT, dpi=170, facecolor="white"); print("saved", OUT)
