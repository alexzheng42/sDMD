#!/usr/bin/env python3
"""Plot the REMD demo.

Run from this folder after analysis.bash. Needs numpy and matplotlib.
Writes
  FT.pdf          heat capacity and mean potential energy vs temperature (WHAM)
  HB_T.pdf        mean alpha-helix and total hydrogen bonds per temperature
  FreeEnergy.pdf  free energy -RT ln P vs RMSD to the alpha-helix at each temperature (WHAM)
  PESurface.pdf   3D surface of the mean potential energy over (RMSD alpha-helix, RMSD beta-hairpin)
"""
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

OUTPUT = "../output/"     # where out_REMD<n>.txt (temperature index vs time) is
CONFIG = "../REMDConfig.txt"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 2,
})


def load(name):
    """Return (column names, data) of an sAnalysis text file."""
    header = []
    with open(name) as f:
        for line in f:
            if line.startswith("#"):
                header = line[1:].split()
            elif line.strip() and not line.startswith("@"):
                break
    return header, np.loadtxt(name, comments=("#", "@"), ndmin=2)


def temperatures():
    with open(CONFIG) as f:
        for line in f:
            if line.split("=")[0].strip() == "Temperatures":
                return [float(x) for x in line.split("=")[1].split()]
    raise SystemExit("no Temperatures line in " + CONFIG)


def style(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.8)


temps = temperatures()

# --- Heat capacity and potential energy vs T --------------------------------
_, d = load("CvPotT.txt")
T, U, Cv = d[:, 0], d[:, 1], d[:, 2]
peak = np.argmax(Cv)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True, constrained_layout=True)
a1.plot(T, Cv, color=BLUE)
a1.plot(T[peak], Cv[peak], "o", color=BLUE, markersize=8, markeredgecolor="white", markeredgewidth=2)
a1.annotate("%.0f K" % T[peak], (T[peak], Cv[peak]), xytext=(8, 0), textcoords="offset points",
            va="center", color=INK)
a1.set_ylabel("C$_v$* (heat capacity)")
a1.set_ylim(bottom=0)
a2.plot(T, U, color=BLUE)
a2.set_ylabel("<E$_{pot}$> (kcal/mol)")
a2.set_xlabel("T (K)")
for ax in (a1, a2):
    style(ax)
fig.savefig("FT.pdf")

# --- Hydrogen bonds per temperature -------------------------------------------
# HBInfo<n>.txt follows replica n; out_REMD<n>.txt says which temperature (1..N)
# replica n had at each time. Average the second half of the run.
hbFiles = sorted(glob.glob("HBInfo*.txt"), key=lambda s: int(s[6:-4]))
if hbFiles and os.path.exists(OUTPUT + "out_REMD0.txt"):
    acc = {k: [[] for _ in temps] for k in ("alpha", "total")}
    for name in hbFiles:
        n = int(name[6:-4])
        cols, h = load(name)
        _, tr = load(OUTPUT + "out_REMD%d.txt" % n)
        h = h[h[:, 0] >= h[:, 0].max() / 2]
        idx = np.clip(np.searchsorted(tr[:, 0], h[:, 0], side="right") - 1, 0, len(tr) - 1)
        tIndex = tr[idx, 2].astype(int) - 1
        for k in acc:
            for i, v in zip(tIndex, h[:, cols.index(k)]):
                acc[k][i].append(v)
    fig, ax = plt.subplots(figsize=(6, 3.8), constrained_layout=True)
    for k, color, label in (("total", AQUA, "total"), ("alpha", BLUE, "α-helix")):
        y = [np.mean(v) if v else np.nan for v in acc[k]]
        ax.plot(temps, y, color=color, label=label, marker="o", markersize=7,
                markeredgecolor="white", markeredgewidth=2)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("Mean number of hydrogen bonds")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right")
    style(ax)
    fig.savefig("HB_T.pdf")

# --- Free energy vs RMSD to the helix -----------------------------------------
# F(RC) = -RT ln P(RC), from the WHAM probability of each RMSD bin at each temperature.
# Adjacent bins are merged (MERGE) to smooth the curves.
R, MERGE = 0.0019872, 3                             # kcal/(mol K)
name = "Probability_RC.txt" if os.path.exists("Probability_RC.txt") else "Probability.txt"
_, d = load(name)
d = d[: len(d) // MERGE * MERGE]
rc = d[:, 0].reshape(-1, MERGE).mean(1)
P = d[:, 1:].reshape(-1, MERGE, d.shape[1] - 1).sum(1)
cmap, norm = plt.get_cmap("Blues"), Normalize(temps[0] - 0.4 * (temps[-1] - temps[0]), temps[-1])
fig, ax = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
for j, t in enumerate(temps[:P.shape[1]]):
    with np.errstate(divide="ignore"):
        f = np.where(P[:, j] > 0, -R * t * np.log(P[:, j]), np.nan)
    ax.plot(rc, f - np.nanmin(f), color=cmap(norm(t)), linewidth=1.5)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cb = fig.colorbar(sm, ax=ax, pad=0.02)
cb.set_label("T (K)")
cb.ax.set_ylim(temps[0], temps[-1])
ax.set_xlabel("RMSD to α-helix (nm)")
ax.set_ylabel("Free energy (kcal/mol)")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
style(ax)
fig.savefig("FreeEnergy.pdf")

# --- Potential-energy surface -------------------------------------------------
# 3D surface of the mean potential energy over the two RMSDs, as in the original
# gnuplot figure: the 20 x 20 map is smoothed with a Gaussian weighted average
# (gnuplot "dgrid3d gauss 0.1,0.15"); bins that were never visited count as 0.
_, d = load("CombinedPotMap.txt")
x, y, pe = d[:, 0], d[:, 1], d[:, 2]
xs, ys = np.unique(x), np.unique(y)
gx, gy = np.meshgrid(xs, ys)
w = np.exp(-((gx[..., None] - x) / 0.1) ** 2 - ((gy[..., None] - y) / 0.15) ** 2)
gz = (w * pe).sum(-1) / w.sum(-1)
fig = plt.figure(figsize=(7, 5.5), constrained_layout=True)
ax = fig.add_subplot(projection="3d")
floor = gz.min() - 0.6 * (gz.max() - gz.min())
surf = ax.plot_surface(gx, gy, gz, cmap="plasma", rstride=1, cstride=1, linewidth=0.3,
                edgecolor=(1, 1, 1, 0.35), antialiased=True)
ax.contour(gx, gy, gz, levels=10, zdir="z", offset=floor, cmap="plasma", linewidths=1)
ax.set_zlim(floor, max(0, gz.max()))
ax.view_init(elev=30, azim=225)
ax.set_xlabel("RMSD to α-helix (nm)")
ax.set_ylabel("RMSD to β-hairpin (nm)")
ax.set_zlabel("")
cb = fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08)
cb.set_label("Mean potential energy (kcal/mol)")
for a in (ax.xaxis, ax.yaxis, ax.zaxis):
    a.pane.set_facecolor("white")
    a._axinfo["grid"]["color"] = GRID
fig.savefig("PESurface.pdf")

print("Wrote FT.pdf, HB_T.pdf, FreeEnergy.pdf, PESurface.pdf")
