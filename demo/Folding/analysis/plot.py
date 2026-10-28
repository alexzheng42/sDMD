#!/usr/bin/env python3
"""Plot the Folding demo: hydrogen bonds, potential energy, Rg and RMSD vs time.

Run from this folder after analysis.bash. Needs numpy and matplotlib.
Writes HBNum.pdf, PE.pdf, RG.pdf and RMSD.pdf.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Cooling schedule of fold.bash: each run ends at pTime with temperature pTemp.
pTime = [20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000, 180000, 200000, 220000]
pTemp = [0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40]

SMOOTH = 1000            # running-average window (time units)
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


def smooth(t, y):
    n = max(1, int(round(SMOOTH / np.median(np.diff(t)))))
    if n == 1:
        return y
    # Divide by the number of points actually in the window so the ends are not pulled to 0.
    w = np.ones(n)
    return np.convolve(y, w, mode="same") / np.convolve(np.ones_like(y), w, mode="same")


def figure(ylabel):
    fig, ax = plt.subplots(figsize=(7, 3.8), constrained_layout=True)
    ax.set_xlim(0, pTime[-1])
    ax.set_xlabel("t (time units, 1 t.u. = 48.88 fs)")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    # Mark the temperature stages; label each stage with its T* on top.
    starts = [0] + pTime[:-1]
    for s in starts[1:]:
        ax.axvline(s, color=GRID, linewidth=0.8, zorder=0)
    top = ax.secondary_xaxis("top")
    top.set_xticks([(a + b) / 2 for a, b in zip(starts, pTime)][1::2])
    top.set_xticklabels(["%.2f" % T for T in pTemp[1::2]])
    top.tick_params(length=0, labelsize=9, colors=MUTED)
    top.set_xlabel("T* of each stage", color=MUTED, fontsize=9)
    top.spines["top"].set_visible(False)
    return fig, ax


def series(ax, t, y, color, label=None):
    ax.plot(t, y, color=color, alpha=0.18, linewidth=0.6)
    ax.plot(t, smooth(t, y), color=color, label=label)


# Hydrogen bonds by type
cols, d = load("HBInfo.txt")
fig, ax = figure("Number of hydrogen bonds")
for key, color, label in (("total", AQUA, "total"), ("alpha", BLUE, "α-helix"), ("beta", ORANGE, "β")):
    series(ax, d[:, 0], d[:, cols.index(key)], color, label)
ax.set_ylim(0, 10)
ax.legend(loc="upper left", ncol=3)
fig.savefig("HBNum.pdf")

# Potential energy
cols, d = load("Energy.txt")
fig, ax = figure("Potential energy (kcal/mol)")
series(ax, d[:, 0], d[:, 3], BLUE)
fig.savefig("PE.pdf")

# Radius of gyration
_, d = load("RG.txt")
fig, ax = figure("Rg (nm)")
series(ax, d[:, 0], d[:, 1], BLUE)
fig.savefig("RG.pdf")

# RMSD to the two reference structures (-ref ref_alpha.gro ref_beta.gro)
_, d = load("RMSD.txt")
fig, ax = figure("RMSD (nm)")
series(ax, d[:, 0], d[:, 1], BLUE, "to α-helix")
if d.shape[1] > 2:
    series(ax, d[:, 0], d[:, 2], ORANGE, "to β-hairpin")
ax.set_ylim(bottom=0)
ax.legend(loc="upper right", ncol=2)
fig.savefig("RMSD.pdf")

print("Wrote HBNum.pdf, PE.pdf, RG.pdf, RMSD.pdf")
