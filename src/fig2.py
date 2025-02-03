import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import omniplate as om
import pandas as pd
import seaborn as sns
import sys
import os

basedir = os.path.expanduser("~/huo2025/")

sys.path.append(basedir + "src/utils/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = basedir + "data/fig2/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

# Framework
fig, axes = plt.subplots(2, 2, figsize=(8, 8), dpi=300)
axes = axes.flat

# Fig 2A left
ax = axes[0]

cmap = plt.get_cmap("tab10")
cmap2 = plt.get_cmap("Paired")
cmap3 = plt.get_cmap("Accent")
colors = [cmap2(6), cmap(0), cmap3(7), cmap2(2)]
from matplotlib.colors import ListedColormap

cm = ListedColormap(colors, name="gr")
try:
    import matplotlib as mpl

    mpl.colormaps.register(cmap=cm)
except:
    pass
cmap = plt.get_cmap("gr")

p = om.platereader(ls=False, wdir=datadir)
p.importdf(
    [
        "Experiment_20210811_4_copy",
        "Experiment_20211001_4",
        "Experiment_20201023_1",
        "Experiment_20210922_1",
    ]
)
add_genotype(p)
sugars = ["Glu", "Fru", "Suc", "Gal"]
cond = [f"0.5% {s}, 0.5% Pal" for s in sugars]
q = "genotype =='FY4' and condition.isin(@cond).values and time < 40.0"
data = p.s.query(q)
# plot
lineplot_with_error(data=data, x="time", y="OD_mean", ax=ax, cmap="gr", alpha=0.6)
# axis-Labels
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_xlim(0, 44)
ax.set_yscale("log", base=2)
ax.set_ylim(0.25, 1.2)

# legend
legend_names = sorted(cond)
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize="small",
)
ax.tick_params(reset=True, direction="in")

# axin = ax.inset_axes([0.18, 0.6, 0.28, 0.38], alpha=1)
axin = axes[1]

axin.legend(
    handles=legend_handles,
    loc="upper right",
    fontsize="small",
)

lineplot_with_error(data=data, x="time", y="gr", ax=axin, cmap="gr", alpha=0.6)
axin.set_xlabel("time (h)")
axin.set_ylabel("growth rate (h$^{-1}$)")
axin.arrow(
    18,
    -0.05,
    0,
    0.1,
    color=cmap(2),
    width=0.4,
    head_length=0.03,
    length_includes_head=True,
    zorder=11,
)
axin.arrow(
    27,
    -0.05,
    0,
    0.1,
    color=cmap(1),
    width=0.4,
    head_length=0.03,
    length_includes_head=True,
    zorder=10,
)
axin.set_xlim(0, 44)
axin.tick_params(reset=True, direction="in")

######
p = om.platereader(ls=False, wdir=datadir)
p.importdf(
    [
        "Experiment_20200913",
        "Experiment_20201009_3",
        "Experiment_20201009_4",
    ]
)
add_maxOD_to_sc(p)
results = get_peak_valley(p, "Gal", "Pal")
# Constants
alpha = 0.9  # set unrelated alphas
filt1 = (p.s["strain"] == "FY4_1") & (p.s["condition"] == "0.2% Gal, 0.8% Pal")
ymax = 0.34
ymin = -0.02
tmin = -0.02
tmax = 32
filt_blue = (results["strain"] == "FY4_1") & (
    results["condition"] == "0.2% Gal, 0.8% Pal"
)
od_yield1 = results.loc[filt_blue, "OD at gr valley"].to_numpy()[0]

deltaOD = results.loc[filt_blue, "second OD yield"].to_numpy()[0]
od_max = results.loc[filt_blue, "max OD"].to_numpy()[0]


# First plot
ax = axes[2]
sns.lineplot(
    data=p.s.loc[filt1], x="OD_mean", y="gr", hue=None, legend=False, alpha=alpha, ax=ax
)
ax.vlines(x=od_yield1, ymin=ymin, ymax=ymax, color="k")


import matplotlib.patches as patches

rect = patches.Rectangle(
    (0, ymin),
    od_yield1,
    ymax - ymin,
    linewidth=1,
    edgecolor="none",
    facecolor="yellow",
    alpha=0.2,
)
ax.add_patch(rect)
ax.arrow(
    x=0,
    y=0.02,
    dx=od_yield1,
    dy=0,
    length_includes_head=True,
    alpha=1,
    width=0.005,
    color="k",
)
ax.arrow(
    x=od_yield1,
    y=0.02,
    dx=-od_yield1,
    dy=0,
    length_includes_head=True,
    alpha=1,
    width=0.005,
    color="k",
)
ax.text(od_yield1 / 3, 0.03, "$OD_1$", size=15)

rect2 = patches.Rectangle(
    (od_yield1, ymin),
    deltaOD,
    ymax - ymin,
    linewidth=1,
    edgecolor="none",
    facecolor="red",
    alpha=0.2,
)
ax.add_patch(rect2)
ax.arrow(
    x=od_yield1,
    y=0.02,
    dx=deltaOD,
    dy=0,
    length_includes_head=True,
    alpha=1,
    width=0.005,
    color="k",
)
ax.arrow(
    x=od_max,
    y=0.02,
    dx=-deltaOD,
    dy=0,
    length_includes_head=True,
    alpha=1,
    width=0.005,
    color="k",
)
ax.text(od_yield1 + deltaOD / 3, 0.03, "$OD_2$", size=15)
ax.tick_params(reset=True, direction="in")

ax.set_ylim(-0.02, ymax)
ax.set_xlabel("OD (A.U.)")
ax.set_ylabel("growth rate (h$^{-1}$)")

cmap = plt.get_cmap("tab10")
legend_names = [
    "0.1% Gal, 0.4% Pal",
]
legend_handles = [
    mlines.Line2D([], [], marker="", markersize=5, color=cmap(i), label=j, alpha=0.7)
    for i, j in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    loc="lower left",
    fontsize="small",
    bbox_to_anchor=(0.0, -0.01),
)

############
# Second plot in inset
axin = ax.inset_axes([0.55, 0.64, 0.45, 0.45])
sns.lineplot(
    data=p.s.loc[filt1],
    x="time",
    y="OD_mean",
    hue=None,
    legend=False,
    alpha=alpha,
    ax=axin,
)
axin.hlines(y=od_yield1, xmin=tmin, xmax=tmax, color="k")

rect = patches.Rectangle(
    (0, 0),
    tmax - tmin,
    od_yield1,
    linewidth=1,
    edgecolor="none",
    facecolor="yellow",
    alpha=0.2,
)
axin.add_patch(rect)
axin.arrow(
    y=0,
    x=3,
    dy=od_yield1,
    dx=0,
    length_includes_head=True,
    alpha=1,
    width=0.5,
    head_length=0.04,
    head_width=2,
    color="k",
)
axin.arrow(
    y=od_yield1,
    x=3,
    dy=-od_yield1,
    dx=0,
    length_includes_head=True,
    alpha=1,
    width=0.5,
    head_length=0.04,
    head_width=2,
    color="k",
)
axin.text(4.0, od_yield1 / 3, "$OD_1$", size=15)

rect2 = patches.Rectangle(
    (0, od_yield1),
    tmax - tmin,
    deltaOD,
    linewidth=1,
    edgecolor="none",
    facecolor="red",
    alpha=0.2,
)
axin.add_patch(rect2)
axin.arrow(
    y=od_yield1,
    x=3,
    dy=deltaOD,
    dx=0,
    length_includes_head=True,
    alpha=1,
    width=0.5,
    head_length=0.04,
    head_width=2,
    color="k",
)
axin.arrow(
    y=od_max,
    x=3,
    dy=-deltaOD,
    dx=0,
    length_includes_head=True,
    alpha=1,
    width=0.5,
    head_length=0.04,
    head_width=2,
    color="k",
)
axin.text(4.0, od_yield1 + deltaOD / 3, "$OD_2$", size=15)

axin.set_ylabel("OD (A.U.)")
axin.set_xlabel("time (h)")
axin.set_yticks([0, 1])
axin.tick_params(reset=True, direction="in")

ax = axes[3]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(
    [
        "Experiment_20200913",
        "Experiment_20201009_3",
        "Experiment_20201009_4",
    ]
)

add_sugar_column(p, ["Gal", "Pal", "Glu"], asstr=False)
add_maxOD_to_sc(p)
results2 = get_peak_valley(p, "Gal", "Pal")
plot_OD_vs_conc(p, results2, "Gal", "Pal", ax)
ax.tick_params(reset=True, direction="in")

plt.subplots_adjust(hspace=0.3, wspace=0.3)

# Add subplot labels
for n, ax in enumerate(axes):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )
plt.savefig(figdir + "fig2.png", bbox_inches="tight")
