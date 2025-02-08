import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import omniplate as om
import pandas as pd
import seaborn as sns
import os

scriptdir = os.path.dirname(os.path.realpath(__file__))
basedir = os.path.abspath(os.path.join(scriptdir, ".."))

from utils.om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = os.path.join(basedir, "data", "supp5")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flat

p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20210825_2")
cond = ["0.2% Suc", "0.2% Suc, 0.8% Pal"]
q = "condition.isin(@cond) and strain.str.contains('FY4_2')"
data = p.s.query(q)
# linestyle = ["--", "-"]
for j, y in enumerate(["OD_mean", "gr"]):
    ax = axes[j]
    # filt_c = data["condition"] == c
    lineplot_with_error(
        data=data, x="time", y=y, ax=ax, units="strain", hue="condition", alpha=0.5
    )
    ax.set_xlabel("time (h)")
    if j == 0:
        ax.set_ylabel("OD (A.U.)")
        ax.set_yscale("log", base=2)
        ax.set_ylim(0.05, 1.5)
    else:
        ax.set_ylabel("growth rate (h$^{-1}$)")
        ax.set_ylim(
            0,
        )
    ax.tick_params(reset=True, direction="in")

cmap = plt.get_cmap("tab10")
legend_names = sorted(cond)
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    loc="upper right",
    # fontsize="small",
)
ax.tick_params(reset=True, direction="in")

ax = axes[2]
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20210417_1")
cond = ["0.1% Suc, 0.4% Pal"]
strains = ["1514_1", "1514_2"]
q = "condition.isin(@cond) and strain.isin(@strains)"
data = p.s.query(q)
lineplot_with_error(
    data=data,
    x="OD_mean",
    y="c-GFPperOD",
    ax=ax,
    units="strain",
    hue="condition",
    alpha=0.5,
)
ax.set_ylabel("Ima5-GFP per OD (A.U.)")
ax.set_xlabel("OD (A.U.)")
legend_names = sorted(cond)
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    # fontsize="small",
)
ax.tick_params(reset=True, direction="in")


p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20230527_2")
cond = ["0.1% Fru", "0.1% Fru, 0.4% Pal"]
q = "condition.isin(@cond) and strain.str.contains('BY4741')"
data = p.s.query(q)
# linestyle = ["--", "-"]
for j, y in enumerate(["OD_mean", "gr"], start=3):
    ax = axes[j]
    # filt_c = data["condition"] == c
    lineplot_with_error(
        data=data, x="time", y=y, ax=ax, units="strain", hue="condition", alpha=0.5
    )
    ax.set_xlabel("time (h)")
    ax.set_xticks([10, 20, 30, 40])
    if j == 3:
        ax.set_ylabel("OD (A.U.)")
        ax.set_yscale("log", base=2)
        ax.set_ylim(0.05, 1.5)
    else:
        ax.set_ylabel("growth rate (h$^{-1}$)")
        ax.set_ylim(
            0,
        )
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.tick_params(reset=True, direction="in")

cmap = plt.get_cmap("tab10")
legend_names = sorted(cond)
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    loc="upper right",
    # fontsize="small",
)
ax.tick_params(reset=True, direction="in")

plt.subplots_adjust(wspace=0.3)

axes[-1].axis("off")

# Add subplot labels
for n, ax in enumerate(axes[:-1]):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )

plt.savefig(os.path.join(figdir, "supp5.png"), bbox_inches="tight", dpi=300)
