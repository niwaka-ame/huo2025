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
sys.path.append(basedir + "src/intermediate/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = basedir + "data/fig3/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

# Framework
fig, axes = plt.subplots(2, 1, figsize=(4, 9), dpi=300)

import metabolomics as mb

ax = axes[1]
sns.lineplot(
    data=mb.met,
    x="Time",
    y="Area",
    hue="Sugar",
    marker="o",
    legend=None,
    markersize=8,
    ax=ax,
    alpha=0.7,
)
ax2 = ax.twinx()
cmap = plt.get_cmap("tab10")
sns.lineplot(
    data=mb.gdf,
    x="Time",
    y="CC OD",
    style="Condition",
    style_order=["0.1% Gal, 0.4% Pal", "0.1% Gal"],
    color=cmap(2),
    marker="o",
    ax=ax2,
    legend=None,
    markersize=8,
    alpha=0.7,
)
ax.set_ylabel("normalised sugar concentration")
ax.set_xlabel("time (h)")
ax2.set_ylabel("OD (A.U.)", color=cmap(2))
ax2.tick_params(axis="y", labelcolor=cmap(2))
# ax2.grid(False)
legend_names = [
    "Normalised [Gal] in media",
    "Normalised [Pal] in media",
    "OD in 0.1% Gal, 0.4% Pal",
]
legend_handles = [
    mlines.Line2D([], [], marker="o", markersize=5, color=cmap(i), label=j, alpha=0.7)
    for i, j in enumerate(legend_names)
]
legend_handles.append(
    mlines.Line2D(
        [],
        [],
        marker="o",
        markersize=5,
        color=cmap(len(legend_names) - 1),
        linestyle="dashed",
        label="OD in 0.1% Gal",
        alpha=0.7,
    )
)
ax.vlines(x=20, ymin=0, ymax=1, color="k", linestyle="dotted")
ax.legend(
    handles=legend_handles,
    fontsize="medium",
    loc="upper left",
    bbox_to_anchor=(1.1, 1),
)
ax.tick_params(reset=True, direction="in", right=False)
ax2.tick_params(direction="in")

ax = axes[0]
# load data
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20230527_2")
add_genotype(p, rename_dict={"1514": "IMA5-GFP", "BY4741": "WT"})
# choose data
cond_list = ["0.1% Gal, 0.4% Pal", "0.1% Fru, 0.4% Pal"]
strain_list = ["IMA5-GFP"]
linestyle = ["-"]
filt = p.s["condition"].isin(cond_list)
for s, ls in zip(strain_list, linestyle):
    filt_s = p.s["genotype"] == s
    data = p.s.loc[filt & filt_s]
    lineplot_with_error(
        data,
        x="OD_mean",
        y="c-GFPperOD",
        ax=ax,
        hue="condition",
        units="strain",
        ls=ls,
        reverse=True,
        alpha=0.5,
    )
# vline
ax.vlines(x=0.35, ymin=0.0, ymax=5000, color="k", linestyle="dotted")
# legend
cmap = plt.get_cmap("tab10")
legend_names = cond_list
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
ax.legend(
    handles=legend_handles,
    bbox_to_anchor=(1.04, 0),
    loc="lower left",
    fontsize="medium",
)
# labels
ax.set_ylabel("Ima5-GFP per OD (A.U.)")
ax.set_xlabel("OD (A.U.)")
ax.tick_params(reset=True, direction="in")

axin = ax.inset_axes([1.2, 0.5, 0.5, 0.5])
for s, ls in zip(strain_list, linestyle):
    filt_s = p.s["genotype"] == s
    data = p.s.loc[filt & filt_s]
    lineplot_with_error(
        data,
        x="time",
        y="OD_mean",
        ax=axin,
        hue="condition",
        units="strain",
        ls=ls,
        reverse=True,
        alpha=0.5,
    )
axin.hlines(y=0.35, xmin=0, xmax=40, linestyle="dotted", color="k")
axin.set_ylabel("OD (A.U.)")
axin.set_xlabel("time (h)")
axin.set_yscale("log", base=2)
axin.set_ylim(
    0.05,
)
axin.tick_params(reset=True, direction="in")

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
plt.savefig(figdir + "fig3.png", bbox_inches="tight")
