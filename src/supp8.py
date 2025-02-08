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

datadir = os.path.join(basedir, "data", "supp8")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")

fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
axes = axes.flat

# Fig S8A
ax = axes[0]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20230415_2"])
rename_dict = {"1624": "gal80Δ MAL11-OE", "FY4": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["gal80Δ MAL11-OE", "WT"]
cond = ["0.1% Gal", "0.1% Gal, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data, x="time", y="OD_mean", ax=ax, hue="genotype", units="strain", ls=ls
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(0.05, 1)

# Legend
strain_list.sort()
legend_names = strain_list
cmap = plt.get_cmap("tab10")
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
legend_handles = legend_handles + [
    mlines.Line2D([], [], linestyle=ls, color="k", label=n)
    for n, ls in zip(cond, linestyle)
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize="small",
)

# Fig S8B
ax = axes[1]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20210611_1"])
rename_dict = {"1537": "gal80Δ ima1Δ", "BY4741": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["gal80Δ ima1Δ", "WT"]
cond = ["0.1% Gal", "0.1% Gal, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data, x="time", y="OD_mean", ax=ax, hue="genotype", units="strain", ls=ls
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(0.05, 1)

# Legend
strain_list.sort()
legend_names = strain_list
cmap = plt.get_cmap("tab10")
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
legend_handles = legend_handles + [
    mlines.Line2D([], [], linestyle=ls, color="k", label=n)
    for n, ls in zip(cond, linestyle)
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize="small",
)

# Fig S8C
ax = axes[2]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20230531_2"])
rename_dict = {"1360": "ima5Δ", "BY4741": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["ima5Δ", "WT"]
cond = ["0.1% Gal", "0.1% Gal, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data, x="time", y="OD_mean", ax=ax, hue="genotype", units="strain", ls=ls
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(0.05, 1)

# Legend
strain_list.sort()
legend_names = strain_list
cmap = plt.get_cmap("tab10")
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
legend_handles = legend_handles + [
    mlines.Line2D([], [], linestyle=ls, color="k", label=n)
    for n, ls in zip(cond, linestyle)
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize="small",
)

plt.subplots_adjust(wspace=0.25)

for n, ax in enumerate(axes):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )
    ax.tick_params(reset=True, direction="in")

plt.savefig(os.path.join(figdir, "supp8.png"), bbox_inches="tight")
