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

datadir = os.path.join(basedir, "data", "fig4")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")

# Framework
fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
fig.subplots_adjust(wspace=0.25)
axes = axes.flat

# Fig 4A: gal80
ax = axes[0]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20210227_3"])
rename_dict = {"1503": "gal80Δ", "BY4741": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["gal80Δ", "WT"]
cond = ["0.1% Gal", "0.1% Gal, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    filt_t = True
    data = p.s.loc[filt & filt_c & filt_t]
    lineplot_with_error(
        data,
        x="time",
        y="OD_mean",
        ax=ax,
        hue="genotype",
        units="strain",
        ls=ls,
        reverse=False,
        alpha=0.5,
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
# Legend
cmap = plt.get_cmap("tab10")
strain_list.sort()
legend_names = strain_list
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
ax.set_title("gal80$\\Delta$ in Gal+Pal")
ax.set_yscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_ylim(
    0.04,
)
ax.set_xlim(
    0,
)
ax.tick_params(reset=True, direction="in")

# Fig 4B: gal80 gal2
ax = axes[1]
# data
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20220318_1")
rename_dict = {"1556": "gal80Δ gal2Δ", "1489": "gal80Δ", "FY4": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
cond = ["0.1% Fru", "0.1% Fru, 0.4% Pal"]
linestyle = ["--", "-"]
strain_list = ["WT", "gal80Δ", "gal80Δ gal2Δ"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data,
        x="time",
        y="OD_mean",
        ax=ax,
        hue="genotype",
        units="strain",
        ls=ls,
        alpha=0.5,
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.tick_params(reset=True, direction="in")

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
ax.set_title("gal80$\\Delta$ & gal80$\\Delta$ gal2$\\Delta$ in Fru+Pal")
ax.set_yscale("log", base=2)
ax.set_ylim(
    0.04,
)
ax.set_xlim(
    0,
)
ax.tick_params(reset=True, direction="in")

# Fig 4C: GAL2-OE
ax = axes[2]
# data
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20220608_1")
rename_dict = {"1566": "GAL2-OE", "FY4": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["GAL2-OE", "WT"]
cond = ["0.1% Fru", "0.1% Fru, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data,
        x="time",
        y="OD_mean",
        ax=ax,
        hue="genotype",
        units="strain",
        ls=ls,
        reverse=True,
        alpha=0.5,
    )
    ax.set_xlabel("time (h)")
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(
    0.04,
)
ax.set_xlim(
    0,
)
ax.tick_params(reset=True, direction="in")

# Legend
legend_names = sorted(strain_list, reverse=True)
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
ax.set_title("GAL2-OE in Fru+Pal")

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
plt.savefig(os.path.join(figdir, "fig4.png"), bbox_inches="tight")
