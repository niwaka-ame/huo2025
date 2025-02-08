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

datadir = os.path.join(basedir, "data", "fig5")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")

# Framework
fig, axes = plt.subplots(1, 4, figsize=(13, 4), dpi=300)
fig.subplots_adjust(wspace=0.25)
# Fig 5ABC
df = pd.read_csv(os.path.join(datadir, "../rnaseq/processed_reads.csv"))
df["gc"] = df["genotype"] + " in " + df["condition"]
tp_name = ["mid-log", "+10h", "+16h"]
d = {i + 1: tp_name[i] for i in range(3)}
df["timepoint"].replace(d, inplace=True)
for j, g in enumerate(["MAL11", "IMA1", "IMA5"], start=0):
    ax = axes[j]
    ax.set_title(g)
    q = "Gene.str.contains(@g).values and cpm > 0"
    sns.stripplot(
        data=df.query(q),
        x="timepoint",
        y="cpm",
        palette="tab20",
        ax=ax,
        hue="gc",
        alpha=0.5,
    )
    sns.pointplot(
        data=df.query(q),
        x="timepoint",
        y="cpm",
        palette="tab20",
        ax=ax,
        hue="gc",
        errorbar="sd",
        marker="_",
        markersize=20,
        linestyle=":",
        legend=False,
    )
    ax.legend().remove()
    ax.set_ylim(3, 5000)
    ax.set_xlabel("")
    ax.set_ylabel("CPM")
    ax.set_yscale("log")
    ax.tick_params(reset=True, direction="in")
axes[0].legend(fontsize="small")
# Fig 5D
ax = axes[3]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20230408_2"])
rename_dict = {"FY4": "WT", "1622": "MAL11-OE"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["WT", "MAL11-OE"]
cond = ["0.1% Gal", "0.1% Gal, 0.4% Pal"]
linestyle = ["--", "-"]
filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    filt_t = p.s["time"] <= 40
    data = p.s.loc[filt & filt_c & filt_t]
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
ax.set_ylabel("OD (A.U.)")
# ax.set_ylim(-0.02, 1.)
ax.set_ylim(0.05)
ax.set_xlim(0, 42)
ax.set_yscale("log", base=2)
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

plt.subplots_adjust(wspace=0.3)

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
plt.savefig(os.path.join(figdir, "fig5.png"), bbox_inches="tight")
