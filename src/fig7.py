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

datadir = basedir + "data/fig7/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

# Framework
fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)
axes = axes.flat
fig.subplots_adjust(wspace=0.25, hspace=0.25)

# Fig 7A
ax = axes[0]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20230527_2"])
rename_dict = {"1538": "ima1Δ IMA5-GFP", "1514": "WT IMA5-GFP", "BY4741": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["WT IMA5-GFP", "ima1Δ IMA5-GFP"]
cond = ["0.1% Gal, 0.4% Pal"]
linestyle = ["-"]
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
        reverse=False,
        alpha=0.5,
    )
ax.hlines(y=0.35, xmin=0, xmax=40, linestyle="dotted", color="k")
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
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
ax.set_yscale("log", base=2)
ax.set_ylim(
    0.05,
)
ax.set_xlim(-1, 45)
ax.tick_params(reset=True, direction="in")

# Fig 7B
ax = axes[1]

ax.legend(
    handles=legend_handles,
    loc="upper left",
    fontsize="small",
)

filt = p.s["genotype"].isin(strain_list)
for c, ls in zip(cond, linestyle):
    filt_c = p.s["condition"] == c
    data = p.s.loc[filt & filt_c]
    lineplot_with_error(
        data,
        x="OD_mean",
        y="c-GFPperOD",
        ax=ax,
        hue="genotype",
        units="strain",
        ls=ls,
        reverse=False,
        alpha=0.5,
    )
ax.vlines(x=0.35, ymin=0, ymax=6500, linestyle="dotted", color="k")
ax.set_xlabel("OD (A.U.)")
ax.set_ylabel("GFP per OD (A.U.)")
ax.tick_params(reset=True, direction="in")

# Fig 7C
ax = axes[2]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20201230_4"])
rename_dict = {"1425": "ima1Δ", "1360": "ima5Δ", "BY4741": "WT"}
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["WT", "ima1Δ"]
filt = p.s["genotype"].isin(strain_list)
filt2 = p.s["condition"] == "2% Pal"
data = p.s.loc[filt & filt2]
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

legend_names = strain_list
cmap = plt.get_cmap("tab10")
cond = ["2% Pal"]
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
    loc="upper left",
    fontsize="small",
)
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(0.05)
ax.tick_params(reset=True, direction="in")

axin = ax.inset_axes([0.7, 0.1, 0.28, 0.38])
df1 = pd.read_csv(datadir + "../fig1/gr_data.csv")
df1["genotype"] = "WT"
df2 = pd.read_csv(datadir + "inset/gr_data_ima1.csv")
df2["genotype"] = "ima1$\Delta$"
q1 = "condition == '2% Pal'"
df1 = df1.query(q1)
df = pd.concat([df1, df2])
df.dropna(axis=0, how="any", inplace=True, subset=["genotype", "gradient"])
sns.stripplot(
    data=df,
    x="genotype",
    y="gradient",
    hue="genotype",
    ax=axin,
    alpha=0.3,
)
sns.pointplot(
    data=df,
    x="genotype",
    y="gradient",
    hue="genotype",
    ax=axin,
    linestyle="none",
    marker="_",
    markersize=20,
    errorbar="sd",
)
axin.set_ylabel("growth rate (h$^{-1}$)")
axin.set_xlabel("")
axin.set_ylim(0.2, 0.4)
ax.set_xlim(0, 45)
ax.tick_params(reset=True, direction="in")
axin.tick_params(reset=True, direction="in")

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
plt.savefig(figdir + "fig7.png", bbox_inches="tight")
