import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import omniplate as om
import pandas as pd
import seaborn as sns
import string
import sys
import os

basedir = os.path.expanduser("~/huo2025/")

sys.path.append(basedir + "src/utils/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = basedir + "data/supp6/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Fig S6A: gal80 gal4
ax = axes[0]
# data
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20220219_2")
d = {
    "1489": "gal80Δ",
    "1558": "gal80Δ gal4Δ",
    "FY4": "WT",
}
p.rename(d)
strain_list = ["gal80Δ", "gal80Δ gal4Δ", "WT"]
p.addcolumn("genotype", "strain", strain_list)
cond = ["0.1% Fru", "0.1% Fru, 0.4% Pal"]
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

###################
# Fig S6B
ax = axes[1]
p = om.platereader(ls=False, wdir=datadir)
p.importdf("Experiment_20230602_2")
d = {"FY4": "WT", "1536": "gal80Δ gal1-10-7Δ"}
p.rename(d)
strain_list = ["gal80Δ gal1-10-7Δ", "WT"]
p.addcolumn("genotype", "strain", strain_list)
cond = ["0.1% Fru", "0.1% Fru, 0.4% Pal"]
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

###################
# Fig S6C
ax = axes[2]
ax.axis("off")
data = plt.imread(figdir + "supp6c.png")
ax.imshow(data)

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
plt.savefig(figdir + "supp6.png", bbox_inches="tight", dpi=300)
