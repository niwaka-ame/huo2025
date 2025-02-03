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

datadir = basedir + "data/supp4/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

# Framework
fig, axes = plt.subplots(2, 2, figsize=(8, 8))
fig.subplots_adjust(wspace=0.25, hspace=0.25)
axes = axes.flat

########################################

# Fig S4ABC - history

p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20210409_3", "Experiment_20210419_3"])
add_sugar_column(p, ["Pal", "Gal"])
d = {"FY4_1": "FY4(Pal)_1", "FY4_2": "FY4(Pal)_2"}
p.rename(d)
p.addcolumn("history", "strain", ["Raf", "Glu", "Pal"])

q = "strain != 'Null' and condition != '0% Gal'"
# q3 = "not (history == 'Pal' and strain == 'FY4(Pal)_2' and condition == '0.05% Gal, 0.45% Pal')"
data = p.s.query(q)  # .query(q3)
history_d = {"Raf": "raffinose", "Glu": "glucose", "Pal": "palatinose"}

# Legend
cmap = plt.get_cmap("tab10")
legend_names = p.s["condition"].value_counts().index.to_list()
legend_names.remove("0% Gal")
legend_handles = [
    mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
    for i, n in enumerate(legend_names)
]
for i, h in enumerate(["Glu", "Raf", "Pal"]):
    ax = axes[i]
    q2 = "history == @h"
    filt_h = p.s["history"] == h
    lineplot_with_error(data.query(q2), "time", "OD_mean", ax)
    ax.set_xlabel("time (h)")
    ax.set_title(history_d[h] + " history")
    ax.set_ylim(0.05, 1)
    ax.set_yscale("log", base=2)
    ax.set_ylabel("OD (A.U.)")
    ax.set_xlabel("time (h)")
    if i == 0:
        ax.legend(handles=legend_handles, loc="upper left", fontsize="small")
    else:
        ax.legend(handles=legend_handles, loc="lower right", fontsize="small")
    ax.tick_params(reset=True, direction="in")


##############
# Fig S4D: ethanol
ax = axes[3]
p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20230531_2"])
rename_dict = {"1626": "icl1Δ", "BY4741": "WT"}
split = lambda s: s.split("_")[0]
p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
strain_list = ["icl1Δ", "WT"]
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
    )
ax.set_xlabel("time (h)")
ax.set_ylabel("OD (A.U.)")
ax.set_yscale("log", base=2)
ax.set_ylim(0.05, 1)
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
ax.legend(handles=legend_handles, loc="lower right", fontsize="small")
ax.tick_params(reset=True, direction="in")

for n, ax in enumerate(axes):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )
plt.savefig(figdir + "supp4.png", bbox_inches="tight")
