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

datadir = os.path.join(basedir, "data", "supp3")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")


fig, axes = plt.subplots(2, 2, figsize=(8, 8))
fig.subplots_adjust(wspace=0.25, hspace=0.25)
axes = axes.flat

########################
# Fig S3A
dfs = [
    "Experiment_20210813_4",
    "Experiment_20210811_4",
    "Experiment_20210922_1",
    "Experiment_20210816_4",
]
sugars = ["Glu", "Gal", "Fru", "Suc"]
q = "genotype =='FY4' and time < 34.0 and total == 1.0"
# plot
for j, (s, x) in enumerate(zip(sugars, [0, 1, 2, 3])):
    ax = axes[x]
    p = om.platereader(ls=False, wdir=datadir)
    p.importdf(dfs[j])
    add_sugar_column(p, [s, "Pal"], asstr=False)
    add_genotype(p)
    data = p.s.query(q)
    lineplot_with_error(data=data, x="time", y="OD_mean", ax=ax)
    # axis-Labels
    ax.set_xlabel("time (h)")
    ax.set_ylabel("OD (A.U.)")
    ax.set_xlim(-0.5, 36)
    # title
    ax.set_title(f"{s} + Pal")
    ax.set_ylim(0.05, 1.5)
    ax.set_yscale("log", base=2)
    # legend
    legend_names = [
        f"0.2% {s}, 0.8% Pal",
        f"0.4% {s}, 0.6% Pal",
        f"0.5% {s}, 0.5% Pal",
        f"0.6% {s}, 0.4% Pal",
        f"0.8% {s}, 0.2% Pal",
    ]
    cmap = plt.get_cmap("tab10")
    legend_handles = [
        mlines.Line2D([], [], linestyle="-", color=cmap(i), label=n)
        for i, n in enumerate(legend_names)
    ]
    if x < 3:
        ax.legend(
            handles=legend_handles,
            loc="lower right",
            fontsize="small",
        )
    else:
        ax.legend(
            handles=legend_handles,
            loc="upper left",
            fontsize="small",
        )
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
plt.savefig(os.path.join(figdir, "supp3.png"), bbox_inches="tight", dpi=300)
