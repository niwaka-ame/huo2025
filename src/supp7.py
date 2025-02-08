import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import omniplate as om
import pandas as pd
import seaborn as sns
import string
import os

scriptdir = os.path.dirname(os.path.realpath(__file__))
basedir = os.path.abspath(os.path.join(scriptdir, ".."))

from utils.om_extra import *
from utils.diff_exp import *

sns.set_theme(context="paper", style="white")

datadir = os.path.join(basedir, "data/rnaseq")
figdir = os.path.join(basedir, "fig")
genefile = os.path.join(datadir, "yeast.txt")

sns.set_theme(context="paper", style="white")

# fig, axes = plt.subplots(4, 4, figsize=(12, 12))
fig, axes = plt.subplot_mosaic(
    """
    AABBCC
    DDDEEE
    FFGGHH
""",
    figsize=(12, 16),
    dpi=300,
)
axes = np.asarray(list(axes.values()))
fig.subplots_adjust(hspace=0.25)

#####
# Fig S7A
ax = axes[0]
ax.axis("off")
data = plt.imread(os.path.join(figdir, "supp7a.png"))
ax.imshow(data)

#####
# Fig S7B
ax = axes[1]
time = ["mid-log", "+10h", "+16h"] * 12
cond = (["0.1% Fru"] * 9 + ["0.1% Fru, 0.9% Pal"] * 9) * 2
geno = ["WT"] * 18 + ["gal80Δ"] * 18
biorep = (["A"] * 3 + ["B"] * 3 + ["C"] * 3) * 4

wt_fru = (
    [0.30, 0.31 * 5, 0.38 * 5] + [0.27, 0.3 * 5, 0.37 * 5] + [0.35, 0.3 * 5, 0.42 * 5]
)
wt_fp = (
    [0.29, 0.36 * 5, 0.28 * 10]
    + [0.26, 0.33 * 5, 0.26 * 10]
    + [0.29, 0.35 * 5, 0.28 * 10]
)
gal80_fru = (
    [0.30, 0.30 * 5, 0.37 * 5] + [0.30, 0.29 * 5, 0.36 * 5] + [0.29, 0.29 * 5, 0.37 * 5]
)
gal80_fp = (
    [0.29, 0.31 * 5, 0.42 * 5] + [0.29, 0.31 * 5, 0.39 * 5] + [0.29, 0.33 * 5, 0.38 * 5]
)

od = wt_fru + wt_fp + gal80_fru + gal80_fp
d = {
    "time point": time,
    "condition": cond,
    "genotype": geno,
    "biological replicate": biorep,
    "OD": od,
}
df = pd.DataFrame(data=d)
df["genotype & condition"] = df["genotype"] + " in " + df["condition"]
sns.stripplot(
    data=df,
    x="time point",
    y="OD",
    hue="genotype & condition",
    palette="tab20",
    ax=ax,
    alpha=0.5,
)
sns.pointplot(
    data=df,
    x="time point",
    y="OD",
    hue="genotype & condition",
    palette="tab20",
    errorbar="sd",
    ax=ax,
    marker="_",
    markersize=20,
    linestyle=":",
    legend=False,
)
ax.set_ylim(0, 3.5)
ax.set_ylabel("OD (A.U.)")
ax.tick_params(reset=True, direction="in")

######
# Fig S7C - PCA
ax = axes[2]
ax.axis("off")
data = plt.imread(os.path.join(figdir, "supp7c.png"))
ax.imshow(data)

######
# Fig S7DE
expt_list = [
    "group_gal80.Fru.0h_vs_WT.Fru.0h.csv",
    "group_WT.FruPal.0h_vs_WT.Fru.0h.csv",
]
title_list = [
    "gal80$\Delta$-Fru-midlog vs WT-Fru-midlog",
    "WT-FruPal-midlog vs WT-Fru-midlog",
]
for i, e in enumerate(expt_list, start=3):
    ax = axes[i]
    palette = ["grey", "red", "blue"]
    hue_order = ["insignificant", "upregulated", "downregulated"]
    df = pd.read_csv(os.path.join(datadir, e))
    process_deseq2_output(df, gene_file=genefile)
    df["log2baseMean"] = df["baseMean"].apply(np.log2)
    length = df["change"].value_counts().shape[0]
    if length == 2:
        palette = palette[:2]
        hue_order = hue_order[:2]
    if i == 3:
        legend_flag = True
    else:
        legend_flag = False
    sns.scatterplot(
        data=df,
        y="log2FoldChange",
        x="log2baseMean",
        hue="change",
        ax=ax,
        palette=palette,
        alpha=0.3,
        hue_order=hue_order,
        legend=legend_flag,
    )
    axes[i].set_title(title_list[i - 3])
    # Gene labels
    if i == 4:
        q = "change != 'insignificant'"
        dfchanged = df.query(q)
        ax.set_ylabel("")
        ax.set_xlabel("log$_2$ base mean")
    else:
        goi = ["GAL1", "GAL2", "GAL7", "GAL4", "GAL10", "GAL80", "PGM1", "PGM2", "GAL3"]
        filter_fun = lambda x: (
            bool(set(x.split(";")).intersection(set(goi))) if x else False
        )
        filt = df["Gene"].apply(filter_fun)
        dfchanged = df.loc[filt]
        ax.set_ylabel("log$_2$ fold change")
        ax.set_xlabel("log$_2$ base mean")
    xx = dfchanged["log2baseMean"].to_list()
    yy = dfchanged["log2FoldChange"].to_list()
    label = dfchanged["Gene"].to_list()
    points = zip(xx, yy, label)
    for a, b, l in points:
        ax.text(a, b, l.split(";")[0], rotation=0, rotation_mode="anchor", size=9)
    ax.tick_params(reset=True, direction="in")

#####
# Fig S7F
df = pd.read_csv(os.path.join(datadir, "processed_reads.csv"))
df["gc"] = df["genotype"] + " in " + df["condition"]
tp_name = ["mid-log", "+10h", "+16h"]
d = {i + 1: tp_name[i] for i in range(3)}
df["timepoint"].replace(d, inplace=True)
for j, g in enumerate(["MAL12", "MAL13", "ZNF1"]):
    ax = axes[j + 5]
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
    ax.set_ylim(3, 2000)
    ax.set_xlabel("")
    if j == 0:
        ax.set_ylabel("CPM")
    else:
        ax.set_ylabel("")
    ax.set_yscale("log")
    if j == 2:
        ax.legend(loc="upper center")
    ax.tick_params(reset=True, direction="in")
plt.subplots_adjust(hspace=0.3, wspace=0.25)

for n, ax in enumerate(axes[:6]):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )

plt.subplots_adjust(wspace=0.25)
plt.savefig(os.path.join(figdir, "supp7.png"), bbox_inches="tight")
