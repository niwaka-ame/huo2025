import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import omniplate as om
import pandas as pd
import seaborn as sns
import os

scriptdir = os.path.dirname(os.path.realpath(__file__))
basedir = os.path.abspath(os.path.join(scriptdir, ".."))

from utils.diff_exp import *
from utils.yeastEnrichR import *
import string

sns.set_theme(context="paper", style="white")

datadir = os.path.join(basedir, "data", "rnaseq")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")
genefile = os.path.join(datadir, "yeast.txt")

DFGENE = pd.read_csv(
    genefile,
    delimiter="\s+",
    skiprows=58,
    usecols=[0, 1],
    names=["gene", "ORF"],
    index_col=None,
)

fig, axes = plt.subplot_mosaic(
    """
    AB
    XC
    XD
""",
    figsize=(8, 12),
    dpi=300,
)

axes = np.asarray(list(axes.values()))
axes[2].axis("off")

#####
# Fig S10AB

set_dict = {}
for t in [0, 10, 16]:
    df = pd.read_csv(os.path.join(datadir, f"group_gal80.Fru.{t}h_vs_WT.Fru.{t}h.csv"))
    process_deseq2_output(df, p_thr=0.05, gene_file=genefile)
    q = "change == 'upregulated' or change == 'downregulated'"
    set_temp = set(translate_gene_list(df.query(q)["SysName"].to_list(), df=DFGENE))
    set_temp.remove("GAL80")
    set_dict[t] = set_temp

set_dict_p01 = {}
for t in [0, 10, 16]:
    df = pd.read_csv(os.path.join(datadir, f"group_gal80.Fru.{t}h_vs_WT.Fru.{t}h.csv"))
    process_deseq2_output(df, p_thr=0.01, gene_file=genefile)
    q = "change == 'upregulated' or change == 'downregulated'"
    set_temp = set(translate_gene_list(df.query(q)["SysName"].to_list(), df=DFGENE))
    set_temp.remove("GAL80")
    set_dict_p01[t] = set_temp

from matplotlib_venn import venn3

tps = [0, 10, 16]
set_labels = [f"mid-log + {t} h" if t > 0 else "mid-log" for t in tps]
titles = [f"p < {p}" for p in ["0.05", "0.01"]]
axs = [axes[0], axes[1]]
for sd, ax, tt in zip([set_dict, set_dict_p01], axs, titles):
    v = venn3([sd[k] for k in tps], set_labels=set_labels, ax=ax)
    ax.set_title(tt)

ax = axes[3]

for t in [0, 10, 16]:
    df = pd.read_csv(os.path.join(datadir, f"group_gal80.Fru.{t}h_vs_WT.Fru.{t}h.csv"))
    process_deseq2_output(df, p_thr=0.05, gene_file=genefile)
    qup = "change == 'upregulated'"
    qdown = "change == 'downregulated'"
    gene_list_up = translate_gene_list(df.query(qup)["SysName"].to_list(), df=DFGENE)
    gene_list_down = translate_gene_list(
        df.query(qdown)["SysName"].to_list(), df=DFGENE
    )
    if t == 0:
        gene_set_up = set(gene_list_up)
        gene_set_down = set(gene_list_down)
    else:
        gene_set_up = gene_set_up.intersection(set(gene_list_up))
        gene_set_down = gene_set_down.intersection(set(gene_list_down))
gene_set_down.remove("GAL80")
gene_set = gene_set_up.union(gene_set_down)

# WARNING: the next line will query the enrichR server - do not abuse!
df_dict = run_enrichR(list(gene_set))
df_dict = remove_rows(df_dict, p_thr=0.01)  # Make results concise
for k in df_dict:
    df_dict[k]["Group size"] = df_dict[k]["Overlapping genes"].apply(len)
# rename
d = {
    "maturation of LSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA) (GO:0000463)": "maturation of LSU-rRNA from tricistronic rRNA transcript (GO:0000463)"
}
df_dict["GO_Biological_Process_2018"]["Term name"].replace(d, inplace=True)
# KEGG
order = (
    df_dict["KEGG_2018"]
    .sort_values(by="Adjusted p-value", ascending=True)["Term name"]
    .to_list()
)
sns.barplot(
    data=df_dict["KEGG_2018"],
    y="Term name",
    x="Group size",
    order=order,
    palette="husl",
    ax=ax,
)
join_genes = lambda l: ", ".join(l)
# GO
ax = axes[4]
order = (
    df_dict["GO_Biological_Process_2018"]
    .iloc[0:20]
    .sort_values(by="Adjusted p-value", ascending=True)["Term name"]
    .to_list()
)
sns.barplot(
    data=df_dict["GO_Biological_Process_2018"].iloc[0:20],
    y="Term name",
    x="Group size",
    order=order,
    palette="husl",
    ax=ax,
)
plt.subplots_adjust(wspace=0.25)

for n, ax in enumerate(axes[[0, 1, 3, 4]]):
    if n < 2:
        indent = -0.05
    else:
        indent = -1.3
    ax.text(
        indent,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )

plt.savefig(os.path.join(figdir, "supp10.png"), bbox_inches="tight")
