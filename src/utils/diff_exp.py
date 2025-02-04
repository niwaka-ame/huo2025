import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def process_deseq2_output(df, gene_file, log2fc_thr=0.5, p_thr=0.05, p_adjusted=True):
    """
    Add gene information and differential expression labels.
    """
    df.rename(columns={df.columns[0]: "feature"}, inplace=True)
    _add_gene_columns(df, gene_file)
    if p_adjusted:
        pcol = "padj"
    else:
        pcol = "pvalue"
    filt_p = df[pcol] < p_thr
    filt_up = df["log2FoldChange"] > log2fc_thr
    filt_down = df["log2FoldChange"] <= -log2fc_thr
    df["change"] = "insignificant"
    df.loc[filt_p & filt_up, "change"] = "upregulated"
    df.loc[filt_p & filt_down, "change"] = "downregulated"
    df["-log10 pvalue"] = -np.log10(df[pcol])


def _add_gene_columns(df, gene_file):
    """
    internal: add systematic name and gene names from feature.
    """
    split = lambda x: x.split("_")[0]
    df["SysName"] = df["feature"].apply(split)
    dfgene = pd.read_csv(
        gene_file,
        delimiter="\s+",
        skiprows=58,
        usecols=[0, 1],
        names=["gene", "ORF"],
        index_col=None,
    )
    df["Gene"] = df["SysName"].apply(_find_gene_name, args=(dfgene,))


def _find_gene_name(sysname, df):
    """
    internal: find the name of genes by systematic names. Inherit if the gene is not named yet.
    """
    gene = df.loc[df["ORF"] == sysname[0:7], "gene"].to_list()
    if gene:
        return gene[0]
    else:  # In case that locus has no name yet
        return sysname


def _find_goi(entry, goi):
    return bool(set(entry.split(";")).intersection(set(goi)))


def find_goi(df, goi, col_oi=["Gene", "padj", "log2FoldChange", "change"]):
    return df.loc[df["Gene"].apply(_find_goi, args=(goi,)), col_oi]
