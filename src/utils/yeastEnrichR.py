import json
import requests
import pandas as pd


def run_enrichR(gene_list):
    """
    run YeastEnrichr. Return a dictionary of dataframes.
    """
    userid = post_enrichR(gene_list)
    data = get_enrichR(userid)
    return parse_enrichR(data)


def post_enrichR(gene_list):
    ENRICHR_URL = "http://maayanlab.cloud/YeastEnrichr/addList"
    genes_str = "\n".join(gene_list)
    description = "Example gene list"
    payload = {"list": (None, genes_str), "description": (None, description)}

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception("Error analyzing gene list")

    data = json.loads(response.text)
    print(data)
    return data["userListId"]


def get_enrichR(userid, lib_list=None):
    if lib_list is None:
        lib_list = [
            "Gene_Interaction_Hubs_BioGRID_2018",
            "GO_Biological_Process_2018",
            "GO_Cellular_Component_2018",
            "GO_Molecular_Function_2018",
            "InterPro_Domains_2019",
            "KEGG_2018",
            "KEGG_2019",
            "Pfam_Domains_2019",
            "PPI_Hubs_BioGRID_2018",
            "TF2DNA_2018",
        ]
    data = dict()
    for l in lib_list:
        ENRICHR_URL = "http://maayanlab.cloud/YeastEnrichr/enrich"
        query_string = "?userListId=%s&backgroundType=%s"
        user_list_id = userid
        gene_set_library = l
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception("Error fetching enrichment results")

        data[l] = json.loads(response.text)
    return data


def parse_enrichR(data):
    df_dict = dict()
    col_names = "Rank, Term name, P-value, Z-score, Combined score, Overlapping genes, Adjusted p-value, Old p-value, Old adjusted p-value".split(
        ", "
    )
    for lib in data:
        df = pd.DataFrame(data[lib][lib], columns=col_names)
        df_dict[lib] = df
    return df_dict


def _remove_rows(df, min_group_size, p_thr):
    big_group = lambda g: len(g) >= min_group_size
    low_p = lambda p: p < p_thr
    keep = df["Overlapping genes"].apply(big_group)
    keep2 = df["Adjusted p-value"].apply(low_p)
    df = df.loc[keep & keep2]
    return df


def remove_rows(df_dict, min_group_size=2, p_thr=0.01):
    """
    Remove rows of df_dict. Return a new dict.
    """
    df_dict_new = dict()
    for k in df_dict:
        df_dict_new[k] = _remove_rows(df_dict[k], min_group_size, p_thr)
    return df_dict_new


def find_gene_name(sysname, df):
    """
    internal: find the name of genes by systematic names. Inherit if the gene is not named yet.
    """
    gene = df.loc[df["ORF"] == sysname[0:7], "gene"].to_list()
    if gene:
        return gene[0]
    else:  # In case that locus has no name yet
        return sysname


def translate_gene_list(gene_list, df):
    """
    Translate systematic names into Entrez-like gene names.
    """
    rm_extra_names = lambda x: x.split(";")[0] if ";" in x else x
    all_names = [find_gene_name(g, df) for g in gene_list]
    return list(map(rm_extra_names, all_names))
