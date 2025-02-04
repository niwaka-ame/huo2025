import omniplate as om
from nunchaku import Nunchaku
import numpy as np
import sys
import os

basedir = os.path.expanduser("~/huo2025/")

sys.path.append(basedir + "src/utils/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = basedir + "data/fig7/inset/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"

p = om.platereader(wdir=datadir, ls=False)
expts = [
    "Experiment_20210116_1.xlsx",
    "Experiment_20201230_4.xlsx",
]
contents = [
    "Experiment_20210116_1_Content.xlsx",
    "Experiment_20201230_4_Content.xlsx",
]
p.load(expts, contents)
p.correctmedia()
p.correctOD()
q = "strain.str.contains('1425') and condition == '2% Pal'"
data = p.sc.query(q)


def run_nc(row, p):
    e = row["experiment"]
    s = row["strain"]
    c = row["condition"]
    q2 = "experiment == @e and strain == @s and condition == @c and time >= 2 and time <= 40"
    df = p.r.query(q2)
    df["log OD"] = df["OD"].apply(np.log)
    X = np.unique(df["time"].to_numpy())
    logOD = [
        df.query(f"well=='{well}'")["log OD"].to_numpy()
        for well in sorted(list(set(df["well"].to_list())))
    ]
    Y = np.stack(logOD)
    nc = Nunchaku(X, Y, estimate_err=False, prior=[0, 5], minlen=5)
    M_max = 10
    numseg, evi = nc.get_number(M_max)
    avg, _ = nc.get_iboundaries(numseg)
    info_df = (
        nc.get_info(avg)
        .sort_values(by="gradient", ascending=False)
        .reset_index(drop=True)
    )
    row["gradient"] = info_df.iloc[0].gradient
    return row


data = data.apply(run_nc, args=(p,), axis=1)
data.to_csv(datadir + "gr_data_ima1.csv")
