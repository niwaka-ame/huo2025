import omniplate as om
import scipy.signal as ss
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import sys
import os

basedir = os.path.expanduser("~/huo2025/")

sys.path.append(basedir + "src/utils/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

datadir = basedir + "data/supp1/"
figdir = basedir + "fig/"
svgdir = basedir + "svg/"


def most_prom_peak(x, num):
    """
    Find the position of the num-th most prominent peak of growth rate.

    num: positive integer. num=1 is the first highest peak; num=2 is the second highest, etc.
    """
    peaks, properties = ss.find_peaks(x, prominence=1e-6)
    properties["position"] = peaks
    df = pd.DataFrame(properties)

    # sort
    df = df.sort_values(by="prominences", ascending=False)
    # Drop any small peaks on the left of the max growth rate
    df = df.drop(df.index[df["position"] < (df["position"].to_numpy()[0])])

    return (df["position"].to_numpy())[0:num]


def most_prom_valley(x, num):
    peaks, properties = ss.find_peaks(-x, prominence=1e-6)
    properties["position"] = peaks
    df = pd.DataFrame(properties)

    # Drop the small valleys on the left of the max growth rate
    # peak_pos = most_prom_peak(x, 1)
    # df = df.drop(df.index[df["position"] < peak_pos[0]])
    df = df.drop(df.index[df["position"] < 90])  # more than 15 hrs

    # sort
    df = df.sort_values(by="prominences", ascending=False)
    return (df["position"].to_numpy())[0:num]


def get_peak_valley(p, s1, s2):
    """
    Get the time and OD at the growth rate valley and second peak.

    s1, s2: the first and second sugar in the NAME of the mixture.
    """
    key_list = [
        # "2nd gr peak",
        "gr valley",
        # "time at 2nd gr peak",
        "time at gr valley",
        # "OD at 2nd gr peak",
        "OD at gr valley",
        "strain",
        "condition",
        s1,
        s2,
        "second OD yield",
        "max OD",
        "experiment",
    ]
    results = {key: [] for key in key_list}
    all_conditions = reduce(lambda x, y: x + y, p.allconditions.values())
    for e in p.allexperiments:
        for s in p.allstrains[e]:
            for c in p.allconditions[e]:
                if (s != "Null") and (s1 in c) and (s2 in c) and (", " in c):
                    # Get gr series for each condition and strain
                    filt = (
                        (p.s["condition"] == c)
                        & (p.s["strain"] == s)
                        & (p.s["experiment"] == e)
                    )
                    filt_sc = (
                        (p.sc["condition"] == c)
                        & (p.sc["strain"] == s)
                        & (p.sc["experiment"] == e)
                    )
                    gr_series = p.s.loc[filt, "gr"].to_numpy()

                    # Find peaks and valleys
                    # peak_pos = most_prom_peak(gr_series, 2)
                    valley_pos = most_prom_valley(gr_series, 1)

                    if len(valley_pos) == 0:
                        print("no valley found")
                        continue

                    # Find corresponding OD, time and gr
                    # if len(peak_pos) <= 1:
                    #     continue
                    # t_2nd_peak = p.s.loc[filt, "time"].to_numpy()[peak_pos[1]]
                    t_valley = p.s.loc[filt, "time"].to_numpy()[valley_pos[0]]

                    # od_2nd_peak = p.s.loc[filt, "OD_mean"].to_numpy()[peak_pos[1]]
                    od_valley = p.s.loc[filt, "OD_mean"].to_numpy()[valley_pos[0]]
                    second_OD_yield = (
                        p.sc.loc[filt_sc, "max OD"].to_numpy()[0] - od_valley
                    )
                    if second_OD_yield < 0.01:
                        print("wrong valley")
                        continue
                    od_max = p.sc.loc[filt_sc, "max OD"].to_numpy()[0]

                    # second_peak = p.s.loc[filt, "gr"].to_numpy()[peak_pos[1]]
                    valley = p.s.loc[filt, "gr"].to_numpy()[valley_pos[0]]

                    sugar1 = float(c.split(f"% {s1}, ")[0])
                    sugar2 = float(c.split(f"% {s1}, ")[1].split("%")[0])

                    # Append to results
                    res_list = [
                        # second_peak,
                        valley,
                        # t_2nd_peak,
                        t_valley,
                        # od_2nd_peak,
                        od_valley,
                        s,
                        c,
                        sugar1,
                        sugar2,
                        second_OD_yield,
                        od_max,
                        e,
                    ]
                    for i, k in enumerate(key_list):
                        results[k].append(res_list[i])

    return pd.DataFrame(results)


p = om.platereader(ls=False, wdir=datadir)
p.importdf(["Experiment_20210813_4"])

fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
add_sugar_column(p, ["Gal", "Pal", "Glu"], asstr=False)
add_maxOD_to_sc(p)
results2 = get_peak_valley(p, "Glu", "Pal")
plot_OD_vs_conc(p, results2, "Glu", "Pal", ax)
ax.tick_params(reset=True, direction="in")
plt.savefig(figdir + "supp1.png", bbox_inches="tight")
