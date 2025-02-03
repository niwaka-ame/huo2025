import scipy.signal as ss
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
import string
import matplotlib.lines as mlines
import seaborn as sns


def add_sugar_column(p, sugars, asstr=False):
    """
    Add concentrations of each sugar and total conc in a mixture.

    p: platereader object;
    sugars: list of sugars of interest;
    asstr: concentration as string.
    """
    for s in sugars:
        p.addnumericcolumn(
            newcolumnname=s, oldcolumn="condition", leftsplitstr=s, asstr=False
        )
        p.s[s] = p.s[s].fillna(0.0)
    for i, s in enumerate(sugars):
        if i == 0:
            p.s["total"] = p.s[s]
            p.sc["total"] = p.sc[s]
        else:
            p.s["total"] = p.s["total"] + p.s[s]
            p.sc["total"] = p.sc["total"] + p.sc[s]
    # Convert sugar concentration to string
    if asstr:
        for s in sugars:
            p.addnumericcolumn(
                newcolumnname=s, oldcolumn="condition", leftsplitstr=s, asstr=True
            )


def add_subplot_labels(axes, flat=True, size=20):
    if flat:
        axes = axes.flat
    for n, ax in enumerate(axes):
        ax.text(
            -0.05,
            1.05,
            string.ascii_uppercase[n],
            transform=ax.transAxes,
            size=size,
            weight="bold",
        )


def add_genotype(p, rename_dict=None):
    if not rename_dict:
        rename_dict = {}
    split = lambda s: s.split("_")[0]
    p.s["genotype"] = p.s["strain"].apply(split).replace(rename_dict)
    p.sc["genotype"] = p.sc["strain"].apply(split).replace(rename_dict)


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
    peak_pos = most_prom_peak(x, 1)
    df = df.drop(df.index[df["position"] < peak_pos[0]])

    # sort
    df = df.sort_values(by="prominences", ascending=False)
    return (df["position"].to_numpy())[0:num]


def get_peak_valley(p, s1, s2):
    """
    Get the time and OD at the growth rate valley and second peak.

    s1, s2: the first and second sugar in the NAME of the mixture.
    """
    key_list = [
        "2nd gr peak",
        "gr valley",
        "time at 2nd gr peak",
        "time at gr valley",
        "OD at 2nd gr peak",
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
                    peak_pos = most_prom_peak(gr_series, 2)
                    valley_pos = most_prom_valley(gr_series, 1)

                    # Find corresponding OD, time and gr
                    if len(peak_pos) <= 1:
                        continue
                    t_2nd_peak = p.s.loc[filt, "time"].to_numpy()[peak_pos[1]]
                    t_valley = p.s.loc[filt, "time"].to_numpy()[valley_pos[0]]

                    od_2nd_peak = p.s.loc[filt, "OD_mean"].to_numpy()[peak_pos[1]]
                    od_valley = p.s.loc[filt, "OD_mean"].to_numpy()[valley_pos[0]]
                    second_OD_yield = (
                        p.sc.loc[filt_sc, "max OD"].to_numpy()[0] - od_valley
                    )
                    od_max = p.sc.loc[filt_sc, "max OD"].to_numpy()[0]

                    second_peak = p.s.loc[filt, "gr"].to_numpy()[peak_pos[1]]
                    valley = p.s.loc[filt, "gr"].to_numpy()[valley_pos[0]]

                    sugar1 = float(c.split(f"% {s1}, ")[0])
                    sugar2 = float(c.split(f"% {s1}, ")[1].split("%")[0])

                    # Append to results
                    res_list = [
                        second_peak,
                        valley,
                        t_2nd_peak,
                        t_valley,
                        od_2nd_peak,
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


def lineplot_with_error(data, x, y, ax, cmap="tab10", hue="condition", units="strain", reverse=False, **kwargs):
    hues = list(set(data[hue].to_list()))
    expts = list(set(data["experiment"].to_list()))
    expts.sort()
    hues.sort(reverse=reverse)
    units_list = set(data[units].to_list())
    cm = plt.get_cmap(cmap)
    if y == "OD_mean":
        y_err = "OD_err"
    else:
        y_err = y + "_err"
    for i, c in enumerate(hues):
        for e in expts:
            for s in units_list:
                filt_cs = (data[hue] == c) & (data[units] == s) & (data["experiment"] == e)
                data_cs = data.loc[filt_cs].sort_values(by="time")
                if data_cs.shape[0] < 1:
                    continue
                X = data_cs[x].to_numpy()
                Y = data_cs[y].to_numpy()
                Y_err = data_cs[y_err].to_numpy()
                Y_upper = Y + Y_err
                Y_lower = Y - Y_err
                ax.plot(X, Y, color=cm(i), label=c, **kwargs)
                ax.fill_between(X, Y_upper, Y_lower, alpha=0.1, color=cm(i))


def add_maxOD_to_sc(p):
    g = p.s.groupby(by=["experiment", "strain", "condition"])["OD_mean"].max()
    g = g.reset_index()

    def query_g(row, g):
        filt1 = g["experiment"] == row["experiment"]
        filt2 = g["strain"] == row["strain"]
        filt3 = g["condition"] == row["condition"]
        filt = filt1 & filt2 & filt3
        data = g.loc[filt, "OD_mean"].to_numpy()
        if len(data) > 0:
            return data[0]
        else:
            print(row)
            return np.nan

    p.sc["max OD"] = p.sc.apply(query_g, args=(g,), axis=1)


def find_stat_at_max_gr(row, df, dtype):
    """
    Query time of local max gr from ROW of p.sc and get corresponding DTYPE from DF.
    Usage: p.sc.apply(find_stat_at_max_gr, axis=1, args=(p.s, field,))
    """
    filt = (
        (df["experiment"] == row["experiment"])
        & (df["strain"] == row["strain"])
        & (df["condition"] == row["condition"])
    )
    filt2 = ((df.loc[filt, "time"] - row["time of local max gr"]).apply(abs)) < 0.1
    if df.loc[filt & filt2].shape[0] > 0:
        return df.loc[filt & filt2, dtype].to_numpy()[0]
    return np.nan


def find_stat_at_OD(row, df, dtype, od=0.1):
    filt = (
        (df["experiment"] == row["experiment"])
        & (df["strain"] == row["strain"])
        & (df["condition"] == row["condition"])
    )
    ind = (
        (df.loc[filt, "OD_mean"] - od)
        .apply(abs)
        .sort_values(ascending=True)
        .index.to_list()
    )
    if ind:
        return df.loc[ind[0], dtype]
    return np.nan


def lineplot_with_maxgr_point(
    data,
    data_sc,
    x,
    y,
    ax,
    cmap="tab10",
    hue="condition",
    units="strain",
    size=100,
    err=False,
    reverse=False,
    hues_order=None,
    key="at max specific gr",
    show_all=1,
):
    if not hues_order:
        hues = list(set(data[hue].to_list()))
        hues.sort(reverse=reverse)
    else:
        hues = hues_order
    expts = list(set(data["experiment"].to_list()))
    expts.sort()
    units_list = set(data[units].to_list())
    if isinstance(cmap, str):
        cm = plt.get_cmap(cmap)
    else:
        cm = lambda i: cmap[i]
    if y == "OD_mean":
        y_err = "OD_err"
    else:
        y_err = y + "_err"
    for i, c in enumerate(hues):
        for e in expts:
            for s in units_list:
                filt_cs = (data[hue] == c) & (data[units] == s) & (data["experiment"] == e)
                data_cs = data.loc[filt_cs].sort_values(by="time")
                if data_cs.shape[0] < 1:
                    continue
                filt = (data_sc[hue] == c) & (data_sc[units] == s) & (data_sc["experiment"] == e)
                x_at_maxgr = data_sc.loc[filt, f"{x} {key}"].to_numpy()[0]
                y_at_maxgr = data_sc.loc[filt, f"{y} {key}"].to_numpy()[0]
                if show_all == 0:
                    filt_all = data[x] < x_at_maxgr
                    data_cs = data_cs.loc[filt_all]
                elif show_all == 2:
                    filt_all = data[x] > x_at_maxgr
                    data_cs = data_cs.loc[filt_all]
                X = data_cs[x].to_numpy()
                Y = data_cs[y].to_numpy()
                if err:
                    Y_err = data_cs[y_err].to_numpy()
                    Y_upper = Y + Y_err
                    Y_lower = Y - Y_err
                    ax.fill_between(X, Y_upper, Y_lower, alpha=0.05, color=cm(i))
                ax.plot(X, Y, color=cm(i), alpha=0.6, label=c)
                ax.scatter(x_at_maxgr, y_at_maxgr, s=size, color=cm(i), alpha=0.8, edgecolors="k")

def add_data_at_maxgr_to_sc(
    p, dtypes=["OD_mean", "c-GFP120", "c-GFP120perOD", "d/dt OD", "flog(OD)"]
):
    for field in dtypes:
        p.sc[field + " at max specific gr"] = p.sc.apply(
            find_stat_at_max_gr, axis=1, args=(p.s, field,)
        )


def add_data_at_OD_to_sc(
    p,
    dtypes=["OD_mean", "c-GFP120", "c-GFP120perOD", "d/dt OD", "flog(OD)", "gr"],
    od=0.1,
):
    for field in dtypes:
        p.sc[field + f" at OD {od}"] = p.sc.apply(
            find_stat_at_OD, axis=1, args=(p.s, field, od,)
        )


def regress_slope_acceleration(row, p, phases=False, x="c-GFP120", y="d/dt OD"):
    e = row["experiment"]
    s = row["strain"]
    c = row["condition"]
    filt = (p.s["experiment"] == e) & (p.s["strain"] == s) & (p.s["condition"] == c)
    if not phases:
        tmax = row["time of local max gr"]
        filt_t = p.s["time"] <= tmax
    else:
        filt_t = p.s["phase"].isin(phases)
    filt_na = (~p.s[x].isna()) & (~p.s[y].isna())
    data = p.s.loc[filt & filt_t & filt_na]
    X = data[x].to_numpy()
    Y = data[y].to_numpy()
    r = linregress(X, Y)
    row["slope"] = r.slope
    row["intercept"] = r.intercept
    row["rsquare"] = r.rvalue**2
    return row


def add_tr_eff(p, gfp_channel="c-GFP120"):
    p.s["tr eff"] = p.s["d/dt OD"] / p.s[gfp_channel]
    # error
    def tr_eff_err(row):
        f = row["tr eff"]
        x = row[gfp_channel]
        y = row["d/dt OD"]
        sigma_x = row[f"{gfp_channel} err"]
        sigma_y = row["d/dt OD err"]
        if (y != 0) and (x != 0):
            return np.abs(f) * np.sqrt((sigma_x / x) ** 2 + (sigma_y / y) ** 2)
        else:
            return np.nan

    p.s["tr eff err"] = p.s.apply(tr_eff_err, axis=1)

def add_growth_phase(p, data, quantile=0.9):
    conds = list(set(data["condition"].to_list()))
    strs = list(set(data["strain"].to_list()))
    expts = list(set(data["experiment"].to_list()))
    for e in expts:
        for c in conds:
            for s in strs:
                q = "experiment == @e and condition == @c and strain == @s"
                if data.query(q).shape[0] == 0:
                    continue
                maxgr = p.sc.query(q)["local max gr"].to_numpy()[0]
                maxgr_thr = maxgr * quantile
                q2 = "gr >= @maxgr_thr"
                gr_indices = data.query(q).query(q2).index.to_numpy()[::-1]
                max_exp = gr_indices.max()
                indices_diff = np.diff(gr_indices)
                mask = np.where(indices_diff < -1, True, False)
                index_region = gr_indices[:-1][mask]
                if len(index_region) > 0:
                    min_exp = index_region[-1]
                else:
                    min_exp = gr_indices.min()
                min_acc = data.query(q).index.to_numpy().min()
                max_ret = data.query(q).index.to_numpy().max()
                p.s.loc[min_acc:min_exp-1, "phase"] = "acc"
                p.s.loc[min_exp:max_exp, "phase"] = "exp"
                p.s.loc[max_exp+1:max_ret, "phase"] = "ret"

def plot_OD_vs_conc(p, results, s1, s2, ax, marker="o"):
    """
    results: dataframe returned by get_peak_valley();
    s1: first sugar to be consumed;
    s2: second sugar to be consumed.
    """
    line1, line2 = linregress_OD_conc(results, s1, s2)
    filt = (results[s1] > 0) & (results[s2] > 0)
    cmap = plt.get_cmap("tab10")
    sns.scatterplot(
        data=results.loc[filt], x=s1, y="OD at gr valley", ax=ax, color=cmap(0), s=50, marker=marker
    )
    sugar1 = results.loc[filt, s1].to_numpy()
    ax.plot(sugar1, line1.slope * sugar1 + line1.intercept, color=cmap(0))
    sns.scatterplot(
        data=results.loc[filt], x=s2, y="second OD yield", ax=ax, color=cmap(1), s=50, marker=marker
    )
    sugar2 = results.loc[filt, s2].to_numpy()
    ax.plot(sugar2, line2.slope * sugar2 + line2.intercept, color=cmap(1))
    ax.set_ylabel("OD yield (A.U.)")
    ax.set_xlabel("sugar concentration (%)")
    # legend
    legend_names = [
        # f"[{s1}] vs $OD_1$: $y={line1.slope:.2f} x + {line1.intercept:.2f}$, $ R^2 = {line1.rvalue**2:.2f}$",
        # f"[{s2}] vs $OD_2$: $y={line2.slope:.2f} x + {line2.intercept:.2f}$, $ R^2 = {line2.rvalue**2:.2f}$",
        f"[{s1}] vs $OD_1$: $ R^2 = {line1.rvalue**2:.2f}$",
        f"[{s2}] vs $OD_2$: $ R^2 = {line2.rvalue**2:.2f}$",
    ]
    legend_handles = [
        mlines.Line2D([], [], linestyle="-", marker=marker, color=cmap(i), label=n)
        for i, n in enumerate(legend_names)
    ]

    ax.legend(
        handles=legend_handles,
        fontsize="small",
    )
    ax.set_xlim(-0.02, 0.85)
    ax.set_ylim(-0.02, 1.5)

def linregress_OD_conc(df, s1, s2):
    """
    df: dataframe returned by get_peak_valley().
    """
    line1 = linregress(df.loc[:, [s1, "OD at gr valley"]])
    line2 = linregress(df.loc[:, [s2, "second OD yield"]])
    return line1, line2
