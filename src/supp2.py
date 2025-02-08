import os
import seaborn as sns
import matplotlib.pyplot as plt

# from mpl_toolkits.mplot3d import Axes3D

import numpy as np

basedir = os.path.expanduser("~/huo2025/")

from wela.dataloader import dataloader
from wela.plotting import kymograph, plot_lineage
import string

sns.set_theme(context="paper", style="white")

datadir = os.path.join(basedir, "data", "supp2")
tsvdir = os.path.join(datadir, "tsv")
figdir = os.path.join(basedir, "fig")

# from genutils import figs2pdf
import matplotlib.cm

datasets = [
    "1713_2023_07_12_ss_switch_2pc_pyr_to_1pc_pal_cy5_gfp_00",
    "1515_2023_03_31_ss_switch_2pc_pyr_to_1pc_pal_0_2_pc_gal_cy5_gfp_00",
    "1498_2023_03_30_ss_switch_2pc_pyr_to_1pc_pal_0p2pc_suc_cy5_gfp_00",
]
fig = plt.figure(figsize=(8, 8))
titles = ["1% pal", "0.2% gal + 1% pal", "0.2% suc + 1% pal"]
cmap_offsets = [16, 12, 8]
axes = []
for i, dataname in enumerate(datasets):
    dl = dataloader(datadir, tsvdir)
    dl.load(dataname, use_tsv=True)
    t_and_ts = dl.get_time_series(signal="median_GFP", group=1121)
    ts = t_and_ts[1]
    ts_oi = np.log10(ts[:, [0, 72, 144, 216]])
    ax = fig.add_subplot(2, 2, i + 1, projection="3d")
    axes.append(ax)
    bins = np.linspace(2, 3.5, 15)
    xs = (bins[:-1] + bins[1:]) / 2
    times = [0, 6, 12, 18]
    cmap = plt.get_cmap("tab20c_r")
    for j, t in enumerate(times):
        data = ts_oi[:, j]
        dens, _ = np.histogram(data, bins=bins, density=True)
        dens_temp = dens[dens > 0]
        xs_temp = xs[dens > 0]
        c = cmap(j + cmap_offsets[i])
        # ax.hist(ts_oi, bins, density=True, histtype='bar', label=str(j))
        ax.bar(xs_temp, dens_temp, zs=t, zdir="y", color=c, ec=c, alpha=0.7, width=0.08)
    ax.set_title("2% pyruvate \n $\\to$ " + titles[i])
    ax.set_xlim(2, 3.5)
    ax.set_ylim(0, 20)
    ax.set_yticks(times)
    ax.set_xticks([2.0, 2.5, 3.0, 3.5])
    # ax.set_zticks([0, 0.01, 0.02, 0.03])
    ax.set_zlim(0, 8)
    ax.set_xlabel(r"$\log_{10}$ Ima1-GFP")
    xlims = [2, 3.5]
    ax.plot(xlims, [4, 4], [0, 0], color="red", alpha=0.5)
    ax.set_zlabel("density")
    ax.set_ylabel("time (h)")
fig.subplots_adjust(wspace=0.2, hspace=0.3)
for n, ax in enumerate(axes):
    ax.text2D(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )
plt.savefig(os.path.join(figdir, "supp2.png"), dpi=300)
