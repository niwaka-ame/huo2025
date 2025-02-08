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

datadir = os.path.join(basedir, "data", "fig1")
figdir = os.path.join(basedir, "fig")
svgdir = os.path.join(basedir, "svg")

# Framework
layout = """
AAB
"""
fig, axes = plt.subplot_mosaic(layout, figsize=(12, 4), dpi=300)
fig.subplots_adjust(wspace=0.25)
axes = np.asarray(list(axes.values()))

# Cartoon - method
ax = axes[0]
ax.axis("off")
cartoon = plt.imread(os.path.join(svgdir, "fig1a.png"))
ax.imshow(cartoon)

cmap = plt.get_cmap("tab10")
cmap2 = plt.get_cmap("Paired")
colors = [cmap(7), cmap2(6), cmap2(2), cmap(0), cmap(3)]
from matplotlib.colors import ListedColormap

cm = ListedColormap(colors, name="gr0")
try:
    import matplotlib as mpl

    mpl.colormaps.register(cmap=cm)
except:
    pass

ax = axes[1]
df = pd.read_csv(os.path.join(datadir, "gr_data.csv"))
order = ["2% Glu", "2% Fru", "2% Suc", "2% Gal", "2% Pal"]
q = "experiment != '20220218_fru'"
sns.stripplot(
    data=df.query(q),
    x="condition",
    y="gradient",
    hue="condition",
    order=order,
    ax=ax,
    palette="gr0",
    alpha=0.5,
)
sns.pointplot(
    data=df.query(q),
    x="condition",
    y="gradient",
    hue="condition",
    order=order,
    ax=ax,
    palette="gr0",
    linestyle="none",
    marker="_",
    markersize=20,
    errorbar="sd",
)
ax.set_ylabel("growth rate (h$^{-1}$)")
ax.set_xlabel("")
ax.set_ylim(0.2)
ax.tick_params(reset=True, direction="in")

# Add subplot labels
for n, ax in enumerate(axes):
    ax.text(
        -0.05,
        1.05,
        string.ascii_uppercase[n],
        transform=ax.transAxes,
        size=20,
        weight="bold",
    )
plt.subplots_adjust(wspace=0.3)
plt.savefig(os.path.join(figdir, "fig1.png"), bbox_inches="tight")
