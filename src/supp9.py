import numpy as np
import seaborn as sns
from scipy.optimize import fsolve
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from sympy import symbols, Symbol, simplify, numer, Poly, nroots
from matplotlib.colors import LogNorm
import string
import os

basedir = os.path.expanduser("~/huo2025/")

sys.path.append(basedir + "src/utils/")
from om_extra import *
import string

sns.set_theme(context="paper", style="white")

figdir = basedir + "fig/"
svgdir = basedir + "svg/"


def uT(p, bT, uTmax, KT, n):
    return (bT + uTmax * (p / KT) ** n) / (1 + (p / KT) ** n)


def uI(p, alpha, n):
    return (alpha * p**n) / (1 + p**n)


def vI(p, KmI, vImax):
    return vImax * p / (KmI + p)


def eqn(p, KT, alpha, n, params):
    dT, uTmax, bT, KmI, vImax = params
    return vI(p, KmI, vImax) * uI(p, alpha, n) / uT(p, bT, uTmax, KT, n) * dT - 1


# bifur
KT1 = 1.2
alpha1 = 1.2
n1 = 2

# params
dT = 0.37
uTmax = 3.7e-5
bT = 3.7e-7
KmI = 18
vImax = 1.7e5

params = (dT, uTmax, bT, KmI, vImax)

p, KT, alpha = symbols("p, K_T, \\alpha")
n = Symbol("n", integer=True)
poly = numer(simplify(eqn(p, KT, alpha, n, params)))


def solve_ss_min(KT1, alpha1, n1, poly):
    poly1 = Poly(poly.subs({KT: KT1, alpha: alpha1, n: n1}), p)
    roots = nroots(poly1)
    roots = [r for r in roots if r.is_real and r > 0]
    return min(roots)


def solve_ss_max(KT1, alpha1, n1, poly):
    poly1 = Poly(poly.subs({KT: KT1, alpha: alpha1, n: n1}), p)
    roots = nroots(poly1)
    roots = [r for r in roots if r.is_real and r > 0]
    return max(roots)


sss_min_v = np.vectorize(solve_ss_min, excluded=["n1", "poly"])
sss_max_v = np.vectorize(solve_ss_max, excluded=["n1", "poly"])

# make grid
npoints = 30
KTs = np.linspace(1, 10, npoints)
alphas = np.linspace(6e-10, 5e-9, npoints)
KTv, alphav = np.meshgrid(KTs, alphas)
n1s = (2, 4)

from skimage import measure

fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=300)
fig.subplots_adjust(wspace=0.25)
for j, n1 in enumerate(n1s):
    # solve
    low_matrix = sss_min_v(KTv, alphav, n1=n1, poly=poly)
    high_matrix = sss_max_v(KTv, alphav, n1=n1, poly=poly)
    high_matrix2 = high_matrix.copy()
    high_matrix2[high_matrix == low_matrix] = np.nan
    contours = measure.find_contours(high_matrix == low_matrix)
    # plot
    ax = axes[j]
    cax1 = ax.imshow(
        low_matrix.astype(float) / KTv,
        origin="lower",
        aspect="auto",
        cmap="bwr",
        norm=LogNorm(vmin=0.1, vmax=10),
        # vmax=2.3,
        # vmin=-2.3,
    )
    for contour in contours:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=1, color="yellow")

    ax.set_xlabel("$K_T$: repression strength")
    ax.set_ylabel(r"$u_{I,\max} / v_T$: isomaltase expression ($\times 10^{-9}$)")
    ax.set_xticks(np.linspace(0, npoints - 1, 5), np.round(np.linspace(1, 10, 5), 2))
    ax.set_yticks(np.linspace(0, npoints - 1, 5), np.round(np.linspace(0.6, 5, 5), 2))
    ax.set_title(f"$n = {n1}$")
    cbar1 = fig.colorbar(cax1, ax=ax)

    axin = ax.inset_axes([0.45, 0.55, 0.5, 0.4])
    axin.imshow(
        (high_matrix2.astype(float) / KTv)[: int(npoints * 4 / 5), :],
        origin="lower",
        aspect="auto",
        cmap="bwr",
        norm=LogNorm(vmin=0.1, vmax=10),
    )
    axin.set_xticks(np.linspace(0, npoints - 1, 3), np.round(np.linspace(1, 10, 3), 2))
    axin.set_yticks(
        np.linspace(0, npoints * 4 / 5 - 1, 4), np.round(np.linspace(0.6, 3.9, 4), 2)
    )
    axin.tick_params(colors="white")

fig_dir = "/home/yu/papers/gal-mal-paper/figs/"

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

plt.savefig(figdir + "supp9.png", bbox_inches="tight")
