"""Figure 1 of paper/paper.qmd: three ways of looking at the same data.

The observed reaction times sit on the left. Three arrows, each carrying a
magnifying glass, lead to the three approaches. Every row shows what that
model is made of and the distribution that follows from it, drawn back over
the same data.

The trials are the lexical decision data of Wagenmakers et al. (2008) used in
Examples 1 and 2: correct responses from participants 1-3, 0.2 s < RT < 2 s.
The ex-Gaussian curves use the posterior means reported in Example 1. The
shifted Wald curves are maximum-likelihood fits to the same trials, since the
paper reports drift and boundary but not the non-decision time.

    Rscript make_figures.R          # writes paper/wagenmakers2008.csv
    python make_fig_approaches.py

Writes paper/figures/fig_approaches.png.
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch
from scipy.optimize import minimize
from scipy.stats import exponnorm, norm

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "wagenmakers2008.csv")
OUT = os.path.join(HERE, "figures", "fig_approaches.png")

ACC, SPD = "#2E7D32", "#B03A2E"          # matches the posterior predictive figure
INK, MUTED, PALE = "#22252E", "#6B7280", "#C9CDD6"
ACCENT = "#7048C4"

W_IN, H_IN, DPI = 8.6, 6.19, 400         # drawn oversize; the PDF scales it down
XLIM = (0.2, 1.6)
GRID = np.linspace(*XLIM, 600)

# labels and arrows sit above the curves (which are drawn at zorder 4), with a
# white halo so they stay readable where they cross one
Z_ARROW, Z_TEXT = 7, 8
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]

# figure-fraction layout
COL_DATA = (0.035, 0.215)
COL_MECH = (0.358, 0.528)
COL_DENS = (0.560, 0.988)
ARROW_X = (0.221, 0.352)
PANEL_H = 0.165
ROW_Y0 = (0.665, 0.395, 0.125)           # top, middle, bottom
DATA_Y = (0.275, 0.680)

TITLE = "From means to mechanisms"
SUBTITLE = "Three ways of modelling the same reaction times"


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def load():
    rt = {"Accuracy": [], "Speed": []}
    with open(CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            # `Error` is written as 0/1 by write.csv(); older versions of the
            # dataset stored it as a logical, hence "TRUE" as well.
            if r["Error"] in ("TRUE", "1") or int(r["Participant"]) > 3:
                continue
            # RT <= 2 s is the only exclusion, matching the models exactly.
            # A lower bound here would drop four fast trials the fitted models
            # keep, and put this figure's n out from every other n in the paper.
            t = float(r["RT"])
            if t <= 2.0:
                rt[r["Condition"]].append(t)
    return {k: np.array(v) for k, v in rt.items()}


# --------------------------------------------------------------------------
# shifted Wald (first passage of a single accumulator, diffusion sd = 1)
# --------------------------------------------------------------------------
def wald_pdf(t, drift, bound, shift):
    d = np.asarray(t, dtype=float) - shift
    out = np.zeros_like(d)
    ok = d > 1e-6
    d = d[ok]
    out[ok] = (bound / np.sqrt(2 * np.pi * d ** 3)
               * np.exp(-((drift * d - bound) ** 2) / (2 * d)))
    return out


def wald_fit(x):
    def nll(p):
        drift, bound, shift = np.exp(p[0]), np.exp(p[1]), p[2]
        if shift <= 0 or shift >= x.min() - 1e-4:
            return 1e10
        dens = wald_pdf(x, drift, bound, shift)
        if np.any(dens <= 0):
            return 1e10
        return -np.log(dens).sum()

    res = minimize(nll, [np.log(3.0), np.log(1.2), x.min() * 0.6],
                   method="Nelder-Mead",
                   options=dict(maxiter=4000, xatol=1e-6, fatol=1e-6))
    return np.exp(res.x[0]), np.exp(res.x[1]), res.x[2]


# --------------------------------------------------------------------------
# panel helpers
# --------------------------------------------------------------------------
def bare(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def rt_axis(ax, rt, labelled, alpha=0.20, bins=70):
    """Observed RTs for both conditions, as a backdrop."""
    edges = np.linspace(*XLIM, bins)
    for cond, col in (("Accuracy", ACC), ("Speed", SPD)):
        ax.hist(rt[cond], bins=edges, density=True, color=col, alpha=alpha,
                edgecolor="none", zorder=1)
    ax.set_xlim(*XLIM)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(PALE)
    ax.tick_params(labelsize=8, colors=MUTED, length=2.5)
    if labelled:
        ax.set_xlabel("Reaction time (s)", fontsize=8.5, color=INK)
    else:
        ax.set_xticklabels([])


def note(ax, text, y=0.99):
    ax.text(0.985, y, text, transform=ax.transAxes, fontsize=7.6,
            style="italic", color=MUTED, ha="right", va="top", zorder=Z_TEXT,
            path_effects=HALO)


def heading(fig, x, panel_top, title, subtitle):
    fig.text(x, panel_top + 0.050, title, fontsize=10, weight="bold",
             color=INK, ha="left", va="baseline")
    fig.text(x, panel_top + 0.018, subtitle, fontsize=7.4, style="italic",
             color=MUTED, ha="left", va="baseline")


def magnifier(fig, cx, cy):
    """A magnifying glass, drawn in its own square axes."""
    w, h = 0.052, 0.052 * W_IN / H_IN
    ax = fig.add_axes([cx - w / 2, cy - h / 2, w, h], zorder=8)
    bare(ax)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect("equal")
    ax.patch.set_visible(False)
    ax.plot([0.26, 0.80], [-0.26, -0.80], color=INK, lw=2.4, zorder=3,
            solid_capstyle="round")
    ax.add_patch(Circle((-0.12, 0.12), 0.62, facecolor="white",
                        edgecolor=INK, lw=1.7, zorder=4))
    ax.add_patch(Circle((-0.12, 0.12), 0.62, facecolor=ACCENT, alpha=0.10,
                        edgecolor="none", zorder=5))
    ax.plot([-0.42, -0.16], [0.36, 0.44], color=ACCENT, lw=1.2, alpha=0.7,
            zorder=6, solid_capstyle="round")


# --------------------------------------------------------------------------
# left column: what the model is made of
# --------------------------------------------------------------------------
def sketch_location(ax):
    """A single point on a line."""
    bare(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.plot([0.8, 9.2], [2.4, 2.4], color=PALE, lw=1.4, solid_capstyle="round")
    ax.plot([5.0], [2.4], marker="o", ms=8, color=ACCENT, zorder=4)
    ax.annotate("", xy=(5.0, 2.4), xytext=(5.0, 4.9),
                arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=1.2,
                                shrinkA=0, shrinkB=2))
    ax.text(5.0, 5.2, r"$\mu$", fontsize=11, color=ACCENT, ha="center",
            va="bottom")
    ax.text(5.0, 1.0, "One number per condition", fontsize=7.4, color=MUTED,
            ha="center", va="center", style="italic")


def sketch_shape(ax):
    """One skewed density, with its three parameters marked on it."""
    bare(ax)
    mu, sg, ta = 3.2, 0.72, 1.9
    x = np.linspace(0, 10, 500)
    dens = exponnorm.pdf(x, ta / sg, loc=mu, scale=sg)
    peak = dens.max()
    ax.fill_between(x, 0, dens, color=ACCENT, alpha=0.16, zorder=2)
    ax.plot(x, dens, color=ACCENT, lw=1.5, zorder=3)
    ax.set_xlim(0, 10)
    ax.set_ylim(-peak * 0.34, peak * 1.62)

    ax.plot([mu, mu], [0, peak * 1.14], color=MUTED, lw=0.7, ls=":", zorder=1)
    ax.text(mu, peak * 1.18, r"$\mu$", fontsize=10, color=INK, ha="center",
            va="bottom", zorder=Z_TEXT, path_effects=HALO)
    ysig = peak * 0.52
    ax.annotate("", xy=(mu - sg, ysig), xytext=(mu + sg, ysig), zorder=Z_ARROW,
                arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=0.9,
                                shrinkA=0, shrinkB=0))
    ax.text(mu + sg + 0.25, ysig, r"$\sigma$", fontsize=9.5, color=INK,
            ha="left", va="center", zorder=Z_TEXT, path_effects=HALO)
    ytau = peak * 0.24
    ax.annotate("", xy=(9.5, ytau), xytext=(mu + 1.5, ytau), zorder=Z_ARROW,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9,
                                shrinkA=0, shrinkB=0))
    ax.text((mu + 1.5 + 9.5) / 2, ytau * 1.30, r"$\tau$", fontsize=10,
            color=INK, ha="center", va="bottom", zorder=Z_TEXT,
            path_effects=HALO)
    ax.text(5.0, -peak * 0.20, "A location, a width and a tail", fontsize=7.4,
            color=MUTED, ha="center", va="center", style="italic")


def sketch_process(ax, drift, bound, shift, seed=4):
    """Noisy accumulation to a boundary; the crossing times are the data."""
    bare(ax)
    rng = np.random.default_rng(seed)
    tmax, dt = 0.95, 0.001
    n = int((tmax - shift) / dt)

    paths = []                               # keep trials that finish in frame
    while len(paths) < 7:
        path = np.cumsum(drift * dt + rng.normal(0, np.sqrt(dt), n))
        idx = int(np.argmax(path >= bound))
        if path[idx] < bound:
            continue
        paths.append(path[:idx + 1])
    rug = min(p.min() for p in paths) - 0.26

    ax.axvspan(0, shift, color=PALE, alpha=0.28, lw=0, zorder=1)
    ax.plot([shift, tmax], [bound, bound], color=INK, lw=1.3, zorder=4)
    ax.plot([shift, tmax], [0, 0], color=PALE, lw=1.0, zorder=2)

    hits = []
    for path in paths:
        t = shift + np.arange(len(path)) * dt
        ax.plot(t, path, color=ACCENT, lw=0.7, alpha=0.6, zorder=3)
        hits.append(t[-1])

    ax.plot([shift, tmax], [rug, rug], color=PALE, lw=0.9, zorder=2)
    for h in hits:
        ax.plot([h], [bound], marker="v", ms=3.0, color=ACCENT, zorder=5)
        ax.plot([h, h], [rug - 0.08, rug + 0.08], color=ACCENT, lw=1.0, zorder=5)
        ax.plot([h, h], [rug + 0.08, bound], color=ACCENT, lw=0.5, ls=":",
                alpha=0.3, zorder=2)

    ax.set_xlim(0, tmax)
    ax.set_ylim(rug - 0.44, bound * 1.34)
    ax.text(tmax * 0.99, bound * 1.04, "boundary", fontsize=7.4, color=INK,
            ha="right", va="bottom", zorder=Z_TEXT, path_effects=HALO)
    ax.annotate("", xy=(shift + 0.19, drift * 0.19), xytext=(shift, 0),
                zorder=Z_ARROW,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8,
                                mutation_scale=6, shrinkA=0, shrinkB=0))
    ax.text(shift + 0.20, drift * 0.19 + 0.06, "drift", fontsize=7.4,
            color=INK, ha="left", va="bottom", zorder=Z_TEXT,
            path_effects=HALO)
    ax.text(shift / 2, (rug + bound) / 2, "non-decision time", fontsize=7.0,
            color=MUTED, ha="center", va="center", rotation=90, zorder=Z_TEXT,
            path_effects=HALO)
    ax.text(tmax * 0.99, rug - 0.14, "response times", fontsize=7.2,
            color=MUTED, ha="right", va="top", style="italic")


# --------------------------------------------------------------------------
def main():
    rt = load()
    n = sum(len(v) for v in rt.values())
    fig = plt.figure(figsize=(W_IN, H_IN))
    fig.patch.set_facecolor("white")
    fig.text(COL_DATA[0], 0.955, TITLE, fontsize=13.5, weight="bold",
             color=INK, ha="left", va="baseline")
    fig.text(COL_DATA[0], 0.915, SUBTITLE, fontsize=8.6, style="italic",
             color=MUTED, ha="left", va="baseline")

    # -- the data, once, on the left ---------------------------------------
    ax_data = fig.add_axes([COL_DATA[0], DATA_Y[0],
                            COL_DATA[1] - COL_DATA[0], DATA_Y[1] - DATA_Y[0]])
    rt_axis(ax_data, rt, labelled=True, alpha=0.62, bins=44)
    heading(fig, COL_DATA[0], DATA_Y[1], "Observed data",
            f"{n:,} decisions")

    # -- three rows --------------------------------------------------------
    rows = []
    for y0 in ROW_Y0:
        axm = fig.add_axes([COL_MECH[0], y0, COL_MECH[1] - COL_MECH[0], PANEL_H])
        axd = fig.add_axes([COL_DENS[0], y0, COL_DENS[1] - COL_DENS[0], PANEL_H])
        rows.append((axm, axd))
    for i, (_, axd) in enumerate(rows):
        rt_axis(axd, rt, labelled=i == 2)

    # -- arrows and magnifying glasses -------------------------------------
    src = (ARROW_X[0], (DATA_Y[0] + DATA_Y[1]) / 2)
    for y0 in ROW_Y0:
        dst = (ARROW_X[1], y0 + PANEL_H / 2)
        fig.add_artist(FancyArrowPatch(
            src, dst, transform=fig.transFigure, arrowstyle="-|>",
            mutation_scale=13, lw=1.4, color=MUTED, alpha=0.85,
            shrinkA=2, shrinkB=3, zorder=6,
            connectionstyle="arc3,rad=0.0",
        ))
        magnifier(fig, (src[0] + dst[0]) / 2, (src[1] + dst[1]) / 2)

    # -- row 1: summary statistics -----------------------------------------
    axm, ax = rows[0]
    sketch_location(axm)
    heading(fig, COL_MECH[0], ROW_Y0[0] + PANEL_H,
            "Summary statistics approach", "The model is a location")

    resid = np.concatenate([v - v.mean() for v in rt.values()])
    sd = np.sqrt((resid ** 2).sum() / (n - 2))
    for cond, col in (("Accuracy", ACC), ("Speed", SPD)):
        m = rt[cond].mean()
        ax.plot(GRID, norm.pdf(GRID, m, sd), color=col, lw=1.7, zorder=4)
        ax.axvline(m, color=col, lw=0.9, ls=(0, (3, 2)), zorder=3)
    top = ax.get_ylim()[1]
    ma, ms = rt["Accuracy"].mean(), rt["Speed"].mean()
    ax.annotate("", xy=(ma, top * 0.70), xytext=(ms, top * 0.70),
                zorder=Z_ARROW,
                arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.0,
                                shrinkA=0, shrinkB=0))
    ax.text((ma + ms) / 2, top * 0.80, "141 ms", fontsize=8, color=INK,
            ha="center", va="bottom", weight="bold", zorder=Z_TEXT,
            path_effects=HALO)
    note(ax, "Symmetric, equally wide,\nwith mass where none occurred")

    # -- row 2: distributional ---------------------------------------------
    axm, ax = rows[1]
    sketch_shape(axm)
    heading(fig, COL_MECH[0], ROW_Y0[1] + PANEL_H,
            "Distributional approach", "The model is a shape")

    # Posterior means from the ex-Gaussian of Example 2 (mu, sigma, tau), in
    # seconds. `mu` is on an identity link, so it is read off the posterior
    # directly rather than through softplus. Regenerate with make_tables.R if
    # the model is refitted.
    pars = {"Accuracy": (0.476, 0.058, 0.213), "Speed": (0.429, 0.051, 0.120)}
    for cond, col in (("Accuracy", ACC), ("Speed", SPD)):
        mu, sg, ta = pars[cond]
        ax.plot(GRID, exponnorm.pdf(GRID, ta / sg, loc=mu, scale=sg),
                color=col, lw=1.7, zorder=4)
    for (x0, x1), y, lab in (
        ((0.245, 0.175), 0.52, r"bulk: 476 $\rightarrow$ 429 ms"),
        ((0.560, 0.430), 0.20, r"tail: 213 $\rightarrow$ 120 ms"),
    ):
        ax.annotate("", xy=(x1, y), xytext=(x0, y), xycoords="axes fraction",
                    textcoords="axes fraction", zorder=Z_ARROW,
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0))
        ax.text(x0 + 0.015, y, lab, fontsize=7.6, transform=ax.transAxes,
                color=INK, ha="left", va="center", zorder=Z_TEXT,
                path_effects=HALO)
    note(ax, "The same 141 ms, split into\na shift and a shorter tail")

    # -- row 3: computational ----------------------------------------------
    #
    # Posterior means from the shifted Wald of Example 2 (drift, boundary,
    # ndt), in seconds, on the response scale. Regenerate with make_tables.R.
    #
    # These are the *fitted* values, not a maximum-likelihood fit to the same
    # trials, and the difference matters. Four of these 4,285 responses are
    # faster than 200 ms; the fitted model hands all four to its outlier
    # component (p_outlier = 1.00 for each), while an MLE shifted Wald has no
    # such component and is forced to drop `ndt` to zero to give them positive
    # density. `wald_fit()` below is kept for reference but is no longer what
    # the figure draws.
    axm, ax = rows[2]
    fits = {"Accuracy": (2.732, 1.003, 0.323), "Speed": (3.853, 0.965, 0.299)}
    sketch_process(axm, *fits["Accuracy"])
    heading(fig, COL_MECH[0], ROW_Y0[2] + PANEL_H,
            "Computational approach", "The model is a process")

    for cond, col in (("Accuracy", ACC), ("Speed", SPD)):
        ax.plot(GRID, wald_pdf(GRID, *fits[cond]), color=col, lw=1.7, zorder=4)
    note(ax, "The same 141 ms, expressed as\nevidence, threshold and encoding")

    # -- legend, tucked under the data it describes -------------------------
    handles = [plt.Line2D([], [], color=c, lw=2.4, label=l)
               for l, c in (("Accuracy instructions", ACC),
                            ("Speed instructions", SPD))]
    ax_data.legend(handles=handles, loc="upper center", frameon=False,
                   bbox_to_anchor=(0.5, -0.17), fontsize=8, handlelength=1.5,
                   handletextpad=0.6, labelspacing=0.5, borderpad=0)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=DPI, facecolor="white")
    print("wrote", OUT, "| Wald fits:",
          {k: tuple(round(v, 3) for v in p) for k, p in fits.items()})


if __name__ == "__main__":
    main()
