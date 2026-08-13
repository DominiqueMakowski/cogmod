import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import exponnorm

from _anim_utils import (
    build_path,
    build_softplus_param_panels,
    softplus_inv,
)

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
N_FRAMES = 150
FPS = 20
TRAIL_LEN = 25
LINK_LABEL = "Softplus link - log(1+e\u1dbb)"
GRID_N = 400  # resolution for all curves/grids
TAU_SIZE_RANGE = (6, 26)  # marker diameter (pt) mapped from tau

# ----------------------------------------------------------------------
# Priors (linked space): all three ExGaussian parameters (mu, sigma, tau)
# are given a Normal prior on an unconstrained ("linked") scale and reach
# the world through a softplus link, keeping them positive as required for
# reaction-time modelling.
# ----------------------------------------------------------------------
MU_PRIOR_LOC, MU_PRIOR_SCALE = 0, 1
SIGMA_PRIOR_LOC, SIGMA_PRIOR_SCALE = 0, 1
TAU_PRIOR_LOC, TAU_PRIOR_SCALE = 0, 1

# ----------------------------------------------------------------------
# Path the (mu, sigma, tau) point takes through parameter space: a loop
# through 5 landmark stops, closing back on the start. All three stay
# strictly positive, as required by the softplus link.
# ----------------------------------------------------------------------
LANDMARKS = [
    (0.6, 0.5, 0.3),
    (2.2, 0.5, 0.3),
    (2.2, 1.6, 2.2),
    (0.5, 1.6, 2.2),
    (0.5, 0.35, 0.8),
]
KEYFRAME_PATH = LANDMARKS + [LANDMARKS[0]]

mu_path, sigma_path, tau_path = build_path(KEYFRAME_PATH, N_FRAMES)

# ----------------------------------------------------------------------
# Automatic axis ranges
# ----------------------------------------------------------------------
mu_landmarks = [lm[0] for lm in LANDMARKS]
sigma_landmarks = [lm[1] for lm in LANDMARKS]
tau_landmarks = [lm[2] for lm in LANDMARKS]

# parameter-space plane: mu and sigma are both softplus-linked, so both
# have an invalid (<= 0) strip
MU_MAX = max(mu_landmarks)
MU_RANGE = (-0.15 * MU_MAX, MU_MAX * 1.15)
SIGMA_MAX = max(sigma_landmarks)
SIGMA_RANGE = (-0.15 * SIGMA_MAX, SIGMA_MAX * 1.15)

# tau is shown as marker size rather than as a plane axis
TAU_LO, TAU_HI = min(tau_landmarks), max(tau_landmarks)


def size_for_tau(tau):
    frac = (tau - TAU_LO) / (TAU_HI - TAU_LO)
    return TAU_SIZE_RANGE[0] + np.clip(frac, 0, 1) * (
        TAU_SIZE_RANGE[1] - TAU_SIZE_RANGE[0]
    )


# density plot: wide enough for the fullest ExGaussian on the path
means = [mu + tau for mu, _, tau in LANDMARKS]
sds = [np.sqrt(sigma**2 + tau**2) for _, sigma, tau in LANDMARKS]
x_hi = max(m + 4.5 * s for m, s in zip(means, sds))
x_lo = min(mu - 4.5 * sigma for mu, sigma, _ in LANDMARKS)
X_RANGE = np.linspace(min(x_lo, 0), x_hi, GRID_N)
peak_density = max(
    exponnorm.pdf(X_RANGE, tau / sigma, loc=mu, scale=sigma).max()
    for mu, sigma, tau in LANDMARKS
)

# ----------------------------------------------------------------------
# Figure & layout
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(11, 13))
outer = gridspec.GridSpec(
    2,
    1,
    height_ratios=[0.75, 1.6],
    hspace=0.45,
    top=0.93,
    bottom=0.05,
    left=0.09,
    right=0.95,
)
top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=0.3)
bottom_gs = gridspec.GridSpecFromSubplotSpec(
    3, 2, subplot_spec=outer[1], wspace=0.65, hspace=0.9
)

ax_density = fig.add_subplot(top_gs[0, 0])
ax_plane = fig.add_subplot(top_gs[0, 1])
ax_mu_linked = fig.add_subplot(bottom_gs[0, 0])
ax_mu_unlinked = fig.add_subplot(bottom_gs[0, 1])
ax_sigma_linked = fig.add_subplot(bottom_gs[1, 0])
ax_sigma_unlinked = fig.add_subplot(bottom_gs[1, 1])
ax_tau_linked = fig.add_subplot(bottom_gs[2, 0])
ax_tau_unlinked = fig.add_subplot(bottom_gs[2, 1])

fig.suptitle("The ExGaussian Model", fontsize=14, fontweight="bold")

# ---- Top-left: density -------------------------------------------------
ax_density.set_xlim(X_RANGE.min(), X_RANGE.max())
ax_density.set_ylim(0, peak_density * 1.15)
ax_density.set_title("The ExGaussian Distribution")
ax_density.set_xlabel("x (RT)")
ax_density.set_ylabel("density")
ax_density.set_yticks([])
(line_density,) = ax_density.plot([], [], lw=2.5, color="#1f77b4")
text_params = ax_density.text(
    0.97,
    0.93,
    "",
    transform=ax_density.transAxes,
    fontsize=10,
    va="top",
    ha="right",
    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85),
)
for lm_mu, lm_sigma, lm_tau in LANDMARKS:
    ax_density.plot(
        X_RANGE,
        exponnorm.pdf(X_RANGE, lm_tau / lm_sigma, loc=lm_mu, scale=lm_sigma),
        color="#1f77b4",
        lw=1.2,
        alpha=0.25,
        zorder=1,
    )

# ---- Top-right: parameter plane (mu, sigma), tau shown as marker size --
ax_plane.set_xlim(*MU_RANGE)
ax_plane.set_ylim(*SIGMA_RANGE)
ax_plane.set_facecolor("#8fd694")  # valid region (mu >= 0 and sigma >= 0)
ax_plane.axhspan(SIGMA_RANGE[0], 0, facecolor="#d9534f", alpha=0.6, zorder=0)
ax_plane.axvspan(MU_RANGE[0], 0, facecolor="#d9534f", alpha=0.6, zorder=0)
ax_plane.set_title("Parameter space (\u03bc, \u03c3) \u2014 size = \u03c4")
ax_plane.set_xlabel("location (\u03bc)")
ax_plane.set_ylabel("SD (\u03c3)")
(trail_plane,) = ax_plane.plot([], [], "-", color="black", lw=1, alpha=0.35, zorder=4)
(point_plane,) = ax_plane.plot(
    [], [], "o", color="black", mec="white", mew=1.5, zorder=5
)

ax_plane.plot(
    [p[0] for p in KEYFRAME_PATH],
    [p[1] for p in KEYFRAME_PATH],
    "--",
    color="black",
    lw=1,
    alpha=0.25,
    zorder=2,
)
ax_plane.plot(
    [p[0] for p in LANDMARKS],
    [p[1] for p in LANDMARKS],
    "o",
    color="black",
    ms=6,
    alpha=0.35,
    zorder=3,
)

legend_taus = [TAU_LO, (TAU_LO + TAU_HI) / 2, TAU_HI]
legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="black",
        markeredgecolor="white",
        markersize=size_for_tau(t),
        label=f"\u03c4={t:.2g}",
    )
    for t in legend_taus
]
ax_plane.legend(
    handles=legend_handles,
    loc="upper right",
    fontsize=7,
    title="size = \u03c4",
    title_fontsize=7,
    framealpha=0.85,
)

# ---- Bottom: priors, one linked/unlinked row per parameter --------------
mu_panels = build_softplus_param_panels(
    fig,
    ax_mu_linked,
    ax_mu_unlinked,
    "\u03bc",
    MU_PRIOR_LOC,
    MU_PRIOR_SCALE,
    mu_landmarks,
    "#9467bd",
    LINK_LABEL,
    grid_n=GRID_N,
)
sigma_panels = build_softplus_param_panels(
    fig,
    ax_sigma_linked,
    ax_sigma_unlinked,
    "\u03c3",
    SIGMA_PRIOR_LOC,
    SIGMA_PRIOR_SCALE,
    sigma_landmarks,
    "#d62728",
    LINK_LABEL,
    grid_n=GRID_N,
)
tau_panels = build_softplus_param_panels(
    fig,
    ax_tau_linked,
    ax_tau_unlinked,
    "\u03c4",
    TAU_PRIOR_LOC,
    TAU_PRIOR_SCALE,
    tau_landmarks,
    "#2ca02c",
    LINK_LABEL,
    grid_n=GRID_N,
)

PARAMS = [
    (mu_path, mu_panels),
    (sigma_path, sigma_panels),
    (tau_path, tau_panels),
]

# ----------------------------------------------------------------------
# Animation
# ----------------------------------------------------------------------
ALL_MARKERS = [point_plane] + [
    a
    for panels in (mu_panels, sigma_panels, tau_panels)
    for a in (panels["marker_linked"], panels["marker_unlinked"], panels["marker_link"])
]
ALL_VLINES = [
    a
    for panels in (mu_panels, sigma_panels, tau_panels)
    for a in (panels["vline_linked"], panels["vline_unlinked"])
]
ALL_ARTISTS = [line_density, trail_plane, text_params] + ALL_MARKERS + ALL_VLINES


def init():
    line_density.set_data([], [])
    trail_plane.set_data([], [])
    for marker in ALL_MARKERS:
        marker.set_data([], [])
    text_params.set_text("")
    return ALL_ARTISTS


def update(frame):
    mu, sigma, tau = mu_path[frame], sigma_path[frame], tau_path[frame]

    line_density.set_data(
        X_RANGE, exponnorm.pdf(X_RANGE, tau / sigma, loc=mu, scale=sigma)
    )
    text_params.set_text(f"\u03bc = {mu:.2f}\n\u03c3 = {sigma:.2f}\n\u03c4 = {tau:.2f}")

    point_plane.set_data([mu], [sigma])
    point_plane.set_markersize(size_for_tau(tau))
    idx = np.arange(frame - TRAIL_LEN, frame + 1) % N_FRAMES
    trail_plane.set_data(mu_path[idx], sigma_path[idx])

    for path, panels in PARAMS:
        value = path[frame]
        linked_value = float(softplus_inv(value))
        panels["marker_linked"].set_data(
            [linked_value], [panels["prior_link"].pdf(linked_value)]
        )
        panels["vline_linked"].set_xdata([linked_value, linked_value])
        panels["marker_unlinked"].set_data(
            [value], [panels["prior_natural"].pdf(value)]
        )
        panels["vline_unlinked"].set_xdata([value, value])
        panels["marker_link"].set_data([linked_value], [value])

    return ALL_ARTISTS


anim = FuncAnimation(
    fig, update, frames=N_FRAMES, init_func=init, blit=False, interval=1000 / FPS
)
anim.save("anim_priors_exgaussian.gif", writer=PillowWriter(fps=FPS))
print("Saved figure")
