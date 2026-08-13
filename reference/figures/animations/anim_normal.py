import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.stats import norm

from _anim_utils import (
    HalfStudentT,
    build_path,
    padded_range,
    prior_support,
    setup_prior_panel,
    add_link_panel,
)

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
N_FRAMES = 150
FPS = 20
TRAIL_LEN = 25
MU_LINK_LABEL = "Identity link - 1:1"
SIGMA_LINK_LABEL = "Identity link - 1:1"
GRID_N = 400  # resolution for all curves/grids

# ----------------------------------------------------------------------
# Priors (these define the curves in the bottom row)
# ----------------------------------------------------------------------
MU_PRIOR_LOC, MU_PRIOR_SCALE = 0, 1
prior_mu = norm(loc=MU_PRIOR_LOC, scale=MU_PRIOR_SCALE)

# sigma is treated the way `brms` treats scale parameters by default: no
# link at all (identity), the parameter being declared directly on its
# natural scale with a lower bound of 0, and positivity enforced by a
# *truncated* prior - the half Student-t student_t(3, 0, 2.5). The
# "linked" and "unlinked" spaces are therefore one and the same.
SIGMA_PRIOR_DF, SIGMA_PRIOR_LOC, SIGMA_PRIOR_SCALE = 3, 0, 2.5
prior_sigma = HalfStudentT(
    df=SIGMA_PRIOR_DF, loc=SIGMA_PRIOR_LOC, scale=SIGMA_PRIOR_SCALE
)

# ----------------------------------------------------------------------
# Path the (mu, sigma) point takes through parameter space: a loop
# through 5 landmark stops, closing back on the centre.
# ----------------------------------------------------------------------
LANDMARKS = [(0.0, 1.15), (3.0, 1.15), (3.0, 2.0), (-3.0, 2.0), (-3.0, 0.5)]
KEYFRAME_PATH = LANDMARKS + [LANDMARKS[0]]

mu_path, sigma_path = build_path(KEYFRAME_PATH, N_FRAMES)

# ----------------------------------------------------------------------
# Automatic axis ranges,
# ----------------------------------------------------------------------
mu_landmarks = [lm[0] for lm in LANDMARKS]
sigma_landmarks = [lm[1] for lm in LANDMARKS]

# parameter-space plane: padded to the landmarks' own range
MU_RANGE = padded_range(
    min(mu_landmarks), max(mu_landmarks), pad_frac=0.25, min_pad=0.5
)
SIGMA_MAX = max(sigma_landmarks)
SIGMA_RANGE = (
    -0.15 * SIGMA_MAX,
    SIGMA_MAX * 1.15,
)  # small negative strip for the invalid zone

# prior panels: wide enough for each prior's bulk and for every landmark
mu_grid = np.linspace(*prior_support(prior_mu, mu_landmarks), GRID_N)

# fixed window rather than the prior quantiles: the half-t is heavy-tailed,
# and reaching well below 0 makes its truncation visible (the density is
# exactly zero to the left of 0). Identity link, so one grid serves both
# panels; 0 itself is forced into the grid to keep that edge crisp.
SIGMA_VIEW = (-3.0, 6.0)
sigma_grid = np.unique(np.concatenate([np.linspace(*SIGMA_VIEW, GRID_N), [0.0]]))

# density plot: wide enough for the fullest bell curve on the path (~4.5 SDs)
x_half_width = max(abs(v) for v in mu_landmarks) + SIGMA_MAX * 4.5
X_RANGE = np.linspace(-x_half_width, x_half_width, GRID_N)

# ----------------------------------------------------------------------
# Figure & layout
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(11, 10))
outer = gridspec.GridSpec(
    2,
    1,
    height_ratios=[1.05, 1],
    hspace=0.55,
    top=0.90,
    bottom=0.07,
    left=0.09,
    right=0.95,
)
top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=0.3)
bottom_gs = gridspec.GridSpecFromSubplotSpec(
    2, 2, subplot_spec=outer[1], wspace=0.65, hspace=0.85
)

ax_density = fig.add_subplot(top_gs[0, 0])
ax_plane = fig.add_subplot(top_gs[0, 1])
ax_mu_linked = fig.add_subplot(bottom_gs[0, 0])
ax_mu_unlinked = fig.add_subplot(bottom_gs[0, 1])
ax_sigma_linked = fig.add_subplot(bottom_gs[1, 0])
ax_sigma_unlinked = fig.add_subplot(bottom_gs[1, 1])

fig.suptitle("The Gaussian Model", fontsize=14, fontweight="bold")

# ---- Top-left: density -------------------------------------------------
ax_density.set_xlim(X_RANGE.min(), X_RANGE.max())
peak_density = norm.pdf(
    0, scale=sigma_path.min()
)  # tallest curve occurs at the smallest sigma
ax_density.set_ylim(0, peak_density * 1.15)
ax_density.set_title("The Normal Distribution")
ax_density.set_xlabel("x")
ax_density.set_ylabel("density")
ax_density.set_yticks([])
(line_density,) = ax_density.plot([], [], lw=2.5, color="#1f77b4")
text_params = ax_density.text(
    0.03,
    0.93,
    "",
    transform=ax_density.transAxes,
    fontsize=10,
    va="top",
    ha="left",
    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85),
)
for lm_mu, lm_sigma in LANDMARKS:
    ax_density.plot(
        X_RANGE,
        norm.pdf(X_RANGE, lm_mu, lm_sigma),
        color="#1f77b4",
        lw=1.2,
        alpha=0.25,
        zorder=1,
    )

# ---- Top-right: parameter plane ----------------------------------------
ax_plane.set_xlim(*MU_RANGE)
ax_plane.set_ylim(*SIGMA_RANGE)
ax_plane.set_facecolor("#8fd694")  # valid region (sigma >= 0)
ax_plane.axhspan(
    SIGMA_RANGE[0], 0, facecolor="#d9534f", alpha=0.6, zorder=0
)  # invalid: sigma < 0
ax_plane.set_title("Parameter space (\u03bc, \u03c3)")
ax_plane.set_xlabel("location (\u03bc)")
ax_plane.set_ylabel("SD (\u03c3)")
(trail_plane,) = ax_plane.plot([], [], "-", color="black", lw=1, alpha=0.35, zorder=4)
(point_plane,) = ax_plane.plot(
    [], [], "o", color="black", ms=9, mec="white", mew=1.5, zorder=5
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


# ---- Bottom: priors -----------------------------------------------------
MU_DIST_LABEL = f"N({MU_PRIOR_LOC:g}, {MU_PRIOR_SCALE:g})"
SIGMA_DIST_LABEL = (
    f"half-t({SIGMA_PRIOR_DF:g}, {SIGMA_PRIOR_LOC:g}, {SIGMA_PRIOR_SCALE:g})"
)

marker_mu_linked, vline_mu_linked = setup_prior_panel(
    ax_mu_linked,
    "Prior for \u03bc — linked space",
    mu_grid,
    prior_mu,
    "#9467bd",
    "\u03bc",
    mu_landmarks,
    MU_DIST_LABEL,
)
marker_mu_unlinked, vline_mu_unlinked = setup_prior_panel(
    ax_mu_unlinked,
    "Prior for \u03bc — unlinked space",
    mu_grid,
    prior_mu,
    "#9467bd",
    "\u03bc",
    mu_landmarks,
    MU_DIST_LABEL,
)
marker_sigma_linked, vline_sigma_linked = setup_prior_panel(
    ax_sigma_linked,
    "Prior for \u03c3 — linked space",
    sigma_grid,
    prior_sigma,
    "#d62728",
    "\u03c3",
    sigma_landmarks,
    SIGMA_DIST_LABEL,
)
marker_sigma_unlinked, vline_sigma_unlinked = setup_prior_panel(
    ax_sigma_unlinked,
    "Prior for \u03c3 — unlinked space",
    sigma_grid,
    prior_sigma,
    "#d62728",
    "\u03c3",
    sigma_landmarks,
    SIGMA_DIST_LABEL,
)


# ---- Connectors between linked <-> unlinked panels ----------------------
marker_mu_link = add_link_panel(
    fig,
    ax_mu_linked,
    ax_mu_unlinked,
    MU_LINK_LABEL,
    lambda z: z,
    mu_grid,
    (mu_grid.min(), mu_grid.max()),
    "#9467bd",
)
marker_sigma_link = add_link_panel(
    fig,
    ax_sigma_linked,
    ax_sigma_unlinked,
    SIGMA_LINK_LABEL,
    lambda z: z,
    sigma_grid,
    (sigma_grid.min(), sigma_grid.max()),
    "#d62728",
)

# ----------------------------------------------------------------------
# Animation
# ----------------------------------------------------------------------
MU_PAIRS = [
    (marker_mu_linked, vline_mu_linked),
    (marker_mu_unlinked, vline_mu_unlinked),
]
ALL_MARKERS = [m for m, _ in MU_PAIRS] + [
    marker_sigma_linked,
    marker_sigma_unlinked,
    marker_mu_link,
    marker_sigma_link,
]

ALL_ARTISTS = (
    [line_density, point_plane, trail_plane, text_params]
    + [a for pair in MU_PAIRS for a in pair]
    + [
        marker_sigma_linked,
        vline_sigma_linked,
        marker_sigma_unlinked,
        vline_sigma_unlinked,
        marker_mu_link,
        marker_sigma_link,
    ]
)


def init():
    line_density.set_data([], [])
    point_plane.set_data([], [])
    trail_plane.set_data([], [])
    for marker in ALL_MARKERS:
        marker.set_data([], [])
    text_params.set_text("")
    return ALL_ARTISTS


def update(frame):
    mu, sigma = mu_path[frame], sigma_path[frame]

    line_density.set_data(X_RANGE, norm.pdf(X_RANGE, mu, sigma))
    text_params.set_text(f"\u03bc = {mu:.2f}\n\u03c3 = {sigma:.2f}")

    point_plane.set_data([mu], [sigma])
    idx = np.arange(frame - TRAIL_LEN, frame + 1) % N_FRAMES
    trail_plane.set_data(mu_path[idx], sigma_path[idx])

    mu_pdf_val = prior_mu.pdf(mu)
    for marker, vline in MU_PAIRS:
        marker.set_data([mu], [mu_pdf_val])
        vline.set_xdata([mu, mu])
    marker_mu_link.set_data([mu], [mu])  # identity link: natural == linked

    sigma_pdf_val = prior_sigma.pdf(sigma)
    marker_sigma_linked.set_data([sigma], [sigma_pdf_val])
    vline_sigma_linked.set_xdata([sigma, sigma])
    marker_sigma_unlinked.set_data([sigma], [sigma_pdf_val])
    vline_sigma_unlinked.set_xdata([sigma, sigma])
    marker_sigma_link.set_data([sigma], [sigma])  # identity link: natural == linked

    return ALL_ARTISTS


anim = FuncAnimation(
    fig, update, frames=N_FRAMES, init_func=init, blit=False, interval=1000 / FPS
)
anim.save("anim_priors_normal.gif", writer=PillowWriter(fps=FPS))
print("Saved figure")
