"""Where the Shifted Wald comes from.

The Wald is the first-passage-time distribution of a single Wiener
accumulator: evidence starts at 0, drifts upwards at a constant drift rate
while being buffeted by Gaussian noise, and the response is emitted the moment
it touches `boundary`. Shift the whole thing right by a non-decision time
`ndt` and you have `cogmod_invgaussian()`.

The animation makes that literal, and in a single frame of reference. A dot
walks through parameter space; as it walks, the arena re-draws itself - the
boundary rises and falls, the non-decision band widens and narrows, the drift
ray tilts. When the dot reaches a landmark it stops and *fires*: two hundred
accumulators race upwards at once, and the moments at which they touch the
threshold are the response times. The analytic Wald density those hits are
drawn from is traced out **along the threshold line itself**, left to right,
keeping pace with the responses arriving underneath it - which is the whole
claim of the model: this shape is not assumed, it is what the process does.

Each stop keeps its own colour, and its density stays on screen once revealed.
Every curve rides the boundary line as it currently stands, so the three share
a baseline and can be read straight off against each other.

Fixed drift throughout (`sigmadrift = 0`, the classic Wald) and no outlier
component, so the curve really is the distribution of the hits that produced
it.

Within-trial noise is the usual s = 1 convention, which is the one
`cogmod_invgaussian()` uses: the decision time is Inverse Gaussian with mean
`boundary / drift` and shape `boundary^2`. The drift rate is `mu` in the brms
family, since brms requires a parameter of that name, and it is called that in
the code below - but it is labelled "drift rate" everywhere the viewer sees it.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection
from matplotlib.patches import FancyArrowPatch
from matplotlib.animation import FuncAnimation, PillowWriter

from _anim_utils import build_path, padded_range

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
SEED = 42

TRAVEL_FRAMES = 20  # moving from one landmark to the next
FIRE_FRAMES = 26  # the race itself, cursor sweeping ndt -> T_MAX and tracing
#                   the density out behind it
HOLD_FRAMES = 14  # a pause on the finished curve before the parameters move on
FPS = 20

N_PATHS = 200  # accumulators fired per landmark
SIM_DT = 0.005  # Euler step for the Wiener paths, in seconds
T_MAX = 1.6  # right edge of the time axis, in seconds
GRID_N = 500
GHOST_ALPHA = 0.45  # a density that has already been revealed, kept for comparison

# Fixed annotations: the boundary and the drift ray mean the same thing at
# every stop, so they keep the same colour throughout. Everything that belongs
# to one particular stop - its paths and its density - takes that stop's
# colour instead.
COL_BOUNDARY = "#d62728"
COL_DRIFT = "#2ca02c"
COL_NDT = "#7f7f7f"

# ----------------------------------------------------------------------
# The tour through (drift rate, boundary, ndt). Three stops, so three legs,
# and each leg is meant to be read on its own:
#   A -> B  boundary alone
#   B -> C  drift and ndt together, which stays legible because their effects
#           are visually orthogonal - one narrows the density, the other
#           slides it bodily to the right
#   C -> A  the reset
# ----------------------------------------------------------------------
LANDMARKS = [
    (3.0, 0.60, 0.20),  # reference point
    (3.0, 1.20, 0.20),  # boundary up: slower *and* more spread out (caution)
    (6.0, 1.20, 0.45),  # drift up, ndt up: much tighter, and slid to the right
]
COLOURS = ["#1f77b4", "#9467bd", "#e08214"]
assert len(COLOURS) == len(LANDMARKS)
N_LANDMARKS = len(LANDMARKS)

mu_lm = [lm[0] for lm in LANDMARKS]
boundary_lm = [lm[1] for lm in LANDMARKS]
ndt_lm = [lm[2] for lm in LANDMARKS]

rng = np.random.default_rng(SEED)


# ----------------------------------------------------------------------
# The distribution itself
# ----------------------------------------------------------------------
def wald_pdf(t, mu, boundary, ndt):
    """Shifted Wald density: the first passage of a Wiener process with drift
    rate `mu` to `boundary`, displaced by `ndt`. Same parameterisation as
    `dcogmod_invgaussian(..., sigmadrift = 0)`.
    """
    t = np.asarray(t, dtype=float)
    d = t - ndt
    out = np.zeros_like(t)
    ok = d > 0
    dd = d[ok]
    out[ok] = (
        boundary
        / np.sqrt(2 * np.pi * dd**3)
        * np.exp(-((boundary - mu * dd) ** 2) / (2 * dd))
    )
    return out


def simulate_paths(mu, boundary, ndt, n_paths, dt, t_max, generator):
    """Euler-Maruyama on dX = mu*dt + dW, absorbed at `boundary`.

    Returns one (time, evidence) pair per accumulator, already truncated at
    its own hit, the hit times, and the indices of three paths spread over the
    quantiles of those times - the ones worth drawing boldly, since a fast, a
    typical and a slow trial say more together than three arbitrary ones.

    Hit times are NaN for the paths still short of the boundary when the
    window runs out. Those are censored rather than thrown away, which is why
    the response count reported on screen can fall short of `n_paths`.
    """
    n_steps = int(np.ceil((t_max - ndt) / dt))
    steps = mu * dt + np.sqrt(dt) * generator.standard_normal((n_paths, n_steps))
    evidence = np.column_stack([np.zeros(n_paths), np.cumsum(steps, axis=1)])
    times = ndt + np.arange(n_steps + 1) * dt

    crossed = evidence >= boundary
    ever = crossed.any(axis=1)
    # argmax on a boolean row gives the first True, i.e. the absorbing step
    first = np.where(ever, crossed.argmax(axis=1), n_steps)

    xs, ys, hits = [], [], np.full(n_paths, np.nan)
    for i in range(n_paths):
        k = first[i]
        x, y = times[: k + 1].copy(), evidence[i, : k + 1].copy()
        if ever[i]:
            # a whole Euler step can carry the path well past the threshold,
            # which then draws as a spray of lines above a boundary nothing is
            # supposed to cross. Land the last point exactly on it instead, and
            # take the sub-step hit time from the same interpolation.
            gap = y[k] - y[k - 1]
            frac = (boundary - y[k - 1]) / gap if gap > 0 else 1.0
            x[k] = x[k - 1] + frac * dt
            y[k] = boundary
            hits[i] = x[k]
        xs.append(x)
        ys.append(y)

    finite = hits[~np.isnan(hits)]
    if finite.size:
        targets = np.quantile(finite, [0.15, 0.5, 0.85])
        highlight = {int(np.nanargmin(np.abs(hits - q))) for q in targets}
    else:
        highlight = set(range(min(3, n_paths)))
    return xs, ys, hits, highlight


# ----------------------------------------------------------------------
# Frame schedule: fire at a landmark, reveal the density it produced, travel
# to the next, fire again. Each entry carries everything `update()` needs, so
# `update()` stays a pure "draw this state" function.
# ----------------------------------------------------------------------
schedule = []
for i, lm in enumerate(LANDMARKS):
    for j in range(FIRE_FRAMES + HOLD_FRAMES):
        schedule.append(
            dict(
                phase="fire",
                lm=i,
                params=lm,
                # the cursor sweeps ndt -> T_MAX during the race and then sits
                # at the end of the axis for the hold
                race=min(j / (FIRE_FRAMES - 1), 1.0),
            )
        )
    nxt = LANDMARKS[(i + 1) % N_LANDMARKS]
    leg = build_path([lm, nxt], TRAVEL_FRAMES)
    for j in range(TRAVEL_FRAMES):
        schedule.append(
            dict(
                phase="travel",
                lm=i,
                params=(leg[0][j], leg[1][j], leg[2][j]),
                race=0.0,
            )
        )
N_FRAMES = len(schedule)

# One pre-simulated batch per landmark, drawn up front so the frame loop does
# no random number generation and the same GIF comes out of every run.
BATCHES = [simulate_paths(*lm, N_PATHS, SIM_DT, T_MAX, rng) for lm in LANDMARKS]

# ----------------------------------------------------------------------
# Axis ranges. The arena holds two things stacked: the evidence space, which
# reaches up to the tallest boundary, and above it the densities that ride on
# those boundaries. DENSITY_SCALE converts a density into evidence units and
# is shared by every landmark, so the curves stay comparable across stops - a
# tighter distribution really does draw taller.
# ----------------------------------------------------------------------
T_GRID = np.linspace(0, T_MAX, GRID_N)
PEAK_DENSITY = max(wald_pdf(T_GRID, *lm).max() for lm in LANDMARKS)
DENSITY_HEIGHT = 0.80 * max(boundary_lm)
DENSITY_SCALE = DENSITY_HEIGHT / PEAK_DENSITY

EVIDENCE_RANGE = (
    -0.10 * max(boundary_lm),
    max(boundary_lm) + DENSITY_HEIGHT * 1.15,
)
MU_RANGE = padded_range(min(mu_lm), max(mu_lm), pad_frac=0.25, min_pad=0.3)
MU_RANGE = (min(MU_RANGE[0], -0.4), MU_RANGE[1])  # room for the drift <= 0 strip
BOUNDARY_RANGE = (-0.12 * max(boundary_lm), 1.25 * max(boundary_lm))
NDT_RANGE = padded_range(min(ndt_lm), max(ndt_lm), pad_frac=0.35, min_pad=0.05)

# ----------------------------------------------------------------------
# Figure & layout
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(11.5, 7.2))
outer = gridspec.GridSpec(
    1,
    2,
    width_ratios=[2.15, 1],
    wspace=0.24,
    top=0.87,
    bottom=0.10,
    left=0.06,
    right=0.97,
)
right_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=outer[1], height_ratios=[3.2, 1], hspace=0.55
)

ax_arena = fig.add_subplot(outer[0])
ax_plane = fig.add_subplot(right_gs[0])
ax_ndt = fig.add_subplot(right_gs[1])

fig.suptitle("The Shifted Wald Model", fontsize=14, fontweight="bold")

# ---- Left: the arena, which is also the density panel -------------------
ax_arena.set_xlim(0, T_MAX)
ax_arena.set_ylim(*EVIDENCE_RANGE)
ax_arena.set_title(
    "Evidence accumulation, and the response times it produces", fontsize=11
)
ax_arena.set_xlabel("time (s)")
ax_arena.set_ylabel("evidence")
ax_arena.set_yticks([0])

# Nothing accumulates during the non-decision time, so it draws as exactly
# that: evidence pinned flat at zero from 0 to ndt, ending in the start point
# every path leaves from.
(line_ndt,) = ax_arena.plot([], [], "-", color=COL_NDT, lw=3.5, zorder=6)
(point_start,) = ax_arena.plot(
    [], [], "o", color="black", ms=7, mec="white", mew=1.4, zorder=7
)
text_ndt = ax_arena.text(
    0,
    EVIDENCE_RANGE[0] * 0.55,
    "ndt",
    fontsize=9,
    color="#555555",
    va="center",
    ha="center",
    style="italic",
    zorder=6,
)

(line_base,) = ax_arena.plot([], [], "-", color="black", lw=1, alpha=0.5, zorder=2)
(line_boundary,) = ax_arena.plot([], [], "-", color=COL_BOUNDARY, lw=2.5, zorder=8)
text_boundary = ax_arena.text(
    T_MAX * 0.995,
    0,
    "boundary",
    fontsize=9,
    color=COL_BOUNDARY,
    va="top",
    ha="right",
    zorder=8,
)

# The drift ray, and the arc that measures its angle. The arc is a
# FancyArrowPatch rather than a plain Arc so that it ends in a head landing on
# the ray - without one the curve reads as a stray mark rather than as "this
# angle, here".
(ray_drift,) = ax_arena.plot([], [], "--", color=COL_DRIFT, lw=1.8, alpha=0.9, zorder=7)
arc_drift = FancyArrowPatch(
    (0, 0),
    (0, 0),
    connectionstyle="arc3,rad=0",
    arrowstyle="-|>",
    mutation_scale=11,
    shrinkA=0,
    shrinkB=0,
    color=COL_DRIFT,
    lw=1.3,
    zorder=7,
)
ax_arena.add_patch(arc_drift)
text_drift = ax_arena.text(
    0, 0, "", fontsize=10, color=COL_DRIFT, va="center", ha="left", zorder=7
)
# the drift geometry has to stay readable on top of two hundred overlapping
# paths, hence the white outline on all three pieces of it
for artist in (ray_drift, arc_drift, text_drift):
    artist.set_path_effects([pe.withStroke(linewidth=3.2, foreground="white")])

lc_paths = LineCollection([], linewidths=0.6, alpha=0.13, zorder=3)
ax_arena.add_collection(lc_paths)
lc_highlight = LineCollection([], linewidths=1.5, alpha=0.9, zorder=4)
ax_arena.add_collection(lc_highlight)

# One density line per stop. Its shape is fixed from the start - that is the
# stop's own distribution, ndt shift included - but its baseline is not: every
# curve is re-seated on whatever the boundary currently is, so the three share
# a line and compare directly. Only that offset and the alpha are animated.
lines_pdf, PDF_SHAPES = [], []
for lm, colour in zip(LANDMARKS, COLOURS):
    grid = T_GRID[T_GRID >= lm[2]]  # the density is exactly zero before ndt
    PDF_SHAPES.append((grid, wald_pdf(grid, *lm) * DENSITY_SCALE))
    (line,) = ax_arena.plot([], [], lw=1.6, color=colour, alpha=0.0, zorder=6)
    lines_pdf.append(line)

legend = ax_arena.legend(
    lines_pdf,
    [f"drift {m:g}, boundary {b:g}, ndt {n:g}" for m, b, n in LANDMARKS],
    loc="upper right",
    bbox_to_anchor=(1.0, 0.90),
    fontsize=8.5,
    frameon=True,
    framealpha=0.85,
    edgecolor="0.75",
)
legend.set_zorder(10)

(cursor_arena,) = ax_arena.plot([], [], "-", color="black", lw=1, alpha=0.3, zorder=1)

text_params = ax_arena.text(
    0.015,
    0.97,
    "",
    transform=ax_arena.transAxes,
    fontsize=10,
    va="top",
    ha="left",
    family="monospace",
    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.88),
    zorder=10,
)
text_phase = ax_arena.text(
    0.985,
    0.97,
    "",
    transform=ax_arena.transAxes,
    fontsize=10,
    va="top",
    ha="right",
    style="italic",
    color="#444444",
    zorder=10,
)

# ---- Top-right: parameter space -----------------------------------------
ax_plane.set_xlim(*MU_RANGE)
ax_plane.set_ylim(*BOUNDARY_RANGE)
ax_plane.set_facecolor("#8fd694")  # valid region: drift > 0 and boundary > 0
ax_plane.axhspan(BOUNDARY_RANGE[0], 0, facecolor="#d9534f", alpha=0.6, zorder=0)
ax_plane.axvspan(MU_RANGE[0], 0, facecolor="#d9534f", alpha=0.6, zorder=1)
ax_plane.set_title("Parameter space (drift rate, boundary)", fontsize=11)
ax_plane.set_xlabel("drift rate")
ax_plane.set_ylabel("threshold (boundary)")

ax_plane.plot(
    mu_lm + mu_lm[:1],
    boundary_lm + boundary_lm[:1],
    "--",
    color="black",
    lw=1,
    alpha=0.25,
    zorder=2,
)
# the landmarks wear their stop's colour here and on the ndt track, which is
# what ties a dot in parameter space to a curve in the arena
ax_plane.scatter(
    mu_lm, boundary_lm, s=70, c=COLOURS, edgecolors="white", linewidths=1.2, zorder=3
)
(trail_plane,) = ax_plane.plot([], [], "-", color="black", lw=1, alpha=0.35, zorder=4)
(halo_plane,) = ax_plane.plot(
    [], [], "o", color=COL_BOUNDARY, ms=11, alpha=0.0, mfc="none", mew=2, zorder=5
)
(point_plane,) = ax_plane.plot(
    [], [], "o", ms=13, mfc="none", mec="black", mew=2, zorder=6
)

# ---- Bottom-right: the ndt track ----------------------------------------
ax_ndt.set_xlim(*NDT_RANGE)
ax_ndt.set_ylim(0, 1)
ax_ndt.set_yticks([])
ax_ndt.set_xlabel("non-decision time (ndt, s)")
ax_ndt.set_title("Parameter space (ndt)", fontsize=10)
for spine in ("left", "right", "top"):
    ax_ndt.spines[spine].set_visible(False)
ax_ndt.axhline(0.5, color=COL_NDT, lw=2, alpha=0.35, zorder=1)
ndt_rows = []
for i, value in enumerate(ndt_lm):
    rank = sum(1 for j in range(i) if ndt_lm[j] == value)
    tied = ndt_lm.count(value)
    ndt_rows.append(0.5 + 0.13 * (rank - (tied - 1) / 2))
ax_ndt.scatter(
    ndt_lm, ndt_rows, s=70, c=COLOURS, edgecolors="white", linewidths=1.2, zorder=2
)
(point_ndt,) = ax_ndt.plot(
    [], [], "o", ms=13, mfc="none", mec="black", mew=2, zorder=3
)

# ----------------------------------------------------------------------
# Aspect correction for the drift-angle arc. The arena is not equal-aspect
# (seconds against evidence), so a slope in data units draws as some other
# angle on screen; ANGLE_K converts between the two, and ARC_RY is the
# y-radius spanning the same number of pixels as ARC_RX does in x - which is
# what makes the arc a circle on screen, and so makes its parametric angle the
# angle actually seen. Depends on the final axes size, hence computed only
# once the layout is settled.
# ----------------------------------------------------------------------
fig.canvas.draw()
_bbox = ax_arena.get_window_extent()
X_SPAN = T_MAX
Y_SPAN = EVIDENCE_RANGE[1] - EVIDENCE_RANGE[0]
ANGLE_K = (_bbox.height / Y_SPAN) / (_bbox.width / X_SPAN)
ARC_RX = 0.11 * X_SPAN
ARC_RY = ARC_RX / ANGLE_K

ALL_ARTISTS = [
    line_ndt,
    point_start,
    text_ndt,
    line_base,
    line_boundary,
    text_boundary,
    ray_drift,
    arc_drift,
    text_drift,
    lc_paths,
    lc_highlight,
    cursor_arena,
    text_params,
    text_phase,
    trail_plane,
    halo_plane,
    point_plane,
    point_ndt,
    *lines_pdf,
]


def clear_race():
    lc_paths.set_segments([])
    lc_highlight.set_segments([])
    cursor_arena.set_data([], [])


def set_densities(boundary, current, t_upto):
    """Seat every density on `boundary`, and set how visible each one is.

    A stop's curve is traced out left to right as its own race runs - `t_upto`
    is how far the responses have got - and once complete it stays on at
    GHOST_ALPHA for the rest of the animation, so the three can be read against
    each other. Pass `t_upto = None` for a stop whose race has not started. The
    legend entries follow their curve's alpha, so the key only ever names what
    is actually on screen.
    """
    for j, (line, (grid, shape)) in enumerate(zip(lines_pdf, PDF_SHAPES)):
        if j < current:
            drawn, alpha, width = slice(None), GHOST_ALPHA, 1.6
        elif j == current and t_upto is not None:
            drawn = grid <= t_upto
            # a curve begins only once the first response has landed, which is
            # also when its legend entry earns its place
            alpha, width = (1.0, 2.6) if drawn.any() else (0.0, 1.6)
        else:
            drawn, alpha, width = slice(0, 0), 0.0, 1.6
        line.set_data(grid[drawn], boundary + shape[drawn])
        line.set_alpha(alpha)
        line.set_linewidth(width)
        legend.get_texts()[j].set_alpha(alpha)
        legend.legend_handles[j].set_alpha(alpha)
    # and the key's own box fades in with the first curve, rather than sitting
    # there as an empty white rectangle through the opening race
    shown = max(line.get_alpha() for line in lines_pdf)
    legend.get_frame().set_alpha(0.85 * min(1.0, shown / GHOST_ALPHA))


def init():
    clear_race()
    set_densities(LANDMARKS[0][1], 0, None)
    line_ndt.set_data([], [])
    point_start.set_data([], [])
    line_base.set_data([], [])
    line_boundary.set_data([], [])
    ray_drift.set_data([], [])
    trail_plane.set_data([], [])
    point_plane.set_data([], [])
    point_ndt.set_data([], [])
    text_params.set_text("")
    text_phase.set_text("")
    return ALL_ARTISTS


def update(frame):
    state = schedule[frame]
    mu, boundary, ndt = state["params"]
    colour = COLOURS[state["lm"]]

    # ---- the arena's geometry *is* the parameters, drawn -----------------
    line_ndt.set_data([0, ndt], [0, 0])
    point_start.set_data([ndt], [0])
    text_ndt.set_x(ndt * 0.5)
    line_base.set_data([ndt, T_MAX], [0, 0])
    line_boundary.set_data([0, T_MAX], [boundary, boundary])
    text_boundary.set_position((T_MAX * 0.995, boundary - 0.012 * Y_SPAN))

    # the drift ray runs from the start point to where a noiseless
    # accumulator would cross, i.e. at the mean decision time boundary / drift
    t_cross = min(ndt + boundary / mu, T_MAX)
    ray_drift.set_data([ndt, t_cross], [0, (t_cross - ndt) * mu])

    theta = np.arctan(mu * ANGLE_K)  # the angle as it appears on screen
    tip = (ndt + ARC_RX * np.cos(theta), ARC_RY * np.sin(theta))  # sits on the ray
    arc_drift.set_positions((ndt + ARC_RX, 0.0), tip)
    # tan(theta/4) is the arc3 curvature at which a quadratic Bezier matches a
    # circular arc subtending `theta`; the sign puts the bulge on the inside,
    # so the arc sweeps towards the corner rather than away from it
    arc_drift.set_connectionstyle(f"arc3,rad={np.tan(theta / 4):.4f}")
    text_drift.set_position(
        (
            ndt + ARC_RX * 1.15 * np.cos(theta / 2),
            ARC_RY * 1.15 * np.sin(theta / 2),
        )
    )
    text_drift.set_text("drift rate")

    text_params.set_text(
        f"drift rate = {mu:.2f}\n"
        f"boundary   = {boundary:.2f}\n"
        f"ndt        = {ndt:.2f}\n"
        f"E[RT]      = {ndt + boundary / mu:.2f} s"
    )

    # ---- parameter space -------------------------------------------------
    point_plane.set_data([mu], [boundary])
    point_ndt.set_data([ndt], [0.5])
    trail = [s["params"] for s in schedule[max(0, frame - TRAVEL_FRAMES) : frame + 1]]
    trail_plane.set_data([p[0] for p in trail], [p[1] for p in trail])

    # ---- the race --------------------------------------------------------
    if state["phase"] == "travel":
        clear_race()
        # its own race is over, so this stop's density joins the ghosts
        set_densities(boundary, state["lm"] + 1, None)
        halo_plane.set_alpha(0.0)
        text_phase.set_text("tuning the parameters ...")
        return ALL_ARTISTS

    xs, ys, hits, highlight = BATCHES[state["lm"]]
    t_cursor = ndt + state["race"] * (T_MAX - ndt)
    set_densities(boundary, state["lm"], t_cursor)
    k = int(round((t_cursor - ndt) / SIM_DT))  # step the race has reached

    segs, high = [], []
    for i in range(N_PATHS):
        m = min(k + 1, len(xs[i]))
        if m < 2:
            continue
        seg = np.column_stack((xs[i][:m], ys[i][:m]))
        (high if i in highlight else segs).append(seg)
    lc_paths.set_segments(segs)
    lc_paths.set_color(colour)
    lc_highlight.set_segments(high)
    lc_highlight.set_color(colour)

    done = hits[~np.isnan(hits)]
    done = done[done <= t_cursor]

    cursor_arena.set_data([t_cursor, t_cursor], list(EVIDENCE_RANGE))
    halo_plane.set_data([mu], [boundary])
    halo_plane.set_alpha(max(0.0, 1.0 - state["race"]))  # a pulse on firing
    text_phase.set_text(f"firing {N_PATHS} accumulators - {done.size} responded")

    return ALL_ARTISTS


anim = FuncAnimation(
    fig, update, frames=N_FRAMES, init_func=init, blit=False, interval=1000 / FPS
)
anim.save("anim_wald.gif", writer=PillowWriter(fps=FPS))
print(f"Saved figure ({N_FRAMES} frames)")
