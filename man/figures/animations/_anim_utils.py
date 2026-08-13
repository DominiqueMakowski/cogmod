"""Shared building blocks for the `animations_priors*.py` scripts.

These are generic helpers for animating a point moving through parameter
space while showing, for each parameter, how a Normal prior placed on a
"linked" (unconstrained) scale looks once pushed through a link function
onto the "unlinked"/natural scale.
"""

import numpy as np
from matplotlib.patches import ConnectionPatch
from scipy.stats import norm, t as student_t

# ----------------------------------------------------------------------
# Path building
# ----------------------------------------------------------------------


def build_path(keyframes, n_frames):
    """Piecewise-linear path through `keyframes` (a sequence of equal-length
    tuples, e.g. (mu, sigma) or (mu, sigma, tau)), eased in/out per segment
    (smoothstep) so the point doesn't snap or jerk at each corner. Returns
    one 1D array per coordinate, e.g. `mu_path, sigma_path = build_path(...)`.
    """
    keyframes = np.asarray(keyframes, dtype=float)
    n_segments = len(keyframes) - 1
    frames_per_segment = n_frames // n_segments
    n_dims = keyframes.shape[1]
    paths = [[] for _ in range(n_dims)]
    for p0, p1 in zip(keyframes[:-1], keyframes[1:]):
        t = np.linspace(0, 1, frames_per_segment, endpoint=False)
        ease = t * t * (3 - 2 * t)
        for d in range(n_dims):
            paths[d].append(p0[d] + (p1[d] - p0[d]) * ease)
    return tuple(np.concatenate(p) for p in paths)


# ----------------------------------------------------------------------
# Axis-range helpers
# ----------------------------------------------------------------------


def padded_range(lo, hi, pad_frac=0.2, min_pad=0.3):
    """(lo, hi) expanded by a fraction of their span, with a floor so a
    zero-span input still gets visible padding."""
    pad = max((hi - lo) * pad_frac, min_pad)
    return lo - pad, hi + pad


def prior_support(dist, extra_vals, q=0.001, pad_frac=0.05):
    """Range that shows a prior's bulk (via its quantiles) while still
    fitting any extra values (e.g. landmarks) that must stay visible."""
    lo, hi = dist.ppf(q), dist.ppf(1 - q)
    lo, hi = min(lo, min(extra_vals)), max(hi, max(extra_vals))
    return padded_range(lo, hi, pad_frac=pad_frac, min_pad=0)


# ----------------------------------------------------------------------
# Softplus link (for parameters that must stay positive, e.g. RT params)
# ----------------------------------------------------------------------


def softplus(z):
    """log(1 + exp(z)), computed in a numerically stable way."""
    return np.logaddexp(0.0, np.asarray(z, dtype=float))


def softplus_inv(y):
    """Inverse of `softplus`, defined for y > 0."""
    y = np.asarray(y, dtype=float)
    return y + np.log(-np.expm1(-y))


class SoftplusNormal:
    """Distribution induced by pushing Z ~ N(loc, scale) through the
    softplus link, i.e. the law of Y = softplus(Z). Mimics the small
    subset of the scipy.stats distribution interface (`.pdf`, `.ppf`)
    used elsewhere in these scripts.
    """

    def __init__(self, loc=0.0, scale=1.0):
        self.base = norm(loc=loc, scale=scale)

    def pdf(self, y):
        y = np.asarray(y, dtype=float)
        z = softplus_inv(y)
        jacobian = 1.0 / -np.expm1(-y)  # dz/dy = 1 / (1 - exp(-y))
        return self.base.pdf(z) * jacobian

    def ppf(self, q):
        # softplus is monotonically increasing, so quantiles transfer directly.
        return softplus(self.base.ppf(q))


class HalfStudentT:
    """Student-t truncated to values above `lower` (a "half-t"). This is the
    prior `brms` places by default on scale parameters such as `sigma`,
    which are declared directly on the natural scale (identity link) with a
    lower bound of 0. Mimics the small subset of the scipy.stats
    distribution interface (`.pdf`, `.ppf`) used elsewhere in these scripts.
    """

    def __init__(self, df=3, loc=0.0, scale=2.5, lower=0.0):
        self.base = student_t(df=df, loc=loc, scale=scale)
        self.lower = lower
        self._lower_cdf = self.base.cdf(lower)
        self._mass = 1.0 - self._lower_cdf  # renormalising constant

    def pdf(self, y):
        y = np.asarray(y, dtype=float)
        return np.where(y < self.lower, 0.0, self.base.pdf(y) / self._mass)

    def ppf(self, q):
        q = np.asarray(q, dtype=float)
        return self.base.ppf(self._lower_cdf + q * self._mass)


# ----------------------------------------------------------------------
# Prior panels (linked space <-> unlinked/natural space)
# ----------------------------------------------------------------------


def setup_prior_panel(ax, title, grid, dist, color, xlabel, landmark_vals, dist_label):
    """Draw a prior curve + faint landmark markers, return the animated
    marker/vline."""
    pdf = dist.pdf(grid)
    ax.set_title(title, fontsize=10)
    ax.set_yticks([])
    ax.set_xlabel(xlabel)
    ax.set_xlim(grid.min(), grid.max())
    ax.set_ylim(0, pdf.max() * 1.15)
    ax.plot(grid, pdf, color=color, lw=2)
    ax.fill_between(grid, pdf, color=color, alpha=0.15)
    ax.text(
        0.97,
        0.93,
        dist_label,
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        ha="right",
        style="italic",
        color=color,
    )
    for v in landmark_vals:
        ax.plot([v], [dist.pdf(v)], "o", color="black", ms=5, alpha=0.3, zorder=3)
        ax.axvline(v, color="black", lw=0.8, ls=":", alpha=0.2, zorder=2)
    (marker,) = ax.plot([], [], "o", color="black", ms=7, zorder=5)
    vline = ax.axvline(0, color="black", lw=1, ls="--", alpha=0.6)
    return marker, vline


def add_link_panel(
    fig, ax_left, ax_right, label, link_fn, linked_grid, natural_range, color
):
    """Arrow + label between `ax_left` (linked space) and `ax_right`
    (natural/unlinked space), plus a small inset curve underneath showing
    the link function itself (linked value -> natural value). Returns the
    animated marker that tracks the current parameter value on that curve.
    """
    con = ConnectionPatch(
        xyA=(1.02, 0.5),
        coordsA=ax_left.transAxes,
        xyB=(-0.02, 0.5),
        coordsB=ax_right.transAxes,
        arrowstyle="-|>",
        mutation_scale=18,
        lw=1.6,
        color="#333333",
    )
    fig.add_artist(con)
    pos_left, pos_right = ax_left.get_position(), ax_right.get_position()
    xmid = (pos_left.x1 + pos_right.x0) / 2
    ymid = (pos_left.y0 + pos_left.y1) / 2
    fig.text(
        xmid, ymid + 0.018, label, ha="center", va="bottom", fontsize=9, style="italic"
    )

    # small inset axes sitting in the same horizontal gap, below the arrow
    gap = pos_right.x0 - pos_left.x1
    mini_w = max(gap * 0.62, 0.05)
    mini_h = (pos_left.y1 - pos_left.y0) * 0.40
    mini_ax = fig.add_axes([xmid - mini_w / 2, pos_left.y0, mini_w, mini_h])
    mini_ax.plot(linked_grid, link_fn(linked_grid), color=color, lw=1.4)
    mini_ax.set_xlim(linked_grid.min(), linked_grid.max())
    mini_ax.set_ylim(*natural_range)
    mini_ax.tick_params(labelsize=5, length=2)
    (marker,) = mini_ax.plot([], [], "o", color="black", ms=5, zorder=5)
    return marker


def build_softplus_param_panels(
    fig,
    ax_linked,
    ax_unlinked,
    symbol,
    loc,
    scale,
    landmarks,
    color,
    link_label,
    grid_n=400,
):
    """Set up a linked/unlinked prior-panel pair for a parameter that uses a
    softplus link: a Normal(loc, scale) prior lives on the linked
    (unconstrained) scale, and its induced natural-scale distribution is a
    `SoftplusNormal`. Wires up the marker/vline pairs on both panels plus
    the connecting link inset, and returns everything `update()` needs.
    """
    prior_link = norm(loc=loc, scale=scale)
    prior_natural = SoftplusNormal(loc=loc, scale=scale)

    link_landmarks = [float(softplus_inv(v)) for v in landmarks]
    linked_grid = np.linspace(*prior_support(prior_link, link_landmarks), grid_n)

    natural_lo, natural_hi = prior_support(prior_natural, landmarks, q=0.05)
    natural_grid = np.linspace(max(natural_lo, 1e-3), natural_hi, grid_n)

    marker_linked, vline_linked = setup_prior_panel(
        ax_linked,
        f"Prior for {symbol} \u2014 linked space",
        linked_grid,
        prior_link,
        color,
        f"{symbol}_link",
        link_landmarks,
        f"N({loc:g}, {scale:g})",
    )
    marker_unlinked, vline_unlinked = setup_prior_panel(
        ax_unlinked,
        f"Prior for {symbol} \u2014 unlinked space",
        natural_grid,
        prior_natural,
        color,
        symbol,
        landmarks,
        "Softplus-Normal",
    )
    marker_link = add_link_panel(
        fig,
        ax_linked,
        ax_unlinked,
        link_label,
        softplus,
        linked_grid,
        (natural_grid.min(), natural_grid.max()),
        color,
    )

    return {
        "prior_link": prior_link,
        "prior_natural": prior_natural,
        "link_inv": softplus_inv,
        "marker_linked": marker_linked,
        "vline_linked": vline_linked,
        "marker_unlinked": marker_unlinked,
        "vline_unlinked": vline_unlinked,
        "marker_link": marker_link,
    }


def build_identity_param_panels(
    fig,
    ax_linked,
    ax_unlinked,
    symbol,
    dist,
    landmarks,
    color,
    link_label,
    dist_label,
    grid_n=400,
    q=0.05,
):
    """Set up a linked/unlinked prior-panel pair for a parameter that uses an
    identity link, i.e. one that is sampled directly on its natural scale.
    Both panels therefore show the *same* distribution, and the link inset is
    the 1:1 line. Used for parameters whose positivity is enforced by a lower
    bound on the parameter itself (plus a truncated prior) rather than by a
    link function.
    """
    lo, hi = prior_support(dist, landmarks, q=q)
    lo = max(lo, getattr(dist, "lower", -np.inf))
    grid = np.linspace(lo, hi, grid_n)

    marker_linked, vline_linked = setup_prior_panel(
        ax_linked,
        f"Prior for {symbol} — linked space",
        grid,
        dist,
        color,
        symbol,
        landmarks,
        dist_label,
    )
    marker_unlinked, vline_unlinked = setup_prior_panel(
        ax_unlinked,
        f"Prior for {symbol} — unlinked space",
        grid,
        dist,
        color,
        symbol,
        landmarks,
        dist_label,
    )
    marker_link = add_link_panel(
        fig,
        ax_linked,
        ax_unlinked,
        link_label,
        lambda z: z,
        grid,
        (grid.min(), grid.max()),
        color,
    )

    return {
        "prior_link": dist,
        "prior_natural": dist,
        "link_inv": lambda y: y,
        "marker_linked": marker_linked,
        "vline_linked": vline_linked,
        "marker_unlinked": marker_unlinked,
        "vline_unlinked": vline_unlinked,
        "marker_link": marker_link,
    }
