"""Overview figure for paper/paper.qmd (Figure 1).

The package logo sits in the middle, with an arrow into each of four regions.
Each region pairs a small drawing of a response format that `cogmod` covers
with the names of the families available for it.

    python man/figures/make_fig_overview.py

Writes paper/figures/fig_overview.png and a copy in man/figures/, which is
where the README picks it up.
"""

import math
import os
import re
import shutil

import matplotlib

matplotlib.use("Agg")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

HERE = os.path.dirname(os.path.abspath(__file__))  # man/figures
ROOT = os.path.dirname(os.path.dirname(HERE))
LOGO = os.path.join(ROOT, "man", "figures", "logo.png")
OUTDIR = os.path.join(ROOT, "paper", "figures")
CLICK_SVG = os.path.join(OUTDIR, "asset_click.svg")

# Logo palette, plus a fourth accent taken from the hex border gradient.
PINK, BLUE, AMBER, VIOLET = "#DB1F58", "#2477C8", "#C57A00", "#7048C4"
INK, MUTED, PALE = "#22252E", "#6B7280", "#C9CDD6"

W_IN, H_IN, DPI = 7.5, 5.65, 400
XMAX = 100.0
YMAX = XMAX * H_IN / W_IN          # 1 data unit = 0.075 in

LOGO_H = 19.0
LOGO_CX, LOGO_CY = 50.0, 41.0
ICON_W, ICON_H = 17.0, 10.0
PAD = 3.0                          # icon-to-label gap


# --------------------------------------------------------------------------
# task icons (centred on cx, cy)
# --------------------------------------------------------------------------
def rounded(ax, cx, cy, w, h, r=0.9, **kw):
    kw.setdefault("linewidth", 1.1)
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2 + r, cy - h / 2 + r), w - 2 * r, h - 2 * r,
        boxstyle=f"round,pad={r},rounding_size={r}", **kw,
    ))


def icon_likert(ax, cx, cy, color):
    """Five-point discrete scale, one option selected."""
    n, step = 5, 3.5
    x0 = cx - step * (n - 1) / 2
    ax.plot([x0 - 1.2, x0 + step * (n - 1) + 1.2], [cy + 0.6, cy + 0.6],
            color=PALE, lw=1.3, zorder=3, solid_capstyle="round")
    for i in range(n):
        x, sel = x0 + i * step, i == 3
        ax.add_patch(Circle((x, cy + 0.6), 1.15, zorder=4, lw=1.2,
                            facecolor=color if sel else "white",
                            edgecolor=color if sel else PALE))
        if sel:
            ax.add_patch(Circle((x, cy + 0.6), 2.1, facecolor="none", lw=0.9,
                                edgecolor=color, alpha=0.45, zorder=4))
    for i, lab in ((0, "1"), (n - 1, "5")):
        ax.text(x0 + i * step, cy - 2.4, lab, fontsize=5.8, color=MUTED,
                ha="center", va="top", zorder=4)


def icon_slider(ax, cx, cy, color):
    """Continuous analog scale with a draggable handle."""
    half, frac = 7.8, 0.70
    xh = cx - half + 2 * half * frac
    ax.plot([cx - half, cx + half], [cy + 0.6, cy + 0.6], color=PALE, lw=2.8,
            zorder=3, solid_capstyle="round")
    ax.plot([cx - half, xh], [cy + 0.6, cy + 0.6], color=color, lw=2.8,
            zorder=4, solid_capstyle="round", alpha=0.8)
    ax.add_patch(Circle((xh, cy + 0.6), 1.6, facecolor="white",
                        edgecolor=color, lw=1.5, zorder=5))
    for x, lab, ha in ((cx - half, "0", "left"), (cx + half, "1", "right")):
        ax.text(x, cy - 2.4, lab, fontsize=5.8, color=MUTED, ha=ha,
                va="top", zorder=4)


def icon_choice(ax, cx, cy, color):
    """Two-alternative speeded decision: a stopwatch and two response keys."""
    sx, sy, r = cx - 6.4, cy + 0.1, 2.2
    ax.add_patch(Circle((sx, sy), r, facecolor="white", edgecolor=INK,
                        lw=1.2, zorder=4))
    ax.plot([sx - 0.85, sx + 0.85], [sy + r + 0.65, sy + r + 0.65], color=INK,
            lw=1.4, zorder=4, solid_capstyle="round")
    ax.plot([sx, sx], [sy + r, sy + r + 0.65], color=INK, lw=1.1, zorder=4)
    ax.plot([sx, sx], [sy, sy + 1.3], color=color, lw=1.1, zorder=5,
            solid_capstyle="round")
    ax.plot([sx, sx + 1.1], [sy, sy + 0.4], color=color, lw=1.1, zorder=5,
            solid_capstyle="round")
    for kx, kwid, lab, on in ((cx + 0.7, 5.6, "TRUE", True),
                              (cx + 7.2, 6.0, "FALSE", False)):
        rounded(ax, kx, cy + 0.1, kwid, 3.6, r=0.85, zorder=4,
                facecolor=color if on else "white",
                edgecolor=color if on else PALE)
        ax.text(kx, cy + 0.1, lab, fontsize=5.2, weight="bold", zorder=5,
                color="white" if on else MUTED, ha="center", va="center")


_SVG_TOKEN = re.compile(r"[MmLlHhVvCcZz]|-?\d*\.?\d+(?:[eE][-+]?\d+)?")
_click_cache = None


def _parse_svg_path(d):
    """Minimal SVG path reader: the M/L/H/V/C/Z subset the asset uses."""
    tok = _SVG_TOKEN.findall(d)
    verts, codes, i, cmd = [], [], 0, None
    cur = start = (0.0, 0.0)
    while i < len(tok):
        if tok[i].isalpha():
            cmd, i = tok[i], i + 1
        rel, c = cmd.islower(), cmd.upper()

        def num():
            nonlocal i
            v = float(tok[i])
            i += 1
            return v

        if c == "M":
            x, y = num(), num()
            if rel:
                x, y = cur[0] + x, cur[1] + y
            cur = start = (x, y)
            verts.append(cur)
            codes.append(Path.MOVETO)
            cmd = "l" if rel else "L"        # further pairs are lineto
        elif c in "LHV":
            if c == "L":
                x, y = num(), num()
                if rel:
                    x, y = cur[0] + x, cur[1] + y
            elif c == "H":
                x, y = num(), cur[1]
                if rel:
                    x += cur[0]
            else:
                x, y = cur[0], num()
                if rel:
                    y += cur[1]
            cur = (x, y)
            verts.append(cur)
            codes.append(Path.LINETO)
        elif c == "C":
            pts = [(num(), num()) for _ in range(3)]
            if rel:
                pts = [(cur[0] + a, cur[1] + b) for a, b in pts]
            verts.extend(pts)
            codes.extend([Path.CURVE4] * 3)
            cur = pts[-1]
        else:                                # Z
            verts.append(start)
            codes.append(Path.CLOSEPOLY)
            cur = start
    return verts, codes


def click_path():
    """The click asset as one matplotlib Path, in SVG units (y pointing down)."""
    global _click_cache
    if _click_cache is None:
        with open(CLICK_SVG) as fh:
            svg = fh.read()
        verts, codes = [], []
        for d in re.findall(r'<path[^>]*\sd="([^"]+)"', svg, re.S):
            v, c = _parse_svg_path(d)
            verts += v
            codes += c
        _click_cache = Path(verts, codes)
    return _click_cache


def draw_click(ax, cx, cy, height, color):
    """Place the click asset, centred on (cx, cy) and flipped upright."""
    path = click_path()
    x0, y0, x1, y1 = path.get_extents().extents
    s = height / (y1 - y0)
    width = (x1 - x0) * s
    t = (Affine2D().translate(-x0, -y0).scale(s, -s)
         .translate(cx - width / 2, cy + height / 2))
    ax.add_patch(PathPatch(t.transform_path(path), facecolor=color,
                           edgecolor="none", zorder=5))
    return width


def icon_rt(ax, cx, cy, color):
    """A stimulus, a single response, and the latency measured between them."""
    sx, hx, base = cx - 5.0, cx + 5.6, cy - 1.2
    rounded(ax, sx, base, 7.2, 5.8, r=0.9, facecolor="white", edgecolor=PALE,
            zorder=3)
    ax.text(sx, base, "event", fontsize=6.2, color=INK, ha="center",
            va="center", zorder=4)
    width = draw_click(ax, hx, base - 0.1, 7.4, color)
    # dimension line spanning stimulus to response: <- RT ->
    y = base + 4.4
    ax.plot([sx, sx], [base + 2.9, y + 0.5], color=PALE, lw=0.7, zorder=3)
    mid = (sx + hx) / 2
    ax.text(mid, y, "RT", fontsize=6.0, color=INK, ha="center", va="center",
            zorder=5)
    for x0, x1 in ((mid - 1.4, sx), (mid + 1.4, hx - width / 2 - 0.3)):
        ax.add_patch(FancyArrowPatch(
            (x0, y), (x1, y), arrowstyle="-|>", mutation_scale=6.0, lw=0.8,
            color=MUTED, shrinkA=0, shrinkB=0, zorder=4,
        ))


# --------------------------------------------------------------------------
# regions
# --------------------------------------------------------------------------
LEFT, RIGHT = (1.0, 43.0), (57.0, 99.0)

REGIONS = [
    dict(
        title="Ordinal ratings", subtitle="Likert-type scale",
        icon=icon_likert, color=AMBER, xbox=LEFT, side="top",
        families=[("Discrete Beta", 10.5, None)],
    ),
    dict(
        title="Analog ratings", subtitle="Slider, visual analog scale",
        icon=icon_slider, color=PINK, xbox=RIGHT, side="top",
        # NOTE: choco() belongs in this inventory, since the package ships it,
        # but it is the subject of a separate paper and is deliberately not
        # discussed anywhere in the manuscript - not in the text, and not in
        # this figure's caption. Leave it listed here and unmentioned there.
        families=[("CHOCO", 12.0, "Choice-confidence"),
                  ("Beta-Gate", 10.5, None)],
    ),
    dict(
        title="Choices and RTs", subtitle="Speeded two-choice decision",
        icon=icon_choice, color=BLUE, xbox=LEFT, side="bottom",
        families=[("DDM", 12.5, "Drift diffusion model"),
                  ("LBA", 12.5, "Linear ballistic accumulator"),
                  ("LNR", 12.5, "Lognormal race"),
                  ("RDM", 12.5, "Racing diffusion model")],
    ),
    dict(
        title="Response times", subtitle="Detection, simple RT",
        icon=icon_rt, color=VIOLET, xbox=RIGHT, side="bottom",
        families=[("ExGaussian", 10.5, None),
                  ("Shifted Wald", 10.5, None),
                  ("Shifted LogNormal", 10.5, None),
                  ("Single-accumulator LBA", 9.0, None)],
    ),
]

GEOM = dict(top=dict(icon_cy=60.5, word_edge=66.5),
            bottom=dict(icon_cy=23.0, word_edge=16.5))


# --------------------------------------------------------------------------
def measure(ax, fig, s, size, **kw):
    t = ax.text(0, 0, s, fontsize=size, **kw)
    fig.canvas.draw()
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    inv = ax.transData.inverted()
    (x0, _), (x1, _) = inv.transform((bb.x0, bb.y0)), inv.transform((bb.x1, bb.y1))
    t.remove()
    return x1 - x0


def pack(items, width, gap):
    """Greedy best-fit packing of items into rows of at most `width`."""
    rows, left = [], list(items)
    while left:
        row, used = [], 0.0
        while True:
            pick = next((it for it in left
                         if used + it["w"] + (gap if row else 0) <= width), None)
            if pick is None:
                break
            used += pick["w"] + (gap if row else 0)
            row.append(pick)
            left.remove(pick)
        if not row:                              # safety net
            it = left.pop(0)
            row, used = [it], it["w"]
        rows.append((row, used))
    return rows


def draw_header(ax, fig, reg):
    """Icon on the side facing the logo, title and subtitle beside it."""
    x0, x1 = reg["xbox"]
    cy = GEOM[reg["side"]]["icon_cy"]
    tw = max(measure(ax, fig, reg["title"], 9.5, weight="bold"),
             measure(ax, fig, reg["subtitle"], 7.0, style="italic"))
    start = (x0 + x1) / 2 - (ICON_W + PAD + tw) / 2
    if x0 < LOGO_CX:                             # left region: icon on the right
        tx, ix = start + tw / 2, start + tw + PAD + ICON_W / 2
    else:
        ix, tx = start + ICON_W / 2, start + ICON_W + PAD + tw / 2
    reg["icon"](ax, ix, cy, reg["color"])
    ax.text(tx, cy + 1.9, reg["title"], fontsize=9.5, weight="bold", color=INK,
            ha="center", va="center", zorder=4)
    ax.text(tx, cy - 1.7, reg["subtitle"], fontsize=7.0, style="italic",
            color=MUTED, ha="center", va="center", zorder=4)
    return ix


def draw_words(ax, fig, reg, gap=2.0):
    """Each family shows its short/main name big and in colour; only families
    with a genuine full-name expansion (e.g. RDM -> racing diffusion model)
    get a small italic subtitle underneath - plain descriptive names (e.g.
    ExGaussian) don't need one."""
    x0, x1 = reg["xbox"]
    color = reg["color"]
    expanded = any(f[2] is not None for f in reg["families"])
    items = []
    for name, size, exp in reg["families"]:
        w = measure(ax, fig, name, size, weight="bold")
        if exp:
            w = max(w, measure(ax, fig, exp, 6.5, style="italic"))
        items.append(dict(text=name, size=size, exp=exp, w=w))

    if reg.get("sort"):                          # widest first packs tighter
        items.sort(key=lambda it: -it["w"])
    rows = pack(items, x1 - x0, gap)
    row_h = 5.8 if expanded else 4.0
    edge = GEOM[reg["side"]]["word_edge"]
    top = edge + row_h * len(rows) if reg["side"] == "top" else edge

    for i, (row, used) in enumerate(rows):
        y = top - (i + 0.5) * row_h
        x = (x0 + x1) / 2 - used / 2
        for j, it in enumerate(row):
            if j:
                x += gap
            xc = x + it["w"] / 2
            ax.text(xc, y + (0.9 if it["exp"] else 0), it["text"],
                    fontsize=it["size"], weight="bold", color=color,
                    ha="center", va="center", zorder=4,
                    alpha=0.92)
            if it["exp"]:
                ax.text(xc, y - 1.6, it["exp"], fontsize=6.5, style="italic",
                        color=MUTED, ha="center", va="center", zorder=4)
            x += it["w"]


# --------------------------------------------------------------------------
def main():
    fig = plt.figure(figsize=(W_IN, H_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, YMAX)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    targets = []
    for reg in REGIONS:
        targets.append((reg, draw_header(ax, fig, reg)))
        draw_words(ax, fig, reg)

    # arrows are drawn first so that they emerge from behind the hex
    for reg, ix in targets:
        ty = GEOM[reg["side"]]["icon_cy"]
        dx, dy = ix - LOGO_CX, ty - LOGO_CY
        n = math.hypot(dx, dy)
        ux, uy = dx / n, dy / n
        ax.add_patch(FancyArrowPatch(
            (LOGO_CX + ux * 7.5, LOGO_CY + uy * 7.5),
            (ix - ux * 9.6, ty - uy * 9.6),
            arrowstyle="-|>", mutation_scale=11,
            connectionstyle="arc3,rad=0.2", lw=1.4, color=reg["color"],
            alpha=0.75, shrinkA=0, shrinkB=0, zorder=2,
        ))

    img = mpimg.imread(LOGO)
    lw = LOGO_H * img.shape[1] / img.shape[0]
    ax.imshow(img, extent=[LOGO_CX - lw / 2, LOGO_CX + lw / 2,
                           LOGO_CY - LOGO_H / 2, LOGO_CY + LOGO_H / 2],
              zorder=6, interpolation="antialiased")

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "fig_overview.png")
    fig.savefig(out, dpi=DPI, facecolor="white")
    print("wrote", out)

    # The README shows the same figure; man/figures/ is what GitHub and
    # pkgdown serve it from.
    readme_out = os.path.join(ROOT, "man", "figures", "fig_overview.png")
    shutil.copyfile(out, readme_out)
    print("wrote", readme_out)


if __name__ == "__main__":
    main()
