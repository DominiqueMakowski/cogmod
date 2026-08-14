"""Generate the cogmod hex sticker SVG.

Geometry is real: the trajectories are simulated diffusion processes and the
curves are Wiener first-passage-time densities (large-time expansion).

    python man/figures/make_logo.py man/figures/logo.svg

logo.png (520 x 600, transparent outside the hex) is that SVG rasterised with
headless Chrome/Edge:

    msedge --headless=new --disable-gpu --hide-scrollbars \
           --default-background-color=00000000 --window-size=520,600 \
           --screenshot=man/figures/logo.png man/figures/logo.svg

--window-size MUST match the width/height on the <svg> element. A smaller
window CROPS the image rather than scaling it, and the result still has
plausible dimensions and transparent corners -- so look at the png, don't
just check its size. To rasterise larger, add --force-device-scale-factor=N
(integer) and keep the window at 520,600.
"""
import math
import os
import random
import re

PINK, BLUE, AMBER = "#F0286B", "#31A8FF", "#FFB020"

W, H = 520.0, 600.0
CX, CY = 260.0, 300.0
HW = 257.2          # half width (center -> side vertex)
HH = 297.0          # half height (center -> top/bottom vertex)

# ---------------------------------------------------------------- plot frame
X0, X1 = 66.0, 454.0        # time axis span (x0 = non-decision time offset)
TOP, BOT = 214.0, 368.0     # upper / lower decision boundary (svg y)
W_START = 0.57              # start point, biased toward the upper boundary
START_Y = BOT - W_START * (BOT - TOP)

TMAX = 2.35                 # seconds of decision time spanned by X0..X1


def t2x(t):
    return X0 + (t / TMAX) * (X1 - X0)


def e2y(e):
    """evidence in [0, 1] (0 = lower bound, 1 = upper bound) -> svg y"""
    return BOT - e * (BOT - TOP)


# ------------------------------------------------------- diffusion sampling
def simulate(v, w, seed, dt=0.002, tmax=TMAX, s=0.47):
    """Random walk between 0 and 1 with drift v, start w. Returns pts, bound."""
    rng = random.Random(seed)
    x, t = w, 0.0
    pts = [(t, x)]
    sq = s * math.sqrt(dt)
    while t < tmax:
        x += v * dt + rng.gauss(0.0, 1.0) * sq
        t += dt
        if x >= 1.0:
            pts.append((t, 1.0))
            return pts, "upper"
        if x <= 0.0:
            pts.append((t, 0.0))
            return pts, "lower"
        pts.append((t, x))
    return pts, None


def thin(pts, every):
    out = [pts[i] for i in range(0, len(pts), every)]
    if out[-1] != pts[-1]:
        out.append(pts[-1])
    return out


def polyline(pts):
    return " ".join("%.2f,%.2f" % (t2x(t), e2y(e)) for t, e in pts)


# ------------------------------------------------- Wiener FPT density (WFPT)
def wfpt_lower(t, v, a, w, K=40):
    """Density of hitting the LOWER boundary at time t (Ratcliff parameterisation)."""
    if t <= 0:
        return 0.0
    tt = t / (a * a)
    s = sum(k * math.exp(-(k * k) * (math.pi ** 2) * tt / 2.0) * math.sin(k * math.pi * w)
            for k in range(1, K + 1))
    p = math.pi * s
    return max(0.0, p * math.exp(-v * a * w - (v * v) * t / 2.0) / (a * a))


def density_curve(v, a, w, upper, n=260):
    """Sampled (t, density) for one boundary."""
    if upper:                       # flip drift & start for the upper bound
        f = lambda t: wfpt_lower(t, -v, a, 1.0 - w)
    else:
        f = lambda t: wfpt_lower(t, v, a, w)
    ts = [TMAX * (i / (n - 1.0)) for i in range(n)]
    return [(t, f(t)) for t in ts]


def trim(curve, level):
    """Drop the far tail so the curve returns to the boundary instead of
    drawing a long flat line across the whole sticker."""
    last = max(i for i, (_, d) in enumerate(curve) if d > level)
    return curve[:last + 1]


def area_path(curve, base_y, height, scale, up=True):
    """Filled area anchored on a boundary line, extending away from it."""
    d = ["M %.2f,%.2f" % (t2x(curve[0][0]), base_y)]
    for t, dens in curve:
        y = base_y - dens * scale * height if up else base_y + dens * scale * height
        d.append("L %.2f,%.2f" % (t2x(t), y))
    d.append("L %.2f,%.2f" % (t2x(curve[-1][0]), base_y))
    d.append("Z")
    return " ".join(d)


# ------------------------------------------- inner-ring clearance accounting
RING = 0.955                # the thin grey hexagon drawn inside the border
CLEARANCE = 13.0            # keep artwork this far off that ring
_checked = []


def ring_clearance(x, y):
    """How far a point sits inside the inner ring (negative = outside)."""
    dx, dy = abs(x - CX), abs(y - CY)
    lo, hi = -80.0, 200.0
    for _ in range(40):                     # bisect the inset that just contains it
        m = (lo + hi) / 2.0
        hw, hh = HW * RING - m, HH * RING - m
        inside = hw > 0 and dx <= hw and dy <= hh / 2.0 * (2.0 - dx / hw)
        lo, hi = (m, hi) if inside else (lo, m)
    return lo


def check(label, pts):
    worst = min(ring_clearance(x, y) for x, y in pts)
    _checked.append((label, worst))
    return worst


# same process, rescaled to the s = 1 convention the analytic density assumes
DRIFT, SEP, S = 0.27, 1.0, 0.47   # ~18% of trials terminate at the lower boundary
w_start = W_START
up_raw = density_curve(DRIFT / S, SEP / S, w_start, upper=True)
lo_raw = density_curve(DRIFT / S, SEP / S, w_start, upper=False)

# The start point sits above the midline, so lower-boundary ("error") responses
# come out slower than upper ones -- that shape difference is real. Their mass
# is not: the lower curve is drawn LO_BOOST times taller so it reads at sticker
# size, since honest scaling leaves it a barely-visible sliver.
PEAK = max(d for _, d in up_raw)
LO_BOOST = 1.8
up_curve = trim(up_raw, 0.018 * PEAK)
lo_curve = trim(lo_raw, 0.018 * PEAK / LO_BOOST)

UP_H = 86.0
up_area = area_path(up_curve, TOP - 5, UP_H, 1.0 / PEAK, up=True)
lo_area = area_path(lo_curve, BOT + 5, UP_H * LO_BOOST, 1.0 / PEAK, up=False)


def mean_t(curve):
    return sum(t * d for t, d in curve) / sum(d for _, d in curve)


print("lower peak %.0f%% of upper (drawn %.0f%%); mean RT %.2f up / %.2f low"
      % (100 * max(d for _, d in lo_raw) / PEAK,
         100 * LO_BOOST * max(d for _, d in lo_raw) / PEAK,
         mean_t(up_raw), mean_t(lo_raw)))


def curve_pts(curve, base_y, h, up):
    return [(t2x(t), base_y - d * h / PEAK if up else base_y + d * h / PEAK)
            for t, d in curve]


up_pts = curve_pts(up_curve, TOP - 5, UP_H, True)
lo_pts = curve_pts(lo_curve, BOT + 5, UP_H * LO_BOOST, False)
check("upper density", up_pts)
check("lower density", lo_pts)
check("boundaries", [(X0 - 24, TOP), (X1 + 16, TOP), (X0 - 24, BOT), (X1 + 16, BOT)])

# ------------------------------------------------------------- trajectories
hero = None
for seed in range(1, 40000):                    # pick a good-looking hero path
    pts, b = simulate(DRIFT, w_start, seed)
    if b != "upper" or not (1.40 < pts[-1][0] < 1.62):
        continue
    if min(e for _, e in pts) > 0.32:           # wants a dip toward the wrong bound
        continue
    if max(e for t, e in pts if t < 0.95) > 0.82:   # but no early flirt with the top
        continue
    hero = pts
    break
hero_pts = thin(hero, 9)
hero_end = (t2x(hero[-1][0]), e2y(1.0))

# Same process as the hero, no drift tricks -- but ask for one error among four,
# so the trace set matches the response proportions the densities describe.
ghosts, seed = [], 8000
for target in ("upper", "upper", "lower", "upper"):
    while seed < 400000:
        seed += 1
        pts, b = simulate(DRIFT, w_start, seed)
        if b != target or not (0.45 < pts[-1][0] < 2.2):
            continue
        if abs(pts[-1][0] - hero[-1][0]) < 0.22:
            continue
        if any(abs(pts[-1][0] - g[0][-1][0]) < 0.34 for g in ghosts):
            continue
        ghosts.append((thin(pts, 13), b))
        break
ghosts.sort(key=lambda g: g[0][-1][0])

# ---------------------------------------------------- head silhouette + cogs
# Profile taken from head_silhouette.svg (SVG Repo) sitting next to this script.
# That file is one path holding the head plus two gear-shaped holes, and two
# more paths for the gear hubs. Split them apart so the head can be a soft
# translucent shape and the cogs can be gold on top of it.
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "head_silhouette.svg")
SRC_BOX = (16.05, 0.0, 139.73, 155.74)      # artwork extent within its viewBox

with open(SRC, encoding="utf-8", errors="replace") as fh:
    _paths = re.findall(r'\sd="([^"]+)"', fh.read())
if len(_paths) != 3:
    raise SystemExit("head_silhouette.svg: expected 3 paths, found %d" % len(_paths))

hub_a, hub_b, _combined = _paths
_subs = [p for p in re.split(r"(?=M)", _combined) if p.strip()]
if len(_subs) != 3:
    raise SystemExit("head_silhouette.svg: expected head + 2 cogs, found %d" % len(_subs))
HEAD, cog_a, cog_b = _subs
# evenodd turns each hub into a hole punched back out of its gold cog
COGS = " ".join([cog_a, hub_b, cog_b, hub_a])

HEAD_H = 88.0                               # drawn height of the silhouette
HEAD_CX, HEAD_BOTTOM = 306.0, 170.0         # centred over the threshold crossing
sx0, sy0, sx1, sy1 = SRC_BOX
HEAD_S = HEAD_H / (sy1 - sy0)
HEAD_X = HEAD_CX - (sx0 + sx1) / 2.0 * HEAD_S
HEAD_Y = HEAD_BOTTOM - sy1 * HEAD_S

# the source icon faces left; the trailing flip turns it to face right
HEAD_T = ("translate(%.2f,%.2f) scale(%.4f) translate(%.2f,0) scale(-1,1)"
          % (HEAD_X, HEAD_Y, HEAD_S, sx0 + sx1))
HEAD_TMPL = ('<g transform="%s">\n'
             '      <path d="%s" fill="#C7D8F7" fill-opacity=".16" stroke="#C7D8F7" '
             'stroke-opacity=".70" stroke-width="%.2f" stroke-linejoin="round"/>\n'
             '      <path d="%s" fill-rule="evenodd" fill="%s" fill-opacity=".92"/>\n'
             '    </g>')
head_svg = HEAD_TMPL % (HEAD_T, HEAD, 2.3 / HEAD_S, COGS, AMBER)

head_corners = [(HEAD_X + x * HEAD_S, HEAD_Y + y * HEAD_S)
                for x in (sx0, sx1) for y in (sy0, sy1)]
check("head", head_corners)

# ----------------------------------------------------------------------- type
WORD_Y, WORD_W, WORD_SIZE = 481.0, 216.0, 65.0   # baseline / locked width / size
SUB_Y, SUB_W, SUB_SIZE = 514.0, 176.0, 12.2


def text_box(w, y, size, descends):
    top, bot = y - 0.72 * size, y + (0.23 * size if descends else 2.0)
    return [(260 - w / 2, top), (260 + w / 2, top),
            (260 - w / 2, bot), (260 + w / 2, bot)]


check("wordmark", text_box(WORD_W, WORD_Y, WORD_SIZE, True))
check("subtitle", text_box(SUB_W, SUB_Y, SUB_SIZE, False))

# ----------------------------------------------------------------- assemble

hexpath = ("M %.1f,%.1f L %.1f,%.1f L %.1f,%.1f L %.1f,%.1f L %.1f,%.1f L %.1f,%.1f Z"
           % (CX, CY - HH, CX + HW, CY - HH / 2, CX + HW, CY + HH / 2,
              CX, CY + HH, CX - HW, CY + HH / 2, CX - HW, CY - HH / 2))

ghost_svg = []
for pts, b in ghosts:
    col = PINK if b == "upper" else BLUE
    ghost_svg.append(
        '    <polyline points="%s" stroke="%s" stroke-width="2" opacity=".34"/>'
        % (polyline(pts), col))
    ex, ey = t2x(pts[-1][0]), e2y(1.0 if b == "upper" else 0.0)
    ghost_svg.append('    <circle cx="%.2f" cy="%.2f" r="3.2" fill="%s" stroke="none" opacity=".6"/>'
                     % (ex, ey, col))

svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="520" height="600" viewBox="0 0 520 600" role="img"
     aria-label="cogmod: an R package for cognitive models">
  <title>cogmod</title>

  <defs>
    <!-- hexagon used both as fill shape and as clip -->
    <path id="hex" d="{hexpath}"/>
    <clipPath id="hexclip"><use xlink:href="#hex"/></clipPath>

    <linearGradient id="bg" x1="0" y1="0" x2="0.35" y2="1">
      <stop offset="0%" stop-color="#1B2A4A"/>
      <stop offset="45%" stop-color="#111B31"/>
      <stop offset="100%" stop-color="#080C16"/>
    </linearGradient>

    <radialGradient id="glowbg" cx="0.5" cy="0.42" r="0.62">
      <stop offset="0%" stop-color="{PINK}" stop-opacity=".20"/>
      <stop offset="60%" stop-color="{PINK}" stop-opacity=".05"/>
      <stop offset="100%" stop-color="{PINK}" stop-opacity="0"/>
    </radialGradient>

    <linearGradient id="rim" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="{AMBER}"/>
      <stop offset="42%" stop-color="{PINK}"/>
      <stop offset="100%" stop-color="{BLUE}"/>
    </linearGradient>

    <linearGradient id="trace" x1="0" y1="1" x2="1" y2="0">
      <stop offset="0%" stop-color="{BLUE}"/>
      <stop offset="45%" stop-color="{PINK}"/>
      <stop offset="100%" stop-color="{AMBER}"/>
    </linearGradient>

    <linearGradient id="updens" x1="0" y1="1" x2="0" y2="0">
      <stop offset="0%" stop-color="{PINK}" stop-opacity=".55"/>
      <stop offset="100%" stop-color="{PINK}" stop-opacity=".05"/>
    </linearGradient>
    <linearGradient id="lodens" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="{BLUE}" stop-opacity=".50"/>
      <stop offset="100%" stop-color="{BLUE}" stop-opacity=".04"/>
    </linearGradient>

    <filter id="soft" x="-30%" y="-30%" width="160%" height="160%">
      <feGaussianBlur stdDeviation="5" result="b"/>
      <feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
    </filter>
    <filter id="tight" x="-40%" y="-40%" width="180%" height="180%">
      <feGaussianBlur stdDeviation="3.2" result="b"/>
      <feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
    </filter>
  </defs>

  <!-- ============================== body ============================== -->
  <use xlink:href="#hex" fill="url(#bg)"/>

  <g clip-path="url(#hexclip)">
    <rect x="0" y="0" width="520" height="600" fill="url(#glowbg)"/>

    <!-- the mind the decision lands in: profile with meshing cogs -->
    {head_svg}

    <!-- decision boundaries -->
    <g stroke="#8FA3C8" stroke-width="2" opacity=".55" stroke-linecap="round">
      <line x1="{X0 - 24:.1f}" y1="{TOP}" x2="{X1 + 16:.1f}" y2="{TOP}" stroke-dasharray="1 9"/>
      <line x1="{X0 - 24:.1f}" y1="{BOT}" x2="{X1 + 16:.1f}" y2="{BOT}" stroke-dasharray="1 9"/>
    </g>

    <!-- first-passage densities -->
    <g filter="url(#soft)">
      <path d="{up_area}" fill="url(#updens)"/>
      <path d="{up_area}" fill="none" stroke="{PINK}" stroke-width="3" stroke-linejoin="round"/>
      <path d="{lo_area}" fill="url(#lodens)"/>
      <path d="{lo_area}" fill="none" stroke="{BLUE}" stroke-width="2.6" stroke-linejoin="round"/>
    </g>

    <!-- sampled trajectories -->
    <g fill="none" stroke-linejoin="round" stroke-linecap="round">
{chr(10).join(ghost_svg)}
    </g>

    <!-- hero trajectory -->
    <polyline points="{polyline(hero_pts)}" fill="none" stroke="url(#trace)"
              stroke-width="4.6" stroke-linejoin="round" stroke-linecap="round"
              filter="url(#tight)"/>
    <circle cx="{hero_end[0]:.2f}" cy="{hero_end[1]:.2f}" r="8" fill="{AMBER}" filter="url(#tight)"/>
    <circle cx="{hero_end[0]:.2f}" cy="{hero_end[1]:.2f}" r="3.2" fill="#FFF7E6"/>

    <!-- starting point z -->
    <circle cx="{t2x(0.0):.2f}" cy="{START_Y}" r="6.5" fill="#0B1120" stroke="#E8EEF9" stroke-width="2.6"/>
    <line x1="{X0 - 24:.1f}" y1="{START_Y}" x2="{t2x(0.0) - 9:.1f}" y2="{START_Y}"
          stroke="#E8EEF9" stroke-width="2" opacity=".45" stroke-linecap="round"/>
  </g>

  <!-- =============================== type =============================== -->
  <g font-family="'Fira Sans','Segoe UI','Helvetica Neue',Helvetica,Arial,sans-serif"
     text-anchor="middle">
    <text x="260" y="{WORD_Y}" font-size="{WORD_SIZE}" font-weight="700" letter-spacing="-1.3"
          fill="#F4F7FC" textLength="{WORD_W}" lengthAdjust="spacingAndGlyphs">cog<tspan fill="{PINK}">mod</tspan></text>
    <text x="260" y="{SUB_Y}" font-size="{SUB_SIZE}" font-weight="600" letter-spacing="2.8"
          fill="#8497BC" textLength="{SUB_W}" lengthAdjust="spacing">COGNITIVE MODELS IN R</text>
  </g>

  <!-- =============================== rim =============================== -->
  <use xlink:href="#hex" fill="none" stroke="#080C16" stroke-width="16"/>
  <use xlink:href="#hex" fill="none" stroke="url(#rim)" stroke-width="9"/>
  <use xlink:href="#hex" fill="none" stroke="#FFFFFF" stroke-width="1.6" opacity=".28"
       transform="translate(260,300) scale(0.955) translate(-260,-300)"/>
</svg>
'''

import sys
out = sys.argv[1] if len(sys.argv) > 1 else "logo.svg"
with open(out, "w", encoding="utf-8") as fh:
    fh.write(svg)
print("wrote", out)
print("inner-ring clearance (need >= %.0f):" % CLEARANCE)
for label, gap in sorted(_checked, key=lambda c: c[1]):
    print("  %-15s %6.1f %s" % (label, gap, "" if gap >= CLEARANCE else "<-- TOO CLOSE"))
print("hero decision time %.2f s at x=%.1f" % (hero[-1][0], hero_end[0]))
print("ghosts:", [(round(p[-1][0], 2), b) for p, b in ghosts])
