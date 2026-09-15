"""Record the interactive figure of the RT and decision-making articles as a GIF.

The other scripts in this directory draw their frames with matplotlib. This one
cannot: what it is recording is a live web figure - `vignettes/articles/_widget.qmd`,
an ojs cell that draws itself in the browser and answers to tabs, sliders and
drag handles on the figure. So the frames come out of a real browser, and the
script's job is to render a page holding the widget, serve it, drive it through
a fixed list of actions, and encode what it saw.

    python anim_widget.py                        # man/figures/anim_widget.gif
    python anim_widget.py --fps 8 --width 700    # smaller
    python anim_widget.py --mp4                  # an mp4 beside the gif

Needs `quarto` on PATH, and the python side needs `playwright` (with chromium
installed: `python -m playwright install chromium`), `pillow` and
`imageio-ffmpeg`, which is where the ffmpeg binary comes from.

## What is worth knowing before editing this

The widget animates on its own - a fresh set of trials every 7 s, taking 6.2 s
to sweep in - and three of its habits decide how the capture has to work:

- the volley is off entirely under `prefers-reduced-motion`, so the browser
  context is opened at `no-preference` rather than left to the machine's
  setting;
- it is paused by an IntersectionObserver whenever the figure is off screen, so
  the viewport is made tall enough to hold the whole widget and nothing ever
  scrolls;
- the traces are Web Animations driven by the compositor while the volley is a
  plain `setInterval`, so there is no one clock to step and no way to render
  this frame by frame. The capture is therefore in real time: the choreography
  advances by however long the last screenshot actually took, and each frame is
  written down with the moment it was taken, so a slow screenshot costs a frame
  rather than bending the timeline. ffmpeg resamples to a constant rate at the
  end from those timestamps.

Nothing is aimed at by pixel. Tabs, sliders and segments are found by selector,
and the figure's own drag handles - the drift arrow, the boundary line - are
found by *asking the widget*: it sets `plot.style.cursor` when the pointer comes
within reach of one, so sweeping synthetic pointermoves over the figure and
watching the cursor change maps the handles wherever the figure has since moved
them. That is the part that makes this script survive a redesign, which is the
whole reason for having a script rather than a screen recording.

## Size

GIF is a poor container for thirty wiggling traces, and the knobs that matter
are `--fps`, `--width` and `--colors`, in that order. The script prints what it
wrote. `--mp4` is there for when the page embedding this can take a video.
"""

from __future__ import annotations

import argparse
import functools
import http.server
import shutil
import socketserver
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import imageio_ffmpeg
from PIL import Image, ImageDraw
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ARTICLES = ROOT / "vignettes" / "articles"

# The wrapper page. The widget is an include partial and cannot be rendered on
# its own, and the two articles that do include it carry a few hundred lines of
# prose and fitted models with them - so the capture gets a page of its own,
# next to `_widget.qmd` because quarto resolves an include relative to the file
# that asks for it. Underscored, and both it and everything quarto writes beside
# it are already in .gitignore. `cogmod-start` is the assumption set the figure
# opens on, exactly as in an article's front matter.
PAGE_QMD = """---
title: "widget"
cogmod-start: @START@
format:
  html:
    toc: false
---

{{< include _widget.qmd >}}
"""
PAGE_STEM = "_anim_widget_page"

# Tall enough that the whole widget is on screen at once and nothing scrolls -
# which is what keeps the volley running, the figure being watched by an
# IntersectionObserver. The width is a desktop article's; the widget lays itself
# out in one row above ~900px and stacks below it.
VIEWPORT = {"width": 1100, "height": 1400}


# ----------------------------------------------------------------------
# The choreography
# ----------------------------------------------------------------------
# This is the part meant to be edited. A step is a target to move to, a click,
# a drag, or a pause, and each carries how long it takes in seconds. Targets are
# resolved at the moment the step begins rather than up front, because half of
# them do not exist yet when it does: switching to a ballistic model brings a
# `sigmabias` slider into the column, and the figure is rebuilt - handles and
# all - every time any parameter moves.


@dataclass(frozen=True)
class Sel:
    """The centre of the first element matching a CSS selector."""

    css: str


@dataclass(frozen=True)
class Seg:
    """One segment of an assumption toggle, by the text on it.

    `group` is the toggle's own class - `.cogmod-assumption-ratedist` and the
    rest - and `text` is the label, matched whole so that "Normal" does not
    find "LogNormal".
    """

    group: str
    text: str


@dataclass(frozen=True)
class Slider:
    """A point on a parameter's slider.

    `par` is the name in the class the markup puts on the row, so "boundary" is
    `.cogmod-par-boundary`. With `value` given the point is where the thumb
    would sit at that value, which is what a drag aims at; without one it is
    where the thumb is now, which is what a drag has to start from.
    """

    par: str
    value: float | None = None


@dataclass(frozen=True)
class Handle:
    """A drag handle on the figure, found by the cursor the widget shows over it.

    `cursor` is one of the shapes `dragHandles` sets - "move" for a drift arrow,
    "ns-resize" for the boundary, the start point and the start-point range,
    "ew-resize" for ndt. Where several answer to the same shape, `pick` chooses:
    "middle" takes the centre of the cluster, "top" and "bottom" its extremes.
    """

    cursor: str
    pick: str = "middle"


@dataclass(frozen=True)
class Offset:
    """A point so many pixels from wherever the pointer is now."""

    dx: float
    dy: float


@dataclass(frozen=True)
class Move:
    target: object
    secs: float = 0.8


@dataclass(frozen=True)
class Drag:
    target: object
    secs: float = 1.0


@dataclass(frozen=True)
class Click:
    secs: float = 0.3


@dataclass(frozen=True)
class Hold:
    secs: float = 1.0


# A pass through the widget: a parameter moved from the column, a parameter
# moved from the figure itself, a model changed from the tabs, an assumption
# changed from the bar, and a parameter that only exists once that assumption
# has been changed. Each of the five is a different way in, which is the thing
# about the widget worth a moving picture.
#
# The pauses are not padding. A volley takes 6.2 s to sweep in and comes round
# every 7 s, so a change with nothing after it is a change the reader never sees
# land - and the beat before each click is what shows the hover state the thing
# under the pointer is answering with.
STORY = [
    Hold(1.6),                                            # a volley fills
    Move(Slider("boundary"), 0.9),
    Hold(0.4),                                            # the ring comes up
    Drag(Slider("boundary", 1.55), 1.0),
    Hold(0.6),
    Drag(Slider("boundary", 0.85), 0.8),
    Hold(0.9),
    Move(Handle("move"), 0.9),                            # the drift arrow
    Hold(0.3),
    Drag(Offset(30, -62), 1.0),
    Hold(0.5),
    Drag(Offset(-22, 40), 0.8),
    Hold(1.0),
    Move(Sel('.cogmod-tab[data-model="lba2"]'), 1.0),     # over the other tabs
    Hold(0.5),
    Click(0.3),
    Hold(2.4),
    Move(Seg(".cogmod-assumption-ratedist", "LogNormal"), 1.0),
    Hold(0.5),
    Click(0.3),
    Hold(2.4),
    Move(Slider("sigmabias"), 0.8),
    Hold(0.3),
    Drag(Slider("sigmabias", 0.45), 1.1),
    Hold(2.0),
    Move(Sel('.cogmod-tab[data-model="ddm"]'), 1.0),      # back where it began
    Hold(0.4),
    Click(0.3),
    Hold(1.8),
]


# ----------------------------------------------------------------------
# The page
# ----------------------------------------------------------------------


def build_page(start: str, reuse: bool = False) -> Path:
    """Render the wrapper page, and hand back the html quarto wrote."""
    src = ARTICLES / f"{PAGE_STEM}.qmd"
    html = ARTICLES / f"{PAGE_STEM}.html"
    if reuse and html.exists():
        print(f"reusing {html.name}")
        return html
    quarto = shutil.which("quarto")
    if quarto is None:
        sys.exit("quarto is not on PATH")
    src.write_text(PAGE_QMD.replace("@START@", start), encoding="utf-8")
    print(f"rendering {src.name} ...")
    subprocess.run([quarto, "render", src.name, "--to", "html"],
                   cwd=ARTICLES, check=True,
                   stdout=subprocess.DEVNULL)
    return html


def drop_page() -> None:
    """Take the wrapper page and quarto's output back out of the tree."""
    for path in (ARTICLES / f"{PAGE_STEM}.qmd",
                 ARTICLES / f"{PAGE_STEM}.html"):
        path.unlink(missing_ok=True)
    shutil.rmtree(ARTICLES / f"{PAGE_STEM}_files", ignore_errors=True)


class _Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):        # one line per request, otherwise
        pass


def serve(directory: Path):
    """Serve `directory` on a free local port.

    The page has to come over http rather than off the disk: ojs loads itself as
    ES modules from the `_files` directory quarto writes beside the html, and a
    file:// origin is not allowed to.
    """
    handler = functools.partial(_Quiet, directory=str(directory))
    server = socketserver.ThreadingTCPServer(("127.0.0.1", 0), handler)
    server.daemon_threads = True
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, server.server_address[1]


# ----------------------------------------------------------------------
# Asking the page where things are
# ----------------------------------------------------------------------

# The union of the four blocks the widget is made of. There is no single element
# wrapping them - the title, the tab strip, the row of sliders and figure, and
# the assumptions bar are siblings in the article's flow - so the crop is
# whatever box holds all four.
_EXTENT_JS = """() => {
  const parts = ['.cogmod-title', '.cogmod-tabs', '.cogmod-layout',
                 '.cogmod-assumptions'];
  let box = null;
  for (const sel of parts) {
    const el = document.querySelector(sel);
    if (!el) continue;
    const r = el.getBoundingClientRect();
    const b = {x0: r.left + scrollX, y0: r.top + scrollY,
               x1: r.right + scrollX, y1: r.bottom + scrollY};
    box = box === null ? b : {x0: Math.min(box.x0, b.x0),
                              y0: Math.min(box.y0, b.y0),
                              x1: Math.max(box.x1, b.x1),
                              y1: Math.max(box.y1, b.y1)};
  }
  return box;
}"""

_BOX_JS = """(css) => {
  const el = document.querySelector(css);
  if (!el) return null;
  const r = el.getBoundingClientRect();
  return {x: r.left + r.width / 2, y: r.top + r.height / 2};
}"""

# A segment is a label with a hidden radio in it, so the text to match is the
# label's own, trimmed and whole.
_SEG_JS = """([group, text]) => {
  const root = document.querySelector(group);
  if (!root) return null;
  for (const label of root.querySelectorAll('label')) {
    if (label.textContent.trim() === text) {
      const r = label.getBoundingClientRect();
      return {x: r.left + r.width / 2, y: r.top + r.height / 2};
    }
  }
  return null;
}"""

# Where the thumb sits, or would sit at `value`. The thumb is 13px wide and its
# centre travels the track less one thumb, which is the same arithmetic the
# widget's own value bubble does to follow it.
_SLIDER_JS = """([par, value]) => {
  const input = document.querySelector('.cogmod-par-' + par +
                                       ' input[type="range"]');
  if (!input || input.offsetParent === null) return null;
  const r = input.getBoundingClientRect();
  const lo = Number(input.min), hi = Number(input.max);
  const v = value === null ? Number(input.value)
                           : Math.max(lo, Math.min(hi, value));
  const THUMB = 13;
  return {x: r.left + THUMB / 2 + (v - lo) / (hi - lo) * (r.width - THUMB),
          y: r.top + r.height / 2};
}"""

# Sweep the figure with synthetic pointermoves and collect every point the
# widget answers with a cursor. `dragHandles` listens on the plot and reads
# clientX/clientY off the event, so a dispatched one tells the truth - and this
# runs entirely inside the page, so a fine grid costs milliseconds rather than a
# round trip each.
_HANDLES_JS = """(step) => {
  let plot = null, area = 0;
  for (const svg of document.querySelectorAll('.cogmod-figure svg')) {
    const r = svg.getBoundingClientRect();
    if (r.width * r.height > area) { area = r.width * r.height; plot = svg; }
  }
  if (plot === null) return [];
  const r = plot.getBoundingClientRect();
  const was = plot.style.cursor;
  const found = [];
  for (let y = r.top + 1; y < r.bottom; y += step) {
    for (let x = r.left + 1; x < r.right; x += step) {
      plot.style.cursor = '';
      plot.dispatchEvent(new PointerEvent('pointermove',
        {clientX: x, clientY: y, bubbles: true}));
      if (plot.style.cursor) found.push({x, y, cursor: plot.style.cursor});
    }
  }
  plot.style.cursor = was;
  return found;
}"""


def resolve(page, target, at):
    """Turn a target into a point, `at` being where the pointer is now."""
    if isinstance(target, Offset):
        return (at[0] + target.dx, at[1] + target.dy)
    if isinstance(target, Sel):
        got = page.evaluate(_BOX_JS, target.css)
        what = target.css
    elif isinstance(target, Seg):
        got = page.evaluate(_SEG_JS, [target.group, target.text])
        what = f"{target.group} / {target.text}"
    elif isinstance(target, Slider):
        got = page.evaluate(_SLIDER_JS, [target.par, target.value])
        what = f"slider {target.par}"
    elif isinstance(target, Handle):
        got = _handle(page, target)
        what = f"handle {target.cursor}"
    else:
        raise TypeError(f"not a target: {target!r}")
    if got is None:
        raise RuntimeError(f"nothing to aim at: {what}")
    return (got["x"], got["y"])


def _handle(page, target: Handle, step: float = 4.0):
    """One of the figure's drag handles, by the cursor it shows."""
    hits = [h for h in page.evaluate(_HANDLES_JS, step)
            if h["cursor"] == target.cursor]
    if not hits:
        return None
    if target.pick == "top":
        return min(hits, key=lambda h: h["y"])
    if target.pick == "bottom":
        return max(hits, key=lambda h: h["y"])
    # The middle of the cluster rather than of its box: an arrow runs diagonally
    # and the centre of what it spans need not be on it.
    hits.sort(key=lambda h: (h["x"], h["y"]))
    return hits[len(hits) // 2]


# ----------------------------------------------------------------------
# The pointer
# ----------------------------------------------------------------------


def cursor_sprite(scale: float, pressed: bool = False):
    """An arrow, and how far into it the tip sits, to paste at the pointer.

    A browser screenshot has no pointer in it, and a recording of somebody using
    a thing that does not show them using it is a recording of a figure changing
    its mind on its own. Drawn four times over and shrunk, there being no
    antialiasing in a polygon fill.
    """
    over = 4
    pad = 7 if pressed else 0
    w, h = 16 + 2 * pad, 24 + 2 * pad
    im = Image.new("RGBA", (int(w * over), int(h * over)), (0, 0, 0, 0))
    draw = ImageDraw.Draw(im)
    if pressed:
        # A ring where the press landed, which is at the tip and not at the
        # middle of the arrow. Thin, unfilled and held for a moment only - see
        # PRESS_TICK. A ring left up for a whole drag sits on the slider's own
        # thumb and reads as a second one.
        r = 5.5 * over
        c = (pad * over, pad * over)
        draw.ellipse([c[0] - r, c[1] - r, c[0] + r, c[1] + r],
                     outline=(55, 71, 79, 150), width=max(1, over // 2))
    arrow = [(0, 0), (0, 17.5), (4.6, 13.4), (7.3, 19.8), (10.2, 18.5),
             (7.5, 12.4), (12.8, 12.4)]
    poly = [((x + pad) * over, (y + pad) * over) for x, y in arrow]
    draw.polygon(poly, fill=(255, 255, 255, 255))
    draw.line(poly + [poly[0]], fill=(28, 34, 38, 255),
              width=int(1.3 * over), joint="curve")
    size = (max(1, int(w * scale)), max(1, int(h * scale)))
    return im.resize(size, Image.Resampling.LANCZOS), int(pad * scale)


# ----------------------------------------------------------------------
# Capture
# ----------------------------------------------------------------------


# How long the ring stays up after a press. Long enough to be seen at any rate
# the capture manages, short enough that a drag is the pointer and the thing it
# is moving and nothing else.
PRESS_TICK = 0.22


def smoothstep(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3 - 2 * t)


@dataclass
class Live:
    """The state a step needs while it is running."""

    index: int = 0
    began: float = 0.0
    origin: tuple = (0.0, 0.0)
    goal: tuple = (0.0, 0.0)
    down: bool = False
    released: bool = False
    frames: list = field(default_factory=list)


def capture(page, clip, scale, fps, out_dir: Path, story=STORY):
    """Run the story, writing one png per frame and noting when it was taken.

    Nothing is drawn into a frame here. Compositing the pointer means decoding
    and re-encoding a png, and a loop that stops to do that between screenshots
    is a loop taking fewer of them - the timeline would still be right, the
    motion would just be coarser. The frames are stamped afterwards, from the
    pointer positions noted alongside them.
    """
    live = Live()
    pressed_at = -99.0
    at = resolve(page, Sel(".cogmod-title"), (0, 0))   # somewhere harmless
    page.mouse.move(*at)
    total = sum(step.secs for step in story)
    period = 1.0 / fps
    t0 = time.perf_counter()
    entered = False

    def begin(step):
        """Do whatever the step does on its way in."""
        nonlocal pressed_at, entered
        live.origin = at
        live.released = False
        if isinstance(step, (Move, Drag)):
            live.goal = resolve(page, step.target, at)
        if isinstance(step, (Drag, Click)):
            page.mouse.move(*at)
            page.mouse.down()
            live.down = True
            pressed_at = time.perf_counter() - t0
        entered = True

    def end(step):
        """And on its way out, wherever the clock got to."""
        nonlocal at, entered
        if isinstance(step, (Move, Drag)):
            at = live.goal
            page.mouse.move(*at)
        if live.down:
            page.mouse.up()
            live.down = False
        entered = False

    while True:
        now = time.perf_counter() - t0
        if now >= total or live.index >= len(story):
            break

        # Which step we are in. A step is entered when the clock reaches it
        # rather than when the one before it was drawn, so a slow frame costs a
        # frame and never shifts the timeline - and a step the clock stepped
        # clean over, a screenshot having taken longer than the step lasts, is
        # still opened and closed rather than skipped. A click that quietly did
        # not happen would leave everything after it in the wrong model.
        while (live.index < len(story)
               and now >= live.began + story[live.index].secs):
            step = story[live.index]
            if not entered:
                begin(step)
            end(step)
            live.began += step.secs
            live.index += 1
        if live.index >= len(story):
            break

        step = story[live.index]
        if not entered:
            begin(step)
        into = (now - live.began) / step.secs if step.secs else 1.0

        if isinstance(step, (Move, Drag)):
            ease = smoothstep(into)
            at = (live.origin[0] + (live.goal[0] - live.origin[0]) * ease,
                  live.origin[1] + (live.goal[1] - live.origin[1]) * ease)
            page.mouse.move(*at)
        if isinstance(step, Click) and into >= 0.45 and not live.released:
            page.mouse.up()
            live.down = False
            live.released = True

        shot = out_dir / f"f{len(live.frames):05d}.png"
        page.screenshot(path=str(shot), clip=clip, animations="allow")
        live.frames.append((shot, time.perf_counter() - t0, at,
                            now - pressed_at < PRESS_TICK))

        slack = period - (time.perf_counter() - t0 - now)
        if slack > 0:
            time.sleep(slack)

    if entered and live.index < len(story):
        end(story[live.index])
    if live.down:
        page.mouse.up()
    took = time.perf_counter() - t0
    print(f"{len(live.frames)} frames over {took:.1f} s "
          f"({len(live.frames) / took:.1f} per second captured, "
          f"resampled to {fps})")
    return live.frames


def stamp(frames, clip, scale):
    """Draw the pointer into every frame, where it was when the frame was taken."""
    sprite = {False: cursor_sprite(scale), True: cursor_sprite(scale, True)}
    for path, _when, at, holding in frames:
        art, hot = sprite[holding]
        with Image.open(path) as shot:
            frame = shot.convert("RGBA")
        x = int(round((at[0] - clip["x"]) * scale)) - hot
        y = int(round((at[1] - clip["y"]) * scale)) - hot
        frame.alpha_composite(art, (x, y))
        frame.convert("RGB").save(path)


# ----------------------------------------------------------------------
# Encoding
# ----------------------------------------------------------------------


def encode(frames, out: Path, fps: int, width: int, colors: int, dither: str,
           mp4: bool):
    """Resample the frames onto a constant rate and write the gif.

    The frames were taken as fast as the browser would give them up, so they are
    not evenly spaced. What makes them so is the concat demuxer: every frame is
    listed with how long it was actually on screen, and `fps` then samples that
    at a constant rate. A dropped frame costs a repeat rather than a jump.

    The palette is built from the whole run in one pass - `stats_mode=diff`
    weights it towards what is changing, which here is the traces and the curve
    rather than the white the figure mostly is - and `diff_mode=rectangle` lets
    each frame carry only the box that moved, which on a figure this static is
    most of the file size.
    """
    if not frames:
        sys.exit("no frames")
    folder = frames[0][0].parent
    # A frame's duration is how long until the next one, and the last one gets
    # a nominal frame. The demuxer ignores the final entry's duration, hence the
    # repeat: without it the last frame is dropped.
    listing = []
    for i, (path, when, _at, _held) in enumerate(frames):
        nxt = frames[i + 1][1] if i + 1 < len(frames) else when + 1.0 / fps
        listing.append(f"file '{path.name}'\nduration {max(nxt - when, 0.001):.4f}")
    listing.append(f"file '{frames[-1][0].name}'")
    (folder / "frames.txt").write_text("\n".join(listing) + "\n", encoding="utf-8")

    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    chain = (f"fps={fps},scale={width}:-1:flags=lanczos,split[a][b];"
             f"[a]palettegen=max_colors={colors}:stats_mode=diff[p];"
             f"[b][p]paletteuse=dither={dither}:diff_mode=rectangle")
    run = [ffmpeg, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0",
           "-i", "frames.txt", "-vf", chain, "-loop", "0", str(out)]
    subprocess.run(run, cwd=folder, check=True)
    print(f"wrote {out}  ({out.stat().st_size / 1e6:.1f} MB, "
          f"{len(frames)} frames)")

    if mp4:
        # yuv420p needs both sides even, and the height came out of `-1`.
        video = out.with_suffix(".mp4")
        run = [ffmpeg, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0",
               "-i", "frames.txt", "-vf",
               f"fps={fps},scale={width}:-2:flags=lanczos",
               "-c:v", "libx264", "-crf", "20", "-pix_fmt", "yuv420p",
               str(video)]
        subprocess.run(run, cwd=folder, check=True)
        print(f"wrote {video}  ({video.stat().st_size / 1e6:.1f} MB)")


# ----------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    # Up in man/figures/ rather than beside this script, where the other
    # animations keep theirs: this one is in the README, and every image the
    # README names sits directly in man/figures/.
    ap.add_argument("--out", type=Path, default=HERE.parent / "anim_widget.gif")
    ap.add_argument("--start", default="ddm",
                    help="the model the figure opens on, as `cogmod-start`")
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--width", type=int, default=820,
                    help="width of the gif in px; the crop is downsampled to it")
    ap.add_argument("--scale", type=int, default=2,
                    help="device pixel ratio to capture at, before downsampling")
    ap.add_argument("--colors", type=int, default=192)
    ap.add_argument("--dither", default="bayer:bayer_scale=5",
                    help="paletteuse dither; `none` is smaller and flatter")
    ap.add_argument("--mp4", action="store_true",
                    help="also write an mp4 beside the gif")
    ap.add_argument("--reuse", action="store_true",
                    help="skip the quarto render if the page is already there")
    ap.add_argument("--keep", action="store_true",
                    help="keep the rendered page and the png frames")
    args = ap.parse_args(argv)

    html = build_page(args.start, reuse=args.reuse)
    server, port = serve(ARTICLES)
    frames_dir = Path(tempfile.mkdtemp(prefix="anim_widget_"))
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_context(
                viewport=VIEWPORT,
                device_scale_factor=args.scale,
                color_scheme="light",
                # Left to the machine's own setting this would be the difference
                # between a figure that animates and one that does not: the
                # volley is not started at all under a reduce preference.
                reduced_motion="no-preference",
            ).new_page()
            page.goto(f"http://127.0.0.1:{port}/{html.name}")
            page.wait_for_selector('.cogmod-tab[aria-pressed="true"]',
                                   timeout=30_000)
            page.wait_for_selector(".cogmod-figure svg", timeout=30_000)
            page.wait_for_timeout(1200)          # ojs settling, first volley

            # The crop has to hold every model, not the one it opens on: a race
            # brings a second rate and its spread into the column and a
            # ballistic model a start-point range, and each of them makes the
            # widget taller. So every tab is visited once, the boxes unioned,
            # and the figure put back where it started - a gif's frames being
            # all one size whatever the widget does mid-record.
            tabs = page.eval_on_selector_all(
                ".cogmod-tab", "els => els.map(e => e.dataset.model)")
            box = None
            for model in tabs:
                page.click(f'.cogmod-tab[data-model="{model}"]')
                page.wait_for_timeout(180)
                seen = page.evaluate(_EXTENT_JS)
                box = seen if box is None else {
                    "x0": min(box["x0"], seen["x0"]),
                    "y0": min(box["y0"], seen["y0"]),
                    "x1": max(box["x1"], seen["x1"]),
                    "y1": max(box["y1"], seen["y1"])}
            page.click(f'.cogmod-tab[data-model="{args.start}"]')
            page.wait_for_timeout(600)

            pad = 14
            clip = {"x": max(box["x0"] - pad, 0), "y": max(box["y0"] - pad, 0),
                    "width": round(box["x1"] - box["x0"] + 2 * pad),
                    "height": round(box["y1"] - box["y0"] + 2 * pad)}
            if clip["y"] + clip["height"] > VIEWPORT["height"]:
                sys.exit("the widget is taller than the viewport - raise "
                         "VIEWPORT['height'], or the volley will pause when "
                         "the figure scrolls out of view")
            print(f"crop {clip['width']}x{clip['height']} css px, "
                  f"captured at {args.scale}x")

            frames = capture(page, clip, args.scale, args.fps, frames_dir)
            browser.close()
        stamp(frames, clip, args.scale)
        encode(frames, args.out, args.fps, args.width, args.colors,
               args.dither, args.mp4)
    finally:
        server.shutdown()
        server.server_close()
        if args.keep:
            print(f"frames in {frames_dir}")
        else:
            shutil.rmtree(frames_dir, ignore_errors=True)
            drop_page()


if __name__ == "__main__":
    main()
