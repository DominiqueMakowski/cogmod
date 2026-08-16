"""Put a bold title on the two vignette figures, so each one stands alone.

Both come from the package vignettes, which need fitted brms models and so
are not rebuilt here. This script only adds a title band on top and writes
the result into paper/figures/. (fig_motivation carries its own title,
written by make_figures.R.)

Title sizes are set from how wide each image ends up in the PDF, so that all
three headings print at the same size regardless of the source resolution.

    python add_figure_titles.py

Writes figures/fig_ppcheck.png and figures/fig_choco.png.
"""

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VIG = os.path.join(ROOT, "man", "figures")
OUTDIR = os.path.join(HERE, "figures")

INK, MUTED = "#22252E", "#6B7280"
DPI = 300
TITLE_PT, SUB_PT = 10.5, 8.0             # sizes as printed, not as drawn

FIGURES = [
    dict(src=os.path.join(VIG, "rt_models1.png"),
         dst=os.path.join(OUTDIR, "fig_ppcheck.png"),
         shown_in=6.5,
         title="Nine reaction time families fitted to the same trials",
         subtitle="Posterior predictive draws over the observed data"),
    dict(src=os.path.join(VIG, "subjective_ratings1.png"),
         dst=os.path.join(OUTDIR, "fig_choco.png"),
         shown_in=6.5, crop_top=90,      # drops the plot's own ggplot title
         title="Four models of bimodal slider ratings",
         subtitle="Only CHOCO reproduces both modes and the endpoint masses"),
]


def add_title(src, dst, shown_in, title, subtitle, crop_top=0):
    img = mpimg.imread(src)[crop_top:]
    h_px, w_px = img.shape[0], img.shape[1]
    w_in, h_in = w_px / DPI, h_px / DPI
    scale = shown_in / w_in              # shrinkage applied by the PDF
    band = 0.44 / scale                  # so the band prints at 0.44 in
    total = h_in + band
    frac = h_in / total

    fig = plt.figure(figsize=(w_in, total), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, frac])
    ax.imshow(img, interpolation="none")
    ax.axis("off")

    x = 0.012
    fig.text(x, 1 - 0.42 * (1 - frac), title, fontsize=TITLE_PT / scale,
             weight="bold", color=INK, ha="left", va="center")
    fig.text(x, 1 - 0.80 * (1 - frac), subtitle, fontsize=SUB_PT / scale,
             style="italic", color=MUTED, ha="left", va="center")

    fig.savefig(dst, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"wrote {os.path.relpath(dst, HERE)}  ({w_px}x{h_px} -> shown at "
          f"{shown_in} in, title {TITLE_PT / scale:.1f} pt)")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    for f in FIGURES:
        if not os.path.exists(f["src"]):
            raise SystemExit(f"missing source: {f['src']}")
        add_title(f["src"], f["dst"], f["shown_in"], f["title"], f["subtitle"],
                  f.get("crop_top", 0))


if __name__ == "__main__":
    main()
