"""Checks of a 30 m metrics file made by process_region.py: grid and tags, value
ranges against the England reference (docs/methods.md), seams at processing tiles and download
blocks, and a map.

    source dev.env
    python scripts/check_metrics.py "$BELGIUM_LIDAR_WORK_DIR/flanders/metrics_30m_flanders.tif"

Output: printed checks and figures/<file name>.png next to the file.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import rasterio
from scipy import ndimage

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.canopy import GAP_FILL_SIZE  # noqa: E402
from belgium_lidar_metrics.figures import BACKGROUND, BLUES, EDGE, INK, MUTED  # noqa: E402

CELL = 30
SEAMS = {"processing tiles (6 km)": 6000, "download blocks (10 km)": 10_000}
FOREST_COVER = 0.3   # England forest threshold


def seam_ratio(band, transform, spacing, axis):
    """Mean |difference| between neighbouring cells across the seams, divided by the same elsewhere. About 1 means
    the seams are invisible; well above 1 means an edge effect. NaN if no seam line with data lies in the band.
    Only seams on cell edges can be tested: 10 km block edges meet the 30 m cell edges every 30 km."""
    diff = np.abs(np.diff(band, axis=axis))
    if axis == 1:   # between column j and j + 1: the boundary lies at x of column j + 1
        coords = transform.c + CELL * np.arange(1, band.shape[1])
    else:           # between row i and i + 1: the boundary lies at y of the bottom of row i
        coords = transform.f - CELL * np.arange(1, band.shape[0])
    at_seam = np.round(coords) % spacing == 0
    take = (lambda m: diff[:, m]) if axis == 1 else (lambda m: diff[m, :])
    if np.isnan(take(at_seam)).all():
        return np.nan
    return np.nanmean(take(at_seam)) / np.nanmean(take(~at_seam))


def draw(bands, heading, out_path):
    # a wide area (Flanders) gets the two maps above each other, any other shape side by side
    wide = bands[0].shape[1] > 1.6 * bands[0].shape[0]
    fig, axes = plt.subplots(*((2, 1) if wide else (1, 2)), figsize=(12, 10) if wide else (16, 8),
                             constrained_layout=True, facecolor=BACKGROUND)
    for ax, (band, vmax, title, label) in zip(axes, [
            (bands[0], 1, "Canopy cover (share of 1 m pixels > 2 m)", "cover"),
            (bands[1], 35, "90th percentile height", "height (m)")]):
        im = ax.imshow(band, cmap=BLUES, vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_title(title, color=INK, fontsize=11, loc="left")
        ax.set_xticks([]), ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(EDGE)
        bar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.01)
        bar.set_label(label, color=MUTED)
        bar.ax.tick_params(colors=MUTED)
    fig.suptitle(heading, color=INK, fontsize=13, x=0.01, ha="left")
    fig.text(0.01, 0.003, "EPSG:31370, cells on multiples of 30 m. Gray = NoData (outside the area).", color=MUTED,
             fontsize=9)
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", type=Path)
    args = parser.parse_args()
    with rasterio.open(args.file) as src:
        bands, t = src.read(), src.transform
        print(f"grid: {src.width} x {src.height} cells, {src.crs}, cell {t.a} m, origin ({t.c:.0f}, {t.f:.0f}), "
              f"on multiples of 30: {t.c % CELL == 0 and t.f % CELL == 0}")
        print("bands:", src.descriptions)
        print("tags:", {k: v for k, v in src.tags().items() if k != "AREA_OR_POINT"})
        source = src.tags().get("source", "")
    cover, p90, sd, mean = bands

    # Cells away from the Flanders border (all 8 neighbours have data): border cells are partly NoData and noisier
    valid = ~np.isnan(cover)
    interior = ndimage.binary_erosion(valid, structure=np.ones((3, 3)))
    print(f"cells with data: {valid.sum():,} ({valid.sum() * CELL * CELL / 1e6:,.0f} km2), away from the border: "
          f"{interior.sum():,}; sd NaN in cells with data: {np.isnan(sd[valid]).sum()}")
    for name, band in zip(["cover", "p90", "sd", "mean"], bands):
        q = np.nanpercentile(band[interior], [5, 50, 95])   # sd is NaN in a cell with a single pixel
        print(f"{name:5s} median {q[1]:.2f} (5-95 %: {q[0]:.2f}-{q[2]:.2f}), max {np.nanmax(band):.2f}")
    forest = interior & (cover >= FOREST_COVER)
    print(f"forest cells (cover >= {FOREST_COVER}): {forest.sum() / interior.sum():.1%}; median there: p90 "
          f"{np.median(p90[forest]):.1f} m, sd {np.nanmedian(sd[forest]):.1f} m, mean {np.median(mean[forest]):.1f} m, "
          f"cover {np.median(cover[forest]):.2f}")
    print(f"cells with p90 > 50 m: {(p90 > 50).sum()}")

    for seam, spacing in SEAMS.items():
        ratios = [r for band in bands for axis in (0, 1) if np.isfinite(r := seam_ratio(band, t, spacing, axis))]
        if ratios:
            print(f"seams at {seam}: neighbour difference across / elsewhere = {min(ratios):.2f}-{max(ratios):.2f} "
                  f"({len(ratios)} of 8 band and direction tests have a seam line with data)")
        else:
            print(f"seams at {seam}: no seam line with data in the area")

    figures = args.file.parent / "figures"
    figures.mkdir(exist_ok=True)
    heading = f"{args.file.stem}: 30 m metrics, {GAP_FILL_SIZE} m gap filling\n" + "\n".join(source.split(" | "))
    draw(bands, heading, figures / f"{args.file.stem}.png")
    print(f"figure: figures/{args.file.stem}.png")


if __name__ == "__main__":
    main()
