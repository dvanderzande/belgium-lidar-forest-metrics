"""Quality layer for the Belgian metrics: which survey each 30 m cell comes from, and where the Walloon 2013-14 survey
misses forest.

In some flight blocks of the Walloon 2013-14 survey the DSM lies on the ground in leaf-off broadleaf forest, so cover
and heights are far too low there (compare_border_strip.py, check_wallonia_quality.py). SPW publishes no flight blocks
for that survey, so the zone is drawn from the data: ZONE holds rectangles whose edges are where the share of SPW
mask forest that the survey finds jumps (per 1 km row and column; see docs/methods.md).

    source dev.env
    python scripts/make_quality_mask.py "$BELGIUM_LIDAR_WORK_DIR/belgium/metrics_30m_belgium.tif" \\
        "$BELGIUM_LIDAR_WORK_DIR/inputs"

Output next to the metrics file: quality_30m_belgium.tif (codes in CODES, same grid as the metrics),
quality_zone_wallonia_2013_2014.gpkg (the zone as a polygon) and figures/quality_30m_belgium.png.
"""
import argparse
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import matplotlib
import rasterio
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from rasterio.transform import array_bounds
from rasterio.windows import from_bounds
from shapely.geometry import box
from shapely.ops import unary_union

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parent))       # the measure of check_wallonia_quality.py
from belgium_lidar_metrics.canopy import polygon_mask  # noqa: E402
from belgium_lidar_metrics.figures import BACKGROUND, EDGE, MUTED  # noqa: E402
from check_wallonia_quality import FOREST_COVER, MASK_FOREST, read_window  # noqa: E402

# Rectangles (xmin, ymin, xmax, ymax in m, EPSG:31370) of the affected zone; only their part in Wallonia counts.
ZONE = [
    (0, 76_000, 162_000, 138_000),          # Hainaut south of Tournai and Ath, west of x 162 km
    (162_000, 76_000, 194_000, 98_000),     # Entre-Sambre-et-Meuse west of the Meuse
    (84_000, 158_000, 118_000, 170_000),    # the strip north of y 158 km around Flobecq, Lessines, Enghien
]
CODES = {0: "no data", 1: "DHMV II (Flanders, Brussels)", 2: "Wallonia MNS - MNT 2013-2014, reliable",
         3: "Wallonia MNS - MNT 2013-2014, zone where leaf-off broadleaf forest is missed (cover and heights too low)"}
MADE_BY = "belgium-lidar-forest-metrics/scripts/make_quality_mask.py"


def draw(quality, bounds, out_path):
    """Map of the quality layer at 150 m (every 5th cell is enough for a page-wide map)."""
    left, bottom, right, top = (v / 1000 for v in bounds)
    shown = np.ma.masked_equal(quality[::5, ::5], 0)
    palette = {1: ("#2a6fdb", "1  DHMV II (Flanders, Brussels)"), 2: ("#7fb77e", "2  Wallonia 2013–14, reliable"),
               3: ("#e07b39", "3  Wallonia 2013–14, broadleaf forest missed")}
    fig, ax = plt.subplots(figsize=(9, 6.2), facecolor=BACKGROUND, constrained_layout=True)
    ax.imshow(shown, cmap=ListedColormap([palette[k][0] for k in (1, 2, 3)]), vmin=0.5, vmax=3.5,
              extent=(left, right, bottom, top), interpolation="nearest")
    ax.legend(handles=[Patch(color=c, label=label) for c, label in palette.values()], loc="lower left",
              frameon=False, fontsize=9)
    ax.set_xlabel("x (km, EPSG:31370)", color=MUTED, fontsize=9)
    ax.set_ylabel("y (km)", color=MUTED, fontsize=9)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_facecolor(BACKGROUND)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)


def check(mask_path, cover, quality, profile):
    """The share of SPW mask forest the survey finds as forest, per code: low in the zone, high elsewhere."""
    with rasterio.open(mask_path) as src:
        mask_bounds = src.bounds
    metrics_bounds = array_bounds(profile["height"], profile["width"], profile["transform"])
    overlap = (max(mask_bounds[0], metrics_bounds[0]), max(mask_bounds[1], metrics_bounds[1]),
               min(mask_bounds[2], metrics_bounds[2]), min(mask_bounds[3], metrics_bounds[3]))
    window = from_bounds(*overlap, transform=profile["transform"]).round_offsets().round_lengths()
    (share,), _ = read_window(mask_path, overlap, [1])
    rows, cols = window.toslices()
    cover, quality = cover[rows, cols], quality[rows, cols]
    for code in (2, 3):
        forest = (quality == code) & (share >= MASK_FOREST)
        print(f"code {code}: survey finds {np.mean(cover[forest] >= FOREST_COVER):.1%} of {forest.sum():,} "
              "mask-forest cells")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", type=Path, help="metrics_30m_belgium.tif")
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders")
    args = parser.parse_args()
    wallonia = gpd.read_file(args.inputs_dir / "wallonia_2013_2014" / "SPW_Region_Wallonia.gpkg").union_all()
    zone = wallonia.intersection(unary_union([box(*r) for r in ZONE]))

    with rasterio.open(args.file) as src:
        cover, profile = src.read(1), src.profile
    # a cell belongs to the region its centre lies in (border cells mix both surveys at 1 m)
    in_wallonia = polygon_mask([wallonia], cover.shape, profile["transform"])
    in_zone = polygon_mask([zone], cover.shape, profile["transform"])
    quality = np.where(in_zone, 3, np.where(in_wallonia, 2, 1)).astype("uint8")
    quality[np.isnan(cover)] = 0

    out = args.file.with_name("quality_30m_belgium.tif")
    profile.update(dtype="uint8", nodata=0, count=1, predictor=1)
    with rasterio.open(out, "w", **profile) as dst:
        dst.write(quality, 1)
        dst.set_band_description(1, "source_quality")
        dst.update_tags(made_by=MADE_BY, **{f"code_{k}": v for k, v in CODES.items()},
                        zone_rectangles=str(ZONE))
    gpd.GeoDataFrame({"name": ["Wallonia 2013-2014: leaf-off broadleaf forest missed"], "made_by": [MADE_BY]},
                     geometry=[zone], crs="EPSG:31370").to_file(out.with_name("quality_zone_wallonia_2013_2014.gpkg"))

    counts = np.bincount(quality.ravel(), minlength=4)
    print(f"zone {zone.area / 1e6:,.0f} km2 ({zone.area / wallonia.area:.1%} of Wallonia)")
    for code in (1, 2, 3):
        print(f"code {code}: {counts[code]:,} cells ({counts[code] * 900 / 1e6:,.0f} km2) - {CODES[code]}")

    figures = out.parent / "figures"
    figures.mkdir(exist_ok=True)
    draw(quality, array_bounds(profile["height"], profile["width"], profile["transform"]),
         figures / "quality_30m_belgium.png")
    check(args.inputs_dir / "spw_forest_mask_2018" / "forest_mask_2018_share_30m.tif", cover, quality, profile)
    print(f"written: {out.name}, quality_zone_wallonia_2013_2014.gpkg")


if __name__ == "__main__":
    main()
