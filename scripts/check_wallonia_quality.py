"""Where does the Walloon 2013-14 survey miss forest? Our metrics in the forest of the SPW forest mask 2018.

Along the Flemish border (compare_border_strip.py) the Walloon DSM lies on the ground in leaf-off broadleaf forest in
some flight blocks. Inside Wallonia there is no second survey, so the SPW forest mask 2018 (download_forest_mask.py;
from orthophotos, independent of both LiDAR surveys) stands in for it: in cells that are forest in the mask, a good
survey finds forest (cover >= 0.3) in most of them. That share is mapped per TILE m tile of process_region.py. The
strip along the border, where DHMV II exists as well, shows what a good survey finds in mask forest.

    source dev.env
    python scripts/check_wallonia_quality.py "$BELGIUM_LIDAR_WORK_DIR/belgium/metrics_30m_belgium.tif" \\
        "$BELGIUM_LIDAR_WORK_DIR/inputs" "$BELGIUM_LIDAR_WORK_DIR/wallonia_quality"

The border strip cells are read from border_strip/border_strip_cells.csv next to the output folder.
Output in the last folder: wallonia_quality_tiles.csv, border_strip_mask_forest.csv and figures/wallonia_quality.png.
"""
import argparse
import math
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
from shapely import affinity

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parent))       # the grid of process_region.py
from belgium_lidar_metrics.canopy import polygon_mask  # noqa: E402
from belgium_lidar_metrics.figures import BACKGROUND, BLUES, EDGE, INK, MUTED  # noqa: E402
from belgium_lidar_metrics.metrics import CELL  # noqa: E402
from process_region import GRID_ORIGIN, TILE  # noqa: E402

FOREST_COVER = 0.3   # a cell is forest, as England
MASK_FOREST = 0.99   # cells entirely forest in the mask (all 36 samples), so edges do not count
MIN_CELLS = 400      # mask-forest cells (36 ha) a tile needs to be shown
SECTION = 20_000     # metres of x per section of the border strip, as compare_border_strip.py


def read_window(path, bounds, bands):
    """Bands of a raster on the 30 m grid over bounds (whole cells), and the window's transform."""
    with rasterio.open(path) as src:
        window = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
        return src.read(bands, window=window), src.window_transform(window)


def cell_frame(rows, cols, transform, **values):
    """Cell centres and values of the given cells as a table."""
    return pd.DataFrame({"x": transform.c + (cols + 0.5) * CELL, "y": transform.f - (rows + 0.5) * CELL} | values)


def tile_table(cover, p90, mask_forest, transform):
    """Per tile, in mask-forest cells: how many, the share the survey finds as forest, median cover and p90."""
    rows, cols = np.nonzero(mask_forest & ~np.isnan(cover))
    cells = cell_frame(rows, cols, transform, cover=cover[rows, cols], p90=p90[rows, cols])
    x0, y0 = GRID_ORIGIN
    cells["tile_x"] = ((cells["x"] - x0) // TILE * TILE + x0).astype(int)
    cells["tile_y"] = ((cells["y"] - y0) // TILE * TILE + y0).astype(int)
    cells["found"] = cells["cover"] >= FOREST_COVER
    table = cells.groupby(["tile_x", "tile_y"]).agg(cells=("cover", "size"), found=("found", "mean"),
                                                    cover=("cover", "median"), p90=("p90", "median"))
    return table[table["cells"] >= MIN_CELLS].reset_index()


def strip_table(strip_cells, share, transform):
    """The border strip (both surveys) in mask-forest cells, per SECTION and overall."""
    cols = ((strip_cells["x"] - transform.c) // CELL).astype(int)
    rows = ((transform.f - strip_cells["y"]) // CELL).astype(int)
    strip = strip_cells[share[rows, cols] >= MASK_FOREST].copy()
    strip["section"] = (strip["x"] // SECTION * SECTION + SECTION / 2).astype(int)

    def summary(group):
        row = {"cells": len(group)}
        for survey, name in [("D", "DHMV II"), ("W", "Wallonia")]:
            row[f"found {name}"] = (group[f"cover_{survey}"] >= FOREST_COVER).mean()
            row[f"cover {name}"] = group[f"cover_{survey}"].median()
            row[f"p90 {name}"] = group[f"p90_{survey}"].median()
        return pd.Series(row)

    table = strip.groupby("section").apply(summary, include_groups=False)
    table.loc["whole strip"] = summary(strip)
    return table


def tile_image(table, column, bounds):
    """A tile column as an image over bounds (rows from the top), NaN where no tile is shown."""
    x0, y0 = GRID_ORIGIN
    left = x0 + math.floor((bounds[0] - x0) / TILE) * TILE
    bottom = y0 + math.floor((bounds[1] - y0) / TILE) * TILE
    nx, ny = math.ceil((bounds[2] - left) / TILE), math.ceil((bounds[3] - bottom) / TILE)
    image = np.full((ny, nx), np.nan)
    image[ny - 1 - (table["tile_y"] - bottom) // TILE, (table["tile_x"] - left) // TILE] = table[column]
    return image, (left / 1000, (left + nx * TILE) / 1000, bottom / 1000, (bottom + ny * TILE) / 1000)


def draw(tiles, wallonia, reference, out_path):
    outline = gpd.GeoSeries([affinity.scale(wallonia.boundary, 1e-3, 1e-3, origin=(0, 0))])
    fig, axes = plt.subplots(2, 1, figsize=(12, 13), constrained_layout=True, facecolor=BACKGROUND)
    panels = [("found", "Share of mask-forest cells the Walloon survey finds as forest (cover >= 0.3)", 0, 1),
              ("p90", "Median p90 height in mask-forest cells (m)", 0, 30)]
    for ax, (column, title, vmin, vmax) in zip(axes, panels):
        image, extent = tile_image(tiles, column, wallonia.bounds)
        shown = ax.imshow(image, extent=extent, cmap=BLUES, vmin=vmin, vmax=vmax, interpolation="nearest")
        outline.plot(ax=ax, color=INK, linewidth=0.8)
        ax.set_title(title, color=INK, fontsize=11, loc="left")
        ax.set_xlabel("x (km, EPSG:31370)", color=MUTED, fontsize=9)
        ax.set_ylabel("y (km)", color=MUTED, fontsize=9)
        ax.tick_params(colors=MUTED)
        for spine in ax.spines.values():
            spine.set_color(EDGE)
        fig.colorbar(shown, ax=ax, shrink=0.7).ax.tick_params(colors=MUTED)
    fig.suptitle(f"Walloon 2013-14 survey in the forest of the SPW forest mask 2018, per {TILE // 1000} km tile "
                 f"(tiles with at least {MIN_CELLS} mask-forest cells)\nAlong the Flemish border, in mask forest: "
                 f"DHMV II finds forest in {reference['found DHMV II']:.0%} of cells, the Walloon survey in "
                 f"{reference['found Wallonia']:.0%}", color=INK, fontsize=12, x=0.01, ha="left")
    fig.savefig(out_path, dpi=110, facecolor=fig.get_facecolor())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", type=Path, help="metrics_30m_belgium.tif")
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders")
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    (args.out_dir / "figures").mkdir(parents=True, exist_ok=True)
    wallonia = gpd.read_file(args.inputs_dir / "wallonia_2013_2014" / "SPW_Region_Wallonia.gpkg").union_all()
    x0, y0 = GRID_ORIGIN
    bounds = (x0 + math.floor((wallonia.bounds[0] - x0) / CELL) * CELL,
              y0 + math.floor((wallonia.bounds[1] - y0) / CELL) * CELL,
              x0 + math.ceil((wallonia.bounds[2] - x0) / CELL) * CELL,
              y0 + math.ceil((wallonia.bounds[3] - y0) / CELL) * CELL)

    (cover, p90), transform = read_window(args.file, bounds, [1, 2])
    (share,), mask_transform = read_window(args.inputs_dir / "spw_forest_mask_2018" /
                                           "forest_mask_2018_share_30m.tif", bounds, [1])
    assert mask_transform == transform and share.shape == cover.shape, "the two rasters are not on the same grid"
    inside = polygon_mask([wallonia], cover.shape, transform)
    mask_forest = inside & (share >= MASK_FOREST)
    any_forest = (inside & (share > 0)).sum()
    print(f"Wallonia: {inside.sum():,} cells; mask forest (entire cell) {mask_forest.sum():,} "
          f"({mask_forest.sum() / inside.sum():.1%}); any mask forest {any_forest / inside.sum():.1%}")
    print(f"all of Wallonia, in mask forest: survey finds forest in "
          f"{np.mean(cover[mask_forest & ~np.isnan(cover)] >= FOREST_COVER):.1%}")

    strip = strip_table(pd.read_csv(args.out_dir.parent / "border_strip" / "border_strip_cells.csv"), share, transform)
    strip.to_csv(args.out_dir / "border_strip_mask_forest.csv")
    print("\nborder strip, mask-forest cells, both surveys (section = x centre):")
    print(strip.round(2).to_string())

    tiles = tile_table(cover, p90, mask_forest, transform)
    tiles.to_csv(args.out_dir / "wallonia_quality_tiles.csv", index=False)
    print(f"\n{len(tiles)} tiles with at least {MIN_CELLS} mask-forest cells; share found per tile:")
    print(tiles["found"].describe(percentiles=[0.1, 0.25, 0.5, 0.75, 0.9]).round(2).to_string())
    draw(tiles, wallonia, strip.loc["whole strip"], args.out_dir / "figures" / "wallonia_quality.png")
    print("figure: figures/wallonia_quality.png")


if __name__ == "__main__":
    main()
