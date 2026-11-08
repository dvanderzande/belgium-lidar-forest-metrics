"""Walloon 2013-14 survey against DHMV II on the same ground, along the whole Flemish-Walloon border.

Uses the strip of Wallonia within 2 km of Flanders that download_border_strip.py fetched from DHMV II. Both surveys
get the same rules (PICC buildings, 2.5 and 60 m, 3 m closing) and only pixels with data in both are used; the 30 m
cells (grid of process_region.py) that are complete in both are compared. The border mostly runs
west-east, so the result is also given per SECTION m of x.

    source dev.env
    python scripts/compare_border_strip.py "$BELGIUM_LIDAR_WORK_DIR/inputs" "$BELGIUM_LIDAR_WORK_DIR/border_strip"

Output in the last folder: border_strip_cells.csv (every compared cell), border_strip_sections.csv and
figures/border_strip.png.
"""
import argparse
import math
import sys
import time
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
from shapely.geometry import box

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parent))       # the readers of process_region.py
from belgium_lidar_metrics.canopy import GAP_FILL_SIZE, canopy_height, polygon_mask  # noqa: E402
from belgium_lidar_metrics.figures import BACKGROUND, EDGE, INK, MUTED  # noqa: E402
from belgium_lidar_metrics.metrics import CELL, METRIC_NAMES, cell_metrics  # noqa: E402
from process_region import (GRID_ORIGIN, HALO, REGIONS, TILE, dhmv_blocks, read_buildings,  # noqa: E402
                            read_heights, tiles_in, walloon_files)

SECTION = 20_000      # metres of x per summary section
FOREST_COVER = 0.3    # forest cells, as England; fixed on DHMV II so both surveys are compared on the same cells
TALL_P90 = 20.0


def snap(bounds):
    """Bounds rounded outwards to whole 30 m cells of the grid."""
    x0, y0 = GRID_ORIGIN
    xmin, ymin, xmax, ymax = bounds
    return (x0 + math.floor((xmin - x0) / CELL) * CELL, y0 + math.floor((ymin - y0) / CELL) * CELL,
            x0 + math.ceil((xmax - x0) / CELL) * CELL, y0 + math.ceil((ymax - y0) / CELL) * CELL)


def read_survey(files, outer):
    touching = [k for k, f in files.items() if f[0].intersects(box(*outer))]
    if not touching:
        return None
    dsm, transform = read_heights([files[k][1] for k in touching], outer)
    dtm, _ = read_heights([files[k][2] for k in touching], outer)
    return dsm, dtm, transform


def compare_piece(area, zone, walloon, dhmv, wallonia_folder):
    """Cells of one piece of the strip (area = snapped bounds of the strip within a tile) for both surveys."""
    outer = (area[0] - HALO, area[1] - HALO, area[2] + HALO, area[3] + HALO)
    w, d = read_survey(walloon, outer), read_survey(dhmv, outer)
    if w is None or d is None:
        return None
    transform = w[2]
    inside = polygon_mask([zone.intersection(box(*outer))], w[0].shape, transform)
    both = inside & ~np.isnan(w[0]) & ~np.isnan(w[1]) & ~np.isnan(d[0]) & ~np.isnan(d[1])
    buildings = polygon_mask(read_buildings(REGIONS["wallonia"], wallonia_folder, outer, None), w[0].shape, transform)
    metrics = {}
    for name, (dsm, dtm, _) in {"W": w, "D": d}.items():
        height = canopy_height(dsm, dtm, buildings, both, gap_fill_size=GAP_FILL_SIZE)[HALO:-HALO, HALO:-HALO]
        metrics[name] = cell_metrics(height)
    # cells whose 900 pixels all have data (the same in both surveys, as both use the pixels in `both`)
    full = ~np.isnan(height.reshape(height.shape[0] // CELL, CELL, -1, CELL)).any(axis=(1, 3))
    rows, cols = np.nonzero(full)
    table = {"x": area[0] + (cols + 0.5) * CELL, "y": area[3] - (rows + 0.5) * CELL}
    for name in metrics:
        for metric in METRIC_NAMES:
            table[f"{metric}_{name}"] = metrics[name][metric][rows, cols]
    return pd.DataFrame(table)


def summary(cells):
    """Medians per survey and the per-cell difference W - D, in forest cells (fixed on DHMV II)."""
    forest = cells[cells["cover_D"] >= FOREST_COVER]
    tall = cells[cells["p90_D"] > TALL_P90]
    row = {"cells": len(cells), "forest cells": len(forest),
           "forest share D": len(forest) / len(cells), "forest share W": (cells["cover_W"] >= FOREST_COVER).mean(),
           "cover tall forest D": tall["cover_D"].median(), "cover tall forest W": tall["cover_W"].median()}
    for metric in METRIC_NAMES:
        diff = forest[f"{metric}_W"] - forest[f"{metric}_D"]
        row[f"{metric} D"] = forest[f"{metric}_D"].median()
        row[f"{metric} W"] = forest[f"{metric}_W"].median()
        row[f"{metric} bias W-D"] = diff.mean()
        row[f"{metric} corr"] = forest[[f"{metric}_W", f"{metric}_D"]].corr().iloc[0, 1]
    return row


def draw(sections, out_path):
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5), constrained_layout=True, facecolor=BACKGROUND)
    x = sections.index / 1000
    for ax, (column, label) in zip(axes, [("cover bias W-D", "cover"), ("p90 bias W-D", "p90 (m)"),
                                          ("mean bias W-D", "mean height (m)"), ("forest cells", "forest cells")]):
        ax.bar(x, sections[column], width=SECTION / 1000 * 0.8, color=INK if column == "forest cells" else "#2a6fdb")
        ax.axhline(0, color=MUTED, linewidth=1)
        ax.set_title(label if column == "forest cells" else f"{label}: Wallonia minus DHMV II", color=INK,
                     fontsize=10, loc="left")
        ax.set_xlabel("x of the section (km, EPSG:31370; west to east)", color=MUTED, fontsize=9)
        ax.tick_params(colors=MUTED)
        for spine in ax.spines.values():
            spine.set_color(EDGE)
    fig.suptitle(f"Same ground, both surveys: forest cells (DHMV II cover >= {FOREST_COVER}) in Wallonia within 2 km "
                 f"of Flanders, mean difference per {SECTION // 1000} km section", color=INK, fontsize=12, x=0.01,
                 ha="left")
    fig.savefig(out_path, dpi=120, facecolor=fig.get_facecolor())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders")
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    (args.out_dir / "figures").mkdir(parents=True, exist_ok=True)
    strip = args.inputs_dir / "dhmv2_wallonia_strip"
    zone = gpd.read_file(strip / "strip_zone.gpkg").union_all()
    wallonia_folder = args.inputs_dir / REGIONS["wallonia"]["folder"]
    walloon, dhmv = walloon_files(wallonia_folder), dhmv_blocks(strip)

    tiles = tiles_in(zone)
    print(f"strip {zone.area / 1e6:.0f} km2 in {len(tiles)} tiles", flush=True)
    parts, start = [], time.time()
    for i, (x, y) in enumerate(tiles, start=1):
        piece = zone.intersection(box(x, y, x + TILE, y + TILE))
        if piece.is_empty:
            continue
        cells = compare_piece(snap(piece.bounds), zone, walloon, dhmv, wallonia_folder)
        if cells is not None:
            parts.append(cells)
        print(f"tile {i}/{len(tiles)} done; elapsed {(time.time() - start) / 60:.1f} min", flush=True)
    cells = pd.concat(parts, ignore_index=True).drop_duplicates(["x", "y"])
    cells.to_csv(args.out_dir / "border_strip_cells.csv", index=False)

    overall = pd.Series(summary(cells))
    print("\nwhole strip:\n" + overall.round(3).to_string())
    cells["section"] = (cells["x"] // SECTION * SECTION + SECTION / 2).astype(int)
    sections = pd.DataFrame({s: summary(group) for s, group in cells.groupby("section")
                             if (group["cover_D"] >= FOREST_COVER).sum() >= 200}).T
    sections.to_csv(args.out_dir / "border_strip_sections.csv")
    print("\nper section (x centre, sections with at least 200 forest cells):")
    print(sections[["forest cells", "cover bias W-D", "p90 bias W-D", "mean bias W-D", "p90 corr"]].round(2)
          .to_string())
    draw(sections, args.out_dir / "figures" / "border_strip.png")
    print("figure: figures/border_strip.png")


if __name__ == "__main__":
    main()
