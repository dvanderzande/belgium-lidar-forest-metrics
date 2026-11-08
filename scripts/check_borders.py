"""Is there a step in the metrics at the regional borders of the Belgian map (made by process_region.py)?

For each border, every 30 m cell gets its distance to the border (negative on the first region's side, positive on the
second's). The metrics are summarised in BAND m bands up to REACH m on either side. A difference between the surveys
shows as a step at 0; a change of landscape shows as a gradual trend. Flanders-Brussels is a control: both sides use
the same survey (DHMV II), only the building data differ.

    source dev.env
    python scripts/check_borders.py "$BELGIUM_LIDAR_WORK_DIR/belgium/metrics_30m_belgium.tif" "$BELGIUM_LIDAR_WORK_DIR/inputs"

Output: printed tables and figures/regional_borders.png next to the metrics file.
"""
import argparse
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.canopy import polygon_mask  # noqa: E402
from belgium_lidar_metrics.figures import BACKGROUND, EDGE, INK, MUTED  # noqa: E402

BAND = 500       # metres per distance band
REACH = 5000     # metres on either side of the border
FOREST_COVER = 0.3
BORDERS = {"Flanders | Wallonia": ("flanders", "wallonia"), "Flanders | Brussels": ("flanders", "brussels")}
BORDER_FILES = {"flanders": "flanders_dhmv2/VRBG_Refgew_Flanders.gpkg",
                "brussels": "brussels_dhmv2/UrbIS_Region_Brussels.gpkg",
                "wallonia": "wallonia_2013_2014/SPW_Region_Wallonia.gpkg"}


def signed_distance(first, second):
    """Distance (m) of each cell to the border between two region masks: negative in the first, positive in the second,
    NaN beyond REACH or outside both."""
    # cells of the first region that touch the second form the border line; the border itself lies half a cell
    # beyond them, so distances are shifted by 15 m (the line cells get -15 m, their neighbours across +15 m)
    line = first & ndimage.binary_dilation(second, structure=np.ones((3, 3)))
    distance = ndimage.distance_transform_edt(~line) * 30.0
    signed = np.where(first, -(distance + 15), np.where(second, distance - 15, np.nan))
    signed[np.abs(signed) > REACH] = np.nan
    return signed


def profile(bands, signed):
    """Per distance band: share of forest cells, and medians in forest cells."""
    cover, p90, sd, mean = bands
    edges = np.arange(-REACH, REACH + 1, BAND)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        cells = (signed >= lo) & (signed < hi) & ~np.isnan(cover)
        forest = cells & (cover >= FOREST_COVER)
        rows.append({"distance (m)": (lo + hi) / 2, "cells": int(cells.sum()),
                     "forest share": forest.sum() / max(cells.sum(), 1),
                     "cover (forest)": np.median(cover[forest]), "p90 (forest, m)": np.median(p90[forest]),
                     "SD (forest, m)": np.nanmedian(sd[forest]), "mean (forest, m)": np.median(mean[forest])})
    return pd.DataFrame(rows)


def step(table):
    """Difference between the two bands next to the border (second region minus first)."""
    inner = table[np.abs(table["distance (m)"]) == BAND / 2].set_index("distance (m)")
    return inner.loc[BAND / 2] - inner.loc[-BAND / 2]


def draw(tables, out_path):
    columns = ["forest share", "cover (forest)", "p90 (forest, m)", "mean (forest, m)"]
    fig, axes = plt.subplots(len(tables), len(columns), figsize=(16, 4 * len(tables)), constrained_layout=True,
                             facecolor=BACKGROUND, squeeze=False)
    for i, (name, table) in enumerate(tables.items()):
        for j, column in enumerate(columns):
            ax = axes[i, j]
            ax.plot(table["distance (m)"] / 1000, table[column], color=INK, marker="o", markersize=4, linewidth=1.5)
            ax.axvline(0, color=MUTED, linewidth=1, linestyle="--")
            ax.set_title(f"{name}: {column}", color=INK, fontsize=10, loc="left")
            ax.set_xlabel(f"distance to the border (km); left = {name.split(' | ')[0]}", color=MUTED, fontsize=9)
            ax.tick_params(colors=MUTED)
            for spine in ax.spines.values():
                spine.set_color(EDGE)
    fig.suptitle(f"Metrics in {BAND} m bands either side of the regional borders (forest: cover >= {FOREST_COVER})",
                 color=INK, fontsize=13, x=0.01, ha="left")
    fig.savefig(out_path, dpi=120, facecolor=fig.get_facecolor())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", type=Path, help="metrics_30m_belgium.tif")
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders (border files)")
    args = parser.parse_args()
    with rasterio.open(args.file) as src:
        bands, transform, shape = src.read(), src.transform, (src.height, src.width)
    masks = {name: polygon_mask(gpd.read_file(args.inputs_dir / path).geometry, shape, transform)
             for name, path in BORDER_FILES.items()}

    tables = {}
    for name, (a, b) in BORDERS.items():
        tables[name] = profile(bands, signed_distance(masks[a], masks[b]))
        print(f"\n{name} (negative distance = {a}):")
        print(tables[name].round(3).to_string(index=False))
        print(f"step at the border ({b} minus {a}, bands next to it):")
        print(step(tables[name]).drop(["distance (m)", "cells"], errors="ignore").round(3).to_string())
    figures = args.file.parent / "figures"
    figures.mkdir(exist_ok=True)
    draw(tables, figures / "regional_borders.png")
    print("figure: figures/regional_borders.png")


if __name__ == "__main__":
    main()
