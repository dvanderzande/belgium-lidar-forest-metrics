"""Download the SPW forest mask 2018 of Wallonia as the forest share of each 30 m cell of our grid.

The mask (Gembloux Agro-Bio Tech / ULiège, from orthophotos 2006-2018, woody vegetation of about 2 m and more, 0/1)
is independent of both LiDAR surveys; check_wallonia_quality.py uses it to find where the Walloon 2013-14 survey
misses forest. The service only renders images (at most 4096 px), in Lambert 2008 (EPSG:3812). So per CHUNK the image
is fetched in EPSG:3812, reprojected here with PROJ (nearest neighbour, SAMPLE m pixels) and averaged to the 30 m
cells. Checked on a 3 km box (2026-11-08): the same as the server's own reprojection to EPSG:31370 (99.6 % of samples,
no shift), so the datum shift between the two Lambert systems is handled.

    source dev.env
    python scripts/download_forest_mask.py "$BELGIUM_LIDAR_WORK_DIR/inputs"

Output: <inputs>/spw_forest_mask_2018/forest_mask_2018_share_30m.tif (0-1; also 0 beyond the extent of the mask, NaN
in chunks that do not touch Wallonia) and download_info.json.
"""
import argparse
import io
import json
import math
import sys
import time
from datetime import date
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import Affine, from_bounds
from rasterio.warp import Resampling, reproject, transform_bounds
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.download import CRS, fetch  # noqa: E402
from belgium_lidar_metrics.metrics import CELL  # noqa: E402
from belgium_lidar_metrics.output import write_raster  # noqa: E402

MASK_URL = "https://geoservices.wallonie.be/arcgis/rest/services/FORET/MASQUES_FORESTIERS/MapServer"
LAYER = 1           # "Masque forestier 2018"; layer 0 is 2021-22, partly derived from the 2021-22 LiDAR
MASK_CRS = "EPSG:3812"
SAMPLE = 5          # metres per sample: 36 samples per 30 m cell
CHUNK = 666 * CELL  # metres per request (19.98 km): 3996 samples, within the 4096 px limit, whole 30 m cells
MARGIN = 100        # metres added around a chunk in EPSG:3812, so the reprojected chunk is covered


def fetch_mask(bounds_3812):
    """The rendered mask over a box in EPSG:3812 (widened to whole SAMPLE m) as (forest bool array, transform)."""
    left, bottom = (math.floor(v / SAMPLE) * SAMPLE for v in bounds_3812[:2])
    right, top = (math.ceil(v / SAMPLE) * SAMPLE for v in bounds_3812[2:])
    width, height = (right - left) // SAMPLE, (top - bottom) // SAMPLE
    # the image is asked for directly (f=image): the links of f=json answers gave "Invalid URL" (2026-11-08). The box
    # is a whole number of SAMPLE m pixels in both directions, so the server need not widen it to fit the image size
    params = {"bbox": f"{left},{bottom},{right},{top}", "bboxSR": 3812, "imageSR": 3812, "size": f"{width},{height}",
              "layers": f"show:{LAYER}", "format": "png32", "transparent": "true", "f": "image"}
    image, _ = fetch(f"{MASK_URL}/export", params)
    with rasterio.open(io.BytesIO(image)) as src:
        alpha = src.read(4)
    assert alpha.shape == (height, width), f"image of {alpha.shape}, asked {height} x {width}"
    return alpha > 0, from_bounds(left, bottom, right, top, width, height)


def chunk_share(x, y):
    """Forest share of the 30 m cells of one chunk (lower-left corner x, y in EPSG:31370)."""
    outer = transform_bounds(CRS, MASK_CRS, x, y, x + CHUNK, y + CHUNK, densify_pts=21)
    forest, source_transform = fetch_mask((outer[0] - MARGIN, outer[1] - MARGIN,
                                           outer[2] + MARGIN, outer[3] + MARGIN))
    n = CHUNK // SAMPLE
    samples = np.zeros((n, n), dtype="uint8")
    reproject(forest.astype("uint8"), samples, src_transform=source_transform, src_crs=MASK_CRS,
              dst_transform=Affine(SAMPLE, 0, x, 0, -SAMPLE, y + CHUNK), dst_crs=CRS, resampling=Resampling.nearest)
    k = CELL // SAMPLE
    return samples.reshape(n // k, k, n // k, k).mean(axis=(1, 3))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders")
    args = parser.parse_args()
    out_dir = args.inputs_dir / "spw_forest_mask_2018"
    out_dir.mkdir(exist_ok=True)
    wallonia = gpd.read_file(args.inputs_dir / "wallonia_2013_2014" / "SPW_Region_Wallonia.gpkg").union_all()

    # chunks on the 30 m grid (cells on multiples of 30 m), covering Wallonia
    xmin, ymin = (math.floor(v / CELL) * CELL for v in wallonia.bounds[:2])
    nx, ny = (math.ceil((wallonia.bounds[i + 2] - v) / CHUNK) for i, v in enumerate((xmin, ymin)))
    cells = CHUNK // CELL
    share = np.full((ny * cells, nx * cells), np.nan, dtype="float32")
    chunks = [(i, j) for i in range(nx) for j in range(ny)
              if wallonia.intersects(box(xmin + i * CHUNK, ymin + j * CHUNK, xmin + (i + 1) * CHUNK,
                                         ymin + (j + 1) * CHUNK))]
    print(f"{len(chunks)} chunks of {CHUNK / 1000:.2f} km", flush=True)
    start = time.time()
    for done, (i, j) in enumerate(chunks, start=1):
        row = (ny - 1 - j) * cells   # rows count from the top
        share[row:row + cells, i * cells:(i + 1) * cells] = chunk_share(xmin + i * CHUNK, ymin + j * CHUNK)
        if done % 10 == 0 or done == len(chunks):
            print(f"chunk {done}/{len(chunks)}; elapsed {(time.time() - start) / 60:.1f} min", flush=True)

    transform = Affine(CELL, 0, xmin, 0, -CELL, ymin + ny * CHUNK)
    write_raster(out_dir / "forest_mask_2018_share_30m.tif", [share], {"crs": CRS, "transform": transform},
                 names=["forest_share"],
                 tags={"source": "Masque forestier 2018 (Gembloux Agro-Bio Tech / ULiege, SPW), " + MASK_URL,
                       "method": f"rendered in EPSG:3812, nearest neighbour to {SAMPLE} m EPSG:31370, "
                                 f"mean per {CELL} m cell",
                       "made_by": "belgium-lidar-forest-metrics/scripts/download_forest_mask.py"})
    (out_dir / "download_info.json").write_text(json.dumps({
        "date": date.today().isoformat(), "service": MASK_URL, "layer": LAYER, "crs": CRS, "sample_m": SAMPLE,
        "script": "belgium-lidar-forest-metrics/scripts/download_forest_mask.py"}, indent=2))
    print("written: spw_forest_mask_2018/forest_mask_2018_share_30m.tif", flush=True)


if __name__ == "__main__":
    main()
