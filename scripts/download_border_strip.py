"""Download DHMV II DSM and DTM (1 m) for the strip of Wallonia along the Flemish border, to compare the two surveys
on the same ground (scripts/compare_border_strip.py).

DHMV II covers Flanders plus a strip of about 5 km of Wallonia. This script takes the part of Wallonia within STRIP m
of Flanders, in the 10 km blocks and 1 km pieces of download_flanders.py. Finished files are skipped; a file only
gets its final name when it is complete.

    source dev.env
    python scripts/download_border_strip.py "$BELGIUM_LIDAR_WORK_DIR/inputs" "$BELGIUM_LIDAR_WORK_DIR/inputs/dhmv2_wallonia_strip"

Output: DSM/DHMVII_DSM_1m_<x>_<y>.tif, DTM/..., download_info.json.
"""
import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.download import CRS, LICENCE, RASTERS, WCS_URL, download_raster  # noqa: E402

BLOCK = 10_000   # metres, as in download_flanders.py
STRIP = 2000     # metres of Wallonia along the Flemish border


def strip_zone(inputs_dir):
    """Wallonia within STRIP m of Flanders."""
    flanders = gpd.read_file(inputs_dir / "flanders_dhmv2" / "VRBG_Refgew_Flanders.gpkg").union_all()
    wallonia = gpd.read_file(inputs_dir / "wallonia_2013_2014" / "SPW_Region_Wallonia.gpkg").union_all()
    return wallonia.intersection(flanders.buffer(STRIP))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders (border files)")
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    for sub in ["", *RASTERS]:
        (args.out_dir / sub).mkdir(parents=True, exist_ok=True)

    zone = strip_zone(args.inputs_dir)
    gpd.GeoSeries([zone], crs=CRS).to_file(args.out_dir / "strip_zone.gpkg")
    xmin, ymin, xmax, ymax = zone.bounds
    blocks = [(x, y) for x in range(int(xmin // BLOCK * BLOCK), int(xmax) + 1, BLOCK)
              for y in range(int(ymin // BLOCK * BLOCK), int(ymax) + 1, BLOCK)
              if zone.intersects(box(x, y, x + BLOCK, y + BLOCK))]
    print(f"strip {zone.area / 1e6:.0f} km2, {len(blocks)} blocks", flush=True)
    # written first: process_region.dhmv_blocks reads the block size from it
    (args.out_dir / "download_info.json").write_text(json.dumps({
        "date": date.today().isoformat(), "crs": CRS, "block_m": BLOCK, "strip_m": STRIP, "wcs": WCS_URL,
        "coverages": RASTERS, "licence": LICENCE,
        "script": "belgium-lidar-forest-metrics/scripts/download_border_strip.py"}, indent=2))

    start = time.time()
    for i, (x, y) in enumerate(blocks, start=1):
        for layer, coverage in RASTERS.items():
            path = args.out_dir / layer / f"DHMVII_{layer}_1m_{x}_{y}.tif"
            if path.exists():
                continue
            pieces, retries = download_raster(coverage, (x, y, x + BLOCK, y + BLOCK), path,
                                              keep_piece=lambda b: zone.intersects(box(*b)))
            print(f"block {i}/{len(blocks)} ({x}, {y}) {layer}: {pieces} pieces, {retries} retries; elapsed "
                  f"{(time.time() - start) / 60:.1f} min", flush=True)

    print("all done", flush=True)


if __name__ == "__main__":
    main()
