"""Download DHMV II DSM and DTM (1 m) and GRB buildings for all of Flanders, in 10 x 10 km blocks.

Blocks follow a 10 km grid in Lambert 72 (EPSG:31370). Only the 1 km pieces that touch Flanders plus MARGIN are
downloaded; the rest of a border block is NoData. The run is resumable: finished files are skipped, and a file only
gets its final name when it is complete. Progress goes to download_log.csv in the output folder.

    source dev.env
    python scripts/download_flanders.py "$BELGIUM_LIDAR_WORK_DIR/inputs/flanders_dhmv2"                  # everything
    python scripts/download_flanders.py <folder> --block 150000 210000 --block 220000 200000           # a subset
    python scripts/download_flanders.py <folder> --list                                                # blocks only

Output: DSM/DHMVII_DSM_1m_<x>_<y>.tif, DTM/..., GRB/GRB_<x>_<y>.gpkg (layers GBG and GBA, block + MARGIN),
VRBG_Refgew_Flanders.gpkg, download_info.json. <x>_<y> is the lower-left corner of the block in metres.
"""
import argparse
import csv
import json
import sys
import time
from datetime import date, datetime
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.download import (CRS, LICENCE, RASTERS, VECTORS, WCS_URL, WFS_URLS,  # noqa: E402
                                            download_raster, download_vector, download_vector_by_area)

BLOCK = 10_000   # metres; one file per block and layer (~200 MB per layer for a full block)
MARGIN = 300     # metres around Flanders, so border cells are complete in any 30 m grid chosen later
MAX_GB = 100     # stop before the next block once the output folder holds more than this (expected total ~62 GB)


def folder_gb(folder):
    """Size of everything in the folder, in GB."""
    return sum(f.stat().st_size for f in Path(folder).rglob("*") if f.is_file()) / 1e9


def flanders_zone(out_dir):
    """Flanders plus MARGIN (the border is downloaded once and kept)."""
    path = out_dir / "VRBG_Refgew_Flanders.gpkg"
    if not path.exists():
        download_vector(WFS_URLS["VRBG"], "VRBG:Refgew").to_file(path)
    return gpd.read_file(path).union_all().buffer(MARGIN)


def blocks_in(zone):
    """Lower-left corners of all BLOCK x BLOCK blocks that touch the zone."""
    xmin, ymin, xmax, ymax = zone.bounds
    return [(x, y)
            for x in range(int(xmin // BLOCK * BLOCK), int(xmax) + 1, BLOCK)
            for y in range(int(ymin // BLOCK * BLOCK), int(ymax) + 1, BLOCK)
            if zone.intersects(box(x, y, x + BLOCK, y + BLOCK))]


def log(out_dir, **row):
    path = out_dir / "download_log.csv"
    new = not path.exists()
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["time", "block", "layer", "pieces", "retries", "seconds", "mb"])
        if new:
            writer.writeheader()
        writer.writerow({"time": datetime.now().isoformat(timespec="seconds")} | row)


def download_block(out_dir, x, y, zone):
    """DSM, DTM and buildings of one block; finished files are skipped."""
    block_box = (x, y, x + BLOCK, y + BLOCK)
    for layer, coverage in RASTERS.items():
        path = out_dir / layer / f"DHMVII_{layer}_1m_{x}_{y}.tif"
        if path.exists():
            continue
        start = time.time()
        pieces, retries = download_raster(coverage, block_box, path, keep_piece=lambda b: zone.intersects(box(*b)))
        mb = path.stat().st_size / 1e6 if path.exists() else 0
        log(out_dir, block=f"{x}_{y}", layer=layer, pieces=pieces, retries=retries,
            seconds=round(time.time() - start, 1), mb=round(mb, 1))
    path = out_dir / "GRB" / f"GRB_{x}_{y}.gpkg"
    if not path.exists():
        start = time.time()
        partial = path.with_name("partial_" + path.name)   # keeps the .gpkg extension GDAL expects
        grb_box = (x - MARGIN, y - MARGIN, x + BLOCK + MARGIN, y + BLOCK + MARGIN)
        for name, (service, layer) in VECTORS.items():
            features = download_vector_by_area(WFS_URLS[service], layer, grb_box)
            features.to_file(partial, layer=name.replace("GRB_", ""), driver="GPKG")
        partial.rename(path)
        log(out_dir, block=f"{x}_{y}", layer="GRB", pieces="", retries="", seconds=round(time.time() - start, 1),
            mb=round(path.stat().st_size / 1e6, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--block", type=int, nargs=2, action="append", metavar=("X", "Y"),
                        help="lower-left corner of one block to download (repeatable); default: all blocks")
    parser.add_argument("--list", action="store_true", help="only print the blocks and the number of pieces")
    parser.add_argument("--max-gb", type=float, default=MAX_GB, help=f"size limit of the output folder (default {MAX_GB})")
    args = parser.parse_args()
    for sub in ["", *RASTERS, "GRB"]:
        (args.out_dir / sub).mkdir(parents=True, exist_ok=True)

    zone = flanders_zone(args.out_dir)
    all_blocks = blocks_in(zone)
    blocks = [tuple(b) for b in args.block] if args.block else all_blocks
    unknown = set(blocks) - set(all_blocks)
    if unknown:
        sys.exit(f"not a block of the 10 km grid touching Flanders: {sorted(unknown)}")
    print(f"{len(blocks)} of {len(all_blocks)} blocks", flush=True)
    if args.list:
        print(blocks)
        return

    (args.out_dir / "download_info.json").write_text(json.dumps({
        "date": date.today().isoformat(), "crs": CRS, "block_m": BLOCK, "margin_m": MARGIN, "wcs": WCS_URL,
        "coverages": RASTERS, "wfs": WFS_URLS, "layers": {k: v[1] for k, v in VECTORS.items()},
        "licence": LICENCE, "script": "belgium-lidar-forest-metrics/scripts/download_flanders.py"}, indent=2))
    start = time.time()
    for i, (x, y) in enumerate(blocks, start=1):
        size = folder_gb(args.out_dir)
        if size > args.max_gb:
            sys.exit(f"STOPPED before block {i}/{len(blocks)} ({x}, {y}): output folder holds {size:.1f} GB, "
                     f"more than the limit of {args.max_gb:g} GB. Rerun with a higher --max-gb to continue.")
        download_block(args.out_dir, x, y, zone)
        elapsed = time.time() - start
        print(f"block {i}/{len(blocks)} ({x}, {y}) done; {folder_gb(args.out_dir):.1f} GB on disk; elapsed "
              f"{elapsed / 60:.1f} min, remaining about {elapsed / i * (len(blocks) - i) / 60:.0f} min", flush=True)
    print("all blocks done", flush=True)


if __name__ == "__main__":
    main()
