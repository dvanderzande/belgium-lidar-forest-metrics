"""Download DHMV II DSM and DTM (1 m) for the Brussels region, plus the UrbIS buildings and region border.

DHMV II also covers Brussels, but download_flanders.py skipped it (only Flanders plus 300 m). This script uses the
same 10 km blocks and file names, keeping only the 1 km pieces that touch Brussels plus MARGIN; the rest of a block is
NoData. Finished files are skipped; a file only gets its final name when it is complete.

    source dev.env
    python scripts/download_brussels.py "$BELGIUM_LIDAR_WORK_DIR/inputs/brussels_dhmv2"

Output: DSM/DHMVII_DSM_1m_<x>_<y>.tif, DTM/..., UrbIS_Buildings.gpkg (Brussels plus MARGIN), UrbIS_Region_Brussels.gpkg,
download_info.json. <x>_<y> is the lower-left corner of the block in metres (EPSG:31370).
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
from belgium_lidar_metrics.download import (CRS, RASTERS, URBIS_WFS, WCS_URL, download_raster,  # noqa: E402
                                            download_vector_by_area, fetch)

BLOCK = 10_000   # metres, as in download_flanders.py
MARGIN = 300     # metres around Brussels, as around Flanders
LICENCE = ("DHMV II: Gratis open data licentie Vlaanderen v1.0, source Digitaal Vlaanderen; "
           "UrbIS (paradigm.brussels): CC0 per the catalogue, to confirm for the buildings")


def brussels_border(out_dir):
    """The Brussels-Capital Region (one polygon, downloaded once)."""
    path = out_dir / "UrbIS_Region_Brussels.gpkg"
    if not path.exists():
        content, _ = fetch(URBIS_WFS, {"service": "WFS", "version": "2.0.0", "request": "GetFeature",
                                       "typeNames": "urbisvector:Region", "srsName": CRS,
                                       "outputFormat": "application/json"})
        gpd.read_file(content.decode("utf-8")).set_crs(CRS, allow_override=True).to_file(path)
    return gpd.read_file(path).union_all()


def blocks_in(zone):
    """Lower-left corners of all BLOCK x BLOCK blocks that touch the zone (the grid of download_flanders.py)."""
    xmin, ymin, xmax, ymax = zone.bounds
    return [(x, y)
            for x in range(int(xmin // BLOCK * BLOCK), int(xmax) + 1, BLOCK)
            for y in range(int(ymin // BLOCK * BLOCK), int(ymax) + 1, BLOCK)
            if zone.intersects(box(x, y, x + BLOCK, y + BLOCK))]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    for sub in ["", *RASTERS]:
        (args.out_dir / sub).mkdir(parents=True, exist_ok=True)

    zone = brussels_border(args.out_dir).buffer(MARGIN)
    blocks = blocks_in(zone)
    print(f"{len(blocks)} blocks: {blocks}", flush=True)

    path = args.out_dir / "UrbIS_Buildings.gpkg"
    if not path.exists():
        start = time.time()
        xmin, ymin, xmax, ymax = (int(v) for v in zone.bounds)
        buildings = download_vector_by_area(URBIS_WFS, "urbisvector:Buildings", (xmin, ymin, xmax + 1, ymax + 1),
                                            id_field="INSPIRE_ID")
        partial = path.with_name("partial_" + path.name)   # keeps the .gpkg extension GDAL expects
        buildings.to_file(partial, driver="GPKG")
        partial.rename(path)
        print(f"buildings: {len(buildings):,} in {(time.time() - start) / 60:.1f} min", flush=True)

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

    (args.out_dir / "download_info.json").write_text(json.dumps({
        "date": date.today().isoformat(), "crs": CRS, "block_m": BLOCK, "margin_m": MARGIN, "wcs": WCS_URL,
        "coverages": RASTERS, "wfs": URBIS_WFS, "layers": ["urbisvector:Buildings", "urbisvector:Region"],
        "licence": LICENCE, "script": "belgium-lidar-forest-metrics/scripts/download_brussels.py"}, indent=2))
    print("all done", flush=True)


if __name__ == "__main__":
    main()
