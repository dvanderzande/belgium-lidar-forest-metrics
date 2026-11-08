"""The four 30 m metrics for Flanders, Brussels and Wallonia (any of them, or all), from the downloaded inputs.

Each region has its own heights, buildings and border (REGIONS). The work goes per TILE x TILE m tile of the 30 m
grid. A tile reads the 1 m heights with a HALO, so gap filling and cells across download blocks come out as in one
large raster. Where a tile crosses a regional border, each region's canopy height is made from its own data and
buildings inside its own border; the 1 m heights are joined before the 30 m cells are computed, so border cells
are complete. A pixel inside two borders (they differ by a few metres) goes to the region listed first. Finished tiles
are skipped (the run is resumable); at the end all tiles are merged into one GeoTIFF.

    source dev.env
    python scripts/process_region.py "$BELGIUM_LIDAR_WORK_DIR/inputs" "$BELGIUM_LIDAR_WORK_DIR/belgium"     # all
    python scripts/process_region.py <inputs> <out> --region brussels                                       # one
    python scripts/process_region.py <inputs> <out> --region flanders --tile 150000 210000                  # a test
    python scripts/process_region.py <inputs> <out> --region wallonia --part 1 6     # 1 of 6 runs side by side

Output: tiles/metrics_<x>_<y>.tif (x, y = lower-left corner) and metrics_30m_<regions>.tif, 4 bands named as the
England layers, NaN outside the regions.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.merge import merge
from rasterio.transform import Affine
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.canopy import GAP_FILL_SIZE, canopy_height, polygon_mask   # noqa: E402
from belgium_lidar_metrics.metrics import CELL, METRIC_NAMES, cell_metrics            # noqa: E402
from belgium_lidar_metrics.output import ENGLAND_NAMES, processing_tags, write_raster  # noqa: E402

# 30 m cells start at multiples of CELL from this point in Lambert 72 (EPSG:31370). Provisional: change it here if
# the species distribution model needs another grid.
GRID_ORIGIN = (0, 0)
TILE = 6000   # metres per processing tile (200 x 200 cells, about 2 GB of memory); a multiple of CELL
HALO = 10     # metres read around a tile: gap filling with a 3 x 3 closing looks 2 pixels beyond the tile
MADE_BY = "belgium-lidar-forest-metrics/scripts/process_region.py"

# Per region: input folder (under the inputs folder), border file, and the data sources (see the readers below).
# The order matters only where borders overlap: the region listed first gets the pixel.
REGIONS = {
    "flanders": {"folder": "flanders_dhmv2", "border": "VRBG_Refgew_Flanders.gpkg", "heights": "dhmv",
                 "buildings": "grb", "source": "DHMV II DSM - DTM (Digitaal Vlaanderen); buildings GRB GBG + GBA"},
    "brussels": {"folder": "brussels_dhmv2", "border": "UrbIS_Region_Brussels.gpkg", "heights": "dhmv",
                 "buildings": "urbis", "source": "DHMV II DSM - DTM (Digitaal Vlaanderen); buildings UrbIS"},
    "wallonia": {"folder": "wallonia_2013_2014", "border": "SPW_Region_Wallonia.gpkg", "heights": "wallonia",
                 "buildings": "picc", "source": "MNS - MNT 2013-2014 (SPW); buildings PICC (SPW)"},
}


def dhmv_blocks(folder):
    """DHMV II download blocks as {(x, y): (block box, DSM path, DTM path)}, from the file names."""
    size = json.loads((folder / "download_info.json").read_text())["block_m"]
    blocks = {}
    for dsm in (folder / "DSM").glob("DHMVII_DSM_1m_*.tif"):
        x, y = map(int, dsm.stem.split("_")[-2:])
        blocks[(x, y)] = (box(x, y, x + size, y + size), dsm, folder / "DTM" / dsm.name.replace("DSM", "DTM"))
    return blocks


def walloon_files(folder):
    """Walloon province GeoTIFFs as {name: (bounds box, MNS path, MNT path)}, without the macOS '._' files."""
    files = {}
    for mns in (folder / "MNS").glob("*/*.tif"):
        if mns.name.startswith("._"):
            continue
        mnt = folder / "MNT" / mns.parent.name / mns.name.replace("MNS", "MNT")
        with rasterio.open(mns) as src:
            files[mns.parent.name] = (box(*src.bounds), mns, mnt)
    return files


def read_heights(paths, bounds):
    """One 1 m raster over bounds from several files, float64 with NaN for NoData (also beyond all files)."""
    sources = [rasterio.open(p) for p in paths]
    nodata = sources[0].nodata
    mosaic, transform = merge(sources, bounds=bounds, res=1, nodata=nodata)
    for s in sources:
        s.close()
    data = mosaic[0].astype("float64")
    data[data == nodata] = np.nan
    return data, transform


def read_buildings(region, folder, outer, touching):
    """Building footprints touching the outer box."""
    if region["buildings"] == "grb":
        # each block's GRB file covers the block plus 300 m, so the touching blocks together cover the tile and halo
        return pd.concat([gpd.read_file(folder / "GRB" / f"GRB_{x}_{y}.gpkg", layer=layer, bbox=outer)
                          for x, y in touching for layer in ("GBG", "GBA")]).geometry
    if region["buildings"] == "urbis":
        return gpd.read_file(folder / "UrbIS_Buildings.gpkg", bbox=outer).geometry
    if region["buildings"] == "picc":
        return pd.concat([gpd.read_file(p, bbox=outer) for p in sorted((folder / "PICC").glob("PICC_buildings_*.gpkg"))
                          if not p.name.startswith("._")]).geometry
    raise ValueError(region["buildings"])


def region_height(region, data, outer, inside):
    """Canopy height over the outer box for one region, NaN outside `inside`; None if no input file touches it."""
    folder, files = data
    touching = [k for k, f in files.items() if f[0].intersects(box(*outer))]
    if not touching:
        return None
    dsm, transform = read_heights([files[k][1] for k in touching], outer)
    dtm, _ = read_heights([files[k][2] for k in touching], outer)
    buildings = polygon_mask(read_buildings(region, folder, outer, touching), dsm.shape, transform)
    return canopy_height(dsm, dtm, buildings, inside, gap_fill_size=GAP_FILL_SIZE)


def process_tile(x, y, regions, out_path):
    """Canopy height of the tile plus HALO (joined over the regions), then the metrics of the tile's cells."""
    outer = (x - HALO, y - HALO, x + TILE + HALO, y + TILE + HALO)
    shape = (TILE + 2 * HALO, TILE + 2 * HALO)
    transform = Affine(1, 0, outer[0], 0, -1, outer[3])   # 1 m pixels, as read_heights returns them
    height = np.full(shape, np.nan, dtype="float32")
    taken = np.zeros(shape, dtype=bool)                   # pixels already given to an earlier region
    for region, border, data in regions.values():
        if not border.intersects(box(*outer)):
            continue
        inside = polygon_mask([border.intersection(box(*outer))], shape, transform) & ~taken
        taken |= inside
        part = region_height(region, data, outer, inside)
        if part is not None:
            height = np.where(inside, part, height)
    metrics = cell_metrics(height[HALO:-HALO, HALO:-HALO])

    profile = {"crs": "EPSG:31370", "transform": Affine(CELL, 0, x, 0, -CELL, y + TILE)}
    partial = out_path.with_name("partial_" + out_path.name)   # final name only when complete
    write_raster(partial, [metrics[n] for n in METRIC_NAMES], profile, names=[ENGLAND_NAMES[n] for n in METRIC_NAMES],
                 tags=tags_for(regions))
    partial.rename(out_path)


def tags_for(regions):
    return processing_tags(MADE_BY, source=" | ".join(f"{name}: {r[0]['source']}" for name, r in regions.items()))


def tiles_in(area):
    """Lower-left corners of the tiles that touch the area."""
    x0, y0 = GRID_ORIGIN
    xmin, ymin, xmax, ymax = area.bounds
    return [(x, y)
            for x in range(int(x0 + (xmin - x0) // TILE * TILE), int(xmax) + 1, TILE)
            for y in range(int(y0 + (ymin - y0) // TILE * TILE), int(ymax) + 1, TILE)
            if area.intersects(box(x, y, x + TILE, y + TILE))]


def merge_tiles(tile_paths, out_path, regions):
    """All tiles in one GeoTIFF (NaN where no tile)."""
    sources = [rasterio.open(p) for p in tile_paths]
    mosaic, transform = merge(sources, nodata=np.nan)
    for s in sources:
        s.close()
    write_raster(out_path, list(mosaic), {"crs": "EPSG:31370", "transform": transform},
                 names=[ENGLAND_NAMES[n] for n in METRIC_NAMES],
                 tags=tags_for(regions) | {"grid_origin": str(GRID_ORIGIN), "cell_m": str(CELL)})


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs_dir", type=Path, help="folder holding the region input folders")
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--region", choices=REGIONS, action="append", help="repeatable; default: all three")
    parser.add_argument("--tile", type=int, nargs=2, action="append", metavar=("X", "Y"),
                        help="lower-left corner of one tile to process (repeatable); default: all, then merge")
    parser.add_argument("--part", type=int, nargs=2, metavar=("K", "N"),
                        help="only every N-th tile, starting at the K-th (K = 1..N), and no merge: run N of these "
                             "side by side (about 3 GB of memory each), then once without --part to merge")
    args = parser.parse_args()
    (args.out_dir / "tiles").mkdir(parents=True, exist_ok=True)

    names = [n for n in REGIONS if n in (args.region or REGIONS)]    # always in the order of REGIONS
    regions = {}
    for name in names:
        region = REGIONS[name]
        folder = args.inputs_dir / region["folder"]
        border = gpd.read_file(folder / region["border"]).union_all()
        files = dhmv_blocks(folder) if region["heights"] == "dhmv" else walloon_files(folder)
        regions[name] = (region, border, (folder, files))
    area = gpd.GeoSeries([r[1] for r in regions.values()]).union_all()
    all_tiles = tiles_in(area)
    tiles = [tuple(t) for t in args.tile] if args.tile else all_tiles
    unknown = set(tiles) - set(all_tiles)
    if unknown:
        sys.exit(f"not a tile of the grid touching {', '.join(names)}: {sorted(unknown)}")
    print(f"{', '.join(names)}: {len(tiles)} of {len(all_tiles)} tiles of {TILE} m", flush=True)

    if args.part:   # every n-th tile, so several runs can work side by side (each tile has its own file)
        k, n = args.part
        tiles = tiles[k - 1::n]
    todo = [(x, y) for x, y in tiles if not (args.out_dir / "tiles" / f"metrics_{x}_{y}.tif").exists()]
    print(f"{len(todo)} tiles to do", flush=True)
    start = time.time()
    for done, (x, y) in enumerate(todo, start=1):
        process_tile(x, y, regions, args.out_dir / "tiles" / f"metrics_{x}_{y}.tif")
        elapsed = time.time() - start
        print(f"tile {done}/{len(todo)} ({x}, {y}) done; elapsed {elapsed / 60:.1f} min, remaining about "
              f"{elapsed / done * (len(todo) - done) / 60:.0f} min", flush=True)

    if not args.tile and not args.part:
        out = args.out_dir / f"metrics_30m_{'_'.join(names) if len(names) < len(REGIONS) else 'belgium'}.tif"
        merge_tiles([args.out_dir / "tiles" / f"metrics_{x}_{y}.tif" for x, y in tiles], out, regions)
        print(f"written: {out.name}", flush=True)


if __name__ == "__main__":
    main()
