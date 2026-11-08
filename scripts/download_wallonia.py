"""Download the Walloon LiDAR MNS and MNT 2013-2014 (1 m, EPSG:31370) per province, the PICC buildings per province
and the border of the Walloon Region.

The Walloon server offers no piece-wise raster download, only one ZIP of GeoTIFFs per province and model. A ZIP is
unpacked into <out>/<MNS|MNT>/<province>/ and kept (the user may delete it). The PICC comes as one ZIP of shapefiles
per province; only the building footprints (CONSTR_BATIEMPRISE) are kept, as PICC/PICC_buildings_<province>.gpkg.
Finished files are skipped.

The server does not say how large a ZIP is, so the download stops as soon as all downloads together would pass
--max-gb (default MAX_GB). Every finished download is listed in <out>/downloaded.csv, which the limit counts, so
deleting ZIPs after unpacking does not reset it. The LiDAR ZIPs barely compress: the unpacked files take about as much
space again.

    source dev.env
    python scripts/download_wallonia.py "$BELGIUM_LIDAR_WORK_DIR/inputs/wallonia_2013_2014" --province BRABANT_WALLON

Sources (Métawal catalogue): MNS dataset 7d23d8ab-962a-493f-8771-2054e06ad36f, MNT 6029e738-f828-438b-b10a-85e67f77af92,
PICC b795de68-726c-4bdf-a62a-a42686aa5b6f; border: ArcGIS LIMITES/LIMITES_ADMINISTRATIVES layer 0 (Région).
Terms of use: Conditions_Acces.pdf inside each ZIP.
"""
import argparse
import csv
import json
import shutil
import sys
import tempfile
import time
import zipfile
from datetime import date
from pathlib import Path

import geopandas as gpd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from anywhere without installing the package
from belgium_lidar_metrics.download import CRS, SizeLimitReached, download_file, fetch  # noqa: E402

BASE = "https://geoservices.wallonie.be/geotraitement/spwdatadownload/results"
DATASETS = {"MNS": "7d23d8ab-962a-493f-8771-2054e06ad36f", "MNT": "6029e738-f828-438b-b10a-85e67f77af92"}
PICC = "b795de68-726c-4bdf-a62a-a42686aa5b6f"
REGION = "https://geoservices.wallonie.be/arcgis/rest/services/LIMITES/LIMITES_ADMINISTRATIVES/MapServer/0/query"
PROVINCES = ["BRABANT_WALLON", "HAINAUT", "LIEGE", "LUXEMBOURG", "NAMUR"]
MAX_GB = 120   # limit for all ZIPs together, set by the user 2026-11-07 (LiDAR 112.9 GB + PICC about 5 GB)


def lidar_url(model, province):
    return f"{BASE}/{DATASETS[model]}/RELIEF_WALLONIE_{model}_2013_2014_GEOTIFF_31370_PROV_{province}.zip"


def picc_url(province):
    return f"{BASE}/{PICC}/PICC_vDIFF_SHAPE_31370_PROV_{province}.zip"


def downloaded_gb(out_dir):
    """Size of all finished downloads so far (downloaded.csv), in GB."""
    path = out_dir / "downloaded.csv"
    if not path.exists():
        return 0.0
    with open(path, newline="") as f:
        return sum(int(row["bytes"]) for row in csv.DictReader(f)) / 1e9


def fetch_zip(url, out_dir, max_gb):
    """Download one ZIP into <out>/zip/ unless it is there; stop the run if all downloads would pass max_gb."""
    zip_path = out_dir / "zip" / url.rsplit("/", 1)[1]
    if not zip_path.exists():
        start = time.time()
        try:
            download_file(url, zip_path, max_bytes=(max_gb - downloaded_gb(out_dir)) * 1e9)
        except SizeLimitReached as stop:
            sys.exit(f"STOPPED: {stop}; all downloads would pass {max_gb:g} GB. Delete the partial_ file; "
                     "rerun with a higher --max-gb to continue.")
        ledger = out_dir / "downloaded.csv"
        new = not ledger.exists()
        with open(ledger, "a", newline="") as f:
            writer = csv.writer(f)
            if new:
                writer.writerow(["file", "bytes", "date"])
            writer.writerow([zip_path.name, zip_path.stat().st_size, date.today().isoformat()])
        print(f"{zip_path.name}: {zip_path.stat().st_size / 1e6:,.0f} MB in {(time.time() - start) / 60:.1f} min; "
              f"all downloads {downloaded_gb(out_dir):.1f} GB", flush=True)
    return zip_path


def unpack(zip_path, folder):
    """Extract the GeoTIFFs (and the terms of use); zipfile refuses paths outside the folder."""
    folder.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if n.lower().endswith((".tif", ".tiff", ".tfw", ".pdf", ".xml"))]
        for name in names:
            z.extract(name, folder)
    return names


def wallonia_border(out_dir):
    """The Walloon Region (one polygon) from the SPW administrative boundaries, downloaded once."""
    path = out_dir / "SPW_Region_Wallonia.gpkg"
    if not path.exists():
        content, _ = fetch(REGION, {"where": "1=1", "outFields": "NAMEFRE,LOCALID", "outSR": 31370, "f": "geojson"})
        # GeoJSON is labelled WGS84 by the reader, but outSR makes the server send Lambert 72 coordinates
        gpd.read_file(content.decode("utf-8")).set_crs(CRS, allow_override=True).to_file(path)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--province", choices=PROVINCES, action="append", help="repeatable; default: all five")
    parser.add_argument("--max-gb", type=float, default=MAX_GB, help=f"limit for all ZIPs (default {MAX_GB})")
    args = parser.parse_args()
    (args.out_dir / "zip").mkdir(parents=True, exist_ok=True)
    (args.out_dir / "PICC").mkdir(exist_ok=True)
    print(f"downloaded so far: {downloaded_gb(args.out_dir):.1f} GB; limit {args.max_gb:g} GB", flush=True)
    wallonia_border(args.out_dir)

    for province in args.province or PROVINCES:
        for model in DATASETS:
            folder = args.out_dir / model / province
            if not (folder / ".unpacked").exists():   # once unpacked, the ZIP is no longer needed
                zip_path = fetch_zip(lidar_url(model, province), args.out_dir, args.max_gb)
                names = unpack(zip_path, folder)
                (folder / ".unpacked").write_text("\n".join(names))
                print(f"  unpacked {len(names)} files to {model}/{province}", flush=True)
        buildings = args.out_dir / "PICC" / f"PICC_buildings_{province}.gpkg"
        if not buildings.exists():
            zip_path = fetch_zip(picc_url(province), args.out_dir, args.max_gb)
            layer = gpd.read_file(f"zip://{zip_path}!CONSTR_BATIEMPRISE.shp", columns=["GEOREF_ID", "NATUR_CODE"])
            partial = buildings.with_name("partial_" + buildings.name)   # keeps the .gpkg extension GDAL expects
            # written on the local disk first: building the spatial index of a large layer directly on the exFAT
            # drive crashes GDAL (Hainaut, 1.6 M buildings, 2026-11-08); the copy is then a plain file copy
            with tempfile.TemporaryDirectory() as tmp:
                local = Path(tmp) / buildings.name
                layer.to_file(local, driver="GPKG")
                shutil.copyfile(local, partial)
            partial.rename(buildings)
            print(f"  {len(layer):,} buildings to PICC/{buildings.name}", flush=True)

    (args.out_dir / "download_info.json").write_text(json.dumps({
        "date": date.today().isoformat(), "crs": "EPSG:31370", "datasets": DATASETS | {"PICC": PICC},
        "base_url": BASE, "region_border": REGION, "provinces": args.province or PROVINCES,
        "licence": "see Conditions_Acces.pdf in each ZIP; source SPW",
        "script": "belgium-lidar-forest-metrics/scripts/download_wallonia.py"}, indent=2))


if __name__ == "__main__":
    main()
