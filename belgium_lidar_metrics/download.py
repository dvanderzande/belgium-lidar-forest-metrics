"""Download Flemish open data from the Digitaal Vlaanderen web services (Gratis open data licentie Vlaanderen):
- DSM, DTM: WCS https://geo.api.vlaanderen.be/DHMV/wcs (coverages DHMVII_DSM_1m, DHMVII_DTM_1m), EPSG:31370
- buildings: WFS https://geo.api.vlaanderen.be/GRB/wfs, layers GBG (building at ground level) and GBA (attachments)
- Flanders border: WFS https://geo.api.vlaanderen.be/VRBG/wfs, layer Refgew
Brussels (paradigm.brussels, UrbIS):
- buildings and the region border: WFS URBIS_WFS, layers urbisvector:Buildings and urbisvector:Region
Wallonia (Service public de Wallonie):
- buildings: PICC, ArcGIS REST layer "Construction - Bâtiment - Emprise" (PICC_BUILDINGS)
"""
import tempfile
import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio
import requests  # not urllib: it fails on some WCS responses that requests and curl read fine (2026-11-07)
from rasterio.merge import merge

CRS = "EPSG:31370"
PIECE = 1000      # metres per WCS request: the server refuses requests of about 5 km and more
PAGE = 5000       # features per WFS request when paging (small layers only)
WFS_CAP = 10_000  # the server returns at most this many features per request
TILE = 2000       # metres: building queries per area piece; sorted deep paging is slow (30 s/page at 40 000)
WCS_URL = "https://geo.api.vlaanderen.be/DHMV/wcs"
WFS_URLS = {"GRB": "https://geo.api.vlaanderen.be/GRB/wfs", "VRBG": "https://geo.api.vlaanderen.be/VRBG/wfs"}
RASTERS = {"DSM": "DHMVII_DSM_1m", "DTM": "DHMVII_DTM_1m"}
VECTORS = {"GRB_GBG": ("GRB", "GRB:GBG"), "GRB_GBA": ("GRB", "GRB:GBA")}
LICENCE = "Gratis open data licentie Vlaanderen v1.0; source Digitaal Vlaanderen"
URBIS_WFS = "https://geoservices-urbis.irisnet.be/geoserver/urbisvector/wfs"
PICC_BUILDINGS = "https://geoservices.wallonie.be/arcgis/rest/services/TOPOGRAPHIE/PICC_VDIFF/MapServer/11/query"
ARCGIS_CAP = 2000   # the Walloon ArcGIS server returns at most this many features per request
ARCGIS_TILE = 1000  # metres per query piece: about 500 buildings per km2 in a town


def fetch(url, params, retries=6):
    """GET with retries, in case a service times out or cuts a response off. Returns (content, retries used)."""
    for attempt in range(retries):
        try:
            # gzip explicitly: the WCS sends uncompressed GeoTIFFs otherwise (3x more data over the network)
            response = requests.get(url, params=params, timeout=300, headers={"Accept-Encoding": "gzip"})
            response.raise_for_status()
            return response.content, attempt
        except requests.RequestException as error:
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {type(error).__name__}", flush=True)
            time.sleep(5 * (attempt + 1))


class SizeLimitReached(Exception):
    """A download stopped because it would pass its size limit."""


def download_file(url, out_path, retries=3, max_bytes=None):
    """Stream a large file to out_path; the final name is only given when complete. Returns the size in MB.

    max_bytes: stop (SizeLimitReached) before the file grows beyond this; needed where the server does not say how
    large a file is. The incomplete partial_<name> file is left for the user to delete.
    """
    out_path = Path(out_path)
    partial = out_path.with_name("partial_" + out_path.name)
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=300) as response:
                response.raise_for_status()
                with open(partial, "wb") as f:
                    written = 0
                    for chunk in response.iter_content(chunk_size=8 * 2**20):
                        if max_bytes is not None and written + len(chunk) > max_bytes:
                            raise SizeLimitReached(f"{out_path.name}: stopped at {written / 1e9:.2f} GB, the limit "
                                                   f"is {max_bytes / 1e9:.2f} GB")
                        f.write(chunk)
                        written += len(chunk)
            break
        except requests.RequestException as error:
            if attempt == retries - 1:
                raise
            print(f"  retry {attempt + 1} after {type(error).__name__}", flush=True)
            time.sleep(30 * (attempt + 1))
    partial.rename(out_path)
    return out_path.stat().st_size / 1e6


def download_raster(coverage, box, out_path, keep_piece=None):
    """Download the box in PIECE x PIECE m requests and merge them into one tiled, compressed GeoTIFF.

    keep_piece: optional function (xmin, ymin, xmax, ymax) -> bool; pieces it rejects are not downloaded and become
    NoData. The file is written under a temporary name and renamed when complete, so an interrupted run never leaves
    a file that looks finished. Returns (pieces downloaded, retries), or (0, 0) if no piece was kept (no file).
    """
    xmin, ymin, xmax, ymax = box
    retries = 0
    with tempfile.TemporaryDirectory() as tmp:
        pieces = []
        for x in range(xmin, xmax, PIECE):
            for y in range(ymin, ymax, PIECE):
                bbox = (x, y, min(x + PIECE, xmax), min(y + PIECE, ymax))
                if keep_piece and not keep_piece(bbox):
                    continue
                data, used = fetch(WCS_URL, {"service": "WCS", "version": "1.0.0", "request": "GetCoverage",
                                             "coverage": coverage, "crs": CRS, "bbox": ",".join(map(str, bbox)),
                                             "resx": 1, "resy": 1, "format": "GeoTIFF"})
                retries += used
                if not data.startswith((b"II*", b"MM\x00*")):
                    raise RuntimeError(f"{coverage} {bbox}: not a GeoTIFF: {data[:200]!r}")
                piece = Path(tmp) / f"{x}_{y}.tif"
                piece.write_bytes(data)
                with rasterio.open(piece) as src:   # a cut-off file would have the wrong size or fail to read
                    expected = (bbox[3] - bbox[1], bbox[2] - bbox[0])
                    if src.shape != expected or src.bounds != rasterio.coords.BoundingBox(*bbox):
                        raise RuntimeError(f"{coverage} {bbox}: got {src.shape} {src.bounds}")
                    src.read(1)
                pieces.append(piece)
        if not pieces:
            return 0, 0
        sources = [rasterio.open(p) for p in pieces]
        mosaic, transform = merge(sources, bounds=box, res=1)   # skipped pieces are filled with the NoData value
        profile = sources[0].profile | {"driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2],
                                        "transform": transform, "tiled": True, "blockxsize": 512,
                                        "blockysize": 512, "compress": "deflate", "predictor": 3}
        for s in sources:
            s.close()
    partial = Path(out_path).with_name(Path(out_path).name + ".partial")
    with rasterio.open(partial, "w", **profile) as dst:
        dst.write(mosaic)
    partial.rename(out_path)
    return len(pieces), retries


def download_vector_by_area(url, layer, box, tile=TILE, id_field="OIDN"):
    """All features of a WFS layer touching the box, fetched per tile x tile m piece without sorting or paging.
    A piece that hits the server cap is split in four. Features touching several pieces are kept once (id_field:
    OIDN in GRB, INSPIRE_ID in UrbIS)."""
    xmin, ymin, xmax, ymax = box
    frames = []

    def fetch_piece(bbox):
        params = {"service": "WFS", "version": "2.0.0", "request": "GetFeature", "typeNames": layer,
                  "srsName": CRS, "outputFormat": "application/json", "count": WFS_CAP,
                  "bbox": ",".join(map(str, bbox)) + "," + CRS}
        content, _ = fetch(url, params)
        piece = gpd.read_file(content.decode("utf-8"))
        if len(piece) >= WFS_CAP:   # possibly incomplete: split
            x0, y0, x1, y1 = bbox
            xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
            for sub in [(x0, y0, xm, ym), (xm, y0, x1, ym), (x0, ym, xm, y1), (xm, ym, x1, y1)]:
                fetch_piece(sub)
        elif len(piece):
            frames.append(piece.set_crs(CRS, allow_override=True))

    for x in range(xmin, xmax, tile):
        for y in range(ymin, ymax, tile):
            fetch_piece((x, y, min(x + tile, xmax), min(y + tile, ymax)))
    if not frames:
        return gpd.GeoDataFrame(geometry=[], crs=CRS)
    result = pd.concat(frames, ignore_index=True).drop_duplicates(id_field, ignore_index=True)
    return gpd.GeoDataFrame(result, geometry="geometry", crs=CRS)


def download_vector(url, layer, box=None):
    """All features of a WFS layer (in the box, if given), paged. Fine for small layers; for buildings over large
    areas use download_vector_by_area (sorted deep paging gets slow)."""
    frames, start = [], 0
    while True:
        params = {"service": "WFS", "version": "2.0.0", "request": "GetFeature", "typeNames": layer,
                  "srsName": CRS, "outputFormat": "application/json", "count": PAGE, "startIndex": start,
                  "sortBy": "OIDN"}
        if box:
            params["bbox"] = ",".join(map(str, box)) + "," + CRS
        content, _ = fetch(url, params)
        page = gpd.read_file(content.decode("utf-8"))
        frames.append(page)
        if len(page) < PAGE:
            break
        start += PAGE
    frames = [f for f in frames if len(f)] or frames[:1]   # an empty last page has no CRS; drop it
    result = pd.concat(frames, ignore_index=True)
    return gpd.GeoDataFrame(result, geometry="geometry", crs=CRS)


def download_arcgis_by_area(url, box, tile=ARCGIS_TILE):
    """All features of an ArcGIS REST layer touching the box, per tile x tile m piece (in EPSG:31370). A piece that
    hits the server cap is split in four. Features touching several pieces are kept once (OBJECTID)."""
    xmin, ymin, xmax, ymax = box
    frames = []

    def fetch_piece(bbox):
        params = {"geometry": ",".join(map(str, bbox)), "geometryType": "esriGeometryEnvelope", "inSR": 31370,
                  "outSR": 31370, "spatialRel": "esriSpatialRelIntersects", "outFields": "OBJECTID,NATUR_CODE",
                  "where": "1=1", "f": "geojson"}
        content, _ = fetch(url, params)
        piece = gpd.read_file(content.decode("utf-8"))
        if len(piece) >= ARCGIS_CAP:   # possibly incomplete: split
            x0, y0, x1, y1 = bbox
            xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
            for sub in [(x0, y0, xm, ym), (xm, y0, x1, ym), (x0, ym, xm, y1), (xm, ym, x1, y1)]:
                fetch_piece(sub)
        elif len(piece):
            # GeoJSON is labelled WGS84 by the reader, but outSR makes the server send Lambert 72 coordinates
            frames.append(piece.set_crs(CRS, allow_override=True))

    for x in range(xmin, xmax, tile):
        for y in range(ymin, ymax, tile):
            fetch_piece((x, y, min(x + tile, xmax), min(y + tile, ymax)))
    if not frames:
        return gpd.GeoDataFrame(geometry=[], crs=CRS)
    result = pd.concat(frames, ignore_index=True).drop_duplicates("OBJECTID", ignore_index=True)
    return gpd.GeoDataFrame(result, geometry="geometry", crs=CRS)
