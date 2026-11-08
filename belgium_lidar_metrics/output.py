"""Writing results: GeoTIFFs with England band names and tags that record how they were made."""
import numpy as np
import rasterio

from belgium_lidar_metrics.canopy import GAP_FILL_SIZE, MAX_VEGETATION_HEIGHT, MIN_VEGETATION_HEIGHT

# Band names as in the England pipeline (gridAggregate_VOM_HPC_V2.R / forest_stack_30m.tif)
ENGLAND_NAMES = {"cover": "cover_VOM", "p90": "perc09_height_VOM", "sd": "sd_height_VOM", "mean": "mean_height_VOM"}


FLANDERS_SOURCE = "DHMV II DSM - DTM (Digitaal Vlaanderen); buildings GRB GBG + GBA"


def processing_tags(made_by, source=FLANDERS_SOURCE):
    """The canopy height rules and data sources, written into every output file."""
    return {"source": source,
            "gap_fill_closing_m": str(GAP_FILL_SIZE), "min_vegetation_height_m": str(MIN_VEGETATION_HEIGHT),
            "max_vegetation_height_m": str(MAX_VEGETATION_HEIGHT), "made_by": made_by}


def write_raster(path, bands, profile, names=None, tags=None):
    """Tiled, compressed float32 GeoTIFF with NaN as NoData; tags record how it was made."""
    profile = profile | {"driver": "GTiff", "dtype": "float32", "nodata": np.nan, "count": len(bands),
                         "height": bands[0].shape[0], "width": bands[0].shape[1], "tiled": True,
                         "blockxsize": 256, "blockysize": 256, "compress": "deflate", "predictor": 3}
    with rasterio.open(path, "w", **profile) as dst:
        dst.update_tags(**(tags or {}))
        for i, band in enumerate(bands, start=1):
            dst.write(band.astype("float32"), i)
            if names:
                dst.set_band_description(i, names[i - 1])
