"""A 1 m canopy height model encoded like the Environment Agency VOM (docs/methods.md):
vegetation height where it is at least 2.5 m, 0 for everything else (ground, low vegetation, buildings, water, tall
structures), NaN only where there are no data (or no building data, e.g. outside Flanders)."""
import numpy as np
from rasterio import features
from scipy import ndimage

MIN_VEGETATION_HEIGHT = 2.5   # metres: the EA VOM has no values between 0 and 2.5 m
# Project choice (2026-11-07): fill canopy gaps with a 3 x 3 m closing. Leaf-off DHMV II (March 2014) lets the laser
# reach the ground through bare crowns; without filling, cover in tall forest is 0.60 against 0.90 in the EA VOM.
# A deviation from England, see docs/methods.md.
GAP_FILL_SIZE = 3
# Project choice (2026-11-07): heights above 60 m are not vegetation (cranes, wind turbines, masts, chimneys; the DSM
# reaches 193 m). The tallest Belgian trees are about 50 m. The EA VOM holds vegetation only, so England needs no such
# rule; see docs/methods.md.
MAX_VEGETATION_HEIGHT = 60.0


def polygon_mask(geometries, shape, transform):
    """True for pixels whose centre lies in any of the polygons (rasterio default, as terra::rasterize)."""
    geometries = list(geometries)
    if not geometries:
        return np.zeros(shape, dtype=bool)
    return features.geometry_mask(geometries, out_shape=shape, transform=transform, invert=True)


def fill_canopy_gaps(height, size):
    """Grey-scale closing (max then min filter, size x size pixels): fills holes narrower than the window, such as
    laser hits on the ground through leaf-off crowns, without widening crowns or forest edges. NaN is not allowed."""
    return ndimage.grey_closing(height, size=(size, size), mode="nearest")


def canopy_height(dsm, dtm, buildings, inside, gap_fill_size=None):
    """Canopy height (float32) from DSM and DTM (NaN = NoData), a building mask and a mask of the area with data.

    buildings:     True where a building footprint is -> 0
    inside:        True where heights and building data both exist; elsewhere NaN
    gap_fill_size: None = EA VOM encoding only; an odd number of pixels = fill gaps first with a grey-scale closing
                   (the project uses GAP_FILL_SIZE; a deviation from England, see docs/methods.md)
    """
    height = dsm.astype("float64") - dtm.astype("float64")
    height[height > MAX_VEGETATION_HEIGHT] = 0.0     # first, so gap filling cannot spread a crane's height
    if gap_fill_size:
        no_data = np.isnan(height)
        height[no_data | buildings] = 0.0    # buildings at 0 first, so their height does not spill into gaps
        height = fill_canopy_gaps(height, gap_fill_size)
        height[no_data] = np.nan
    height[height < MIN_VEGETATION_HEIGHT] = 0.0      # includes the small negative values (DSM a few cm below DTM)
    height[buildings] = 0.0
    height[~inside | np.isnan(dsm) | np.isnan(dtm)] = np.nan
    return height.astype("float32")
