"""The four England metrics per output cell, from a 1 m canopy height raster.

Definitions as in Ant_ENM_England_2025 `source/misc_functions.R` (see docs/methods.md). NoData is NaN and is
ignored (R `na.rm = TRUE`); zeros (no vegetation) count like any other height.
- cover: share of non-NoData pixels with height > 2 m (strictly greater, as `x > threshold` in R)
- p90:   90th percentile, R `quantile` type 7 = numpy 'linear'
- sd:    standard deviation with n - 1 (R `sd`); NaN for a cell with fewer than 2 valid pixels, as R gives NA
- mean:  mean height
A cell without any valid pixel is NaN in all four.
"""
import warnings

import numpy as np

COVER_THRESHOLD = 2.0   # metres
CELL = 30               # pixels per output cell side (1 m pixels -> 30 m cells)
METRIC_NAMES = ("cover", "p90", "sd", "mean")   # band order of the output


def to_cells(heights, cell=CELL):
    """(rows * cell, cols * cell) array -> (rows, cols, cell * cell): the pixels of each output cell."""
    n_rows, n_cols = heights.shape[0] // cell, heights.shape[1] // cell
    if heights.shape != (n_rows * cell, n_cols * cell):
        raise ValueError(f"shape {heights.shape} is not a multiple of the cell size {cell}")
    return heights.reshape(n_rows, cell, n_cols, cell).swapaxes(1, 2).reshape(n_rows, n_cols, cell * cell)


def cell_metrics(heights, cell=CELL):
    """Dict of the four metrics, each (rows, cols), from heights (rows * cell, cols * cell) with NaN as NoData."""
    pixels = to_cells(np.asarray(heights, dtype="float64"), cell)
    valid = ~np.isnan(pixels)
    n_valid = valid.sum(axis=2)
    with warnings.catch_warnings():
        # numpy warns on all-NaN cells and on sd with one value; both correctly give NaN
        warnings.simplefilter("ignore", RuntimeWarning)
        cover = np.where(n_valid > 0, (pixels > COVER_THRESHOLD).sum(axis=2) / np.maximum(n_valid, 1), np.nan)
        p90 = np.nanpercentile(pixels, 90, axis=2, method="linear")
        sd = np.nanstd(pixels, axis=2, ddof=1)
        mean = np.nanmean(pixels, axis=2)
    sd[n_valid < 2] = np.nan
    return {"cover": cover, "p90": p90, "sd": sd, "mean": mean}
