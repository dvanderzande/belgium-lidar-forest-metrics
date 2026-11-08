"""Canopy height rules (EA VOM encoding). Run: python -m unittest discover tests"""
import unittest

import numpy as np
from rasterio.transform import from_origin
from shapely.geometry import box

from belgium_lidar_metrics.canopy import MAX_VEGETATION_HEIGHT, canopy_height, polygon_mask

NAN = np.nan


class TestCanopyHeight(unittest.TestCase):

    def test_rules(self):
        dtm = np.full(8, 100.0)
        dsm = 100.0 + np.array([-0.05, 0.0, 2.49, 2.5, 30.0, 30.0, 30.0, 30.0])
        dsm[7] = NAN                                                   # no DSM data
        buildings = np.array([0, 0, 0, 0, 0, 1, 0, 0], dtype=bool)     # pixel 5 is a building
        inside = np.array([1, 1, 1, 1, 1, 1, 0, 1], dtype=bool)        # pixel 6 outside the area with data
        got = canopy_height(dsm, dtm, buildings, inside)
        expected = [0.0, 0.0, 0.0, 2.5, 30.0, 0.0, NAN, NAN]
        np.testing.assert_allclose(got, expected, atol=1e-4)
        self.assertEqual(got.dtype, np.float32)

    def test_tall_structures_are_not_vegetation(self):
        # exactly 60 m is kept, above is 0; with gap filling a crane must not fill the gap next to it
        self.assertEqual(MAX_VEGETATION_HEIGHT, 60.0)
        dtm = np.zeros((1, 3))
        dsm = np.array([[60.0, 60.01, 193.0]])
        got = canopy_height(dsm, dtm, np.zeros((1, 3), bool), np.ones((1, 3), bool))
        np.testing.assert_array_equal(got, [[60.0, 0.0, 0.0]])
        dsm = np.zeros((5, 5))
        dsm[2, 1], dsm[2, 3] = 150.0, 150.0          # two crane pixels with a 1-pixel gap between them
        got = canopy_height(dsm, np.zeros((5, 5)), np.zeros((5, 5), bool), np.ones((5, 5), bool), gap_fill_size=3)
        np.testing.assert_array_equal(got, 0.0)

    def test_gap_filling(self):
        # 9 x 9 forest of 20 m with a 1-pixel hole (filled by 3 x 3) and a 3 x 3 clearing (kept); field to the right
        dtm = np.zeros((9, 12))
        dsm = np.zeros((9, 12))
        dsm[:, :9] = 20.0
        dsm[1, 1] = 0.0                  # single ground hit inside the crown
        dsm[4:7, 4:7] = 0.0              # real 3 x 3 m clearing
        everywhere, no_buildings = np.ones((9, 12), bool), np.zeros((9, 12), bool)
        got = canopy_height(dsm, dtm, no_buildings, everywhere, gap_fill_size=3)
        self.assertEqual(got[1, 1], 20.0)                    # hole filled
        self.assertEqual(got[5, 5], 0.0)                     # centre of the clearing stays open
        np.testing.assert_array_equal(got[:, 9:], 0.0)       # forest edge does not grow into the field
        np.testing.assert_array_equal(got[:, :9][got[:, :9] > 0], 20.0)   # no new heights invented

    def test_gap_filling_keeps_buildings_and_nodata(self):
        dtm = np.zeros((5, 5))
        dsm = np.full((5, 5), 15.0)
        dsm[0, 0] = np.nan
        buildings = np.zeros((5, 5), bool)
        buildings[2, 2] = True           # a small shed inside trees: filled by closing, then masked again
        got = canopy_height(dsm, dtm, buildings, np.ones((5, 5), bool), gap_fill_size=3)
        self.assertEqual(got[2, 2], 0.0)
        self.assertTrue(np.isnan(got[0, 0]))

    def test_polygon_mask_uses_pixel_centres(self):
        # 4 x 4 pixels of 1 m from (0, 4); the box covers x 0-2.4, y 2-4: centres x 0.5, 1.5 in, 2.5 out
        mask = polygon_mask([box(0, 2, 2.4, 4)], (4, 4), from_origin(0, 4, 1, 1))
        expected = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=bool)
        np.testing.assert_array_equal(mask, expected)

    def test_polygon_mask_empty(self):
        self.assertFalse(polygon_mask([], (2, 2), from_origin(0, 2, 1, 1)).any())


if __name__ == "__main__":
    unittest.main()
