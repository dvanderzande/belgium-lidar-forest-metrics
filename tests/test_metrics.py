"""Metric functions against hand-computed values, on 2 x 2 pixel cells. Run: python -m unittest discover tests"""
import unittest

import numpy as np

from belgium_lidar_metrics.metrics import cell_metrics, to_cells

NAN = np.nan


def metrics_of(values):
    """The four metrics of one 2 x 2 cell given as 4 values."""
    result = cell_metrics(np.array(values, dtype=float).reshape(2, 2), cell=2)
    return {name: float(m[0, 0]) for name, m in result.items()}


class TestCellMetrics(unittest.TestCase):

    def assertMetrics(self, values, expected):
        got = metrics_of(values)
        for name, value in expected.items():
            if np.isnan(value):
                self.assertTrue(np.isnan(got[name]), f"{name}: expected NaN, got {got[name]}")
            else:
                self.assertAlmostEqual(got[name], value, places=4, msg=name)

    def test_mixed_heights(self):
        # sorted 0, 1, 3, 10. cover: 3 and 10 are > 2 -> 2/4. mean 14/4.
        # sd: squared deviations 12.25 + 6.25 + 0.25 + 42.25 = 61, / 3 -> sqrt(20.333) = 4.5092
        # p90 type 7: h = (4 - 1) * 0.9 = 2.7 -> 3 + 0.7 * (10 - 3) = 7.9
        self.assertMetrics([0, 1, 3, 10], {"cover": 0.5, "mean": 3.5, "sd": 4.5092, "p90": 7.9})

    def test_nodata_is_ignored(self):
        # valid 0, 4, 6. cover 2/3. mean 10/3. sd: (11.111 + 0.444 + 7.111) / 2 = 9.333 -> 3.0551
        # p90: h = 2 * 0.9 = 1.8 -> 4 + 0.8 * (6 - 4) = 5.6
        self.assertMetrics([NAN, 0, 4, 6], {"cover": 2 / 3, "mean": 10 / 3, "sd": 3.0551, "p90": 5.6})

    def test_all_zero(self):
        self.assertMetrics([0, 0, 0, 0], {"cover": 0, "mean": 0, "sd": 0, "p90": 0})

    def test_all_nodata(self):
        self.assertMetrics([NAN] * 4, {"cover": NAN, "mean": NAN, "sd": NAN, "p90": NAN})

    def test_single_valid_pixel(self):
        # R: sd of one value is NA
        self.assertMetrics([5, NAN, NAN, NAN], {"cover": 1, "mean": 5, "sd": NAN, "p90": 5})

    def test_threshold_is_strict(self):
        # exactly 2 m is not cover (R: x > threshold)
        self.assertMetrics([2, 2, 2.01, 0], {"cover": 0.25})

    def test_cells_are_separated(self):
        # two cells side by side must not mix: left all 0, right all 10
        heights = np.array([[0, 0, 10, 10], [0, 0, 10, 10]], dtype=float)
        result = cell_metrics(heights, cell=2)
        np.testing.assert_allclose(result["mean"], [[0, 10]])
        np.testing.assert_allclose(result["cover"], [[0, 1]])

    def test_to_cells_layout(self):
        # pixel (row 1, col 3) belongs to cell (0, 1)
        heights = np.arange(16, dtype=float).reshape(4, 4)
        self.assertIn(7.0, to_cells(heights, cell=2)[0, 1])
        with self.assertRaises(ValueError):
            to_cells(np.zeros((5, 4)), cell=2)


if __name__ == "__main__":
    unittest.main()
