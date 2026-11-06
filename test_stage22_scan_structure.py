import unittest

import numpy as np

from scan3d import select_fines
from timemap import classify_period_shape, day_to_step


class ScanStructureTests(unittest.TestCase):
    def test_sand_fines_fraction_and_thinnest_first(self):
        solid = np.ones((2, 3, 4), dtype=bool)
        thickness = np.arange(24, dtype=float).reshape(solid.shape)
        distance = np.ones(solid.shape, dtype=float)
        selected, threshold = select_fines(solid, thickness, distance, 0.25, seed=22)
        self.assertEqual(int(selected.sum()), 6)
        self.assertAlmostEqual(float(selected.mean()), 0.25)
        self.assertEqual(threshold, 5.0)
        self.assertTrue(np.all(thickness[selected] <= threshold))
        self.assertTrue(np.all(thickness[~selected] > threshold))

    def test_registered_time_mapping_and_shape(self):
        np.testing.assert_allclose([day_to_step(15, r) for r in (0.36, 2.6, 9.2)], [108.0, 780.0, 2760.0])
        self.assertEqual(classify_period_shape(np.array([1.0, 1.05, 1.08, 1.09])), "flat")
        self.assertEqual(classify_period_shape(np.array([1.0, 1.0, 1.0, 1.2])), "rise")


if __name__ == "__main__":
    unittest.main()
