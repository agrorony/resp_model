import unittest

import numpy as np

from scan3d import select_fines
from timemap import classify_period_shape, day_to_step, physiological_outputs


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
        self.assertEqual(classify_period_shape(np.array([1.0, 1.2, 1.2, 1.0])), "rise-then-fall")

    def test_primary_mapping_physiology_integrals(self):
        values = physiological_outputs(
            np.full(3000, 2.0),
            np.full(3000, 4.0),
            np.full(3000, 0.1),
            initial_om=3000.0,
            r_lit=2.6,
        )
        self.assertAlmostEqual(values["cumulative_respiration_model_units"], 78.0)
        self.assertAlmostEqual(values["cumulative_respiration_fraction_initial_OM"], 0.026)
        self.assertAlmostEqual(values["qCO2_per_day"], 1.3)
        self.assertAlmostEqual(values["CUE"], 0.0)


if __name__ == "__main__":
    unittest.main()
