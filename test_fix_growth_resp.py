import unittest

import numpy as np

import model3d


class GrowthRespirationTests(unittest.TestCase):
    def setUp(self):
        self.bio = {
            "r_max": 1.0,
            "Ks": 0.2,
            "Y": 0.4,
            "m0": 0.02,
            "m_s": 0.1,
            "k_dis": 0.1,
        }
        self.dt = 0.05
        self.fields = {
            "D": np.zeros((1, 1, 1)),
            "K": np.ones((1, 1, 1)),
            "B0": np.full((1, 1, 1), 0.1),
            "OM_a": np.ones((1, 1, 1)),
            "pore": np.ones((1, 1, 1), dtype=bool),
            "matrix": np.zeros((1, 1, 1), dtype=bool),
            "hab_wet": np.ones((1, 1, 1), dtype=bool),
            "hab_dry": np.zeros((1, 1, 1), dtype=bool),
            "k_ex": 0.0,
        }

    def _run(self, legacy_growth_resp=False):
        return model3d.simulate3d(
            self.fields,
            self.bio,
            2,
            self.dt,
            1.0,
            np,
            solver="jacobi",
            iterations=1,
            legacy_growth_resp=legacy_growth_resp,
        )

    def test_one_voxel_respiration_uses_substrate_uptake(self):
        out = self._run()
        index = 1
        substrate = out["source_dt"][0]
        biomass = out["B"][index]
        growth = (
            self.bio["r_max"]
            * (substrate / (substrate + self.bio["Ks"]))
            * biomass
            * (1.0 - biomass / self.fields["K"].item())
        )
        expected = (
            (1.0 - self.bio["Y"]) / self.bio["Y"] * growth
            + self.bio["m0"] * biomass
        )
        self.assertAlmostEqual(out["R"][index], expected, places=14)

    def test_legacy_flag_reproduces_old_formula_exactly(self):
        out = self._run(legacy_growth_resp=True)
        index = 1
        substrate = out["source_dt"][0]
        biomass = out["B"][index]
        growth = (
            self.bio["r_max"]
            * (substrate / (substrate + self.bio["Ks"]))
            * biomass
            * (1.0 - biomass / self.fields["K"].item())
        )
        expected = (1.0 - self.bio["Y"]) * growth + self.bio["m0"] * biomass
        self.assertEqual(out["R"][index], expected)

    def test_growth_closure_passes_corrected_and_fails_legacy(self):
        corrected = self._run()
        legacy = self._run(legacy_growth_resp=True)
        self.assertLessEqual(corrected["closure_growth_rel"], 1e-10)
        self.assertGreater(legacy["closure_growth_rel"], 1e-10)


if __name__ == "__main__":
    unittest.main()
