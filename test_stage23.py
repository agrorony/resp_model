import unittest

import numpy as np

import chem23
import model3d as M
import packing3d as P
import stage23_analysis as analysis
from model import load_yaml
from scan3d_loess import pore_fill_fines_phase


class Stage23Tests(unittest.TestCase):
    def test_fixed_chemistry_derivations(self):
        chem23.validate_chemistry()
        self.assertAlmostEqual(chem23.CHEMISTRY["S"]["OM_density"], 0.5843, places=4)
        self.assertEqual(format(chem23.CHEMISTRY["V"]["B0_density"], ".4g"), "0.01736")
        self.assertAlmostEqual(
            np.mean([chem23.CHEMISTRY[soil]["OM_density"] for soil in ("S", "V", "L")]),
            chem23.MODEL_REFERENCE_OM_DENSITY,
            places=12,
        )
        self.assertAlmostEqual(
            np.mean([chem23.CHEMISTRY[soil]["B0_density"] for soil in ("S", "V", "L")]),
            chem23.MODEL_REFERENCE_B0_DENSITY,
            places=12,
        )

    def test_per_soil_chemistry_passes_run_total_assertions(self):
        phase = np.full((8, 8, 8), P.SOLID, dtype=np.uint8)
        phase[:4, :, :] = P.PORE
        phase[4:6, :, :] = P.MATRIX
        structure = P.Structure()
        structure.kind = "test"
        structure.n = 8
        structure.vox_um = 15.0
        structure.seed = 22
        structure.phase = phase
        structure.lt = P.local_thickness(phase == P.PORE, structure.vox_um, np)
        structure.info = {}
        bio = load_yaml("configs/biology.yaml")
        for code in ("S", "V", "L"):
            hp = M.load_inherited()
            hp["B0_density"] = chem23.CHEMISTRY[code]["B0_density"]
            out, summary = M.run_structure(
                structure,
                hp,
                bio,
                np,
                0.16445,
                chem23.CHEMISTRY[code]["OM_density"],
                T=2,
            )
            self.assertEqual(len(out["R"]), 2)
            self.assertLessEqual(summary["mass_balance_rel"], 1e-10)

    def test_pore_fill_fines_never_change_original_solid(self):
        crop = np.zeros((8, 8, 8), dtype=np.uint8)
        crop[:4, :, :] = 1
        phase, target = pore_fill_fines_phase(crop, 15.0)
        self.assertEqual(int(np.count_nonzero(phase == P.MATRIX)), target)
        self.assertTrue(np.all(phase[crop == 0] == P.SOLID))
        self.assertTrue(np.all(np.isin(phase[crop == 1], (P.PORE, P.MATRIX))))

    def test_stage23_shape_classifier_order(self):
        self.assertEqual(
            analysis.classify_shape_stage23(np.array([1.0, 1.3, 1.3, 1.0])),
            "rise-then-fall",
        )
        self.assertEqual(
            analysis.classify_shape_stage23(np.array([1.0, 1.2, 1.2, 1.1])),
            "rise",
        )
        self.assertEqual(
            analysis.classify_shape_stage23(np.array([1.2, 1.1, 1.0, 1.0])),
            "fall",
        )

    def test_cross_sum_of_squares_partition(self):
        values = {}
        counts = (8, 8, 12)
        for row_index, row in enumerate(("S", "V", "L")):
            for col_index, column in enumerate(("chem_S", "chem_V", "chem_L")):
                values[(row, column)] = [
                    float(row_index * 2 + col_index + ((row_index * col_index) % 2) * 0.2 + i * 0.01)
                    for i in range(counts[row_index])
                ]
        result = analysis._two_way_decomposition(values)
        self.assertAlmostEqual(
            sum(result["sum_of_squares"].values()),
            result["total_sum_of_squares"],
            places=8,
        )
        self.assertAlmostEqual(sum(result["share_of_total"].values()), 1.0, places=10)
if __name__ == "__main__":
    unittest.main()
