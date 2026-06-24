import copy
import json
import unittest
from pathlib import Path

from . import WorkflowError, authored_to_runtime_case


class SpeedClassRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        with (repo_root / "negaccel-authoring.example.json").open("r", encoding="utf-8") as handle:
            cls.example = json.load(handle)

    def _base_case(self) -> dict:
        return copy.deepcopy(self.example)

    def test_explicit_speed_class_survives_materialization(self) -> None:
        case = self._base_case()
        case["gasInteractions"]["reactions"] = [
            {
                "processId": "probe_fast_slow",
                "projectileKind": "H-",
                "projectileFate": "survive",
                "sourcePath": "Cross-sections/background_gas_ionization.cross",
                "fitDegree": 0,
                "coefficients": [1.0],
                "products": [
                    {"particleKind": "H2+", "speedClass": "slow"},
                    {"particleKind": "e-", "count": 2, "speedClass": "fast"},
                ],
            }
        ]

        runtime = authored_to_runtime_case(case)
        products = runtime["physics"]["reactions"][0]["products"]

        self.assertEqual("slow", products[0]["speedClass"])
        self.assertEqual("fast", products[1]["speedClass"])
        self.assertEqual(2, products[1]["count"])

    def test_legacy_reaction_materializes_explicit_speed_classes(self) -> None:
        case = self._base_case()
        case["gasInteractions"]["reactions"] = [
            {
                "reactionId": "positive_ion_charge_exchange",
                "sourcePath": "Cross-sections/positive_ion_charge_exchange.cross",
                "fitDegree": 0,
                "coefficients": [1.0],
            }
        ]

        runtime = authored_to_runtime_case(case)
        products = runtime["physics"]["reactions"][0]["products"]

        self.assertEqual(
            ["slow", "fast"],
            [str(product["speedClass"]) for product in products],
        )

    def test_missing_speed_class_is_rejected(self) -> None:
        case = self._base_case()
        case["gasInteractions"]["reactions"] = [
            {
                "processId": "missing_speed_class",
                "projectileKind": "H-",
                "projectileFate": "consume",
                "sourcePath": "Cross-sections/negative_ion_single_stripping.cross",
                "fitDegree": 0,
                "coefficients": [1.0],
                "products": [{"particleKind": "H0", "count": 1}],
            }
        ]

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("speedClass", str(ctx.exception))

    def test_magnetic_field_directory_is_resolved_into_runtime_file_path(self) -> None:
        runtime = authored_to_runtime_case(self._base_case())

        magnetic = runtime["externalMagneticField"]
        self.assertTrue(magnetic["enabled"])
        self.assertEqual("Bfield_MTF/EXTfield.fld", magnetic["fields"][0]["filePath"])
        self.assertNotIn("directory", magnetic)

    def test_missing_magnetic_field_file_is_rejected(self) -> None:
        case = self._base_case()
        case["magneticField"]["directory"] = "does-not-exist"
        case["magneticField"]["fields"][0]["filePath"] = "missing.fld"

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("filePath does not exist", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()