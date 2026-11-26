import copy
import json
import unittest
from pathlib import Path

from . import WorkflowError, authored_to_runtime_case


class OutputsRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        with (repo_root / "negaccel-authoring.example.json").open("r", encoding="utf-8") as handle:
            cls.example = json.load(handle)

    def _base_case(self) -> dict:
        case = copy.deepcopy(self.example)
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.12,
                "zStartMeters": 0.0,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.01, 0.008, 0.01],
                },
                {
                    "boundaryId": 8,
                    "name": "EG",
                    "kind": "solid",
                    "zProfileMeters": [0.015, 0.023, 0.026],
                    "rProfileMeters": [0.006, 0.006, 0.007],
                },
            ],
        }
        case["diagnostics"] = {
            "planes": {
                "sampleZPositionsMeters": [0.01, 0.02, 0.03],
                "summaryZPositionMeters": 0.031,
                "emitterExportZPositionMeters": 0.032,
            },
            "species": {
                "writePerSpeciesDiagnostics": True,
                "writePerSpeciesGridPower": True,
                "writePerSpeciesPlots": True,
                "writeNegativeIonSummary": True,
            },
            "gridPower": {
                "ranges": [
                    {
                        "id": 0,
                        "includeInTotal": False,
                    },
                    {
                        "id": 8,
                        "includeInTotal": True,
                    },
                ]
            },
            "summary": {
                "apertureRadiusMeters": 0.006,
                "transmissionPlaneZPositionMeters": 0.033,
            },
            "plots": {
                "meniscus": {
                    "enabled": True,
                    "zMinMeters": 0.0,
                    "zMaxMeters": 0.04,
                    "transverseMinMeters": -0.02,
                    "transverseMaxMeters": 0.02,
                }
            },
        }
        return case

    def test_iteration_output_settings_materialize(self) -> None:
        case = self._base_case()
        case["outputs"]["iteration"] = {
            "enabled": True,
            "everyNIterations": 3,
            "exportPlaneDiagnostics": True,
            "exportSimulationState": True,
            "exportTracedParticles": False,
            "planeZPositionsMeters": [0.009, 0.06, 0.06, 0.11],
        }

        runtime = authored_to_runtime_case(case)
        iteration = runtime["outputs"]["iteration"]

        self.assertTrue(iteration["enabled"])
        self.assertEqual(3, iteration["everyNIterations"])
        self.assertTrue(iteration["exportPlaneDiagnostics"])
        self.assertTrue(iteration["exportSimulationState"])
        self.assertFalse(iteration["exportTracedParticles"])
        self.assertEqual([0.009, 0.06, 0.11], iteration["planeZPositionsMeters"])

    def test_iteration_plane_outside_domain_is_rejected(self) -> None:
        case = self._base_case()
        case["outputs"]["iteration"] = {
            "enabled": True,
            "everyNIterations": 2,
            "exportPlaneDiagnostics": True,
            "exportSimulationState": False,
            "exportTracedParticles": False,
            "planeZPositionsMeters": [0.009, 0.14],
        }

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("outputs.iteration.planeZPositionsMeters[1]", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()