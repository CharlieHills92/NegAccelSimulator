import copy
import json
import unittest
from pathlib import Path

from . import WorkflowError, authored_to_runtime_case


class SolverSettingsRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        with (repo_root / "negaccel-authoring.example.json").open("r", encoding="utf-8") as handle:
            cls.example = json.load(handle)

    def _base_case(self) -> dict:
        return copy.deepcopy(self.example)

    def test_bicgstab_settings_materialize_on_selected_solver(self) -> None:
        case = self._base_case()
        case["run"]["solver"] = {
            "type": "bicgstab",
            "bicgstab": {
                "eps": 5.0e-5,
                "maxIterations": 4321,
                "newtonEps": 2.5e-4,
                "newtonMaxIterations": 7,
                "globallyConvergentNewton": False,
            },
            "multigrid": {
                "levels": 3,
                "mgTolerance": 1.0e-5,
            },
        }

        runtime = authored_to_runtime_case(case)
        solver = runtime["simulation"]["solver"]

        self.assertEqual("bicgstab", solver["type"])
        self.assertEqual(5.0e-5, solver["bicgstab"]["eps"])
        self.assertEqual(4321, solver["bicgstab"]["maxIterations"])
        self.assertEqual(2.5e-4, solver["bicgstab"]["newtonEps"])
        self.assertEqual(7, solver["bicgstab"]["newtonMaxIterations"])
        self.assertFalse(solver["bicgstab"]["globallyConvergentNewton"])
        self.assertNotIn("multigrid", solver)

    def test_multigrid_settings_materialize_on_selected_solver(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "name": "SPIDER",
            "meshSizeMeters": 0.0005,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.07,
                "ySizeMeters": 0.07,
                "zSizeMeters": 0.09,
                "zStartMeters": -0.003,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.012, 0.007, 0.01],
                },
                {
                    "boundaryId": 8,
                    "name": "EG",
                    "kind": "solid",
                    "zProfileMeters": [0.015, 0.023, 0.026],
                    "rProfileMeters": [0.0065, 0.0065, 0.0075],
                },
            ],
        }
        case["diagnostics"]["planes"] = {
            "sampleZPositionsMeters": [0.0, 0.009, 0.02, 0.05, 0.0835],
            "summaryZPositionMeters": 0.0835,
            "emitterExportZPositionMeters": 0.0835,
        }
        case["diagnostics"]["summary"] = {
            "apertureRadiusMeters": 0.007,
            "transmissionPlaneZPositionMeters": 0.0835,
        }
        case["run"]["solver"] = {
            "type": "multigrid",
            "bicgstab": {
                "eps": 5.0e-5,
            },
            "multigrid": {
                "levels": 3,
                "mgTolerance": 1.0e-5,
                "maxCycles": 45,
                "gamma": 2,
                "preSmooth": 4,
                "postSmooth": 6,
                "coarseRelaxation": 1.5,
                "coarseMaxIterations": 3000,
                "localPlasmaMaxIterations": 2,
            },
        }

        runtime = authored_to_runtime_case(case)
        solver = runtime["simulation"]["solver"]

        self.assertEqual("multigrid", solver["type"])
        self.assertEqual(3, solver["multigrid"]["levels"])
        self.assertEqual(1.0e-5, solver["multigrid"]["mgTolerance"])
        self.assertEqual(45, solver["multigrid"]["maxCycles"])
        self.assertEqual(2, solver["multigrid"]["gamma"])
        self.assertEqual(4, solver["multigrid"]["preSmooth"])
        self.assertEqual(6, solver["multigrid"]["postSmooth"])
        self.assertEqual(1.5, solver["multigrid"]["coarseRelaxation"])
        self.assertEqual(3000, solver["multigrid"]["coarseMaxIterations"])
        self.assertEqual(2, solver["multigrid"]["localPlasmaMaxIterations"])
        self.assertNotIn("bicgstab", solver)

    def test_invalid_multigrid_levels_are_rejected_before_runtime(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "name": "SPIDER",
            "meshSizeMeters": 0.0005,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.07,
                "ySizeMeters": 0.07,
                "zSizeMeters": 0.09,
                "zStartMeters": -0.003,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.012, 0.007, 0.01],
                },
                {
                    "boundaryId": 8,
                    "name": "EG",
                    "kind": "solid",
                    "zProfileMeters": [0.015, 0.023, 0.026],
                    "rProfileMeters": [0.0065, 0.0065, 0.0075],
                },
            ],
        }
        case["run"]["solver"] = {
            "type": "multigrid",
            "multigrid": {
                "levels": 4,
            },
        }

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("maximum allowed 3", str(ctx.exception))
        self.assertIn("141 x 141 x 181", str(ctx.exception))

    def test_execution_server_command_round_trips_to_runtime(self) -> None:
        case = self._base_case()
        case["execution"] = {
            "serverCommand": "srun --pty /bin/bash",
        }

        runtime = authored_to_runtime_case(case)

        self.assertIn("execution", runtime)
        self.assertEqual("srun --pty /bin/bash", runtime["execution"]["serverCommand"])
        self.assertEqual(case["run"]["iterations"], runtime["simulation"]["iterations"])

    def test_initial_plasma_max_z_round_trips_to_runtime(self) -> None:
        case = self._base_case()
        case["particles"]["plasma"]["initialPlasmaMaxZMeters"] = 0.0125

        runtime = authored_to_runtime_case(case)

        self.assertEqual(0.0125, runtime["simulation"]["solver"]["initialPlasmaMaxZMeters"])


if __name__ == "__main__":
    unittest.main()