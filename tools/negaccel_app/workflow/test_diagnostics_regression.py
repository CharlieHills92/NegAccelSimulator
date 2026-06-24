import copy
import json
import unittest
from pathlib import Path

from . import WorkflowError, authored_to_runtime_case


class DiagnosticsRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        repo_root = Path(__file__).resolve().parents[3]
        with (repo_root / "negaccel-authoring.example.json").open("r", encoding="utf-8") as handle:
            cls.example = json.load(handle)

    def _base_case(self) -> dict:
        return copy.deepcopy(self.example)

    def test_explicit_diagnostics_survive_materialization(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.12,
                "zStartMeters": -0.01,
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
                "writePerSpeciesDiagnostics": False,
                "writePerSpeciesGridPower": True,
                "writePerSpeciesPlots": False,
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
                    }
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

        runtime = authored_to_runtime_case(case)
        diagnostics = runtime["diagnostics"]

        self.assertEqual([0.01, 0.02, 0.03], diagnostics["planes"]["sampleZPositionsMeters"])
        self.assertEqual(0.031, diagnostics["planes"]["summaryZPositionMeters"])
        self.assertEqual(0.032, diagnostics["planes"]["emitterExportZPositionMeters"])
        self.assertEqual(0.006, diagnostics["summary"]["apertureRadiusMeters"])
        self.assertEqual(0.033, diagnostics["summary"]["transmissionPlaneZPositionMeters"])
        self.assertEqual(0, diagnostics["gridPower"]["ranges"][0]["id"])
        self.assertFalse(diagnostics["gridPower"]["ranges"][0]["includeInTotal"])
        self.assertEqual(8, diagnostics["gridPower"]["ranges"][1]["id"])
        self.assertTrue(diagnostics["gridPower"]["ranges"][1]["includeInTotal"])
        self.assertFalse(diagnostics["species"]["writePerSpeciesDiagnostics"])
        self.assertTrue(diagnostics["species"]["writePerSpeciesGridPower"])

    def test_missing_diagnostics_planes_are_derived_from_geometry(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.12,
                "zStartMeters": -0.01,
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
        diagnostics = copy.deepcopy(case["diagnostics"])
        diagnostics["planes"] = {}
        diagnostics.setdefault("summary", {}).pop("transmissionPlaneZPositionMeters", None)
        case["diagnostics"] = diagnostics

        runtime = authored_to_runtime_case(case)
        resolved = runtime["diagnostics"]

        self.assertEqual([0.0, 0.009, 0.015, 0.026], resolved["planes"]["sampleZPositionsMeters"])
        self.assertEqual(0.11, resolved["planes"]["summaryZPositionMeters"])
        self.assertEqual(0.11, resolved["planes"]["emitterExportZPositionMeters"])
        self.assertEqual(0.009, resolved["summary"]["transmissionPlaneZPositionMeters"])

    def test_explicit_diagnostics_plane_outside_domain_is_rejected(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.085,
                "zStartMeters": 0.0,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.01, 0.008, 0.01],
                }
            ],
        }
        case["diagnostics"] = {
            "planes": {
                "sampleZPositionsMeters": [0.01, 0.09],
                "summaryZPositionMeters": 0.03,
                "emitterExportZPositionMeters": 0.04,
            },
            "species": {
                "writePerSpeciesDiagnostics": False,
                "writePerSpeciesGridPower": True,
                "writePerSpeciesPlots": False,
                "writeNegativeIonSummary": True,
            },
            "gridPower": {
                "ranges": [
                    {
                        "id": 0,
                        "includeInTotal": False,
                    },
                    {
                        "id": 7,
                        "includeInTotal": True,
                    }
                ]
            },
            "summary": {
                "apertureRadiusMeters": 0.006,
                "transmissionPlaneZPositionMeters": 0.03,
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

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("outside geometry.domain z range", str(ctx.exception))

    def test_source_plane_outside_domain_is_rejected_with_runtime_hint(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.085,
                "zStartMeters": 0.0,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.01, 0.008, 0.01],
                }
            ],
        }
        diagnostics = copy.deepcopy(case["diagnostics"])
        diagnostics["planes"] = {}
        diagnostics.setdefault("summary", {}).pop("transmissionPlaneZPositionMeters", None)
        case["diagnostics"] = diagnostics
        particles = case.setdefault("particles", {})
        sources = particles.setdefault("sources", [])
        if sources:
            sources[0].setdefault("uniform", {})["centerMeters"] = [0.0, 0.0, -0.003]
            sources[0]["uniform"]["heightMeters"] = 0.02

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("outside geometry.domain bounds", str(ctx.exception))
        self.assertIn("bad definitions", str(ctx.exception))

    def test_source_on_domain_boundary_is_accepted(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.085,
                "zStartMeters": 0.0,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.01, 0.008, 0.01],
                }
            ],
        }
        diagnostics = copy.deepcopy(case["diagnostics"])
        diagnostics["planes"] = {}
        diagnostics.setdefault("summary", {}).pop("transmissionPlaneZPositionMeters", None)
        case["diagnostics"] = diagnostics
        particles = case.setdefault("particles", {})
        sources = particles.setdefault("sources", [])
        if sources:
            sources[0].setdefault("uniform", {})["centerMeters"] = [0.0, 0.0, 0.0]
            sources[0]["uniform"]["heightMeters"] = 0.02

        runtime = authored_to_runtime_case(case)

        self.assertEqual([0.0, 0.0, 0.0], runtime["particleSources"][0]["uniform"]["centerMeters"])

    def test_diagnostics_plane_on_shifted_domain_boundary_is_accepted(self) -> None:
        case = self._base_case()
        case["geometry"] = {
            "meshSizeMeters": 0.002,
            "exportGeometryVtk": True,
            "domain": {
                "xSizeMeters": 0.08,
                "ySizeMeters": 0.08,
                "zSizeMeters": 0.088,
                "zStartMeters": -0.003,
            },
            "solids": [
                {
                    "boundaryId": 7,
                    "name": "PG",
                    "kind": "solid",
                    "zProfileMeters": [0.0, 0.006, 0.009],
                    "rProfileMeters": [0.01, 0.008, 0.01],
                }
            ],
        }
        case["diagnostics"] = {
            "planes": {
                "sampleZPositionsMeters": [0.0, 0.009],
                "summaryZPositionMeters": 0.085,
                "emitterExportZPositionMeters": 0.085,
            },
            "species": {
                "writePerSpeciesDiagnostics": False,
                "writePerSpeciesGridPower": True,
                "writePerSpeciesPlots": False,
                "writeNegativeIonSummary": True,
            },
            "gridPower": {
                "ranges": [
                    {
                        "id": 0,
                        "includeInTotal": False,
                    },
                    {
                        "id": 7,
                        "includeInTotal": True,
                    }
                ]
            },
            "summary": {
                "apertureRadiusMeters": 0.006,
                "transmissionPlaneZPositionMeters": 0.085,
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

        runtime = authored_to_runtime_case(case)

        self.assertEqual(0.085, runtime["diagnostics"]["planes"]["summaryZPositionMeters"])
        self.assertEqual(0.085, runtime["diagnostics"]["summary"]["transmissionPlaneZPositionMeters"])

    def test_missing_diagnostics_is_rejected(self) -> None:
        case = self._base_case()
        case.pop("diagnostics", None)

        with self.assertRaises(WorkflowError) as ctx:
            authored_to_runtime_case(case)

        self.assertIn("diagnostics", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()