import json
import tempfile
import unittest
from pathlib import Path

from . import (
    aggregate_case_diagnostics,
    parse_diagnostic_summary_txt,
    parse_grid_power_summary_txt,
    synthesize_scan_results,
)


DIAGNOSTIC_TEXT = """# it z[mm] I[mA] <x>[mm] xmax[mm] xmin[mm] <y>[mm] ymax[mm] ymin[mm] <x'>[mrad] <y'>[mrad] Dx[mrad] Dy[mrad] <V>[V] <B>[mT] rho[1/m3] sigma[m2]
0 10.0 12.5 0.1 0.2 -0.1 0.0 0.3 -0.2 1.0 2.0 3.0 4.0 5000.0 25.0 1.0e20 2.0e-19
0 20.0 11.5 0.2 0.4 -0.2 0.1 0.5 -0.3 1.5 2.5 3.5 4.5 5100.0 26.0 1.1e20 2.1e-19
"""

GRID_POWER_TEXT = """# Grid Power Load Analysis
# ID\tPower[W]\tCurrent[A]\tParticles\tIncludeInTotal\tDescription
0\t1.25\t0.5\t10\ttrue\tVolume remainder
8\t2.50\t0.8\t12\tfalse\tExtractor grid
# Total beam power (included rows only): 1.25 W
"""

LEGACY_GRID_POWER_TEXT = """# Grid Power Load Analysis
# Solid_Index\tPower[W]\tCurrent[A]\tParticles\tDescription
0\t3.5\t-0.2\t9\tStripped
8\t1.5\t-0.1\t3\tGrid_1
# Total beam power (grids only): 5.0 W
"""


class PostProcessingTest(unittest.TestCase):
    def test_parse_diagnostic_summary_txt(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "diagnostic.txt"
            path.write_text(DIAGNOSTIC_TEXT, encoding="utf-8")

            parsed = parse_diagnostic_summary_txt(path)

            self.assertEqual("I[mA]", parsed["columns"][2])
            self.assertEqual(2, len(parsed["rows"]))
            self.assertEqual(12.5, parsed["rows"][0]["I[mA]"])
            self.assertEqual(20.0, parsed["rows"][1]["z[mm]"])

    def test_parse_grid_power_summary_txt(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "grid_power.txt"
            path.write_text(GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            self.assertEqual(2, len(parsed["rows"]))
            self.assertEqual(0, parsed["rows"][0]["ID"])
            self.assertFalse(parsed["rows"][1]["IncludeInTotal"])
            self.assertEqual(1.25, parsed["totalIncludedBeamPowerWatts"])

    def test_parse_legacy_grid_power_summary_txt(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "legacy_grid_power.txt"
            path.write_text(LEGACY_GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            self.assertEqual(2, len(parsed["rows"]))
            self.assertEqual(0, parsed["rows"][0]["ID"])
            self.assertIsNone(parsed["rows"][0]["IncludeInTotal"])
            self.assertEqual("Grid_1", parsed["rows"][1]["Description"])
            self.assertEqual(5.0, parsed["totalIncludedBeamPowerWatts"])

    def test_aggregate_case_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            case_dir = Path(temp_dir) / "case_one"
            summary_dir = case_dir / "Summary"
            summary_dir.mkdir(parents=True)

            case_config = {
                "metadata": {"caseTag": "case_one"},
                "outputs": {
                    "rootDirectory": case_dir.as_posix(),
                    "summary": {"directory": "Summary"},
                },
            }
            (case_dir / "case_one.json").write_text(json.dumps(case_config), encoding="utf-8")
            (summary_dir / "case_one_diagnostic_summary.txt").write_text(DIAGNOSTIC_TEXT, encoding="utf-8")
            (summary_dir / "case_one_NEGIONBEAM_diagnostic_summary.txt").write_text(
                DIAGNOSTIC_TEXT,
                encoding="utf-8",
            )
            (summary_dir / "case_one_ALL_grid_power_summary.txt").write_text(
                GRID_POWER_TEXT,
                encoding="utf-8",
            )

            aggregated = aggregate_case_diagnostics(case_dir)

            self.assertEqual("case_one", aggregated["caseTag"])
            self.assertEqual(2, len(aggregated["diagnosticSummary"]["rows"]))
            self.assertIsNotNone(aggregated["negativeIonSummary"])
            self.assertIn("ALL", aggregated["gridPowerSummaries"])

    def test_synthesize_scan_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            manifest_cases = []
            for index, label in enumerate(("6p8kV", "7p2kV")):
                case_tag = f"scan_case_{index}"
                case_dir = root / case_tag
                summary_dir = case_dir / "Summary"
                summary_dir.mkdir(parents=True)

                case_config = {
                    "metadata": {"caseTag": case_tag},
                    "outputs": {
                        "rootDirectory": case_dir.as_posix(),
                        "summary": {"directory": "Summary"},
                    },
                }
                (case_dir / f"{case_tag}.json").write_text(json.dumps(case_config), encoding="utf-8")
                (summary_dir / f"{case_tag}_diagnostic_summary.txt").write_text(
                    DIAGNOSTIC_TEXT,
                    encoding="utf-8",
                )
                (summary_dir / f"{case_tag}_ALL_grid_power_summary.txt").write_text(
                    GRID_POWER_TEXT,
                    encoding="utf-8",
                )

                manifest_cases.append(
                    {
                        "index": index,
                        "caseTag": case_tag,
                        "label": label,
                        "value": 6800 + 400 * index,
                        "configPath": (case_dir / f"{case_tag}.json").as_posix(),
                    }
                )

            manifest_path = root / "scan-manifest.json"
            manifest_path.write_text(
                json.dumps({"metadata": {"scanTag": "demo_scan"}, "cases": manifest_cases}),
                encoding="utf-8",
            )

            synthesized = synthesize_scan_results(manifest_path)

            self.assertEqual(2, len(synthesized["cases"]))
            self.assertEqual(4, len(synthesized["diagnosticRows"]))
            self.assertEqual(4, len(synthesized["gridPowerRows"]))
            self.assertEqual("scan_case_0", synthesized["diagnosticRows"][0]["caseTag"])


if __name__ == "__main__":
    unittest.main()