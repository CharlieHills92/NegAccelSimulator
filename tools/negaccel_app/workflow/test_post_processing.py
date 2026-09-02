import json
import tempfile
import unittest
from pathlib import Path

from . import (
    aggregate_case_diagnostics,
    parse_diagnostic_summary_txt,
    parse_grid_power_breakdown_txt,
    parse_grid_power_summary_txt,
    synthesize_scan_results,
)
from .common import WorkflowError


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

EXCLUDED_GRID_POWER_TEXT = """# Grid Power Load Analysis
# ID\tPower[W]\tCurrent[A]\tParticles\tIncludeInTotal\tDescription
5\t0.00687\t-0.00014256\t66\tfalse\tz-min
7\t0.223\t-0.0532526\t24654\tfalse\tPG
# Total beam power (included rows only): 0.23 W
# Excluded from the rows above: un-extracted primary negative ions (generation 0, moving upstream)
# ExcludedRow\t5\t1751\t-0.00378216
# ExcludedRow\t7\t285\t-0.0006156
# TotalExcluded\t2036\t-0.00439776
"""

NET_GRID_POWER_TEXT = """# Grid Power Load Analysis
# ID\tPower[W]\tCurrent[A]\tParticles\tIncludeInTotal\tDescription
0\t1.25\t0.5\t10\ttrue\tVolume remainder
7\t10.0\t-0.04\t500\ttrue\tPG
# Total beam power (included rows only): 11.25 W
# Excluded from the rows above: un-extracted primary negative ions (generation 0, moving upstream)
# ExcludedRow\t7\t12\t-0.0003
# TotalExcluded\t12\t-0.0003
# Net accounting: energy and charge carried back off each surface by the secondaries it emitted. Net = Gross + Emitted, Emitted <= 0.
# NetRow\t0\t0\t1.25\t0\t0.5
# NetRow\t7\t-2.5\t7.5\t0.01\t-0.03
# TotalNet\t-2.5\t8.75\t0.01\t0.47
"""

GRID_POWER_BREAKDOWN_TEXT = """# Grid Power Load Breakdown by species and generation
# RowID\tSpecies\tGen\tOrigin\tGrossPower[W]\tEmittedPower[W]\tNetPower[W]\tGrossCurrent[A]\tEmittedCurrent[A]\tNetCurrent[A]\tEquivalentCurrent[A]\tParticles\tIncludeInTotal\tDescription
7\tHM\t0\tprimary\t9.0\t0\t9.0\t-0.04\t0\t-0.04\t0.04\t400\ttrue\tPG
7\tH0\t1\tvolume\t1.0\t0\t1.0\t0\t0\t0\t0.006\t90\ttrue\tPG
7\tE\t101\tsurface\t0.0\t-2.5\t-2.5\t0\t0.01\t0.01\t0.0\t0\ttrue\tPG
7\tunclassified\t0\tprimary\t0.5\t0\t0.5\t0\t0\t0\t0.001\t10\ttrue\tPG
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

    def test_parse_grid_power_summary_excluded_unextracted(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "excluded_grid_power.txt"
            path.write_text(EXCLUDED_GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            # Comment lines must not be mistaken for data rows.
            self.assertEqual(2, len(parsed["rows"]))
            self.assertEqual(66, parsed["rows"][0]["Particles"])

            excluded = parsed["excludedUnextractedRows"]
            self.assertEqual({5, 7}, set(excluded))
            self.assertEqual(1751, excluded[5]["Particles"])
            self.assertAlmostEqual(-0.00378216, excluded[5]["Current[A]"])
            self.assertEqual(285, excluded[7]["Particles"])
            self.assertEqual(2036, parsed["totalExcludedUnextracted"]["Particles"])
            self.assertAlmostEqual(-0.00439776,
                                   parsed["totalExcludedUnextracted"]["Current[A]"])

    def test_parse_grid_power_summary_without_excluded_section(self) -> None:
        # Files written before the excluded tally existed must still parse.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "grid_power.txt"
            path.write_text(GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            self.assertEqual({}, parsed["excludedUnextractedRows"])
            self.assertIsNone(parsed["totalExcludedUnextracted"])

    def test_parse_grid_power_summary_merges_net_into_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "net_grid_power.txt"
            path.write_text(NET_GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            self.assertEqual(2, len(parsed["rows"]))
            by_id = {row["ID"]: row for row in parsed["rows"]}

            # Power[W]/Current[A] must keep their GROSS meaning: the GUI plots them.
            self.assertEqual(10.0, by_id[7]["Power[W]"])
            self.assertEqual(-0.04, by_id[7]["Current[A]"])

            # ...with net merged in alongside, from the comment lines.
            self.assertEqual(-2.5, by_id[7]["EmittedPower[W]"])
            self.assertEqual(7.5, by_id[7]["NetPower[W]"])
            self.assertEqual(0.01, by_id[7]["EmittedCurrent[A]"])
            self.assertEqual(-0.03, by_id[7]["NetCurrent[A]"])

            # Excluded tally merged by ID too, defaulting to zero where absent.
            self.assertEqual(12, by_id[7]["ExcludedParticles"])
            self.assertAlmostEqual(-0.0003, by_id[7]["ExcludedCurrent[A]"])
            self.assertEqual(0, by_id[0]["ExcludedParticles"])

            self.assertEqual(8.75, parsed["totalNet"]["NetPower[W]"])
            self.assertIn("NetPower[W]", parsed["columns"])

    def test_parse_grid_power_summary_net_defaults_when_absent(self) -> None:
        # Without a surface-emission ledger nothing was emitted, so net == gross exactly.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "grid_power.txt"
            path.write_text(GRID_POWER_TEXT, encoding="utf-8")

            parsed = parse_grid_power_summary_txt(path)

            for row in parsed["rows"]:
                self.assertEqual(0.0, row["EmittedPower[W]"])
                self.assertEqual(row["Power[W]"], row["NetPower[W]"])
                self.assertEqual(0.0, row["EmittedCurrent[A]"])
                self.assertEqual(row["Current[A]"], row["NetCurrent[A]"])
            self.assertEqual({}, parsed["netRows"])
            self.assertIsNone(parsed["totalNet"])

    def test_parse_grid_power_breakdown_txt(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "case_grid_power_breakdown.txt"
            path.write_text(GRID_POWER_BREAKDOWN_TEXT, encoding="utf-8")

            parsed = parse_grid_power_breakdown_txt(path)

            self.assertEqual(4, len(parsed["rows"]))
            rows = {(r["Species"], r["Gen"]): r for r in parsed["rows"]}

            # Neutral species: real power and a visible equivalent current, but no current.
            neutral = rows[("H0", 1)]
            self.assertEqual(1.0, neutral["GrossPower[W]"])
            self.assertEqual(0.0, neutral["GrossCurrent[A]"])
            self.assertEqual(0.0, neutral["NetCurrent[A]"])
            self.assertGreater(neutral["EquivalentCurrent[A]"], 0.0)
            self.assertEqual("volume", neutral["Origin"])

            # Surface-emitted electrons: pure debit, raw generation preserved.
            emitted = rows[("E", 101)]
            self.assertEqual("surface", emitted["Origin"])
            self.assertEqual(-2.5, emitted["EmittedPower[W]"])
            self.assertEqual(0, emitted["Particles"])

            # PARTICLE_WRONG gets its own row so the species sums reconcile.
            self.assertIn(("unclassified", 0), rows)
            self.assertTrue(rows[("HM", 0)]["IncludeInTotal"])

            # Gross power over all species/generations reconciles with the row total.
            self.assertAlmostEqual(10.5, sum(r["GrossPower[W]"] for r in parsed["rows"]))

    def test_parse_grid_power_breakdown_rejects_wrong_field_count(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "bad_grid_power_breakdown.txt"
            path.write_text("7\tHM\t0\tprimary\t1.0\n", encoding="utf-8")

            with self.assertRaises(WorkflowError):
                parse_grid_power_breakdown_txt(path)

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
            (summary_dir / "case_one_ALL_grid_power_breakdown.txt").write_text(
                GRID_POWER_BREAKDOWN_TEXT,
                encoding="utf-8",
            )

            aggregated = aggregate_case_diagnostics(case_dir)

            self.assertEqual("case_one", aggregated["caseTag"])
            self.assertEqual(2, len(aggregated["diagnosticSummary"]["rows"]))
            self.assertIsNotNone(aggregated["negativeIonSummary"])
            self.assertIn("ALL", aggregated["gridPowerSummaries"])
            self.assertIn("ALL", aggregated["gridPowerBreakdowns"])
            self.assertEqual(4, len(aggregated["gridPowerBreakdowns"]["ALL"]["rows"]))

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