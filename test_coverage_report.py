import unittest
from io import StringIO

from coverage_report import make_report, read_coverage


HEADER = "#chromosome\tStartPosition EndPosition\tGeneSymbol;Accession\tpercentage20\tpercentage30\n"


def report_for(data, panel, threshold=30):
    coverage, skipped = read_coverage(StringIO(HEADER + data), threshold)
    report, summary = make_report(panel, coverage, skipped, "test", threshold)
    return report.splitlines(), summary


class CoverageReportTests(unittest.TestCase):
    def test_short_failing_exon_shared_by_two_transcripts(self):
        data = (
            "2\t100\t10000\tA;TX1\t100\t100\n"
            "2\t11000\t11010\tA;TX1\t100\t90\n"
            "2\t11000\t11010\tA;TX2\t100\t100\n"
        )
        rows, _ = report_for(data, ["A"])
        self.assertIn("A\t2\t1\t90\tchr2:11000-11010\tFAIL", rows)
        self.assertEqual(rows.count("A\tchr2:11000-11010\t90"), 1)

    def test_missing_gene_and_measured_zero_have_different_results(self):
        data = (
            "1\t10\t20\tA;TX1\t0\t0\n"
            "1\t40\t50\tB;TX2\t100\t100\n"
            "1\t60\t70\tOFF_PANEL;TX3\t10\t10\n"
        )
        rows, summary = report_for(data, ["A", "B", "C"])
        self.assertIn("A\t1\t1\t0\tchr1:10-20\tFAIL", rows)
        self.assertIn("C\t0\tNA\tNA\tNA\tMISSING", rows)
        self.assertFalse(any(row.startswith(("B\t", "OFF_PANEL\t")) for row in rows))
        self.assertEqual(summary, "Sample test: 1/3 genes below 30x threshold, 1 missing")

    def test_bad_rows_warn_and_allow_other_genes_to_be_reported(self):
        data = (
            "1\t10\t20\tA;TX1\t100\tNA\n"
            "1\t20\t30\tA;TX1\t100\t100\n"
            "1\t30\t40\t.\t0\t0\n"
            "1\t40\t50\t;TX2\t0\t0\n"
            "1\t50\t60\tB;TX3\t100\tNaN\n"
            "1\t60\t70\tC;TX4\t100\t99.99\n"
        )
        with self.assertLogs(level="WARNING") as messages:
            rows, summary = report_for(data, ["A", "B", "C"])
        self.assertEqual(len(messages.output), 4)
        self.assertIn("Line 2 (A)", messages.output[0])
        self.assertIn("missing gene symbol", messages.output[1])
        self.assertIn("A\t1\tNA\tNA\tNA\tINCOMPLETE", rows)
        self.assertIn("B\t0\tNA\tNA\tNA\tINCOMPLETE", rows)
        self.assertIn("C\t1\t1\t99.99\tchr1:60-70\tFAIL", rows)
        self.assertEqual(summary, "Sample test: 1/3 genes below 30x threshold, 0 missing, 2 incomplete")

    def test_threshold_selects_the_corresponding_column(self):
        data = "1\t10\t20\tA;TX1\t100\t80\n"
        rows, summary = report_for(data, ["A"], threshold=20)
        self.assertFalse(any(row.startswith("A\t") for row in rows))
        self.assertEqual(summary, "Sample test: 0/1 genes below 20x threshold, 0 missing")
        with self.assertRaisesRegex(ValueError, "percentage50"):
            report_for(data, ["A"], threshold=50)


if __name__ == "__main__":
    unittest.main()
