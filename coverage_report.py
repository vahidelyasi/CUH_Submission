"""Report panel genes with incomplete sequencing coverage."""

import argparse
import logging
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path


def read_panel(path):
    """Read unique gene symbols in the order they appear."""
    panel = []
    with path.open(encoding="utf-8") as source:
        for line in source:
            gene = line.strip()
            if gene and gene not in panel:
                panel.append(gene)

    if not panel:
        raise ValueError("panel gene list is empty")

    return panel


def read_percentage(text, column):
    """Parse a finite coverage percentage between 0 and 100."""
    try:
        value = Decimal(text)
    except InvalidOperation:
        raise ValueError(f"invalid {column} value {text!r}") from None
    if not value.is_finite() or not 0 <= value <= 100:
        raise ValueError(f"{column} outside the range 0 to 100: {text!r}")
    return value


def read_interval(row):
    """Read and validate the chromosome, start and end of one interval."""
    chrom = row["chromosome"].strip().removeprefix("chr")
    start, end = int(row["StartPosition"]), int(row["EndPosition"])
    if chrom in ("", ".") or start < 0 or end < start:
        raise ValueError("invalid exon coordinates")
    return chrom, start, end


def read_coverage(source, threshold):
    """Collect measured intervals by gene and record genes with unusable rows."""
    column = f"percentage{threshold}"
    # The supplied header mixes spaces and tabs. The data rows are tab-separated.
    header = source.readline().lstrip("#").split()
    for name in ("chromosome", "StartPosition", "EndPosition", "GeneSymbol;Accession", column):
        if name not in header:
            raise ValueError(f"input has no {name!r} column")

    coverage = defaultdict(dict)
    skipped = set()
    for number, line in enumerate(source, start=2):
        if not line.strip():
            continue
        fields = line.rstrip("\r\n").split("\t")
        row = dict(zip(header, fields))
        gene = row.get("GeneSymbol;Accession", "").partition(";")[0].strip()
        if gene in ("", "."):
            logging.warning("Line %d: missing gene symbol, skipping row", number)
            continue
        try:
            if len(fields) != len(header):
                raise ValueError("unexpected number of tab-separated fields")
            percent = read_percentage(row[column], column)
            interval = read_interval(row)
        except ValueError as error:
            logging.warning("Line %d (%s): %s, skipping row", number, gene, error)
            skipped.add(gene)
            continue

        exons = coverage[gene]
        # Transcripts can share an interval. Count it once and keep its lowest measurement.
        if interval not in exons or percent < exons[interval]:
            exons[interval] = percent

    return coverage, skipped


def read_coverage_file(path, threshold):
    """Open a sambamba file and read its coverage data."""
    with path.open(encoding="utf-8") as source:
        return read_coverage(source, threshold)


def coordinates(region):
    """Format an interval as chr:start-end for the report."""
    chrom, start, end = region
    return f"chr{chrom}:{start}-{end}"


def make_report(panel, coverage, skipped, sample, threshold):
    """Walk the panel in order so a gene with no input rows is still reviewed."""
    rows = ["GeneSymbol\tTotalExons\tFailingExons\tMinCoverage%\tWorstExon(chr:start-end)\tStatus"]
    details = []
    n_failed = n_missing = n_incomplete = 0
    for gene in panel:
        exons = coverage.get(gene, {})
        # Every base must meet the depth threshold. The lowest exon percentage
        # detects a short failing exon even when the rest of the gene is fully covered.
        failing = []
        for interval, percent in exons.items():
            if percent < 100:
                failing.append((percent, interval))

        failing.sort()
        if failing:
            n_failed += 1
            lowest, worst = failing[0]
            rows.append(
                f"{gene}\t{len(exons)}\t{len(failing)}\t{lowest}\t{coordinates(worst)}\tFAIL"
            )
            for pct, region in failing:
                details.append(f"{gene}\t{coordinates(region)}\t{pct}")
        elif gene in skipped:
            n_incomplete += 1
            rows.append(f"{gene}\t{len(exons)}\tNA\tNA\tNA\tINCOMPLETE")
        elif not exons:
            n_missing += 1
            rows.append(f"{gene}\t0\tNA\tNA\tNA\tMISSING")

    summary = (
        f"Sample {sample}: {n_failed}/{len(panel)} genes below {threshold}x threshold, "
        f"{n_missing} missing"
    )
    if n_incomplete:
        summary += f", {n_incomplete} incomplete"
    rows.extend(["", summary])
    if details:
        rows.extend(["", f"Failing exons at {threshold}x"])
        rows.append("GeneSymbol\tExon(chr:start-end)\tCoverage%")
        rows.extend(details)
    return "\n".join(rows) + "\n", summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="single-sample coverage file")
    parser.add_argument("--genelist", type=Path, required=True, help="one panel gene per line")
    parser.add_argument("--threshold", type=int, default=30, help="read depth (default: 30)")
    parser.add_argument("--output", type=Path, required=True, help="output report path")
    args = parser.parse_args()

    try:
        if args.output.resolve() in (args.input.resolve(), args.genelist.resolve()):
            raise ValueError("choose an output path separate from both input files")
        panel = read_panel(args.genelist)
        coverage, skipped = read_coverage_file(args.input, args.threshold)

        # Use the filename to get the sample label shown in the brief, e.g. NGS148_34.
        # The supplied sampleName column contains '1' on every row.
        stem = args.input.stem.removesuffix(".sambamba_output")
        parts = stem.split("_", 2)
        sample = stem
        if len(parts) >= 2:
            run, sample_number = parts[:2]
            if run.startswith("NGS") and run[3:].isdecimal() and sample_number.isdecimal():
                sample = f"{run}_{sample_number}"
        report, summary = make_report(panel, coverage, skipped, sample, args.threshold)
        args.output.write_text(report, encoding="utf-8")
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(summary)


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
    main()
