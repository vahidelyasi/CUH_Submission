# NGS coverage report

Reports panel genes that fall below the required read depth or are missing from the sambamba output. Requires Python 3.9 or later and uses the standard library.

## Run

Use the two input files supplied with the exercise. Replace the paths below with their locations.

```sh
python3 -B coverage_report.py \
  --input /path/to/sambamba_output.txt \
  --genelist /path/to/panel_genes.txt \
  --threshold 30 \
  --output report.txt
```

`--threshold 30` uses the `percentage30` column. A different depth requires the corresponding column in the input. The supplied example produces seven failing genes across 25 intervals and one missing gene.

## Coverage rule

A gene fails if any reported interval has less than 100% of its bases at the chosen depth. The minimum identifies the worst interval. A length-weighted mean can round to 100% even when a short interval has a coverage gap.

`TotalExons` counts distinct measured intervals within each gene. If transcripts share exactly the same coordinates, the interval counts once and the lower percentage is kept. The report lists all failing intervals and flags genes with no rows as `MISSING`. Unusable rows are logged and skipped. A panel gene with a discarded row and no measured failure is flagged as `INCOMPLETE`.

## Run the tests

```sh
python3 -B -m unittest -v
```
