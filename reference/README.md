# Regenerating genuine R outputs

Normal tests need only Python and read the checked-in fixtures. To regenerate,
install R plus satuRn 1.20.0, limma 3.68.4, locfdr 1.1-8, jsonlite 2.0.0 and
satuRn's Bioconductor dependencies in your chosen R libraries. Inspect the output
runtime metadata if installing a newer package; do not silently relabel versions.

From the repository root:

```bash
python reference/make_inputs.py
python reference/generate.py --rscript Rscript --r-library /path/to/R/library
python -m pytest
python scripts/parity_report.py --output reports/latest
```

`--r-library` can be repeated. Omitting it uses R's normal library discovery.
The generator executes Rscript with --vanilla and checks the process exit status;
it then hashes the resulting inputs, scripts and outputs. It does not fabricate
reference values if R is missing. Small-feature and zero-variance warnings are
expected for the deliberately selected reference cases. R's required plotting
side effect is redirected to a disposable PDF in R's temporary directory.

The manifest is rewritten only after R succeeds. Changes to inputs/generators
require regeneration and review of the numerical report. JSON arrays preserve
transcript order, while string "NA" means a deliberately absent empirical result
and "Inf" denotes infinite prior/posterior df in the direct moderation fixtures.
