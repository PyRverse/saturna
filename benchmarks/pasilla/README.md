# Pasilla real-data benchmark

Saturna versus genuine satuRn 1.20.0 on the PasillaTranscriptExpr kallisto
counts (Brooks et al. 2011; 4 control and 3 knockdown samples). Methods, terms,
results, and every material discrepancy are described in
[docs/realdata.md](../../docs/realdata.md).

| File | Role |
| --- | --- |
| source.json | Pinned Bioconductor archive URL, SHA-256 hashes, license, citations |
| prepare.py | Downloads and verifies the archive, applies the shared filter, writes inputs/ |
| inputs/ | Shared counts (with `analysis_set` labels), design, and filter summary |
| satuRn.R | Re-derives the filter from raw files, runs fitDTU/testDTU for both analysis sets |
| reference.py | Runs satuRn.R with `Rscript --vanilla` and rewrites manifest.json after success |
| reference/ | Committed genuine R outputs (17 significant digits) and runtime metadata |
| compare.py | Runs Saturna, applies the pre-registered criteria and regression gate, writes reports |

```bash
python benchmarks/pasilla/prepare.py
python benchmarks/pasilla/reference.py --rscript Rscript --r-library /path/to/R/library
python benchmarks/pasilla/compare.py --output reports/latest/pasilla
```

The raw archive is cached in `.cache/` (ignored by Git). `compare.py` needs only
Python and the committed files, so it can run without R or network access.
