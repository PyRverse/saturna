# Real-data validation: pasilla knockdown

This is Saturna's first real-data benchmark. It compares Saturna with the genuine
satuRn 1.20.0 reference on a public bulk RNA-seq experiment, using identical
transcript filtering, gene mapping, design, and contrast. The generated report is
[reports/validated/pasilla.md](../reports/validated/pasilla.md); the benchmark
code and committed inputs and outputs are in [benchmarks/pasilla](../benchmarks/pasilla).

**Outcome.** All fit and regular-test quantities agree within the synthetic
suite's tolerances (maximum relative error 7.9e-11). Calls are identical at FDR
0.01, 0.05, and 0.10: 75, 293, and 539 transcripts respectively. One
pre-registered criterion failed: 364 of 10,003 empirical p-values differ by more
than the 3e-6 relative tolerance (maximum 1.1e-5). The cause is the
empirical-null estimate, not the empirical-p formula (discrepancy D2). Saturna
rejects 52 uninformative genes that satuRn fits; if those genes are left in,
satuRn's results for all other transcripts change (discrepancy D1).

This benchmark shows numerical equivalence between the two tools on one real
dataset. It is not evidence of statistical calibration, power, or biological
correctness.

## Dataset and terms

| Item | Value |
| --- | --- |
| Source | Bioconductor 3.23 experiment package PasillaTranscriptExpr 1.40.0, pinned by SHA-256 in [source.json](../benchmarks/pasilla/source.json) |
| Experiment | Brooks et al. (2011), GEO GSE18508 / SRA SRP001537 (modENCODE); *Drosophila melanogaster* S2-DRSC cells |
| Samples | 4 untreated controls and 3 pasilla (CG8144) dsRNA knockdowns, GSM461176–GSM461182 |
| Counts | kallisto 0.42.1 expected (fractional) counts per transcript, Ensembl 70 / BDGP5.70 cDNA |
| Gene mapping | transcript-to-gene IDs from the Ensembl 70 GTF, supplied in the package's counts table |
| Data package license | GPL (>= 3), compatible with Saturna's GPL-3.0-or-later; the filtered derivative is redistributed here with attribution in NOTICE.md |
| Raw-data terms | GEO refers users to the ENCODE/modENCODE data release policy. NHGRI states that these data "are available for all to use without restriction" and asks users to cite the consortium. The 2008 policy's nine-month publication embargo ended long ago (data released 2009, published 2011). |

Brooks et al. and the package vignette describe the libraries as "7 biologically
independent samples". These are separately treated cultures of one cell line,
not different organisms. That fits Saturna's independent-bulk-sample contract;
there is no pairing or repeated measure. Library layout is partly confounded
with condition: controls are 2 single-end and 2 paired-end, knockdowns 1
single-end and 2 paired-end. The design therefore includes layout as a
covariate, as the DRIMSeq pasilla workflow does. Sample GSM461179 was sequenced
in several runs with different read lengths, which the package summed.

Estimated counts are fractional. Both tools accept them, and R raised no
warnings.

## Shared analysis

- **Filter** ([prepare.py](../benchmarks/pasilla/prepare.py)), adapted from the
  DRIMSeq pasilla thresholds. Keep transcripts whose gene total over all
  annotated transcripts is at least 10 in all 7 samples, and whose own count is
  at least 10 in at least 3 samples (the smallest group). Then keep genes with
  at least 2 such transcripts. Result: 27,142 annotated → 10,107 transcripts in
  3,248 genes.
- **Uninformative genes.** 52 two-transcript genes have identical count vectors
  for both transcripts in every sample: kallisto could not distinguish them and
  split the reads equally. Their usage is fixed at 0.5 by construction. The
  primary analysis excludes them, leaving **10,003 transcripts in 3,196 genes**.
  This rule was added after Saturna rejected those genes (D1). It is stated in
  terms of the data, applied identically to both tools, and its effect is
  measured below.
- **Design** `~ 0 + intercept + knockdown + paired_end`; **contrast** `knockdown`
  `(0, 1, 0)`, natural-log odds. Residual df is 4, except for 6 transcripts in 3
  genes. In those, one sample has a zero total over the retained transcripts;
  both tools drop that sample, leaving df 3. This exercises the unequal-df
  moderation path.
- **Same inputs, checked by both tools.** Both read the same
  `inputs/counts.csv` and `inputs/design.csv`. The R script independently
  re-derives the filter labels, count values, and design from the raw
  Bioconductor files, and stops if any of them differ.
- **Reference:** R 4.6.1, satuRn 1.20.0, limma 3.68.4, locfdr 1.1.8, run as
  `fitDTU` then `testDTU` with default arguments. All 10,003 models are
  ordinary `glm` fits with no warnings. Saturna needed at most 9 IRLS
  iterations, below R's limit of 25.

## Pre-registered criteria and results

The criteria were fixed in `compare.py` before the first R reference was
generated:

- an identical transcript set;
- every satuRn model a finite `glm` fit;
- every quantity within the synthetic suite's tolerances (regular rtol 2e-8 /
  atol 2e-10; empirical rtol 3e-6 / atol 3e-8);
- identical regular and empirical calls at FDR 0.05;
- Saturna still rejects the unrestricted set.

Five passed; the tolerance criterion failed for `empirical_pvalue` only.

| Quantity | Max abs. error | Max rel. error |
| --- | ---: | ---: |
| Coefficients, estimate | 8.6e-13 | 7.9e-11 |
| Dispersion / posterior variance | 1.8e-10 / 2.9e-9 | 1.3e-12 |
| Standard error, statistic | 1.0e-12 / 9.4e-12 | 5.3e-12 |
| Regular p-value / FDR | 2.9e-13 / 4.6e-13 | 3.4e-12 |
| Empirical p-value | 5.0e-7 | 1.1e-5 (364 outside tolerance) |
| Empirical FDR | 4.2e-7 | 4.2e-7 |

The limma prior df (2.1643) and prior variance (18.608) agree to 1e-11. The
regular FDR closest to 0.05 in either tool is 3.6e-4 away from the threshold,
nine orders of magnitude larger than the largest FDR difference. Call agreement
is therefore not fragile.

## Material discrepancies

**D1 — Uninformative genes: Saturna rejects them; satuRn fits them and they
change other transcripts' results.**

On the unrestricted 10,107 transcripts, Saturna raises `ModelFitError`
("residual dispersion is effectively zero"). Its support contract rejects
perfectly explained models instead of returning p-values for them.

satuRn fits all 104 transcripts with dispersion exactly 0, estimate 0, and
p-value 1, so they produce no calls themselves. Their zero dispersions still
enter limma's prior estimation:

- prior variance falls from 18.608 to 13.920;
- prior df falls from 2.164 to 2.142;
- on the 10,003 shared transcripts, regular calls at FDR 0.05 rise from 293 to
  **373** (80 gained, none lost), and p-values move by up to 0.064.

In other words, when analysed with satuRn, the quantifier's indistinguishable
isoforms make every other transcript look more significant.

Saturna users must remove such genes before fitting. The rule in `prepare.py`
is one documented way to do so. This is a deliberate support boundary, not a
numerical defect, so Saturna is unchanged.

**D2 — Empirical p-values differ by up to 1.1e-5 (relative) because the fitted
null differs slightly.**

Saturna's empirical null (location −0.02111513, scale 1.6490125) differs from
satuRn's (−0.02111486, 1.6490140) by 2.7e-7 and 1.4e-6. The checks show where
the difference comes from:

1. Evaluating Saturna's empirical-p formula with satuRn's fitted null reproduces
   satuRn's empirical p-values to 3.7e-14 relative. The formula is exact.
2. On satuRn's own final truncation window, the two estimates have equal
   negative log-likelihood to within 3.2e-10 (floating-point resolution). Both
   are the same optimum of a flat likelihood.
3. locfdr's `locmle`, called by satuRn, stops damped Newton iterations once the
   natural-parameter step falls below 1e-5. Saturna's independent SciPy
   L-BFGS-B optimizer uses tighter tolerances. Their stopping points differ at
   about 1e-6.
4. The relative effect on a p-value grows roughly with z² times the relative
   scale difference. Every out-of-tolerance value has |z| ≥ 1.66 (empirical p ≤
   0.097), and the largest is about 1e-5. The synthetic fixtures' tails were not
   long enough to show this.

Matching locfdr exactly would mean translating its GPL-2-only routine, which
NOTICE.md excludes. The difference is therefore documented, not removed.
Empirical FDR stays within tolerance (4.2e-7).

The regression gate keeps every other pre-registered criterion and replaces the
failed one with the D2 checks:

- all non-empirical quantities and empirical FDR within the original tolerance;
- the formula replay reproduces satuRn within the regular tolerance;
- null location and scale differ by less than 1e-5;
- both nulls reach the same likelihood within 1e-8.

The report shows both the pre-registered outcome (FAIL) and the gate outcome
(PASS).

**D3 — Empirical-call agreement is untested here.** Neither tool makes any
empirical call. The lowest empirical FDR is 0.9915 in both, because the fitted
null is wide (scale 1.65) with 7 samples and about 2 prior df. The "identical
empirical calls" criterion is therefore met trivially. An exploratory
Saturna-only run with `~ intercept + knockdown` also gave no empirical calls
(lowest 0.61), so it was not added. Testing empirical calls near a threshold
needs a larger real dataset.

## Limitations

This benchmark covers one organism, one cell line, one quantifier, an old
annotation (Ensembl 70), seven samples, and one two-group contrast with one
covariate. It does not assess ground truth, such as known pasilla-regulated
exons, or null calibration, and it does not test robust or trended moderation,
which Saturna does not implement. The filter and the uninformative-gene rule
are benchmark choices, not recommendations for every dataset.

## Reproduce

From the repository root, with R 4.6.1 and the packages listed above:

```bash
python benchmarks/pasilla/prepare.py      # download, verify SHA-256, filter
python benchmarks/pasilla/reference.py --rscript Rscript --r-library /path/to/R/library
python benchmarks/pasilla/compare.py --output reports/latest/pasilla
python -m pytest tests/test_pasilla.py
```

`prepare.py` and `reference.py` regenerate byte-identical inputs and reference
outputs; `manifest.json` records their SHA-256 hashes. `compare.py` and the
tests need only Python and the committed files.
