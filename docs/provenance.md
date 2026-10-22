# Source and reference provenance

The reference runtime was actually executed on 2026-10-22 with R 4.6.1,
satuRn 1.20.0, limma 3.68.4, locfdr 1.1.8, and jsonlite 2.0.0. Exact runtime
metadata is in reference/outputs/runtime.json. Reference outputs are computed
by R, not relabeled Python outputs. All inputs are original synthetic data with
deterministic seeds in reference/make_inputs.py; no patient or private data is used.

Primary source audit:

- [satuRn 1.20.0 fitDTU source](https://github.com/bioc/satuRn/blob/51f8364f4f1a0ca606bd5f7c6ad4caac47d3cfae/R/fitDTU.R)
- [satuRn 1.20.0 testDTU source](https://github.com/bioc/satuRn/blob/51f8364f4f1a0ca606bd5f7c6ad4caac47d3cfae/R/testDTU.R)
- [Reference limma 3.68.4 source](https://github.com/bioc/limma/tree/7c11dae9f0568f17d879280450e557f8851743d5/R)
- [Audited limma 3.68.5 source](https://github.com/bioc/limma/tree/825d1c8330d024adb6c54af66d79b0f4ec258530/R)
- [locfdr upstream package and license declaration](https://cran.r-project.org/package=locfdr)
- [R GLM documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/glm.html)

The limma adaptation was written from the audited 3.68.5 source and validated
against the installed 3.68.4 runtime; both revisions are retained rather than
misidentifying the runtime as 3.68.5. locfdr's GPL-2-only declaration was treated
as a boundary: copied Newton/covariance routines were excluded. The empirical
estimator is a newly written optimization of the mathematical likelihood target;
the differing optimizer and tolerances are disclosed. See NOTICE.md for full
attribution and the accompanying GPL/Artistic licenses.

reference/manifest.json records SHA-256 hashes for inputs, generator scripts,
and each R output. The generated files use strict JSON: R missing values and
infinities use strings such as "NA" and "Inf", never bare non-JSON tokens.
The Python manifest and parity report use allow_nan=False. The test suite checks
every manifest hash, and generated reports identify the manifest by its hash.
Machine-specific executable/library paths are supplied as command arguments and
are not recorded in the portable provenance.

Covered inputs: 24 balanced transcripts, 24 transcripts with gene-specific zero
samples/unequal residual df, 24 transcripts with a continuous covariate and a
multi-coefficient contrast, and 720 transcripts with empirical correction. A
separate 950-statistic shifted/widened null tests empirical correction. Five
direct limma fixtures cover finite/infinite prior df, zeros and two-feature behavior.
