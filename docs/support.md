# Supported inputs and intentional boundaries

| Area | In 0.1.0a1 | Boundary |
| --- | --- | --- |
| Data | Dense transcripts × independent bulk samples | No single-cell claim, sparse matrices, quantifier import, or pseudo-replication |
| Counts | Finite nonnegative numeric values, including fractional estimated counts | Missing, boolean, complex, negative and infinite values rejected |
| Precision | Float64 arithmetic; gene totals at most 2^52 | Unrepresentable integer inputs rejected before casting mixed lists |
| Mapping | One nonmissing gene ID per transcript; at least two transcripts per gene | Single-isoform genes must be filtered by the caller |
| Design | Finite full-rank matrix; numeric covariates and explicit intercept | No formula parser, offsets, weights, random effects, or block correlation |
| Empty gene samples | Removed per gene before fitting, as in satuRn | Retained design must remain full rank with positive residual df |
| Model failure | Actionable exception identifying the transcript | No silent omission or p-values for failed models |
| Degeneracy | Dispersion must exceed 1e-12; condition number at most 1e8 | Perfectly explained/nearly singular models rejected, unlike upstream's permissive behavior |
| Moderation | limma nonrobust, no covariate trend; equal and unequal positive df | Robust, trend and missing-feature moderation not implemented |
| Contrasts | One finite nonzero numeric vector per test_dtu call | Multi-df omnibus tests not implemented; call repeatedly for multiple contrasts |
| Empirical null | Two-stage truncated-normal likelihood, >=500 valid tests | Independent optimizer; degeneracy/nonconvergence raises an error |
| Multiplicity | BH across transcripts within each contrast | No stageR or gene-level confirmatory inference; no correction across contrasts |
| Output | pandas DataFrame in input transcript order | No plotDTU equivalent or R container interoperability |

Empirical correction assumes a usable central null population; 500 is the
upstream eligibility rule, not a guarantee that a dataset satisfies the statistical
assumptions. `empirical=True` requires eligibility and a successful estimate;
`"auto"` skips with a warning only below 500 tests; `False` explicitly disables it.
The optimizer rejects estimates at its numerical bounds (location ±20, standard
deviation 0.02–20) and central samples without spread. It does not silently fall
back to a nominal null when estimation fails.

Input DataFrame sample labels and gene Series indices must match exactly in
order. No automatic sorting or alignment occurs. String conversion must retain
distinct gene IDs, and transcript/coefficient names must remain unique.

Real quantifier output can contain genes whose transcripts received identical
counts in every sample, because the quantifier could not distinguish them.
Their usage is constant, so Saturna rejects them as zero-dispersion models.
Remove such genes first; benchmarks/pasilla/prepare.py shows one rule.
satuRn instead fits them with zero dispersion, and those zeros lower the
moderation prior for every other transcript (see docs/realdata.md, D1).

This alpha has synthetic numerical parity evidence and one independent-bulk
real-data benchmark (pasilla, docs/realdata.md). On that benchmark, regular
outputs and calls match satuRn; empirical p-values differ by up to 1.1e-5
relative. This is not comprehensive biological validation, a calibration
study, or a performance claim. Larger and more diverse real datasets
(especially ones with empirical calls), single-cell validation, broader
diagnostic coverage, and upstream-complete parity remain future work.
