# Saturna real-data benchmark: pasilla knockdown

Pre-registered criteria: **FAIL** (1 of 6 failed: Every quantity within tolerance)

Regression gate (pre-registered criteria with documented discrepancy checks): **PASS**

Reference: actual R 4.6.1 / satuRn 1.20.0 / limma 3.68.4 / locfdr 1.1.8.
Data: PasillaTranscriptExpr 1.40.0 kallisto counts, 10003 transcripts in 3196 genes, 4 control and 3 knockdown samples. Design `~ 0 + intercept + knockdown + paired_end`; contrast `knockdown`.

## Pre-registered criteria

- Identical transcript set and order: yes
- Every satuRn model is a finite glm fit: yes
- Every quantity within tolerance: NO
- Identical regular calls at FDR 0.05: yes
- Identical empirical calls at FDR 0.05: yes
- Saturna still rejects the unrestricted set (documented): yes

## Documented discrepancy checks (post hoc)

- D1: rejection is a zero-dispersion error on an uninformative transcript: yes
- D2: every quantity except empirical_pvalue within pre-registered tolerance: yes
- D2: Saturna's formula with satuRn's null reproduces satuRn's empirical p-values: yes
- D2: null location and scale differ by less than 1e-5: yes
- D2: both nulls reach the same likelihood on satuRn's final window (within 1e-8): yes

## Per-transcript agreement

| Quantity | Pass | Outside tol. | Max abs. error | Max rel. error | Spearman | rtol | atol |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| coef_intercept | True | 0 | 5.19e-13 | 4.13e-13 | 0.999999999988 | 2e-08 | 2e-10 |
| coef_knockdown | True | 0 | 8.63e-13 | 5.38e-12 | 0.999999999982 | 2e-08 | 2e-10 |
| coef_paired_end | True | 0 | 1.81e-13 | 7.91e-11 | 0.999999999985 | 2e-08 | 2e-10 |
| dispersion | True | 0 | 1.84e-10 | 1.2e-12 | 0.999999992782 | 2e-08 | 2e-10 |
| residual_df | True | 0 | 0 | 0 | 1.000000000000 | 2e-08 | 2e-10 |
| variance_posterior | True | 0 | 2.86e-09 | 1.31e-12 | 0.999999993250 | 2e-08 | 2e-10 |
| estimate | True | 0 | 8.63e-13 | 5.38e-12 | 0.999999999982 | 2e-08 | 2e-10 |
| standard_error | True | 0 | 1.01e-12 | 6.48e-13 | 0.999999993004 | 2e-08 | 2e-10 |
| df | True | 0 | 1.61e-12 | 3.12e-13 | 1.000000000000 | 2e-08 | 2e-10 |
| statistic | True | 0 | 9.36e-12 | 5.3e-12 | 0.999999999985 | 2e-08 | 2e-10 |
| pvalue | True | 0 | 2.86e-13 | 3.39e-12 | 0.999999992824 | 2e-08 | 2e-10 |
| fdr | True | 0 | 4.57e-13 | 3.39e-12 | 1.000000000000 | 2e-08 | 2e-10 |
| empirical_pvalue | False | 364 | 5.01e-07 | 1.11e-05 | 0.999999999940 | 3e-06 | 3e-08 |
| empirical_fdr | True | 0 | 4.17e-07 | 4.2e-07 | 1.000000000000 | 3e-06 | 3e-08 |

Maximum absolute log10 p-value difference: regular 1.47e-12, empirical 4.8e-06.

## Calls (FDR below threshold)

| FDR type @ threshold | Saturna | satuRn | Both | Saturna only | satuRn only | Identical |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| regular@0.05 | 293 | 293 | 293 | 0 | 0 | True |
| empirical@0.05 | 0 | 0 | 0 | 0 | 0 | True |
| regular@0.01 | 75 | 75 | 75 | 0 | 0 | True |
| empirical@0.01 | 0 | 0 | 0 | 0 | 0 | True |
| regular@0.1 | 539 | 539 | 539 | 0 | 0 | True |
| empirical@0.1 | 0 | 0 | 0 | 0 | 0 | True |

Closest FDR to 0.05 in either tool: regular 0.000359, empirical 0.941.

## Moderation and empirical null

| Quantity | Saturna | satuRn |
| --- | ---: | ---: |
| prior df | 2.16426908591 | 2.1642690859 |
| prior variance | 18.6080341182 | 18.6080341182 |
| empirical null location | -0.0211151348706 | -0.0211148617444 |
| empirical null scale | 1.64901253985 | 1.64901396278 |
| empirical null observations | 10003 | 10003 |

Empirical-null differences: location 2.73e-07, scale 1.42e-06. Saturna's empirical-p formula evaluated with satuRn's fitted null reproduces satuRn's empirical p-values with maximum relative error 3.75e-14. Negative log-likelihood on satuRn's final truncation window: Saturna 6119.6542755595046, satuRn 6119.6542755598293 (difference 3.25e-10, at floating-point resolution).

## Unrestricted set (documented discrepancy)

Saturna rejected the 10107-transcript set: True (`Transcript 'FBtr0077941': residual dispersion is effectively zero; a perfectly explained model cannot support this inferential workflow`).
satuRn fitted all 104 uninformative transcripts with maximum dispersion 0, maximum absolute estimate 0, minimum p-value 1, and 0 regular / 0 empirical calls at FDR 0.05.
Including them changed satuRn's prior df from 2.16427 to 2.14222 and its prior variance from 18.608 to 13.9195. On the 10003 shared transcripts the maximum p-value change was 0.0639; regular calls at 0.05 went from 293 to 373 (80 gained, 0 lost), empirical calls from 0 to 0.

See docs/realdata.md for interpretation of every discrepancy.
