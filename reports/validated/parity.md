# Saturna numerical parity

Result: **PASS**

Reference: actual R 4.6.1 / satuRn 1.20.0 / limma 3.68.4.

| Quantity | Pass | Max absolute error | rtol | atol |
| --- | --- | ---: | ---: | ---: |
| balanced/coefficients | True | 5.11e-15 | 2e-08 | 2e-10 |
| balanced/covariance_unscaled | True | 4.51e-17 | 2e-08 | 2e-10 |
| balanced/dispersion | True | 6.75e-14 | 2e-08 | 2e-10 |
| balanced/residual_df | True | 0 | 2e-08 | 2e-10 |
| balanced/variance_posterior | True | 3.5e-12 | 2e-08 | 2e-10 |
| balanced/df_posterior | True | 9.54e-13 | 2e-08 | 2e-10 |
| balanced/estimate | True | 5.11e-15 | 2e-08 | 2e-10 |
| balanced/standard_error | True | 2.94e-14 | 2e-08 | 2e-10 |
| balanced/df | True | 9.54e-13 | 2e-08 | 2e-10 |
| balanced/statistic | True | 9.09e-13 | 2e-08 | 2e-10 |
| balanced/pvalue | True | 7.85e-14 | 2e-08 | 2e-10 |
| balanced/fdr | True | 2.48e-13 | 2e-08 | 2e-10 |
| unequal_df/coefficients | True | 4.22e-15 | 2e-08 | 2e-10 |
| unequal_df/covariance_unscaled | True | 5.55e-17 | 2e-08 | 2e-10 |
| unequal_df/dispersion | True | 4.26e-14 | 2e-08 | 2e-10 |
| unequal_df/residual_df | True | 0 | 2e-08 | 2e-10 |
| unequal_df/variance_posterior | True | 9.06e-12 | 2e-08 | 2e-10 |
| unequal_df/df_posterior | True | 8.95e-12 | 2e-08 | 2e-10 |
| unequal_df/estimate | True | 1.19e-15 | 2e-08 | 2e-10 |
| unequal_df/standard_error | True | 1.57e-13 | 2e-08 | 2e-10 |
| unequal_df/df | True | 8.95e-12 | 2e-08 | 2e-10 |
| unequal_df/statistic | True | 3.07e-12 | 2e-08 | 2e-10 |
| unequal_df/pvalue | True | 5.02e-13 | 2e-08 | 2e-10 |
| unequal_df/fdr | True | 2.98e-13 | 2e-08 | 2e-10 |
| empirical/coefficients | True | 6.44e-15 | 2e-08 | 2e-10 |
| empirical/covariance_unscaled | True | 2.5e-16 | 2e-08 | 2e-10 |
| empirical/dispersion | True | 9.59e-14 | 2e-08 | 2e-10 |
| empirical/residual_df | True | 0 | 2e-08 | 2e-10 |
| empirical/variance_posterior | True | 6.36e-13 | 2e-08 | 2e-10 |
| empirical/df_posterior | True | 4.26e-14 | 2e-08 | 2e-10 |
| empirical/estimate | True | 5.33e-15 | 2e-08 | 2e-10 |
| empirical/standard_error | True | 6.74e-14 | 2e-08 | 2e-10 |
| empirical/df | True | 4.26e-14 | 2e-08 | 2e-10 |
| empirical/statistic | True | 1.94e-12 | 2e-08 | 2e-10 |
| empirical/pvalue | True | 1e-13 | 2e-08 | 2e-10 |
| empirical/fdr | True | 3.47e-13 | 2e-08 | 2e-10 |
| empirical/empirical_pvalue | True | 5.96e-09 | 3e-06 | 3e-08 |
| empirical/empirical_fdr | True | 2.98e-08 | 3e-06 | 3e-08 |
| covariate/coefficients | True | 4.88e-15 | 2e-08 | 2e-10 |
| covariate/covariance_unscaled | True | 5.9e-17 | 2e-08 | 2e-10 |
| covariate/dispersion | True | 5.68e-14 | 2e-08 | 2e-10 |
| covariate/residual_df | True | 0 | 2e-08 | 2e-10 |
| covariate/variance_posterior | True | 3.73e-14 | 2e-08 | 2e-10 |
| covariate/df_posterior | True | 6.04e-14 | 2e-08 | 2e-10 |
| covariate/estimate | True | 3.77e-15 | 2e-08 | 2e-10 |
| covariate/standard_error | True | 5e-16 | 2e-08 | 2e-10 |
| covariate/df | True | 6.04e-14 | 2e-08 | 2e-10 |
| covariate/statistic | True | 5.33e-15 | 2e-08 | 2e-10 |
| covariate/pvalue | True | 2.55e-15 | 2e-08 | 2e-10 |
| covariate/fdr | True | 2.55e-15 | 2e-08 | 2e-10 |
| shifted_null/empirical_pvalue | True | 1.15e-08 | 3e-06 | 3e-08 |

Synthetic numerical checks establish only the documented supported subset.
Empirical calibration uses an independently written SciPy optimizer.
