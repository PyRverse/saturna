# Numerical method

For transcript t in gene g and sample j, Saturna models the paired counts
`(count[t,j] + 1, sum_other_transcripts[g,j] + 1)` using a logit-link
quasi-binomial GLM. Samples with original gene total zero are omitted before
the one-count additions. An explicit numeric design defines the linear predictor.
No library-size normalization precedes this within-gene usage model.

Grouped-binomial IRLS uses R's initial mean `(successes + 0.5)/(total + 1)`
and relative deviance convergence rule with tolerance 1e-8. Covariance is the
inverse last weighted design crossproduct, evaluated through a QR factor rather
than normal equations. Contrasts are normalized before inference; triangular
solves avoid cancellation in their variances. Pearson dispersion uses the final
working residuals and last stored IRLS weights, matching R's GLM conventions.
Saturna defaults to 100 iterations rather than R's 25; it reports nonconvergence
as an exception rather than returning a questionable fitted object.

Across transcripts, nonrobust limma variance moderation estimates a scaled
inverse-chi-squared prior. Equal residual df uses the adjusted-log-variance
moment estimator; unequal df uses limma's newer weighted likelihood estimator.
The posterior combines residual and prior sums of squares. Infinite prior df
uses the common prior variance; fewer than three transcripts stay unmoderated.

A numeric contrast c has estimate c′β, standard error
sqrt(c′ covariance_unscaled c × posterior variance), and a t statistic with
residual + prior df. Regular p-values are two-sided and adjusted by BH within
the tested contrast. Coefficients describe natural-log odds, not log2 expression
fold changes or raw usage differences.

The empirical path follows satuRn's sign convention
`z = qnorm(p / 2) * sign(t)` and excludes nonfinite z and |z| >= 10 from fitting.
It first estimates a truncated normal in a median/IQR-based window, then refits
using a window centered/scaled by the first estimate. The width factor is
`4.3 * exp(-0.26 * log10(N))`. This implementation minimizes the conditional
normal negative log likelihood with analytic derivatives using SciPy L-BFGS-B.
It does not implement locfdr's Newton/covariance routines or diagnostic plots.
Empirical p-values use the resulting centered/scaled z values and normal tails.

The checked-in synthetic comparisons permit rtol=2e-8 and atol=2e-10 for fit and
regular-test quantities; empirical p-values/FDR permit rtol=3e-6 and atol=3e-8
because a different optimizer is used. These tolerances are regression gates,
not a proof of agreement on every dataset. The report records actual differences.
