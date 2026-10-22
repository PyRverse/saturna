# Saturna

**Differential transcript usage in native Python.** Saturna fits quasi-binomial
models to transcript count matrices, moderates dispersion across transcripts,
and tests numeric contrasts. The name is a short, pronounceable nod to satuRn.

Version **0.1.0a1** is an early implementation for **independent bulk samples**.
It reproduces a documented subset of satuRn 1.20.0 and runs without R. Genuine
R-generated fixtures cover model fitting, equal and unequal residual degrees of
freedom, covariates, moderation, and empirical-null correction. This is a port
with a bounded support contract, not a drop-in replacement for the entire R package.

```bash
python -m pip install -e ".[dev]"
```

From the repository root, analyze the bundled synthetic dataset:

```python
import pandas as pd
from saturna import fit_dtu, test_dtu

counts = pd.read_csv("datasets/empirical/counts.csv", index_col="transcript_id")
gene_ids = counts.pop("gene_id")
design = pd.read_csv("datasets/empirical/design.csv", index_col="sample_id")

fit = fit_dtu(counts, gene_ids, design)
result = test_dtu(fit, [0, 1])
print(result.sort_values("empirical_fdr").head())
print(result.attrs["empirical_null"])
```

Rows are transcripts and columns are samples. The design has one row per sample;
include an intercept column yourself. DataFrame labels must match in order.
`dtu(counts, gene_ids, design, contrast)` combines fitting and testing.

The result contains the contrast estimate (natural-log odds), standard error,
posterior degrees of freedom, moderated t statistic, regular p-value and BH FDR,
and empirical p-value and FDR. Empirical correction runs automatically with at
least 500 valid tests. Below that threshold it warns and leaves empirical columns
missing; use `empirical=False` to explicitly request regular tests only.

Supported: dense numeric counts, at least two transcripts per gene, full-rank
numeric designs, continuous covariates, one numeric contrast per call, nonrobust
moderation without a trend, and the two-stage central empirical-null likelihood.
Single-cell inference, repeated-measure correlation, robust/trend moderation,
stageR, sparse or delayed arrays, formulas, quantification, and plotting are outside
this release. See the [support matrix](docs/support.md) and [methods](docs/methods.md).

```bash
python -m pytest
python -m ruff check .
python scripts/parity_report.py --output reports/latest
python examples/bulk_dtu.py
python -m build
```

Tests compare against checked-in R outputs without requiring R. The [reference
instructions](reference/README.md) explain regeneration. The [provenance and
license audit](docs/provenance.md) identify exact upstream sources and the separate
empirical optimizer. [NOTICE.md](NOTICE.md) retains author and license notices.

Saturna is distributed under **GPL-3.0-or-later**. Please credit the original
[satuRn method and package](https://doi.org/10.18129/B9.bioc.satuRn) and
[limma](https://doi.org/10.18129/B9.bioc.limma) when using this work.
