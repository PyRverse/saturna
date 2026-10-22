"""Run a complete DTU analysis on the bundled synthetic bulk counts."""

from pathlib import Path

import pandas as pd

from saturna import fit_dtu, test_dtu

root = Path(__file__).resolve().parents[1]
counts = pd.read_csv(root / "datasets/empirical/counts.csv", index_col="transcript_id")
genes = counts.pop("gene_id")
design = pd.read_csv(root / "datasets/empirical/design.csv", index_col="sample_id")
fitted = fit_dtu(counts, genes, design)
results = test_dtu(fitted, [0, 1])
print(results.sort_values("empirical_fdr").head(8).to_string())
print("Empirical null:", results.attrs["empirical_null"])
