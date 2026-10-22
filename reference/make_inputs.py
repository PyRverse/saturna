"""Create original synthetic inputs; no expected answers are produced here."""

from pathlib import Path

import numpy as np
import pandas as pd


def create(name: str, genes: int, samples: int, seed: int, missing: bool = False,
           covariate: bool = False) -> None:
    rng = np.random.default_rng(seed)
    group = np.repeat([0, 1], samples // 2)
    rows, identifiers, mapping = [], [], []
    for gene in range(genes):
        baseline = rng.dirichlet([3, 4, 3])
        concentration = float(np.exp(rng.uniform(np.log(8), np.log(500))))
        block = np.empty((3, samples), dtype=int)
        for sample in range(samples):
            proportions = baseline.copy()
            if gene % 9 == 0 and group[sample]:
                proportions *= [2.2, 0.55, 1]
                proportions /= proportions.sum()
            probabilities = rng.dirichlet(proportions * concentration)
            block[:, sample] = rng.multinomial(int(rng.integers(80, 350)), probabilities)
        if missing and gene % 3:
            block[:, [0, samples // 2][:gene % 3]] = 0
        rows.extend(block)
        identifiers.extend(f"tx{gene:04d}_{i}" for i in range(3))
        mapping.extend([f"gene{gene:04d}"] * 3)
    output = Path(__file__).resolve().parents[1] / "datasets" / name
    output.mkdir(parents=True, exist_ok=True)
    sample_ids = [f"sample{i:02d}" for i in range(samples)]
    counts = pd.DataFrame(rows, index=identifiers, columns=sample_ids)
    counts.index.name = "transcript_id"
    counts.insert(0, "gene_id", mapping)
    counts.to_csv(output / "counts.csv", lineterminator="\n")
    design = pd.DataFrame({"intercept": 1, "condition": group}, index=pd.Index(
        sample_ids, name="sample_id"
    ))
    if covariate:
        design["covariate"] = rng.normal(size=samples)
    design.to_csv(output / "design.csv", lineterminator="\n")


if __name__ == "__main__":
    create("balanced", 8, 8, 20261022)
    create("unequal_df", 8, 10, 27019, missing=True)
    create("empirical", 240, 12, 2026240)
    create("covariate", 8, 12, 912467, covariate=True)
