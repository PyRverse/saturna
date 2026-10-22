"""Produce strict JSON and Markdown evidence against the genuine R fixtures."""

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

import saturna
from saturna import fit_dtu, test_dtu
from saturna._empirical import empirical_pvalues, estimate_empirical_null


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("reports/latest"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest_path = root / "reference/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for relative, digest in manifest["sha256"].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"Reference provenance mismatch: {relative}")
    comparisons = []

    def compare(name, actual, expected, *, empirical=False):
        a, e = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        rtol, atol = (3e-6, 3e-8) if empirical else (2e-8, 2e-10)
        passed = bool(a.shape == e.shape and np.allclose(a, e, rtol=rtol, atol=atol))
        comparisons.append({"name": name, "passed": passed, "elements": int(e.size),
                            "rtol": rtol, "atol": atol,
                            "max_absolute_error": float(np.max(np.abs(a - e)))})

    for scenario in ("balanced", "unequal_df", "empirical", "covariate"):
        reference = json.loads((root / "reference/outputs" / f"{scenario}.json").read_text())
        counts = pd.read_csv(root / "datasets" / scenario / "counts.csv", index_col=0)
        genes = counts.pop("gene_id")
        design = pd.read_csv(root / "datasets" / scenario / "design.csv", index_col=0)
        fit = fit_dtu(counts, genes, design)
        for name in ("coefficients", "covariance_unscaled", "dispersion", "residual_df",
                     "variance_posterior", "df_posterior"):
            compare(f"{scenario}/{name}", getattr(fit, name), reference[name])
        result = test_dtu(fit, reference["contrast"], empirical=scenario == "empirical")
        columns = {"estimate": "estimates", "standard_error": "se", "df": "df",
                   "statistic": "t", "pvalue": "pval", "fdr": "regular_FDR"}
        if scenario == "empirical":
            columns.update(empirical_pvalue="empirical_pval", empirical_fdr="empirical_FDR")
        for name, r_name in columns.items():
            compare(f"{scenario}/{name}", result[name], reference["result"][r_name],
                    empirical=name.startswith("empirical_"))
    reference = json.loads((root / "reference/outputs/empirical_null.json").read_text())
    fit = estimate_empirical_null(reference["pvalues"], reference["statistics"])
    compare("shifted_null/empirical_pvalue",
            empirical_pvalues(reference["pvalues"], reference["statistics"], fit),
            reference["result"]["pval"], empirical=True)
    passed = all(x["passed"] for x in comparisons)
    report = {"schema_version": 1, "passed": passed, "reference_runtime": manifest["runtime"],
              "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
              "python_runtime": {"python": platform.python_version(),
                                 "saturna": saturna.__version__,
                                 "numpy": np.__version__, "scipy": scipy.__version__,
                                 "pandas": pd.__version__}, "comparisons": comparisons}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "parity.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                            encoding="utf-8", newline="\n")
    lines = ["# Saturna numerical parity", "", f"Result: **{'PASS' if passed else 'FAIL'}**",
             "", "Reference: actual R 4.6.1 / satuRn 1.20.0 / limma 3.68.4.", "",
             "| Quantity | Pass | Max absolute error | rtol | atol |",
             "| --- | --- | ---: | ---: | ---: |"]
    for row in comparisons:
        lines.append(f"| {row['name']} | {row['passed']} | {row['max_absolute_error']:.3g} | "
                     f"{row['rtol']:.3g} | {row['atol']:.3g} |")
    lines += ["", "Synthetic numerical checks establish only the documented supported subset.",
              "Empirical calibration uses an independently written SciPy optimizer.", ""]
    (args.output / "parity.md").write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"{'PASS' if passed else 'FAIL'}: {len(comparisons)} comparisons; {args.output}")
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
