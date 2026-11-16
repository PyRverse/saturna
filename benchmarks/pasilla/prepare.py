"""Download, verify, and filter the pasilla transcript counts shared by both tools.

Run from the repository root. The pinned Bioconductor archive is cached under
benchmarks/pasilla/.cache and verified by SHA-256 before anything is read.
"""

import argparse
import hashlib
import io
import json
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
MIN_COUNT = 10
MIN_TRANSCRIPT_SAMPLES = 3
MIN_TRANSCRIPTS_PER_GENE = 2


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_archive(source: dict, cache: Path, archive: Path | None) -> dict[str, bytes]:
    path = archive or cache / Path(source["url"]).name
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(source["url"], headers={"User-Agent": "saturna-benchmark"})
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    data = path.read_bytes()
    if sha256(data) != source["sha256"]:
        raise RuntimeError(f"{path} does not match the pinned SHA-256 {source['sha256']}")
    members = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        for name, expected in source["members"].items():
            content = tar.extractfile(name).read()
            if sha256(content) != expected:
                raise RuntimeError(f"Archive member {name} does not match its pinned SHA-256")
            members[Path(name).name] = content
    return members


def filter_transcripts(counts: pd.DataFrame, samples: list[str]) -> pd.Series:
    """Label each annotated transcript with the step that removed it, or its analysis set.

    1. Gene total over all annotated transcripts >= 10 in every sample.
    2. Transcript count >= 10 in at least 3 samples (the smallest group size).
    3. Genes keep >= 2 transcripts passing steps 1-2.
    4. Genes whose retained transcripts share one identical count vector are
       "uninformative": the quantifier split reads equally, so usage is fixed
       by construction. They form the R-only unrestricted sensitivity set.
    """
    values = counts[samples]
    gene_totals = values.groupby(counts["gene_id"]).transform("sum")
    gene_ok = (gene_totals >= MIN_COUNT).all(axis=1)
    transcript_ok = (values >= MIN_COUNT).sum(axis=1) >= MIN_TRANSCRIPT_SAMPLES
    passing = gene_ok & transcript_ok
    retained_per_gene = passing.groupby(counts["gene_id"]).transform("sum")
    multi = passing & (retained_per_gene >= MIN_TRANSCRIPTS_PER_GENE)
    retained = counts[multi]
    distinct = retained.groupby("gene_id")[samples].apply(lambda x: len(x.drop_duplicates()))
    uninformative_genes = set(distinct.index[distinct == 1])
    label = pd.Series("gene_expression", index=counts.index)
    label[gene_ok & ~transcript_ok] = "transcript_expression"
    label[passing & ~multi] = "single_transcript_gene"
    label[multi] = "primary"
    label[multi & counts["gene_id"].isin(uninformative_genes)] = "uninformative"
    return label


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, help="Use an already downloaded tarball")
    parser.add_argument("--cache", type=Path, default=HERE / ".cache")
    args = parser.parse_args()
    source = json.loads((HERE / "source.json").read_text())
    members = load_archive(source, args.cache, args.archive)
    args.cache.mkdir(parents=True, exist_ok=True)
    for name, content in members.items():
        (args.cache / name).write_bytes(content)

    counts = pd.read_csv(io.BytesIO(members["counts.txt"]), sep="\t", index_col="feature_id")
    metadata = pd.read_csv(io.BytesIO(members["metadata.txt"]), sep="\t", index_col="SampleName")
    samples = list(metadata.index)
    if list(counts.columns) != ["gene_id", *samples]:
        raise RuntimeError("Count columns do not match the metadata sample order")
    invalid = counts.isna().any().any() or (counts[samples] < 0).any().any()
    if counts.index.has_duplicates or invalid:
        raise RuntimeError("Raw counts must be unique, complete and nonnegative")

    label = filter_transcripts(counts, samples)
    kept = label.isin(["primary", "uninformative"])
    table = counts[kept].copy()
    table.insert(1, "analysis_set", label[kept])
    table.index.name = "transcript_id"
    design = pd.DataFrame({
        "intercept": 1,
        "knockdown": (metadata["condition"] == "KD").astype(int),
        "paired_end": (metadata["LibraryLayout"] == "PAIRED").astype(int),
    }, index=pd.Index(samples, name="sample_id"))

    out = HERE / "inputs"
    out.mkdir(exist_ok=True)
    table.to_csv(out / "counts.csv", lineterminator="\n")
    design.to_csv(out / "design.csv", lineterminator="\n")

    def tally(mask):
        return {"transcripts": int(mask.sum()), "genes": int(counts.loc[mask, "gene_id"].nunique())}

    summary = {
        "rule": {
            "min_count": MIN_COUNT,
            "gene_total_min_count_in_all_samples": True,
            "transcript_min_samples": MIN_TRANSCRIPT_SAMPLES,
            "min_transcripts_per_gene": MIN_TRANSCRIPTS_PER_GENE,
            "uninformative_gene": "all retained transcripts share one identical count vector",
        },
        "annotated": tally(label.notna()),
        "removed": {step: tally(label == step) for step in
                    ("gene_expression", "transcript_expression", "single_transcript_gene")},
        "unrestricted": tally(kept),
        "uninformative": tally(label == "uninformative"),
        "primary": tally(label == "primary"),
        "design": {"columns": list(design.columns), "contrast": [0, 1, 0],
                   "samples": metadata.reset_index()[
                       ["SampleName", "LibraryName", "LibraryLayout", "condition"]
                   ].to_dict(orient="records")},
    }
    (out / "filter.json").write_text(json.dumps(summary, indent=2) + "\n",
                                     encoding="utf-8", newline="\n")
    print(json.dumps({k: summary[k] for k in ("unrestricted", "uninformative", "primary")}))


if __name__ == "__main__":
    main()
