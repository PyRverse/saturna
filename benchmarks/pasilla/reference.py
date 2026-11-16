"""Run the genuine satuRn reference on the prepared pasilla inputs and hash provenance.

Run prepare.py first; it caches the verified raw files that satuRn.R re-filters.
"""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HASHED = ["source.json", "prepare.py", "satuRn.R", "reference.py",
          "inputs/counts.csv", "inputs/design.csv", "inputs/filter.json",
          "reference/primary.csv", "reference/unrestricted.csv", "reference/runtime.json"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--r-library", action="append", default=[])
    parser.add_argument("--cache", type=Path, default=HERE / ".cache")
    args = parser.parse_args()
    env = os.environ.copy()
    if args.r_library:
        env["R_LIBS"] = os.pathsep.join(args.r_library)
    env["LC_ALL"] = "C"
    raw = [str((args.cache / name).resolve()) for name in ("counts.txt", "metadata.txt")]
    subprocess.run([args.rscript, "--vanilla", "benchmarks/pasilla/satuRn.R", *raw],
                   cwd=ROOT, env=env, check=True)
    for path in (HERE / "reference").iterdir():
        path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n"))
    runtime = json.loads((HERE / "reference/runtime.json").read_text())
    manifest = {
        "schema_version": 1,
        "generator": "Actual Rscript execution of benchmarks/pasilla/satuRn.R",
        "data_origin": "PasillaTranscriptExpr 1.40.0 kallisto counts; see source.json",
        "runtime": {"R": runtime["R"], "packages": runtime["packages"]},
        "sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                   for name in HASHED},
    }
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n",
                                        encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
