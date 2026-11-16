"""Real-data benchmark: committed genuine satuRn output on the pasilla knockdown data."""

import importlib.util
import json
from pathlib import Path

import pytest

from saturna._glm import ModelFitError

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pasilla_compare",
                                              ROOT / "benchmarks/pasilla/compare.py")
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


@pytest.fixture(scope="module")
def primary():
    reference = compare.load_reference("primary")
    result = compare.run_saturna()
    return result, reference, compare.compare_primary(result, reference)


def test_benchmark_provenance_and_filter():
    manifest, _ = compare.verify_manifest()
    assert manifest["generator"].startswith("Actual Rscript")
    assert manifest["runtime"]["packages"]["satuRn"] == "1.20.0"
    assert manifest["runtime"]["packages"]["limma"] == "3.68.4"
    summary = json.loads((ROOT / "benchmarks/pasilla/inputs/filter.json").read_text())
    assert summary["primary"] == {"transcripts": 10003, "genes": 3196}
    assert summary["uninformative"] == {"transcripts": 104, "genes": 52}


def test_regular_outputs_and_calls_match_satuRn(primary):
    _, _, report = primary
    assert report["same_transcripts_in_order"] and report["gene_ids_match"]
    assert report["all_reference_models_glm"] and report["all_reference_values_finite"]
    for name, metrics in report["quantities"].items():
        if name != "empirical_pvalue":
            assert metrics["passed"], name
    for kind in ("regular", "empirical"):
        assert report["calls"][f"{kind}@0.05"]["identical"]
    assert report["calls"]["regular@0.05"]["saturna"] == 293


def test_empirical_difference_is_the_null_estimate(primary):
    result, reference, _ = primary
    runtime = json.loads((ROOT / "benchmarks/pasilla/reference/runtime.json").read_text())
    diagnosis = compare.empirical_null_diagnosis(
        result, reference, runtime["analyses"]["primary"]["empirical_null"])
    assert diagnosis["replay_with_satuRn_null"]["passed"]
    assert diagnosis["location_difference"] < compare.NULL_PARAMETER_TOLERANCE
    assert diagnosis["scale_difference"] < compare.NULL_PARAMETER_TOLERANCE
    nll = diagnosis["negative_log_likelihood_on_satuRn_window"]
    assert abs(nll["saturna"] - nll["satuRn"]) < 1e-8


def test_unrestricted_set_is_rejected_as_documented():
    with pytest.raises(ModelFitError, match="effectively zero"):
        compare.fit_dtu(*compare.load_inputs())
