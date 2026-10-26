"""Python-Julia Bridge Tests: the Link Through sbi_batch (UT-2.2).

Jonathan H. Morgan, Ph.D.
26 October 2026

Python writes the draws, calls the compiled binary with the SAS network
option, and reads the records. Run from the repository root:
    python sbi_tests/test_bridge.py
"""

import json
import sys
from pathlib import Path

import numpy as np

# Make the python_env package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.julia_bridge import (read_manifest, read_records, records_to_arrays,  # noqa: E402
                                     run_sbi_batch)
from python_env.prior import THETA_FIELDS, sample_prior  # noqa: E402

# Paths on the EC2 container
GRAPHML = "/workspace/ergm_data/sim_test_data/WeakCore1_3.2_2.graphml"
FIXTURE = REPO / "sbi_tests" / "params_batch_ut14.json"
OUT_DIR = REPO / "sbi_tests" / "output" / "python"
BINARY_REFERENCE = Path("/tmp/bin_sas.jsonl")


# Helper Function for the bridge tests: report a dictionary of checks
def report_checks(title, checks):
    """
    Args:
        title (str): heading of the report
        checks (dict): check name -> bool
    Returns:
        bool: True if every check passed
    Notes:
        Prints the checks in name order.
    """

    # Print the report
    print("=" * 60)
    print(title)
    print("=" * 60)
    for name in sorted(checks):
        print(f"  {name:<50} {'PASS' if checks[name] else 'FAIL'}")

    # Assembling Result
    all_passed = all(checks.values())
    print(f"RESULT: {'ALL CHECKS PASSED' if all_passed else 'SOME CHECKS FAILED'}")
    return all_passed


# Helper Function for the bridge tests: theta in every record equals its draw
def theta_round_trip(records, draws):
    """
    Args:
        records (list of dict): from read_records
        draws (list of dict): the draws that produced them
    Returns:
        bool: True if each record's theta equals its input draw exactly
    Notes:
        Compared as Python floats, so a value that changed in its last bit
        anywhere between Python and Julia fails. None (null) must stay None.
    """

    # Match by draw_id and compare all 12 parameters
    by_id = {d["draw_id"]: d for d in draws}
    for rec in records:
        draw = by_id[rec["meta"]["draw_id"]]
        for f in THETA_FIELDS:
            if rec["theta"][f] != draw.get(f):
                return False

    # Assembling Result
    return True


# The Fixture Batch Through Python
def test_fixture_batch():
    """
    Args:
        None
    Returns:
        bool: True if every check passed
    Notes:
        The same draws, seed, network option, and target as the binary run
        written to /tmp/bin_sas.jsonl. If that file exists, the records must
        be identical byte for byte: the same binary, the same inputs, with
        the parameter file written by Python instead of by hand. The test
        network is connected, so seeds from the largest component are the
        same as seeds from any component.
    """

    # Run the batch
    draws = json.loads(FIXTURE.read_text())
    out_file = OUT_DIR / "fixture_sas.jsonl"
    res = run_sbi_batch(GRAPHML, draws, out_file, "WeakCore1_3.2_2", 42,
                        sas_transformation=True, target_region={"K_max": 1.0})
    records = read_records(out_file)
    manifest = read_manifest(res["manifest_file"])
    checks = {}

    # Records, manifest, and parameters
    checks["20 records"] = len(records) == 20
    checks["manifest: SAS network, seed 42"] = manifest["sas_transformation"] is True and manifest["master_seed"] == 42
    checks["manifest: seeds from the largest component"] = manifest["seed_component"] == "largest"
    checks["theta equals the input exactly (UT-2.2)"] = theta_round_trip(records, draws)
    checks["draw 10 keeps permanent immunity (null)"] = all(
        r["theta"]["imm_scale"] is None for r in records if r["meta"]["draw_id"] == 10)

    # in_target: null unless established
    checks["in_target None unless established"] = all(
        (isinstance(r["target"]["in_target"], bool)) == (r["summary"]["class"] == "established") for r in records)

    # Arrays for training
    arr = records_to_arrays(records)
    n_est = sum(r["summary"]["class"] == "established" for r in records)
    checks["arrays: one row per established run"] = arr["thetas"].shape == (n_est, 12) and arr["summaries"].shape[0] == n_est
    checks["arrays: summaries finite"] = bool(np.isfinite(arr["summaries"]).all())
    checks["arrays: every established run in target"] = bool((arr["in_target"] == 1.0).all())

    # Identical to the binary run made from the shell
    if BINARY_REFERENCE.exists():
        checks["identical to /tmp/bin_sas.jsonl"] = out_file.read_bytes() == BINARY_REFERENCE.read_bytes()
    else:
        print("  (no /tmp/bin_sas.jsonl; byte comparison skipped)")

    # Report the batch output
    print("  " + res["stdout"].strip().replace("\n", "\n  "))

    # Assembling Result
    return report_checks("Fixture Batch Through Python (SAS network)", checks)


# Prior Draws Through Python
def test_prior_batch(n=100, seed=2026):
    """
    Args:
        n (int): number of prior draws (default 100)
        seed (int): prior and batch seed (default 2026)
    Returns:
        bool: True if every check passed
    Notes:
        The full-precision round trip of prior draws is the real test of
        UT-2.2: every parameter must come back from Julia unchanged. The run
        classes preview the prior predictive check (UT-1.5).
    """

    # Draw and run
    draws = sample_prior(n, seed)
    out_file = OUT_DIR / "prior_sas.jsonl"
    run_sbi_batch(GRAPHML, draws, out_file, "WeakCore1_3.2_2", seed,
                  sas_transformation=True, target_region={"K_max": 1.0})
    records = read_records(out_file)
    checks = {}

    # Round trip and record count
    checks[f"{n} records"] = len(records) == n
    checks["prior theta round-trips exactly (UT-2.2)"] = theta_round_trip(records, draws)

    # Run classes and established summaries, for the record
    classes = {}
    for r in records:
        classes[r["summary"]["class"]] = classes.get(r["summary"]["class"], 0) + 1
    print("  Classes: " + ", ".join(f"{k}={v}" for k, v in sorted(classes.items())))
    est = records_to_arrays(records, summary_fields=("K", "t_half", "a", "frac_at_inflection"))
    if est["summaries"].shape[0] > 0:
        med = np.median(est["summaries"], axis=0)
        print(f"  Established runs: median K {med[0]:.3f}, t_half {med[1]:.1f}, a {med[2]:.2f}, "
              f"frac_at_inflection {med[3]:.2f}")

    # Assembling Result
    return report_checks(f"Prior Draws Through Python ({n} draws, SAS network)", checks)


# A Failed Batch Raises an Error
def test_failure():
    """
    Args:
        None
    Returns:
        bool: True if every check passed
    Notes:
        A missing network file must raise RuntimeError with Julia's message,
        and an unknown seed_component must raise ValueError before Julia runs.
    """

    # Run with a path that does not exist
    checks = {}
    try:
        run_sbi_batch("/no/such/file.graphml", sample_prior(2, 1), OUT_DIR / "fail.jsonl", "x", 1)
        checks["missing network raises RuntimeError"] = False
    except RuntimeError as e:
        checks["missing network raises RuntimeError"] = "GraphML path not found" in str(e)

    # Unknown seed component
    try:
        run_sbi_batch(GRAPHML, sample_prior(2, 1), OUT_DIR / "fail.jsonl", "x", 1, seed_component="all")
        checks["unknown seed_component raises ValueError"] = False
    except ValueError:
        checks["unknown seed_component raises ValueError"] = True

    # Assembling Result
    return report_checks("A Failed Batch Raises an Error", checks)


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {"Fixture batch": test_fixture_batch(),
               "Prior batch": test_prior_batch(),
               "Failure handling": test_failure()}
    print("=" * 60)
    print("Bridge Test Summary")
    print("=" * 60)
    for name, passed in results.items():
        print(f"  {name:<24} {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if all(results.values()) else 1)
