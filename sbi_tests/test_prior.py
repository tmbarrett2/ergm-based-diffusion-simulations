"""UT-2.1: Prior Sampling, and the Calibration Functions.

Jonathan H. Morgan, Ph.D.
26 October 2026

Revised 26 October 2026: inf_r is logit-normal (median 0.3 by default) and
b_int ~ N(0, 1); the calibration functions are checked on a hand-built scan,
without running Julia. Run from the repository root:
    python sbi_tests/test_prior.py
"""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

# Make the python_env package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.prior import (PRIOR_SPEC, THETA_FIELDS, format_scan_table, make_prior_spec,  # noqa: E402
                              read_calibration, sample_prior, shift_inf_r, suggest_inf_r_median,
                              write_calibration, write_params)


# Helper Function for the prior tests: report a dictionary of checks
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


# Helper Function for test_calibration: check that a call raises ValueError
def _raises_value_error(call):
    """
    Args:
        call (callable): function of no arguments
    Returns:
        bool: True if the call raised ValueError
    """

    # Run the call
    try:
        call()
    except ValueError:
        return True
    return False


# Helper Function for test_calibration: a hand-built scan
def toy_scan():
    """
    Args:
        None
    Returns:
        dict: a scan in the scan_inf_r format, with established shares 0.1,
            0.3, 0.7, and 0.9 at medians 0.1, 0.2, 0.4, and 0.6
    """

    # One cell per median
    cells = []
    for m, est in ((0.1, 0.1), (0.2, 0.3), (0.4, 0.7), (0.6, 0.9)):
        cells.append({"inf_r_median": m, "runs": 200, "nonstarter": 0.1, "extinct": 0.9 - est,
                      "fit_failed": 0.0, "established": est, "fit_success": est / 0.9,
                      "people_median": 10.0, "people_p90": 100.0, "people_max": 1000.0, "K_median": 0.4})

    # Assembling Result
    return {"village": "toy", "graphml": "toy.graphml", "n_nodes": 2500, "n": 200, "seed": 1,
            "inf_r_logit_sd": 0.75, "b_int_mean": 0.0, "cells": cells}


# UT-2.1: Prior Sampling
def test_prior_sampling(n=1000, seed=2026):
    """
    Args:
        n (int): number of draws (default 1000)
        seed (int): prior seed (default 2026)
    Returns:
        bool: True if every check passed
    Notes:
        Checks criteria (a) to (e) of UT-2.1, plus reproducibility. (c) now
        asks for an inf_r median near the default 0.3 with every draw in
        (0, 1). (e) asks that shape samples be concentrated above 1;
        Gamma(4, 0.5) puts 85.7% of its mass above 1, so the check requires
        at least 80%.
    """

    # Draw and arrange by parameter
    draws = sample_prior(n, seed)
    col = {f: np.array([d[f] for d in draws]) for f in THETA_FIELDS}
    checks = {}

    # (a) Physically meaningful ranges
    checks["a: inf_r in (0, 1)"] = bool(((col["inf_r"] > 0) & (col["inf_r"] < 1)).all())
    checks["a: b_close >= 0"] = bool((col["b_close"] >= 0).all())
    checks["a: p_symp in (0, 1)"] = bool(((col["p_symp"] > 0) & (col["p_symp"] < 1)).all())
    checks["a: Weibull scales and shapes > 0"] = all(bool((col[f] > 0).all()) for f in
                                                     ("rec_scale", "rec_shape", "imm_scale", "imm_shape"))
    checks["a: every value finite"] = all(bool(np.isfinite(col[f]).all()) for f in THETA_FIELDS)

    # (b) p_symp centered near 0.5
    checks["b: p_symp mean in [0.45, 0.55]"] = 0.45 <= col["p_symp"].mean() <= 0.55

    # (c) inf_r median near the prior's median
    med = np.median(col["inf_r"])
    target = PRIOR_SPEC["inf_r"][1]["median"]
    checks["c: inf_r median within 10% of the prior's"] = 0.9 * target <= med <= 1.1 * target

    # (d) Behavioral coefficients: mean near 0, sd near 1
    for f in ("b_int", "b_cxn_peers", "b_cxn_total", "b_cxn_symp"):
        checks[f"d: {f} mean ~ 0, sd ~ 1"] = abs(col[f].mean()) < 0.1 and 0.9 <= col[f].std() <= 1.1

    # (e) Weibull shapes concentrated above 1
    for f in ("rec_shape", "imm_shape"):
        checks[f"e: {f} at least 80% above 1"] = (col[f] > 1).mean() >= 0.80

    # Reproducibility and draw identifiers
    checks["same seed gives the same draws"] = sample_prior(n, seed) == draws
    checks["different seed gives different draws"] = sample_prior(n, seed + 1) != draws
    checks["draw_id runs 1 to n"] = [d["draw_id"] for d in draws] == list(range(1, n + 1))

    # The written file reads back to the same values
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "params.json"
        write_params(draws, path)
        checks["written values read back exactly"] = json.loads(path.read_text()) == draws

    # Report the prior's location for the record
    q = np.percentile(col["inf_r"], [5, 95])
    print(f"  inf_r median {med:.3f} (90% in {q[0]:.3f} to {q[1]:.3f}), p_symp mean {col['p_symp'].mean():.3f}, "
          f"rec_scale median {np.median(col['rec_scale']):.2f}, imm_scale median {np.median(col['imm_scale']):.1f}")

    # Assembling Result
    return report_checks(f"UT-2.1: Prior Sampling ({n} draws)", checks)


# Calibration Functions
def test_calibration(n=1000, seed=2026):
    """
    Args:
        n (int): number of draws (default 1000)
        seed (int): prior seed (default 2026)
    Returns:
        bool: True if every check passed
    Notes:
        A spec with another inf_r median must give draws equal to shifting
        the default draws, since both use the same normal deviates; this is
        what lets a scan stand in for the calibrated prior. The suggestion
        is checked on a hand-built scan, where the crossing of 0.5 lies
        halfway between medians 0.2 and 0.4 on the log-odds scale.
    """

    # Specs and the shift
    checks = {}
    base = sample_prior(n, seed)
    spec = make_prior_spec(inf_r_median=0.1)
    shifted = shift_inf_r(base, 0.3, 0.1)
    direct = sample_prior(n, seed, spec=spec)
    checks["make_prior_spec leaves PRIOR_SPEC unchanged"] = PRIOR_SPEC["inf_r"][1]["median"] == 0.3
    checks["spec moves only inf_r (other draws equal)"] = all(
        a[f] == b[f] for a, b in zip(base, direct) for f in THETA_FIELDS if f != "inf_r")
    checks["shift equals drawing from the spec"] = all(
        abs(a["inf_r"] - b["inf_r"]) < 1e-12 for a, b in zip(shifted, direct))
    checks["shift keeps ranks"] = bool((np.argsort([d["inf_r"] for d in base]) ==
                                        np.argsort([d["inf_r"] for d in shifted])).all())
    checks["bad median rejected"] = _raises_value_error(lambda: make_prior_spec(inf_r_median=1.2))

    # Suggestion from a hand-built scan
    scan = toy_scan()
    s = suggest_inf_r_median(scan, target=0.5)
    expected = 1.0 / (1.0 + np.exp(-(np.log(0.2 / 0.8) + np.log(0.4 / 0.6)) / 2.0))
    checks["suggestion interpolates on the log-odds scale"] = abs(s["inf_r_median"] - round(expected, 3)) < 1e-9
    checks["no suggestion below the grid"] = suggest_inf_r_median(scan, target=0.05)["inf_r_median"] is None
    checks["no suggestion above the grid"] = suggest_inf_r_median(scan, target=0.95)["inf_r_median"] is None
    checks["table shows every median and the suggestion"] = (
        format_scan_table(scan).count("\n") == 4 + len(scan["cells"]) and "Suggested" in format_scan_table(scan))

    # Record and read back
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "toy_prior.json"
        written = write_calibration(scan, 0.3, path)
        record = json.loads(path.read_text())
        checks["record keeps the scan, suggestion, and choice"] = (
            record["scan"] == scan and record["suggested"] == s["inf_r_median"] and
            record["chosen_inf_r_median"] == 0.3)
        checks["read_calibration returns the written spec"] = read_calibration(path) == written

    # Assembling Result
    return report_checks("Calibration Functions (no Julia)", checks)


if __name__ == "__main__":
    results = [test_prior_sampling(), test_calibration()]
    sys.exit(0 if all(results) else 1)
