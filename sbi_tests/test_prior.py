"""UT-2.1: Prior Sampling.

Jonathan H. Morgan, Ph.D.
26 October 2026

Run from the repository root:
    python sbi_tests/test_prior.py
"""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

#   Make the sbi_python package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.prior import THETA_FIELDS, sample_prior, write_params  # noqa: E402


#   Helper for the tests: report a dictionary of checks
def report_checks(title, checks):
    """Print each check in name order and return True if all passed."""

    #   Print the report
    print("=" * 60)
    print(title)
    print("=" * 60)
    for name in sorted(checks):
        print(f"  {name:<46} {'PASS' if checks[name] else 'FAIL'}")
    all_passed = all(checks.values())
    print(f"RESULT: {'ALL CHECKS PASSED' if all_passed else 'SOME CHECKS FAILED'}")
    return all_passed


#   UT-2.1: Prior Sampling
def test_prior_sampling(n=1000, seed=2026):
    """Check criteria (a) to (e) of UT-2.1 on n draws, plus reproducibility.

    Notes
    -----
    (e) asks that shape samples be concentrated above 1. Gamma(4, 0.5)
    puts 85.7% of its mass above 1, so the check requires at least 80%.
    """

    #   Draw and arrange by parameter
    draws = sample_prior(n, seed)
    col = {f: np.array([d[f] for d in draws]) for f in THETA_FIELDS}
    checks = {}

    #   (a) Physically meaningful ranges
    checks["a: inf_r > 0"] = bool((col["inf_r"] > 0).all())
    checks["a: b_close >= 0"] = bool((col["b_close"] >= 0).all())
    checks["a: p_symp in (0, 1)"] = bool(((col["p_symp"] > 0) & (col["p_symp"] < 1)).all())
    checks["a: Weibull scales and shapes > 0"] = all(bool((col[f] > 0).all()) for f in
                                                     ("rec_scale", "rec_shape", "imm_scale", "imm_shape"))
    checks["a: every value finite"] = all(bool(np.isfinite(col[f]).all()) for f in THETA_FIELDS)

    #   (b) p_symp centered near 0.5
    checks["b: p_symp mean in [0.45, 0.55]"] = 0.45 <= col["p_symp"].mean() <= 0.55

    #   (c) inf_r right-skewed with median near 0.02
    med = np.median(col["inf_r"])
    checks["c: inf_r median within 10% of 0.02"] = 0.018 <= med <= 0.022
    checks["c: inf_r right-skewed (mean > median)"] = col["inf_r"].mean() > med

    #   (d) Behavioral coefficients: mean near 0, sd near 1
    for f in ("b_cxn_peers", "b_cxn_total", "b_cxn_symp"):
        checks[f"d: {f} mean ~ 0, sd ~ 1"] = abs(col[f].mean()) < 0.1 and 0.9 <= col[f].std() <= 1.1

    #   (e) Weibull shapes concentrated above 1
    for f in ("rec_shape", "imm_shape"):
        checks[f"e: {f} at least 80% above 1"] = (col[f] > 1).mean() >= 0.80

    #   Reproducibility and draw identifiers
    checks["same seed gives the same draws"] = sample_prior(n, seed) == draws
    checks["different seed gives different draws"] = sample_prior(n, seed + 1) != draws
    checks["draw_id runs 1 to n"] = [d["draw_id"] for d in draws] == list(range(1, n + 1))

    #   The written file reads back to the same values
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "params.json"
        write_params(draws, path)
        checks["written values read back exactly"] = json.loads(path.read_text()) == draws

    #   Report the prior's location for the record
    print(f"  inf_r median {med:.4f}, p_symp mean {col['p_symp'].mean():.3f}, "
          f"rec_scale median {np.median(col['rec_scale']):.2f}, imm_scale median {np.median(col['imm_scale']):.1f}")
    return report_checks(f"UT-2.1: Prior Sampling ({n} draws)", checks)


if __name__ == "__main__":
    passed = test_prior_sampling()
    sys.exit(0 if passed else 1)
