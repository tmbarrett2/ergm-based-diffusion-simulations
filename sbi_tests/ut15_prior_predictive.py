"""UT-1.5: Prior Predictive Check on Village Networks.

Jonathan H. Morgan, Ph.D.
26 October 2026

Draws parameter vectors from the prior, runs them on each village network
(SAS-transformed) through the Python bridge, and reports the run classes, fit
success, coverage of K, shape, timing, and persistence. Run from the
repository root, giving the GraphML files:

    python sbi_tests/ut15_prior_predictive.py <net_43.graphml> <net_58.graphml> <net_83.graphml>

Options: --n (draws per village, default 500), --seed (prior and batch seed,
default 2026), --threads (default 8), --seed-component (largest or any,
default largest), --extinction-threshold (default: the binary's 0.02; needs
the reachable-threshold build), --calibration-dir (directory of
<village>_prior.json files from calibrate_prior.py; without it, every village
uses the default prior). Outputs carry the settings in their names, so runs can
be compared.

Revised 26 October 2026: per-village priors from calibration. Every village
uses the same seed, so with calibrated priors the draws share their underlying
random deviates and differ only through the prior. Criterion (c) judges the
upper end of K against the share of the network that is reachable (the largest
component), since K cannot exceed it.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

# Make the python_env package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.julia_bridge import read_manifest, read_records, run_sbi_batch  # noqa: E402
from python_env.prior import PRIOR_SPEC, read_calibration, sample_prior  # noqa: E402

OUT_DIR = REPO / "sbi_tests" / "output" / "ut15"


# Helper Function for village_metrics: median of one summary field
def median_of(records, field):
    """
    Args:
        records (list of dict): records from read_records
        field (str): summary field
    Returns:
        float: median over the non-null values, NaN if there are none
    """

    # Collect the values
    vals = [r["summary"][field] for r in records if r["summary"][field] is not None]

    # Assembling Result
    return float(np.median(vals)) if vals else float("nan")


# Helper Function for main: metrics for one village
def village_metrics(village, records, manifest, elapsed):
    """
    Args:
        village (str): village name
        records (list of dict): records from read_records
        manifest (dict): manifest from read_manifest
        elapsed (float): run time in seconds
    Returns:
        dict: network and seed structure, class counts, fit success rate
            (established share of viable runs), extinct share of viable runs,
            K range, NaN check, shape, timing, persistence, failure statuses,
            and run time
    """

    # Classes
    by_class = {}
    for r in records:
        by_class.setdefault(r["summary"]["class"], []).append(r)
    est = by_class.get("established", [])
    n = len(records)
    n_viable = n - len(by_class.get("nonstarter", []))
    viable = [r for r in records if r["summary"]["class"] != "nonstarter"]

    # Coverage of K and the NaN check for established runs
    K = [r["summary"]["K"] for r in est]
    richards = ("K", "r_growth", "t_i", "a", "frac_at_inflection", "fit_rmse")
    nan_free = all(r["summary"][f] is not None for r in est for f in richards)

    # Shape, tail, and persistence
    a_limit = sum(r["summary"]["a"] <= -0.9 for r in est)
    tails = [r["summary"]["t_90"] - r["summary"]["t_half"] for r in est]
    persistent = sum(not r["summary"]["extinguished"] for r in viable)

    # Statuses of failed fits
    failed = {}
    for r in by_class.get("fit_failed", []):
        failed[r["summary"]["status"]] = failed.get(r["summary"]["status"], 0) + 1

    # Assembling Result
    return {
        "village": village, "n_nodes": manifest["n_nodes"],
        "seed_component": manifest["seed_component"], "n_components": manifest["n_components"],
        "largest_component": manifest["largest_component"], "n_eligible_seeds": manifest["n_eligible_seeds"],
        "reachable_share": manifest["largest_component"] / manifest["n_nodes"],
        "extinction_threshold": manifest["extinction_threshold"], "inf_r_median": None,
        "draws": n,
        "nonstarter": len(by_class.get("nonstarter", [])), "extinct": len(by_class.get("extinct", [])),
        "established": len(est), "fit_failed": len(by_class.get("fit_failed", [])),
        "fit_success_rate": len(est) / n_viable if n_viable else float("nan"),
        "extinct_share_of_viable": len(by_class.get("extinct", [])) / n_viable if n_viable else float("nan"),
        "K_min": min(K) if K else float("nan"), "K_max": max(K) if K else float("nan"),
        "K_median": float(np.median(K)) if K else float("nan"),
        "nan_free": nan_free,
        "a_median": median_of(est, "a"), "frac_median": median_of(est, "frac_at_inflection"),
        "a_at_limit_share": a_limit / len(est) if est else float("nan"),
        "t_half_median": median_of(est, "t_half"),
        "tail_median": float(np.median(tails)) if tails else float("nan"),
        "persistent_share_of_viable": persistent / n_viable if n_viable else float("nan"),
        "failed_statuses": ", ".join(f"{k}={v}" for k, v in sorted(failed.items())) or "none",
        "seconds": elapsed,
    }


# Helper Function for main: print one village's report
def print_village(m):
    """
    Args:
        m (dict): metrics from village_metrics
    Returns:
        None
    """

    # Network and seeds
    print("=" * 60)
    print(f"Village {m['village']} ({m['n_nodes']} nodes, {m['draws']} draws, {m['seconds']:.1f} s)")
    print("=" * 60)
    print(f"  Network: {m['n_components']} components, largest {m['largest_component']} nodes; "
          f"seeds from {m['seed_component']} ({m['n_eligible_seeds']} eligible); "
          f"extinction threshold {m['extinction_threshold']:g}")
    print(f"  Prior: inf_r median {m['inf_r_median']:g}")

    # Classes and fits
    print(f"  Classes: established={m['established']}, extinct={m['extinct']}, "
          f"nonstarter={m['nonstarter']}, fit_failed={m['fit_failed']} ({m['failed_statuses']})")
    print(f"  Fit success rate (established / viable): {m['fit_success_rate']:.3f}")
    print(f"  Extinct share of viable runs:            {m['extinct_share_of_viable']:.3f}")

    # Established runs
    print(f"  K among established runs: min {m['K_min']:.3f}, median {m['K_median']:.3f}, max {m['K_max']:.3f}")
    print(f"  Shape: median a {m['a_median']:.2f}, median frac_at_inflection {m['frac_median']:.2f}, "
          f"a <= -0.9 in {m['a_at_limit_share']:.1%} of established runs")
    print(f"  Timing: median t_half {m['t_half_median']:.1f} days, median t_90 - t_half {m['tail_median']:.1f} days")
    print(f"  Still infected at the horizon: {m['persistent_share_of_viable']:.1%} of viable runs")
    print(f"  No null Richards field in established runs: {m['nan_free']}")


# Helper Function for main: report the UT-1.5 pass criteria
def print_criteria(metrics):
    """
    Args:
        metrics (list of dict): metrics from village_metrics
    Returns:
        None
    Notes:
        Criteria (a) and (b) are matched to villages by name ("58" and "43").
        An established run has K of at least 0.02 by definition (growth below
        0.02 is extinct), so the lower end of (c) is judged at 0.05.
    """

    # (a) and (b): fit success by village
    print("=" * 60)
    print("UT-1.5 Criteria")
    print("=" * 60)
    for m in metrics:
        name = m["village"]
        if "58" in name:
            print(f"  (a) {name}: fit success {m['fit_success_rate']:.3f} >= 0.70  "
                  f"{'PASS' if m['fit_success_rate'] >= 0.70 else 'FAIL'}")
        if "43" in name:
            print(f"  (b) {name}: fit success {m['fit_success_rate']:.3f} >= 0.40  "
                  f"{'PASS' if m['fit_success_rate'] >= 0.40 else 'FAIL'}")

    # (c) Coverage of K
    K_min = min(m["K_min"] for m in metrics)
    K_max = max(m["K_max"] for m in metrics)
    reach = max(m["K_max"] / m["reachable_share"] for m in metrics)
    print(f"  (c) K spans [0.02, 0.95 of reachable]: observed [{K_min:.3f}, {K_max:.3f}], "
          f"largest K = {reach:.2f} of the reachable share  "
          f"{'PASS' if K_min <= 0.05 and reach >= 0.95 else 'FAIL'} (lower end judged at 0.05)")

    # (d) No null Richards field
    print(f"  (d) no null Richards field when fit_ok: "
          f"{'PASS' if all(m['nan_free'] for m in metrics) else 'FAIL'}")


# UT-1.5: Prior Predictive Check
def main():
    """
    Args:
        None (command-line arguments)
    Returns:
        None
    Notes:
        Without --calibration-dir every village runs the same draws from the
        default prior; with it, each village draws from its own calibrated
        prior with the same seed.
    """

    # Arguments
    parser = argparse.ArgumentParser(description="UT-1.5 prior predictive check")
    parser.add_argument("graphml", nargs="+", help="village GraphML files")
    parser.add_argument("--n", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed-component", choices=("largest", "any"), default="largest")
    parser.add_argument("--extinction-threshold", type=float, default=None)
    parser.add_argument("--calibration-dir", default=None)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Settings passed to the batch, and the tag for the output names
    options = {}
    tag = args.seed_component
    if args.extinction_threshold is not None:
        options["extinction_threshold"] = args.extinction_threshold
        tag += f"_x{args.extinction_threshold:g}"
    tag += "_calibrated" if args.calibration_dir else "_default"

    # Run each village
    metrics = []
    for path in args.graphml:
        village = Path(path).stem
        spec = read_calibration(Path(args.calibration_dir) / f"{village}_prior.json") if args.calibration_dir else None
        draws = sample_prior(args.n, args.seed, spec=spec)
        out_file = OUT_DIR / f"{village}_{tag}.jsonl"
        start = time.time()
        res = run_sbi_batch(path, draws, out_file, village, args.seed, sas_transformation=True,
                            seed_component=args.seed_component, target_region={"K_max": 1.0},
                            n_threads=args.threads, **options)
        elapsed = time.time() - start
        m = village_metrics(village, read_records(out_file), read_manifest(res["manifest_file"]), elapsed)
        m["inf_r_median"] = (spec or PRIOR_SPEC)["inf_r"][1]["median"]
        print_village(m)
        metrics.append(m)

    # Pass criteria (roadmap UT-1.5)
    print_criteria(metrics)

    # Save the metrics
    table = OUT_DIR / f"ut15_metrics_{tag}.csv"
    with open(table, "w", newline="") as io:
        writer = csv.DictWriter(io, fieldnames=list(metrics[0].keys()))
        writer.writeheader()
        writer.writerows(metrics)
    print(f"  Metrics: {table}")
    print(f"  Records: {OUT_DIR}")


if __name__ == "__main__":
    main()
