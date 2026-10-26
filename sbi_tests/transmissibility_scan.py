"""Transmissibility Scan: Where Runs Begin to Establish on a Village Network.

Jonathan H. Morgan, Ph.D.
26 October 2026

Shifts the prior's transmissibility on a grid and runs each cell through the
bridge: the mean of b_int moves in steps, and the median of inf_r moves on the
logit scale, so every inf_r stays a probability in (0, 1).
Every cell uses the same base prior draws and the same batch seed, so only the
shift differs between cells, and a draw keeps its seed nodes in every cell.
The other 10 parameters keep their prior draws.

Revised 26 October 2026: inf_r was scaled by the ratio of medians, which put
draws above 1 at a median of 0.5; Julia rightly rejects those. It is now
shifted on the logit scale, which keeps each draw's rank and the cell's median.
Run from the repository root:

    python sbi_tests/transmissibility_scan.py

Options: --graphml (default net_58), --n (draws per cell, default 200),
--seed (default 2026), --b-int-means (default -2 -1 0 1), --inf-r-medians
(default 0.02 0.05 0.1), --threads (default 8).
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import math

import numpy as np

# Make the python_env package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.julia_bridge import read_manifest, read_records, run_sbi_batch  # noqa: E402
from python_env.prior import PRIOR_SPEC, sample_prior  # noqa: E402

DEFAULT_GRAPHML = "/workspace/ergm_data/sim_test_data/synthetic_village_m_graphml/net_58.graphml"
OUT_DIR = REPO / "sbi_tests" / "output" / "scan"
CLASSES = ("nonstarter", "extinct", "fit_failed", "established")


# Helper Function for shift_draws: log-odds of a probability
def _logit(p):
    """
    Args:
        p (float): probability in (0, 1)
    Returns:
        float: log(p / (1 - p))
    """

    # Convert
    return math.log(p / (1.0 - p))


# Helper Function for shift_draws: probability from log-odds
def _expit(x):
    """
    Args:
        x (float): log-odds
    Returns:
        float: 1 / (1 + exp(-x)), in (0, 1)
    """

    # Convert
    return 1.0 / (1.0 + math.exp(-x))


# Helper Function for main: shift the base draws to one cell of the grid
def shift_draws(base, b_int_mean, inf_r_median):
    """
    Args:
        base (list of dict): draws from sample_prior at the prior's settings
        b_int_mean (float): new mean of b_int
        inf_r_median (float): new median of inf_r
    Returns:
        list of dict: copies of the draws with b_int and inf_r shifted
    Notes:
        b_int is Normal, so it is shifted by the change in mean, and its
        spread is unchanged. inf_r is a probability, so it is shifted on the
        logit scale by the change in the logit of the median: every draw
        stays in (0, 1), keeps its rank, and the median moves exactly to
        inf_r_median. For small values this is close to scaling. Every other
        parameter is unchanged.
    """

    # Prior settings
    shift = b_int_mean - PRIOR_SPEC["b_int"][1]["mean"]
    logit_shift = _logit(inf_r_median) - _logit(PRIOR_SPEC["inf_r"][1]["median"])

    # Shifted copies
    draws = []
    for d in base:
        draw = dict(d)
        draw["b_int"] = d["b_int"] + shift
        draw["inf_r"] = _expit(_logit(d["inf_r"]) + logit_shift)
        draws.append(draw)

    # Assembling Result
    return draws


# Helper Function for main: metrics for one cell
def cell_metrics(b_int_mean, inf_r_median, records, n_nodes):
    """
    Args:
        b_int_mean (float): mean of b_int in this cell
        inf_r_median (float): median of inf_r in this cell
        records (list of dict): records from read_records
        n_nodes (int): nodes in the network (final_size is a share of them)
    Returns:
        dict: class shares; final size in people among viable runs (median,
            90th percentile, largest); median K of established runs
    """

    # Class shares
    n = len(records)
    counts = {c: 0 for c in CLASSES}
    for r in records:
        counts[r["summary"]["class"]] = counts.get(r["summary"]["class"], 0) + 1

    # Final size of viable runs, in people
    people = [r["summary"]["final_size"] * n_nodes for r in records
              if r["summary"]["class"] != "nonstarter" and r["summary"]["final_size"] is not None]
    pct = np.percentile(people, [50, 90, 100]) if people else [float("nan")] * 3

    # K of established runs
    K = [r["summary"]["K"] for r in records if r["summary"]["class"] == "established"]

    # Assembling Result
    result = {"b_int_mean": b_int_mean, "inf_r_median": inf_r_median, "runs": n}
    for c in CLASSES:
        result[c] = counts[c] / n
    result.update({"viable_people_median": float(pct[0]), "viable_people_p90": float(pct[1]),
                   "viable_people_max": float(pct[2]),
                   "K_median": float(np.median(K)) if K else float("nan")})
    return result


# Helper Function for main: print the grid as a table
def print_table(rows):
    """
    Args:
        rows (list of dict): metrics from cell_metrics
    Returns:
        None
    """

    # Header
    print("=" * 96)
    print(f"{'b_int':>6} {'inf_r':>6} | {'nonstart':>8} {'extinct':>8} {'failed':>7} {'estab':>7} | "
          f"{'viable people: med':>18} {'p90':>6} {'max':>6} | {'K med':>6}")
    print("=" * 96)

    # One line per cell
    for m in rows:
        print(f"{m['b_int_mean']:>6.1f} {m['inf_r_median']:>6.3f} | {m['nonstarter']:>8.1%} {m['extinct']:>8.1%} "
              f"{m['fit_failed']:>7.1%} {m['established']:>7.1%} | {m['viable_people_median']:>18.0f} "
              f"{m['viable_people_p90']:>6.0f} {m['viable_people_max']:>6.0f} | {m['K_median']:>6.3f}")


# Transmissibility Scan
def main():
    """
    Args:
        None (command-line arguments)
    Returns:
        None
    Notes:
        Each cell is one batch on the SAS-transformed network, with seeds
        from the largest component and the target bounded only by K_max = 1.
    """

    # Arguments
    parser = argparse.ArgumentParser(description="Transmissibility scan")
    parser.add_argument("--graphml", default=DEFAULT_GRAPHML)
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--b-int-means", type=float, nargs="+", default=[-2.0, -1.0, 0.0, 1.0])
    parser.add_argument("--inf-r-medians", type=float, nargs="+", default=[0.02, 0.05, 0.1])
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    # Validation
    if args.n < 1:
        raise ValueError("--n must be at least 1")
    if any(not 0 < m < 1 for m in args.inf_r_medians):
        raise ValueError("--inf-r-medians must lie in (0, 1)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # The same base draws for every cell
    base = sample_prior(args.n, args.seed)
    village = Path(args.graphml).stem

    # Run each cell
    rows = []
    start = time.time()
    for b_int_mean in args.b_int_means:
        for inf_r_median in args.inf_r_medians:
            out_file = OUT_DIR / f"{village}_bint{b_int_mean:g}_infr{inf_r_median:g}.jsonl"
            res = run_sbi_batch(args.graphml, shift_draws(base, b_int_mean, inf_r_median), out_file, village,
                                args.seed, sas_transformation=True, target_region={"K_max": 1.0},
                                n_threads=args.threads)
            n_nodes = read_manifest(res["manifest_file"])["n_nodes"]
            rows.append(cell_metrics(b_int_mean, inf_r_median, read_records(out_file), n_nodes))

    # Report
    print(f"Village {village}: {args.n} draws per cell, {len(rows)} cells, {time.time() - start:.1f} s")
    print_table(rows)

    # Save the grid
    table = OUT_DIR / f"scan_{village}.csv"
    with open(table, "w", newline="") as io:
        writer = csv.DictWriter(io, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Grid: {table}")
    print(f"  Records: {OUT_DIR}")


if __name__ == "__main__":
    main()
