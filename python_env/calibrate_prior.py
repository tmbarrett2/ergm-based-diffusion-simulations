"""Calibrate the Prior's Transmissibility to Each Village Network.

Jonathan H. Morgan, Ph.D.
26 October 2026

For each GraphML file: scans the median of inf_r, prints the table with the
suggested median for the target, asks the user to choose (Enter accepts the
suggestion), and writes <village>_prior.json with the scan, the choice, and the
prior spec. Run from the repository root:

    python sbi_tests/calibrate_prior.py <net_43.graphml> <net_58.graphml> <net_83.graphml>

Options: --target (default 0.5), --metric (established or fit_success,
default established), --medians (scan grid, default 0.02 to 0.8), --n (draws
per median, default 200), --seed (default 2026), --threads (default 8),
--choose (a median applied to every village without asking, for unattended
runs), --out-dir (default sbi_tests/output/calibration).
"""

import argparse
import sys
from pathlib import Path

# Make the python_env package importable from the repository root
REPO = Path("/workspace/ergm-based-diffusion-simulations")
sys.path.insert(0, str(REPO))

from python_env.prior import (SCAN_MEDIANS, format_scan_table, scan_inf_r,  # noqa: E402
                              suggest_inf_r_median, write_calibration)

OUT_DIR = REPO / "sbi_tests" / "output" / "calibration"


# Helper Function for main: the user's choice of median
def choose_median(suggested):
    """
    Args:
        suggested (float or None): the suggested median, if any
    Returns:
        float: the median entered, or the suggestion on Enter
    Notes:
        Asks again until the answer is a number in (0, 1), or Enter with a
        suggestion available.
    """

    # Ask until the answer is valid
    prompt = f"Choose the inf_r median [{suggested:g}]: " if suggested is not None else "Choose the inf_r median: "
    while True:
        answer = input(prompt).strip()
        if answer == "" and suggested is not None:
            return suggested
        try:
            value = float(answer)
        except ValueError:
            print("  Enter a number in (0, 1).")
            continue
        if 0.0 < value < 1.0:
            return value
        print("  Enter a number in (0, 1).")


# Prior Calibration
def main():
    """
    Args:
        None (command-line arguments)
    Returns:
        None
    """

    # Arguments
    parser = argparse.ArgumentParser(description="Calibrate the inf_r median to each village network")
    parser.add_argument("graphml", nargs="+", help="village GraphML files")
    parser.add_argument("--target", type=float, default=0.5)
    parser.add_argument("--metric", choices=("established", "fit_success"), default="established")
    parser.add_argument("--medians", type=float, nargs="+", default=list(SCAN_MEDIANS))
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--choose", type=float, default=None)
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    args = parser.parse_args()

    # Validation
    if args.choose is not None and not 0.0 < args.choose < 1.0:
        raise ValueError("--choose must lie in (0, 1)")
    out_dir = Path(args.out_dir)

    # Scan, present, choose, and record each village
    for path in args.graphml:
        scan = scan_inf_r(path, out_dir, medians=args.medians, n=args.n, seed=args.seed, n_threads=args.threads)
        print(format_scan_table(scan, args.target, args.metric))
        suggested = suggest_inf_r_median(scan, args.target, args.metric)["inf_r_median"]
        chosen = args.choose if args.choose is not None else choose_median(suggested)
        record = out_dir / f"{scan['village']}_prior.json"
        write_calibration(scan, chosen, record, args.target, args.metric)
        print(f"  {scan['village']}: inf_r median {chosen:g} recorded in {record}")


if __name__ == "__main__":
    main()
