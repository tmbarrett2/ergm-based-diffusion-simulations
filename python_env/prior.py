"""Prior over the 12 simulator parameters, with calibration to a network.

Jonathan H. Morgan, Ph.D.
26 October 2026

Draws parameter vectors from the weakly informative prior and writes them as
the JSON array that the Julia sbi_batch mode reads.

Revised 26 October 2026: inf_r, a daily transmission probability, is
logit-normal, so every draw lies in (0, 1); its median is set per network.
Network topology changes how transmissible a disease must be to establish, so
the median of inf_r is chosen from a scan on each network: scan_inf_r runs the
scan, format_scan_table presents it with a suggested median for a settable
target, and the user selects. write_calibration records the scan and the
choice; read_calibration returns the prior spec for that network. b_int ~
N(0, 1). The other parameters are as in roadmap Section 5.1.
"""

import copy
import json
import math
from datetime import date
from pathlib import Path

import numpy as np


# Order of the 12 parameters, matching the Julia record (THETA_FIELDS)
THETA_FIELDS = (
    "inf_r", "b_int", "b_close", "b_cxn_peers", "b_cxn_total", "b_cxn_symp",
    "b_cls_x_smp", "p_symp", "rec_scale", "rec_shape", "imm_scale", "imm_shape",
)

# Default prior specification: name -> (distribution, parameters)
PRIOR_SPEC = {
    "inf_r": ("logitnormal", {"median": 0.3, "sd": 0.75}),
    "b_int": ("normal", {"mean": 0.0, "sd": 1.0}),
    "b_close": ("truncated_normal", {"mean": 1.0, "sd": 1.0, "lower": 0.0}),
    "b_cxn_peers": ("normal", {"mean": 0.0, "sd": 1.0}),
    "b_cxn_total": ("normal", {"mean": 0.0, "sd": 1.0}),
    "b_cxn_symp": ("normal", {"mean": 0.0, "sd": 1.0}),
    "b_cls_x_smp": ("normal", {"mean": 0.0, "sd": 0.5}),
    "p_symp": ("beta", {"a": 2.0, "b": 2.0}),
    "rec_scale": ("lognormal", {"median": 10.0, "sigma": 0.4}),
    "rec_shape": ("gamma", {"shape": 4.0, "scale": 0.5}),
    "imm_scale": ("lognormal", {"median": 180.0, "sigma": 0.4}),
    "imm_shape": ("gamma", {"shape": 4.0, "scale": 0.5}),
}

# Default scan grid of inf_r medians, wide enough for sparse and dense networks
SCAN_MEDIANS = (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# Run classes reported by the scan
SCAN_CLASSES = ("nonstarter", "extinct", "fit_failed", "established")


# Helper Function for the prior: log-odds of a probability
def _logit(p):
    """
    Args:
        p (float or numpy.ndarray): probability in (0, 1)
    Returns:
        float or numpy.ndarray: log(p / (1 - p))
    """

    # Convert
    return np.log(p / (1.0 - p))


# Helper Function for the prior: probability from log-odds
def _expit(x):
    """
    Args:
        x (float or numpy.ndarray): log-odds
    Returns:
        float or numpy.ndarray: 1 / (1 + exp(-x)), in (0, 1)
    """

    # Convert
    return 1.0 / (1.0 + np.exp(-x))


# Helper Function for sample_prior: truncated normal by rejection
def _truncated_normal(rng, mean, sd, lower, size):
    """
    Args:
        rng (numpy.random.Generator): random number generator
        mean (float): mean of the untruncated normal
        sd (float): standard deviation of the untruncated normal
        lower (float): lower truncation bound
        size (int): number of draws
    Returns:
        numpy.ndarray: shape (size,)
    Notes:
        Draws below the bound are redrawn. With the bound one standard
        deviation below the mean, about 84% of draws are accepted, so a
        few passes suffice.
    """

    # Initial draws
    out = rng.normal(mean, sd, size)
    bad = out < lower

    # Redraw the rejected entries until none remain
    while bad.any():
        out[bad] = rng.normal(mean, sd, bad.sum())
        bad = out < lower

    # Assembling Result
    return out


# Helper Function for sample_prior: draws of one parameter
def _draw_parameter(rng, dist, par, size):
    """
    Args:
        rng (numpy.random.Generator): random number generator
        dist (str): 'logitnormal', 'lognormal', 'normal', 'truncated_normal',
            'beta', or 'gamma'
        par (dict): distribution parameters, as in PRIOR_SPEC
        size (int): number of draws
    Returns:
        numpy.ndarray: shape (size,)
    Notes:
        LogitNormal and LogNormal are parameterized by their medians; the
        logit-normal sd is on the log-odds scale. Gamma uses the scale form,
        so Gamma(4, 0.5) has mean 2.
    """

    # Dispatch on the distribution name
    if dist == "logitnormal":
        return _expit(rng.normal(_logit(par["median"]), par["sd"], size))
    if dist == "lognormal":
        return rng.lognormal(math.log(par["median"]), par["sigma"], size)
    if dist == "normal":
        return rng.normal(par["mean"], par["sd"], size)
    if dist == "truncated_normal":
        return _truncated_normal(rng, par["mean"], par["sd"], par["lower"], size)
    if dist == "beta":
        return rng.beta(par["a"], par["b"], size)
    if dist == "gamma":
        return rng.gamma(par["shape"], par["scale"], size)

    # Unknown distribution
    raise ValueError(f"unknown distribution: {dist}")


# Prior Specification with Chosen Transmissibility
def make_prior_spec(inf_r_median=0.3, inf_r_logit_sd=0.75, b_int_mean=0.0):
    """
    **Description**
    Build a prior specification: the default prior with the median of inf_r,
    its spread, and the mean of b_int set as given.

    **Usage**
    make_prior_spec(inf_r_median=0.3, inf_r_logit_sd=0.75, b_int_mean=0.0)

    **Arguments**
    - inf_r_median (float): Median of inf_r, in (0, 1); usually chosen from a
      scan of the network (scan_inf_r).
    - inf_r_logit_sd (float): Standard deviation of inf_r on the log-odds scale
      (default 0.75: at a median of 0.3, about 90% of draws lie in 0.1 to 0.6).
    - b_int_mean (float): Mean of b_int (default 0).

    **Details**
    The other ten parameters keep their PRIOR_SPEC distributions. The result is
    a new dict; PRIOR_SPEC is not changed.

    **Value**
    dict in the PRIOR_SPEC format.

    **Examples**
    >>> spec = make_prior_spec(inf_r_median=0.2)
    >>> spec["inf_r"]
    ('logitnormal', {'median': 0.2, 'sd': 0.75})

    **See Also**
    sample_prior, scan_inf_r, read_calibration
    """

    # Validation
    if not 0.0 < inf_r_median < 1.0:
        raise ValueError("inf_r_median must lie in (0, 1)")
    if inf_r_logit_sd <= 0.0:
        raise ValueError("inf_r_logit_sd must be positive")

    # Copy the default and set the transmissibility
    spec = copy.deepcopy(PRIOR_SPEC)
    spec["inf_r"] = ("logitnormal", {"median": float(inf_r_median), "sd": float(inf_r_logit_sd)})
    spec["b_int"] = ("normal", {"mean": float(b_int_mean), "sd": PRIOR_SPEC["b_int"][1]["sd"]})

    # Assembling Result
    return spec


# Prior Draws
def sample_prior(n, seed, first_draw_id=1, spec=None):
    """
    **Description**
    Draw parameter vectors from the prior over the 12 simulator parameters,
    reproducibly from a seed, ready to be written as an sbi_batch parameter file.

    **Usage**
    sample_prior(n, seed, first_draw_id=1, spec=None)

    **Arguments**
    - n (int): Number of draws (at least 1).
    - seed (int): Seed for numpy.random.default_rng.
    - first_draw_id (int): draw_id of the first draw (default 1); later draws
      count up from it.
    - spec (dict): Prior specification (default PRIOR_SPEC); see
      make_prior_spec and read_calibration.

    **Details**
    The default prior is logit-normal for inf_r (median 0.3), Normal for the
    behavioral coefficients (b_int ~ N(0, 1)), a Normal truncated at 0 for
    b_close, Beta(2, 2) for p_symp, LogNormal for the Weibull scales, and
    Gamma(4, 0.5) for the Weibull shapes. Each parameter is drawn across all n
    draws in THETA_FIELDS order, so the draws for a given seed do not depend on
    how the result is later split into batches. Immunity is always drawn, so
    imm_scale and imm_shape are never null.

    **Value**
    List of n dicts, each with draw_id and the 12 parameters as Python floats,
    in THETA_FIELDS order.

    **Examples**
    >>> draws = sample_prior(3, seed=2026)
    >>> [d["draw_id"] for d in draws]
    [1, 2, 3]
    >>> sample_prior(3, seed=2026, spec=read_calibration("net_58_prior.json"))

    **See Also**
    make_prior_spec, write_params, julia_bridge.run_sbi_batch
    """

    # Validation
    if n < 1:
        raise ValueError("n must be at least 1")
    spec = PRIOR_SPEC if spec is None else spec
    missing = [f for f in THETA_FIELDS if f not in spec]
    if missing:
        raise ValueError(f"prior spec lacks {missing}")

    # Draw each parameter across all draws
    rng = np.random.default_rng(seed)
    columns = {name: _draw_parameter(rng, *spec[name], n) for name in THETA_FIELDS}

    # Assemble one dict per draw
    draws = []
    for k in range(n):
        draw = {"draw_id": first_draw_id + k}
        for name in THETA_FIELDS:
            draw[name] = float(columns[name][k])
        draws.append(draw)

    # Assembling Result
    return draws


# Write a Parameter File for sbi_batch
def write_params(draws, path):
    """
    **Description**
    Write parameter draws as the JSON array that the Julia sbi_batch mode reads.

    **Usage**
    write_params(draws, path)

    **Arguments**
    - draws (list of dict): Draws from sample_prior, or hand-built with the same keys.
    - path (str or pathlib.Path): File to write.

    **Details**
    Every draw must have a unique draw_id and all 12 parameters; imm_scale and
    imm_shape may both be None for permanent immunity. Python's json writes
    floats in their shortest round-trip form, so Julia receives exactly the
    values drawn here (checked by UT-2.2).

    **Value**
    int: the number of draws written.

    **Examples**
    >>> n = write_params(sample_prior(10, seed=1), "params.json")
    >>> n
    10

    **See Also**
    sample_prior
    """

    # Validation: unique draw_id and every parameter present
    ids = [d["draw_id"] for d in draws]
    if len(set(ids)) != len(ids):
        raise ValueError("draw_id values must be unique")
    for d in draws:
        missing = [f for f in THETA_FIELDS if f not in d]
        if missing:
            raise ValueError(f"draw {d['draw_id']} lacks {missing}")

    # Write the file
    with open(path, "w") as io:
        json.dump(draws, io)

    # Assembling Result
    return len(draws)


# Helper Function for scan_inf_r: move every inf_r draw to a new median
def shift_inf_r(draws, from_median, to_median):
    """
    Args:
        draws (list of dict): parameter draws
        from_median (float): median of inf_r in the prior the draws came from
        to_median (float): median of inf_r wanted
    Returns:
        list of dict: copies of the draws with inf_r shifted on the log-odds scale
    Notes:
        Every draw stays in (0, 1) and keeps its rank, the spread on the
        log-odds scale is unchanged, and the median moves exactly. The
        result is a draw from make_prior_spec(to_median, ...) with the same
        underlying normal deviates, so the cells of a scan differ only in
        transmissibility.
    """

    # Shift on the log-odds scale
    shift = _logit(to_median) - _logit(from_median)
    out = []
    for d in draws:
        draw = dict(d)
        draw["inf_r"] = float(_expit(_logit(d["inf_r"]) + shift))
        out.append(draw)

    # Assembling Result
    return out


# Helper Function for scan_inf_r: metrics for one cell
def _scan_cell(inf_r_median, records, n_nodes):
    """
    Args:
        inf_r_median (float): median of inf_r in this cell
        records (list of dict): records from read_records
        n_nodes (int): nodes in the network (final_size is a share of them)
    Returns:
        dict: class shares; fit success (established among runs that spread);
            final size in people among runs that spread (median, 90th
            percentile, largest); median K of established runs
    """

    # Class shares
    n = len(records)
    counts = {c: 0 for c in SCAN_CLASSES}
    for r in records:
        counts[r["summary"]["class"]] = counts.get(r["summary"]["class"], 0) + 1
    viable = n - counts["nonstarter"]

    # Final size of runs that spread, in people
    people = [r["summary"]["final_size"] * n_nodes for r in records
              if r["summary"]["class"] != "nonstarter" and r["summary"]["final_size"] is not None]
    pct = np.percentile(people, [50, 90, 100]) if people else [float("nan")] * 3

    # K of established runs
    K = [r["summary"]["K"] for r in records if r["summary"]["class"] == "established"]

    # Assembling Result
    cell = {"inf_r_median": inf_r_median, "runs": n}
    for c in SCAN_CLASSES:
        cell[c] = counts[c] / n
    cell.update({"fit_success": counts["established"] / viable if viable else float("nan"),
                 "people_median": float(pct[0]), "people_p90": float(pct[1]), "people_max": float(pct[2]),
                 "K_median": float(np.median(K)) if K else float("nan")})
    return cell


# Transmissibility Scan of a Network
def scan_inf_r(graphml, out_dir, medians=SCAN_MEDIANS, n=200, seed=2026, spec=None, n_threads=8,
               sas_transformation=True, seed_component="largest"):
    """
    **Description**
    Scan the median of inf_r on one network: simulate the same prior draws at
    each median and report how often runs establish.

    **Usage**
    scan_inf_r(graphml, out_dir, medians=SCAN_MEDIANS, n=200, seed=2026,
               spec=None, n_threads=8, sas_transformation=True,
               seed_component="largest")

    **Arguments**
    - graphml (str or Path): Village network file.
    - out_dir (str or Path): Directory for each cell's records and manifest.
    - medians (sequence of float): inf_r medians to scan, in (0, 1); the default
      spans 0.02 to 0.8.
    - n (int): Draws per cell (default 200).
    - seed (int): Prior seed and batch master seed (default 2026).
    - spec (dict): Prior for the other parameters and the inf_r spread
      (default PRIOR_SPEC).
    - n_threads (int): Julia threads (default 8).
    - sas_transformation (bool), seed_component (str): Batch settings, as in
      run_sbi_batch (defaults True and "largest").

    **Details**
    The base draws come from spec once; each cell shifts only inf_r, on the
    log-odds scale (shift_inf_r), and uses the same batch seed, so a draw keeps
    its seed nodes in every cell and the cells differ only in transmissibility.
    Network topology sets how transmissible a disease must be to establish, so
    the scan is run on each network before its prior is fixed. Each cell takes
    about 10 seconds at 200 draws on a 2,500-node network with 8 threads.

    **Value**
    dict: village, graphml, n_nodes, n, seed, inf_r_logit_sd, b_int_mean, and
    cells (list of dict, one per median, from smallest to largest).

    **Examples**
    >>> scan = scan_inf_r("net_58.graphml", "output/calibration")
    >>> print(format_scan_table(scan, target=0.5))

    **See Also**
    suggest_inf_r_median, format_scan_table, write_calibration
    """

    # The bridge imports this module, so it is imported here, at call time
    from .julia_bridge import read_manifest, read_records, run_sbi_batch

    # Validation
    medians = sorted(float(m) for m in medians)
    if not medians or any(not 0.0 < m < 1.0 for m in medians):
        raise ValueError("medians must be a non-empty set of values in (0, 1)")
    if n < 1:
        raise ValueError("n must be at least 1")
    spec = PRIOR_SPEC if spec is None else spec
    if spec["inf_r"][0] != "logitnormal":
        raise ValueError("the scan requires a logit-normal prior for inf_r")

    # The same base draws for every cell
    base = sample_prior(n, seed, spec=spec)
    from_median = spec["inf_r"][1]["median"]
    village = Path(graphml).stem
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Run each cell
    cells = []
    n_nodes = None
    for m in medians:
        out_file = out_dir / f"{village}_scan_infr{m:g}.jsonl"
        res = run_sbi_batch(graphml, shift_inf_r(base, from_median, m), out_file, village, seed,
                            sas_transformation=sas_transformation, seed_component=seed_component,
                            target_region={"K_max": 1.0}, n_threads=n_threads)
        n_nodes = read_manifest(res["manifest_file"])["n_nodes"]
        cells.append(_scan_cell(m, read_records(out_file), n_nodes))

    # Assembling Result
    return {"village": village, "graphml": str(graphml), "n_nodes": n_nodes, "n": n, "seed": seed,
            "inf_r_logit_sd": spec["inf_r"][1]["sd"], "b_int_mean": spec["b_int"][1]["mean"],
            "cells": cells}


# Suggested inf_r Median for a Target
def suggest_inf_r_median(scan, target=0.5, metric="established"):
    """
    **Description**
    Suggest the inf_r median at which a scan metric first reaches a target.

    **Usage**
    suggest_inf_r_median(scan, target=0.5, metric="established")

    **Arguments**
    - scan (dict): Result of scan_inf_r.
    - target (float): Value the metric should reach (default 0.5).
    - metric (str): "established" (share of all runs, default) or
      "fit_success" (established among runs that spread).

    **Details**
    Finds the first pair of neighboring cells between which the metric crosses
    the target and interpolates linearly in the log-odds of the median. If the
    smallest median already reaches the target, or no median does, there is no
    suggestion and the reason is given; widen the grid. The metric is noisy at
    200 draws per cell, so the suggestion is a guide for the user's choice.

    **Value**
    dict: inf_r_median (float or None), target, metric, and note (str).

    **Examples**
    >>> suggest_inf_r_median(scan, target=0.5)["inf_r_median"]
    0.3

    **See Also**
    scan_inf_r, format_scan_table
    """

    # Validation
    if metric not in ("established", "fit_success"):
        raise ValueError("metric must be 'established' or 'fit_success'")
    if not 0.0 < target < 1.0:
        raise ValueError("target must lie in (0, 1)")
    cells = scan["cells"]

    # Early return: target reached at the smallest median
    if cells[0][metric] >= target:
        return {"inf_r_median": None, "target": target, "metric": metric,
                "note": f"{metric} is already {cells[0][metric]:.2f} at the smallest median; widen the grid downward"}

    # First crossing, interpolated in the log-odds of the median
    for lo, hi in zip(cells[:-1], cells[1:]):
        if lo[metric] < target <= hi[metric]:
            w = (target - lo[metric]) / (hi[metric] - lo[metric])
            x = _logit(lo["inf_r_median"]) + w * (_logit(hi["inf_r_median"]) - _logit(lo["inf_r_median"]))
            return {"inf_r_median": round(float(_expit(x)), 3), "target": target, "metric": metric,
                    "note": f"between {lo['inf_r_median']:g} and {hi['inf_r_median']:g}"}

    # Assembling Result: target not reached
    return {"inf_r_median": None, "target": target, "metric": metric,
            "note": f"{metric} reaches only {max(c[metric] for c in cells):.2f}; widen the grid upward"}


# Scan Table for the User
def format_scan_table(scan, target=0.5, metric="established"):
    """
    **Description**
    Format a scan as a table, with the suggested inf_r median for a target, so
    the user can choose the median of the network's prior.

    **Usage**
    format_scan_table(scan, target=0.5, metric="established")

    **Arguments**
    - scan (dict): Result of scan_inf_r.
    - target (float), metric (str): As in suggest_inf_r_median.

    **Details**
    One line per median: the share of nonstarter, extinct, fit_failed, and
    established runs; fit success; final size in people among runs that
    spread; and the median K of established runs. The suggestion is printed
    below the table.

    **Value**
    str: the table and the suggestion.

    **Examples**
    >>> print(format_scan_table(scan, target=0.5))

    **See Also**
    scan_inf_r, suggest_inf_r_median
    """

    # Header
    lines = [f"Village {scan['village']}: {scan['n']} draws per median, inf_r logit sd {scan['inf_r_logit_sd']:g}, "
             f"b_int mean {scan['b_int_mean']:g}",
             "=" * 98,
             f"{'inf_r':>6} | {'nonstart':>8} {'extinct':>8} {'failed':>7} {'estab':>7} {'fit succ':>8} | "
             f"{'people: med':>11} {'p90':>6} {'max':>6} | {'K med':>6}",
             "=" * 98]

    # One line per cell
    for c in scan["cells"]:
        lines.append(f"{c['inf_r_median']:>6.3f} | {c['nonstarter']:>8.1%} {c['extinct']:>8.1%} "
                     f"{c['fit_failed']:>7.1%} {c['established']:>7.1%} {c['fit_success']:>8.2f} | "
                     f"{c['people_median']:>11.0f} {c['people_p90']:>6.0f} {c['people_max']:>6.0f} | "
                     f"{c['K_median']:>6.3f}")

    # Suggestion
    s = suggest_inf_r_median(scan, target, metric)
    median = "none" if s["inf_r_median"] is None else f"{s['inf_r_median']:g}"
    lines.append(f"Suggested inf_r median for {metric} = {target:g}: {median} ({s['note']})")

    # Assembling Result
    return "\n".join(lines)


# Record a Calibration
def write_calibration(scan, inf_r_median, path, target=0.5, metric="established"):
    """
    **Description**
    Record a network's scan, the suggestion, and the user's chosen inf_r median,
    with the resulting prior specification.

    **Usage**
    write_calibration(scan, inf_r_median, path, target=0.5, metric="established")

    **Arguments**
    - scan (dict): Result of scan_inf_r.
    - inf_r_median (float): The median the user chose, in (0, 1).
    - path (str or Path): JSON file to write, for example <village>_prior.json.
    - target (float), metric (str): The target the suggestion was made for.

    **Details**
    The record keeps the whole scan, so every village's prior can be traced to
    the evidence for it. The spec uses the scan's inf_r spread and b_int mean.

    **Value**
    dict: the prior specification written (as make_prior_spec returns).

    **Examples**
    >>> spec = write_calibration(scan, 0.3, "net_58_prior.json")

    **See Also**
    read_calibration, scan_inf_r
    """

    # Validation
    if not 0.0 < inf_r_median < 1.0:
        raise ValueError("inf_r_median must lie in (0, 1)")

    # Spec and record
    spec = make_prior_spec(inf_r_median, scan["inf_r_logit_sd"], scan["b_int_mean"])
    record = {"village": scan["village"], "graphml": scan["graphml"], "date": date.today().isoformat(),
              "target": target, "metric": metric,
              "suggested": suggest_inf_r_median(scan, target, metric)["inf_r_median"],
              "chosen_inf_r_median": float(inf_r_median), "spec": spec, "scan": scan}

    # Write the file
    with open(path, "w") as io:
        json.dump(record, io, indent=2)

    # Assembling Result
    return spec


# Read a Calibration
def read_calibration(path):
    """
    **Description**
    Read the prior specification chosen for a network.

    **Usage**
    read_calibration(path)

    **Arguments**
    - path (str or Path): File written by write_calibration.

    **Details**
    JSON stores each (distribution, parameters) pair as a list; it is returned
    as a tuple, as in PRIOR_SPEC.

    **Value**
    dict in the PRIOR_SPEC format, ready for sample_prior(spec=...).

    **Examples**
    >>> draws = sample_prior(500, 2026, spec=read_calibration("net_58_prior.json"))

    **See Also**
    write_calibration, sample_prior
    """

    # Parse the file
    with open(path) as io:
        record = json.load(io)

    # Assembling Result
    return {name: (entry[0], entry[1]) for name, entry in record["spec"].items()}
