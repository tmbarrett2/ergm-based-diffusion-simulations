"""Bridge to the Julia simulator: run sbi_batch and read its records.

Jonathan H. Morgan, Ph.D.
26 October 2026

Python writes the parameter draws, calls the compiled diffusion_sim binary once
per batch as a subprocess, and reads the JSONL records (schema version 3) and
the manifest it writes. The command is configurable, so the same code can call
the binary or Julia directly.
"""

import json
import os
import subprocess
from pathlib import Path

import numpy as np

from .prior import THETA_FIELDS, write_params


# Default binary on the EC2 container
DEFAULT_BINARY = "/workspace/ergm-based-diffusion-simulations/julia_env/build/diffusion_sim_app/bin/diffusion_sim"

# Record layout expected from Julia
SCHEMA_VERSION = 3

# Target bounds: Python name -> sbi_batch option
TARGET_OPTIONS = {
    "I_peak_max": "--target-peak-max",
    "I_surge_max": "--target-surge-max",
    "surge_time": "--target-surge-time",
    "S_max": "--target-slope-max",
    "P_max_thresh": "--target-prev-max",
    "K_max": "--target-K-max",
}


# Helper Function for run_sbi_batch: manifest path beside the output
def manifest_path(out_file):
    """
    Args:
        out_file (str or Path): JSONL output path
    Returns:
        str: results.manifest.json for results.jsonl, as Julia names it
    """

    # Replace the extension
    out_file = str(out_file)
    if out_file.endswith(".jsonl"):
        return out_file[: -len(".jsonl")] + ".manifest.json"
    return out_file + ".manifest.json"


# Helper Function for run_sbi_batch: command line for one call
def build_command(command, graphml, village, params_file, out_file, master_seed,
                  sas_transformation, seed_component, target_region, maxtime,
                  extinct_policy, n_seed_nodes, replicates, seed_file):
    """
    Args:
        command (list of str): program to run, for example [DEFAULT_BINARY]
        other arguments: as in run_sbi_batch
    Returns:
        list of str: the full argument list
    Notes:
        A target bound set to None is left out, so it stays unconstrained.
    """

    # Required options
    args = list(command) + ["sbi_batch", "--graphml", str(graphml), "--village", str(village),
                            "--params", str(params_file), "--out", str(out_file),
                            "--max-time", str(maxtime), "--extinct-policy", extinct_policy,
                            "--n-seed-nodes", str(n_seed_nodes), "--replicates", str(replicates),
                            "--seed-component", seed_component]

    # Seeds: master seed (random mode) or seed file (testing mode)
    if seed_file is not None:
        args += ["--seed-mode", "file", "--seed-file", str(seed_file)]
    else:
        args += ["--seed", str(master_seed)]

    # Network form
    if sas_transformation:
        args.append("--sas-transformation")

    # Target bounds
    for name, value in (target_region or {}).items():
        if name not in TARGET_OPTIONS:
            raise ValueError(f"unknown target bound: {name}")
        if value is not None:
            args += [TARGET_OPTIONS[name], str(value)]

    # Assembling Result
    return args


# Run a Batch of Draws Through sbi_batch
def run_sbi_batch(graphml, draws, out_file, village, master_seed=None, *,
                  sas_transformation=True, seed_component="largest", target_region=None,
                  maxtime=200, extinct_policy="exclude", n_seed_nodes=1, replicates=1,
                  seed_file=None, n_threads=8, command=None, params_file=None):
    """
    **Description**
    Simulate a batch of parameter draws with the Julia simulator: write the
    draws, call sbi_batch once as a subprocess, and return the paths of its
    records and manifest.

    **Usage**
    run_sbi_batch(graphml, draws, out_file, village, master_seed=None, *,
                  sas_transformation=True, seed_component="largest",
                  target_region=None, maxtime=200, extinct_policy="exclude",
                  n_seed_nodes=1, replicates=1, seed_file=None, n_threads=8,
                  command=None, params_file=None)

    **Arguments**
    - graphml (str or Path): Village network file.
    - draws (list of dict): Parameter draws with draw_id (from prior.sample_prior).
    - out_file (str or Path): JSONL output; the manifest is written beside it.
    - village (str): Identifier written to every record.
    - master_seed (int): Batch master seed; required unless seed_file is given.
    - sas_transformation (bool): Load the network with SAS preprocessing
      (default True), as the validated replication did.
    - seed_component (str): "largest" (default) draws seeds from the largest
      weakly connected component; "any" from any node with a tie.
    - target_region (dict): Bounds keyed I_peak_max, I_surge_max, surge_time,
      S_max, P_max_thresh, K_max; a missing or None bound is unconstrained.
    - maxtime, extinct_policy, n_seed_nodes, replicates: As in sbi_batch
      (defaults 200, "exclude", 1, 1).
    - seed_file (str or Path): Seed file for testing mode, in place of master_seed.
    - n_threads (int): JULIA_NUM_THREADS for the subprocess (default 8).
    - command (list of str): Program to run (default [DEFAULT_BINARY]).
    - params_file (str or Path): Where to write the draws (default
      <out>.params.json beside out_file).

    **Details**
    The subprocess boundary is crossed once per batch. To call Julia directly
    instead of the binary, pass for example
    command=["julia", "--project=/workspace/ergm-based-diffusion-simulations/julia_env",
             "-e", "using diffusion_sim; diffusion_sim.CLI.cli_main(ARGS)"];
    Julia passes the remaining arguments to the script as ARGS. A failed batch
    raises RuntimeError carrying Julia's stderr and stdout.

    **Value**
    dict with out_file, manifest_file, params_file, and the subprocess stdout.

    **Examples**
    >>> from python_env.prior import sample_prior
    >>> res = run_sbi_batch("net_58.graphml", sample_prior(100, 1), "out/net_58.jsonl",
    ...                     "net_58", 1, target_region={"K_max": 1.0})
    >>> records = read_records(res["out_file"])

    **See Also**
    read_records, read_manifest, records_to_arrays
    """

    # Validation
    if seed_file is None and master_seed is None:
        raise ValueError("master_seed is required unless seed_file is given")
    if seed_component not in ("largest", "any"):
        raise ValueError("seed_component must be 'largest' or 'any'")

    # Write the draws beside the output
    out_file = Path(out_file)
    if params_file is None:
        params_file = out_file.with_suffix(".params.json")
    write_params(draws, params_file)

    # Run the subprocess with the requested thread count
    args = build_command(command or [DEFAULT_BINARY], graphml, village, params_file, out_file,
                         master_seed, sas_transformation, seed_component, target_region, maxtime,
                         extinct_policy, n_seed_nodes, replicates, seed_file)
    env = dict(os.environ, JULIA_NUM_THREADS=str(n_threads))
    proc = subprocess.run(args, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"sbi_batch failed (exit {proc.returncode}):\n{proc.stderr}\n{proc.stdout}")

    # Assembling Result
    return {"out_file": str(out_file), "manifest_file": manifest_path(out_file),
            "params_file": str(params_file), "stdout": proc.stdout}


# Read the Records of a Batch
def read_records(out_file):
    """
    **Description**
    Read the JSONL output of a batch into one dict per run, checking its layout.

    **Usage**
    read_records(out_file)

    **Arguments**
    - out_file (str or Path): JSONL output of run_sbi_batch.

    **Details**
    JSON null becomes None. Every record must have the theta, summary, target,
    and meta blocks and schema version 3; otherwise ValueError names the record.

    **Value**
    List of dicts in file order (draw order, then replicate order).

    **Examples**
    >>> records = read_records("out/net_58.jsonl")
    >>> records[0]["summary"]["class"]
    'established'

    **See Also**
    read_manifest, records_to_arrays
    """

    # Parse and check each line
    records = []
    with open(out_file) as io:
        for k, line in enumerate(io, start=1):
            rec = json.loads(line)
            for block in ("theta", "summary", "target", "meta"):
                if block not in rec:
                    raise ValueError(f"record {k} lacks the {block} block")
            if rec["meta"]["schema_version"] != SCHEMA_VERSION:
                raise ValueError(f"record {k}: schema version {rec['meta']['schema_version']}, "
                                 f"expected {SCHEMA_VERSION}")
            records.append(rec)

    # Assembling Result
    return records


# Read a Batch Manifest
def read_manifest(manifest_file):
    """
    **Description**
    Read the manifest of batch settings written beside a batch's records.

    **Usage**
    read_manifest(manifest_file)

    **Arguments**
    - manifest_file (str or Path): The <out>.manifest.json file.

    **Details**
    The manifest holds the settings shared by every run: the network and its
    component structure, the horizon, the seed mode and master seed, the
    extinct-run policy, and the target bounds.

    **Value**
    dict of the manifest fields.

    **Examples**
    >>> read_manifest("out/net_58.manifest.json")["sas_transformation"]
    True

    **See Also**
    read_records
    """

    # Parse the file
    with open(manifest_file) as io:
        return json.load(io)


# Helper Function for records_to_arrays: JSON value to float
def _as_float(x):
    """
    Args:
        x: a JSON value (number or None)
    Returns:
        float: NaN for None, otherwise float(x)
    """

    # Convert, with null as NaN
    return np.nan if x is None else float(x)


# Arrays for Training
def records_to_arrays(records, summary_fields=("K", "t_i", "peak_incidence", "frac_at_inflection"),
                      classes=("established",)):
    """
    **Description**
    Turn batch records into the arrays BayesFlow trains on.

    **Usage**
    records_to_arrays(records, summary_fields=("K", "t_i", "peak_incidence",
                      "frac_at_inflection"), classes=("established",))

    **Arguments**
    - records (list of dict): From read_records.
    - summary_fields (tuple of str): Summary quantities forming the conditioning
      vector. The default is a placeholder for the open decision on the
      summary vector (roadmap open decision 4).
    - classes (tuple of str): Run classes to keep (default established runs
      only, which conditions on establishment as the default extinct policy does).

    **Details**
    Null values become NaN, including the immunity parameters of a draw with
    permanent immunity. Rows follow the record order.

    **Value**
    dict of numpy arrays: thetas (N, 12) in THETA_FIELDS order; summaries
    (N, d); in_target (N,) as float; draw_id (N,); replicate (N,).

    **Examples**
    >>> arr = records_to_arrays(read_records("out/net_58.jsonl"))
    >>> arr["thetas"].shape[1]
    12

    **See Also**
    read_records
    """

    # Keep the requested classes
    kept = [r for r in records if r["summary"]["class"] in classes]

    # Build the arrays
    thetas = np.array([[_as_float(r["theta"][f]) for f in THETA_FIELDS] for r in kept],
                      dtype=float).reshape(-1, len(THETA_FIELDS))
    summaries = np.array([[_as_float(r["summary"][f]) for f in summary_fields] for r in kept],
                         dtype=float).reshape(-1, len(summary_fields))
    in_target = np.array([_as_float(r["target"]["in_target"]) for r in kept], dtype=float)
    draw_id = np.array([r["meta"]["draw_id"] for r in kept], dtype=int)
    replicate = np.array([r["meta"]["replicate"] for r in kept], dtype=int)

    # Assembling Result
    return {"thetas": thetas, "summaries": summaries, "in_target": in_target,
            "draw_id": draw_id, "replicate": replicate}
