"""Bridge to the Julia simulator: run sbi_batch and read its records.

Jonathan H. Morgan, Ph.D.
26 October 2026

Python writes the parameter draws, calls the compiled diffusion_sim binary
once per batch as a subprocess, and reads the JSONL records (schema
version 3) and the manifest that it writes. The command is configurable,
so the same code can call the binary or Julia directly.
"""

import json
import os
import subprocess
from pathlib import Path

import numpy as np

from .prior import THETA_FIELDS, write_params


#   Default binary on the EC2 container
DEFAULT_BINARY = "/workspace/ergm-based-diffusion-simulations/julia_env/build/diffusion_sim_app/bin/diffusion_sim"

#   Record layout expected from Julia (schema version 3)
SCHEMA_VERSION = 3

#   Target bounds: Python name -> sbi_batch option
TARGET_OPTIONS = {
    "I_peak_max": "--target-peak-max",
    "I_surge_max": "--target-surge-max",
    "surge_time": "--target-surge-time",
    "S_max": "--target-slope-max",
    "P_max_thresh": "--target-prev-max",
    "K_max": "--target-K-max",
}


# ============================================================
#   RUNNING A BATCH
# ============================================================

#   Helper for run_sbi_batch: the manifest path beside the output
def manifest_path(out_file):
    """Return results.manifest.json for results.jsonl, as Julia names it."""

    #   Replace the extension
    out_file = str(out_file)
    if out_file.endswith(".jsonl"):
        return out_file[: -len(".jsonl")] + ".manifest.json"
    return out_file + ".manifest.json"


#   Helper for run_sbi_batch: build the command line
def build_command(command, graphml, village, params_file, out_file, master_seed,
                  sas_transformation, target_region, maxtime, extinct_policy,
                  n_seed_nodes, replicates, seed_file):
    """Return the argument list for one sbi_batch call.

    Parameters
    ----------
    command : list of str
        The program to run, for example [DEFAULT_BINARY].
    Other parameters are as in run_sbi_batch.

    Returns
    -------
    list of str
    """

    #   Required options
    args = list(command) + ["sbi_batch", "--graphml", str(graphml), "--village", str(village),
                            "--params", str(params_file), "--out", str(out_file),
                            "--max-time", str(maxtime), "--extinct-policy", extinct_policy,
                            "--n-seed-nodes", str(n_seed_nodes), "--replicates", str(replicates)]

    #   Seeds: master seed (random mode) or seed file (testing mode)
    if seed_file is not None:
        args += ["--seed-mode", "file", "--seed-file", str(seed_file)]
    else:
        args += ["--seed", str(master_seed)]

    #   Network form and target bounds
    if sas_transformation:
        args.append("--sas-transformation")
    for name, value in (target_region or {}).items():
        if name not in TARGET_OPTIONS:
            raise ValueError(f"unknown target bound: {name}")
        if value is not None:
            args += [TARGET_OPTIONS[name], str(value)]
    return args


#   Run a Batch of Draws Through sbi_batch
def run_sbi_batch(graphml, draws, out_file, village, master_seed=None, *,
                  sas_transformation=True, target_region=None, maxtime=200,
                  extinct_policy="exclude", n_seed_nodes=1, replicates=1,
                  seed_file=None, n_threads=8, command=None, params_file=None):
    """Simulate a batch of parameter draws with the Julia simulator.

    Parameters
    ----------
    graphml : str or Path
        Village network file.
    draws : list of dict
        Parameter draws with draw_id (from prior.sample_prior).
    out_file : str or Path
        JSONL output; the manifest is written beside it.
    village : str
        Identifier written to every record.
    master_seed : int, optional
        Batch master seed; required unless seed_file is given.
    sas_transformation : bool, default True
        Load the network with SAS preprocessing, as the validated
        replication did. True until the network representation is
        decided after the pipeline is built.
    target_region : dict, optional
        Bounds keyed I_peak_max, I_surge_max, surge_time, S_max,
        P_max_thresh, K_max; None or omitted leaves a bound unset.
    maxtime, extinct_policy, n_seed_nodes, replicates : as in sbi_batch
    seed_file : str or Path, optional
        Seed file for testing mode, in place of master_seed.
    n_threads : int, default 8
        JULIA_NUM_THREADS for the subprocess.
    command : list of str, optional
        Program to run; default [DEFAULT_BINARY]. To call Julia directly,
        pass for example ["julia", "--project=/workspace/ergm-based-diffusion-simulations/julia_env",
        "-e", "using diffusion_sim; diffusion_sim.CLI.cli_main(ARGS)"]; Julia passes
        the remaining arguments to the script as ARGS.
    params_file : str or Path, optional
        Where to write the draws; default beside out_file as
        <out>.params.json.

    Returns
    -------
    dict
        out_file, manifest_file, params_file, and the subprocess stdout.

    Raises
    ------
    RuntimeError
        If the subprocess fails; its stderr is included.
    """

    #   Validation
    if seed_file is None and master_seed is None:
        raise ValueError("master_seed is required unless seed_file is given")

    #   Write the draws beside the output
    out_file = Path(out_file)
    if params_file is None:
        params_file = out_file.with_suffix(".params.json")
    write_params(draws, params_file)

    #   Run the subprocess with the requested thread count
    args = build_command(command or [DEFAULT_BINARY], graphml, village, params_file, out_file,
                         master_seed, sas_transformation, target_region, maxtime,
                         extinct_policy, n_seed_nodes, replicates, seed_file)
    env = dict(os.environ, JULIA_NUM_THREADS=str(n_threads))
    proc = subprocess.run(args, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"sbi_batch failed (exit {proc.returncode}):\n{proc.stderr}\n{proc.stdout}")

    #   Assembling Result
    return {"out_file": str(out_file), "manifest_file": manifest_path(out_file),
            "params_file": str(params_file), "stdout": proc.stdout}


# ============================================================
#   READING RECORDS
# ============================================================

#   Read the Records of a Batch
def read_records(out_file):
    """Read a JSONL output into a list of dicts, one per run.

    Notes
    -----
    JSON null becomes None. Each record is checked for the four blocks
    and the schema version.
    """

    #   Parse each line
    records = []
    with open(out_file) as io:
        for k, line in enumerate(io, start=1):
            rec = json.loads(line)
            for block in ("theta", "summary", "target", "meta"):
                if block not in rec:
                    raise ValueError(f"record {k} lacks the {block} block")
            if rec["meta"]["schema_version"] != SCHEMA_VERSION:
                raise ValueError(f"record {k}: schema version {rec['meta']['schema_version']}, expected {SCHEMA_VERSION}")
            records.append(rec)
    return records


#   Read a Batch Manifest
def read_manifest(manifest_file):
    """Read the manifest of batch settings written beside the records."""

    #   Parse
    with open(manifest_file) as io:
        return json.load(io)


#   Helper for records_to_arrays: JSON value to float
def _as_float(x):
    """Return NaN for None (JSON null), otherwise float(x)."""
    return np.nan if x is None else float(x)


#   Arrays for Training
def records_to_arrays(records, summary_fields=("K", "t_i", "peak_incidence", "frac_at_inflection"),
                      classes=("established",)):
    """Turn records into arrays for BayesFlow.

    Parameters
    ----------
    records : list of dict
        From read_records.
    summary_fields : tuple of str
        Summary quantities that form the conditioning vector. The default
        is a placeholder for the open decision on the summary vector
        (roadmap open decision 4).
    classes : tuple of str
        Run classes to keep; default established runs only, which
        conditions on establishment as the default extinct policy does.

    Returns
    -------
    dict of numpy.ndarray
        thetas (N, 12) in THETA_FIELDS order, with NaN for a null
        immunity parameter; summaries (N, d); in_target (N,) as float,
        NaN where the record has null; draw_id (N,); replicate (N,).
    """

    #   Keep the requested classes
    kept = [r for r in records if r["summary"]["class"] in classes]

    #   Build the arrays
    thetas = np.array([[_as_float(r["theta"][f]) for f in THETA_FIELDS] for r in kept], dtype=float).reshape(-1, len(THETA_FIELDS))
    summaries = np.array([[_as_float(r["summary"][f]) for f in summary_fields] for r in kept], dtype=float).reshape(-1, len(summary_fields))
    in_target = np.array([_as_float(r["target"]["in_target"]) for r in kept], dtype=float)
    draw_id = np.array([r["meta"]["draw_id"] for r in kept], dtype=int)
    replicate = np.array([r["meta"]["replicate"] for r in kept], dtype=int)
    return {"thetas": thetas, "summaries": summaries, "in_target": in_target,
            "draw_id": draw_id, "replicate": replicate}
