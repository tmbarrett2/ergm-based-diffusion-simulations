"""Prior over the 12 simulator parameters (roadmap Section 5.1).

Jonathan H. Morgan, Ph.D.
26 October 2026

Draws parameter vectors from the weakly informative prior and writes them
as the JSON array that the Julia sbi_batch mode reads. The prior is the
initial specification; UT-1.5 (the prior predictive check) decides whether
it is revised.
"""

import json
import math

import numpy as np


# ============================================================
#   PARAMETER LAYOUT
# ============================================================

#   Order of the 12 parameters, matching the Julia record (THETA_FIELDS)
THETA_FIELDS = (
    "inf_r", "b_int", "b_close", "b_cxn_peers", "b_cxn_total", "b_cxn_symp",
    "b_cls_x_smp", "p_symp", "rec_scale", "rec_shape", "imm_scale", "imm_shape",
)

#   Prior specification: name -> (distribution, parameters)
PRIOR_SPEC = {
    "inf_r": ("lognormal", {"median": 0.02, "sigma": 0.5}),
    "b_int": ("normal", {"mean": -2.0, "sd": 1.0}),
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


# ============================================================
#   SAMPLING
# ============================================================

#   Helper for sample_prior: truncated normal by rejection
def _truncated_normal(rng, mean, sd, lower, size):
    """Draw from Normal(mean, sd) truncated below at lower.

    Parameters
    ----------
    rng : numpy.random.Generator
    mean, sd, lower : float
    size : int

    Returns
    -------
    numpy.ndarray of shape (size,)

    Notes
    -----
    Rejection sampling: draws that fall below the bound are redrawn. With
    the bound one standard deviation below the mean, about 84% of draws are
    accepted, so a few passes suffice.
    """

    #   Draw, then redraw the rejected entries until none remain
    out = rng.normal(mean, sd, size)
    bad = out < lower
    while bad.any():
        out[bad] = rng.normal(mean, sd, bad.sum())
        bad = out < lower
    return out


#   Helper for sample_prior: one parameter's draws
def _draw_parameter(rng, dist, par, size):
    """Draw `size` values from one prior distribution.

    Parameters
    ----------
    rng : numpy.random.Generator
    dist : str
        'lognormal', 'normal', 'truncated_normal', 'beta', or 'gamma'.
    par : dict
        Distribution parameters, as in PRIOR_SPEC.
    size : int

    Returns
    -------
    numpy.ndarray of shape (size,)

    Notes
    -----
    LogNormal(ln m, s) is parameterized by its median m. Gamma(shape,
    scale) uses the scale form, so Gamma(4, 0.5) has mean 2.
    """

    #   Dispatch on the distribution name
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
    raise ValueError(f"unknown distribution: {dist}")


#   Prior Draws
def sample_prior(n, seed, first_draw_id=1):
    """Draw n parameter vectors from the prior.

    Parameters
    ----------
    n : int
        Number of draws.
    seed : int
        Seed for numpy.random.default_rng, so a batch of draws is
        reproducible.
    first_draw_id : int, default 1
        draw_id of the first draw; later draws count up from it.

    Returns
    -------
    list of dict
        One dict per draw: draw_id and the 12 parameters as Python floats,
        in THETA_FIELDS order.

    Notes
    -----
    Parameters are drawn one at a time across all n draws, in
    THETA_FIELDS order, so the draws for a given seed do not depend on
    how the result is later split into batches. The prior always draws
    immunity, so imm_scale and imm_shape are never null here.
    """

    #   Validation
    if n < 1:
        raise ValueError("n must be at least 1")

    #   Draw each parameter across all draws
    rng = np.random.default_rng(seed)
    columns = {name: _draw_parameter(rng, *PRIOR_SPEC[name], n) for name in THETA_FIELDS}

    #   Assemble one dict per draw
    draws = []
    for k in range(n):
        draw = {"draw_id": first_draw_id + k}
        for name in THETA_FIELDS:
            draw[name] = float(columns[name][k])
        draws.append(draw)
    return draws


# ============================================================
#   OUTPUT
# ============================================================

#   Write a Parameter File for sbi_batch
def write_params(draws, path):
    """Write draws as the JSON array that sbi_batch reads.

    Parameters
    ----------
    draws : list of dict
        From sample_prior, or hand-built with the same keys.
    path : str or pathlib.Path

    Returns
    -------
    int
        Number of draws written.

    Notes
    -----
    Python's json writes floats with repr, the shortest form that reads
    back to the same double, so Julia receives the exact values drawn
    here (checked by the round-trip test).
    """

    #   Validation: unique draw_id and every parameter present
    ids = [d["draw_id"] for d in draws]
    if len(set(ids)) != len(ids):
        raise ValueError("draw_id values must be unique")
    for d in draws:
        missing = [f for f in THETA_FIELDS if f not in d]
        if missing:
            raise ValueError(f"draw {d['draw_id']} lacks {missing}")

    #   Write
    with open(path, "w") as io:
        json.dump(draws, io)
    return len(draws)
