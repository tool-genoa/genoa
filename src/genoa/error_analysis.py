# ================================================================================
#
#   GENOA v3: the Generator of Optimized Atmospheric chemical mechanisms
#
#   Copyright (C) 2025 CE-CERT (UCR) - ACOM (NCAR) - CEREA (ENPC) - INERIS
#
#   Distributed under the terms of the GNU General Public License v3 (GPLv3).
#
# ================================================================================

"""
This module provides functions to analysize reduction errors.

"""


import math

import numpy as np


ERR_TYPE_ALL = {
    "re",
    "mre",
    "mse",
    "rmse",
    "rmse/ave",
    "R",
    "R2",
    "mngb",
    "mnge",
    "mfb",
    "mfe",
    "nmb",
    "nme",
    "fb",
    "fe",
    "few",
    "mfe_24",
    "few_24",
}

EPS = 1e-30


def cmp_error(val: np.ndarray, ref: np.ndarray, indicator: str = "mfe_24") -> float:
    """Compute the error depending on the input indicator, with optional weighting.
    - mfe_[n1, n2, ...]: max fractional error over sections: [0, n1], [n1, n2], ...
    - few_[n1, n2, ...]: weighted fractional error over sections: [0, n1], [n1, n2], ...
    """

    val, ref = np.array(val), np.array(ref)

    # Return if equal
    if np.array_equal(val, ref):
        return 0.0

    # Check lengths
    if len(val) == 0 or len(val) != len(ref):
        raise ValueError("Lengths of val is 0 or " + f"lengths are not equal: {len(val)} vs {len(ref)}")

    # Get error maping
    indicator_funcs = {
        "re": compute_re,
        "mre": compute_mre,
        "mse": compute_mse,
        "rmse": compute_rmse,
        "rmse/ave": compute_rmse_ave,
        "R": compute_correlation,
        "R2": compute_r_squared,
        "mngb": compute_mngb,
        "mnge": compute_mnge,
        "mfb": compute_mfb,
        "mfe": compute_mfe,
        "mfew": compute_mfew,
        "nmb": compute_nmb,
        "nme": compute_nme,
        "fb": compute_fb,
        "fe": compute_fe,
        "few": compute_few,
    }

    # Compute error
    if indicator in indicator_funcs:
        return indicator_funcs[indicator](val, ref)

    if indicator.startswith("mfe_"):
        ncuts = [int(i) for i in indicator.split("_")[1].split(",") if i.isdigit()]
        ncuts = [0] + ncuts + [len(ref)]
        return compute_mfe_sections(val, ref, ncuts)

    if indicator.startswith("few_"):
        ncuts = [int(i) for i in indicator.split("_")[1].split(",") if i.isdigit()]
        ncuts = [0] + ncuts + [len(ref)]
        return compute_few_sections(val, ref, ncuts)

    raise ValueError(f"Indicator not found: {indicator}")


def compute_errors_between_files(file1: str, file2: str, indicator: str) -> list:
    """Compute errors between two files."""
    # Read files
    val = np.loadtxt(file1)
    ref = np.loadtxt(file2)

    # Compute errors for each column
    errors = []
    for i in range(val.shape[1]):
        errors.append(cmp_error(val[:, i], ref[:, i], indicator))
    return errors


def compute_re(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute relative error (re)."""
    return np.mean(np.abs(np.where(ref > EPS, (val - ref) / ref, 0)))


def compute_mre(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute maximum relative error (mre)."""
    return np.max(np.abs(np.where(ref > EPS, (val - ref) / ref, 0)))


def compute_mse(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean square error (mse)."""
    return np.mean((val - ref) ** 2)


def compute_rmse(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute root mean square error (rmse)."""
    mse = compute_mse(val, ref)
    return math.sqrt(mse)


def compute_rmse_ave(val, ref) -> float:
    """Compute root mean square error divided by average (rmse/ave)."""
    rmse = compute_rmse(val, ref)
    return rmse / np.mean(ref)


def compute_correlation(val, ref) -> float:
    """Compute correlation coefficient (R)."""
    return np.corrcoef(val, ref)[0, 1]


def compute_r_squared(val, ref) -> float:
    """Compute R squared (R2)."""
    sstot = np.sum((val - np.mean(val)) ** 2)
    ssres = np.sum((val - ref) ** 2)
    return 1 - (ssres / sstot) if sstot != 0 else 0


def compute_mngb(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean normalized gross bias (mngb)."""
    return np.mean(np.where(ref > EPS, (val - ref) / ref, 0))


def compute_mnge(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean normalized gross error (mnge)."""
    return np.mean(np.abs(np.where(ref > EPS, (val - ref) / ref, 0)))


def compute_mfb(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean fractional bias (mfb)."""
    return np.mean(np.where(val + ref > EPS, (val - ref) / (val + ref) * 2, 0))


def compute_mfe(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean fractional error (mfe)."""
    return np.mean(np.where(val + ref > EPS, (val - ref) / (val + ref) * 2, 0))


def compute_mfew(val: np.ndarray, ref: np.ndarray) -> float:
    """Compute mean fractional error weighted (mfew)."""
    return np.mean(np.abs(np.where(val + ref > EPS, (val - ref) / (val + ref) * 2 * ref, 0)))


def compute_nmb(val, ref) -> float:
    """Compute normalized mean bias (nmb)."""
    return np.sum(val - ref) / np.sum(ref)


def compute_nme(val, ref) -> float:
    """Compute normalized mean error (nme)."""
    return np.sum(np.abs(val - ref)) / np.sum(ref)


def compute_fb(val, ref) -> float:
    """Compute fractional bias (fb)."""
    numer = np.sum(val - ref)
    denom = np.sum(val + ref)
    return 2 * numer / denom if denom != 0 else 0


def compute_fe(val, ref) -> float:
    """Compute fractional error (fe)."""
    numer = np.sum(np.abs(val - ref))
    denom = np.sum(val + ref)
    return 2 * numer / denom if denom != 0 else 0


def compute_few(val, ref) -> float:
    """Compute weighted fractional error (few)."""
    numer = np.sum(np.abs(val - ref) * ref)
    denom = np.sum((val + ref) * ref)
    return 2 * numer / denom if denom != 0 else 0


def compute_mfe_sections(val, ref, ncuts) -> float:
    """Compute max fractional error over sections (mfe_)."""
    err = 0.0
    for i0, i1 in zip(ncuts[:-1], ncuts[1:]):
        numer = np.sum(np.abs(val[i0:i1] - ref[i0:i1]))
        denom = np.sum(val[i0:i1] + ref[i0:i1])
        section_err = 2 * numer / denom if denom != 0 else 0.0
        err = max(err, section_err)
    return err


def compute_few_sections(val, ref, ncuts) -> float:
    """Compute weighted fractional error over sections (few_)."""
    err = 0.0
    total_ref_sum = np.sum(ref)
    for i0, i1 in zip(ncuts[:-1], ncuts[1:]):
        numer = np.sum(np.abs(val[i0:i1] - ref[i0:i1]) * ref[i0:i1])
        denom = np.sum((val[i0:i1] + ref[i0:i1]) * ref[i0:i1])
        total = np.sum(ref[i0:i1])
        section_err = 2 * numer / denom * total if denom != 0 else 0.0
        err += section_err
    return err / total_ref_sum if total_ref_sum != 0 else 0.0
