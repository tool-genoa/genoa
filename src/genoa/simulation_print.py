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
This module contains functions to record simulation information to screen or a file.
"""

import numpy as np

from .logger import setup_logger
from .simulation_init import RunSmlSetting
from .simulation_tconc import get_genoa_concs


# Logger
logger = setup_logger(__name__)


def record_simulation_info(opts: RunSmlSetting) -> None:
    """Record simulation information to screen or a file if need."""

    if opts.loge and not opts.loge.muted:  # Error log file

        write_einfo_to_loge(opts)

        log_error_statistics(opts)  # Error statistics


def log_error_statistics(opts: RunSmlSetting) -> None:
    """Record error statistics"""

    if opts.nerr < 1:
        return

    if not opts.err_arr:  # Check error array
        raise ValueError("Error array is not available.")

    # Check if containd np.nan
    if opts.err_arr.has_nan:
        logger.error("Error array contain np.nan. Got:\n%s", opts.err_arr.arrs)
        return

    if np.any(opts.err_arr.arrs == -1):
        logger.debug("Error array contain -1. Got:\n%s", opts.err_arr.arrs)

    if opts.err_arr.emax is None:
        opts.err_arr.update_statistics()

    # nchem, nsml, nerr, nsps = opts.err_arr.arrs.shape

    # Write max/ave errors for all simulations
    parts = ["\nError statistics:", "Mech\tErr\tSps\tMax\tAve"]  # Header
    err_max = np.max(opts.err_arr.arrs, axis=1)  # Max errors
    err_ave = np.mean(opts.err_arr.arrs, axis=1)  # Ave errors

    for i, imech in enumerate(opts.mech_names):
        for j in range(opts.nerr):
            for k, sp in enumerate(opts.err_sps):
                parts.append(f"{imech}\t{j}\t{sp}\t{err_max[i, j, k]:.5f}\t{err_ave[i, j, k]:.5f}")

    # Write max & ave errors for all
    parts.append(f"\nTotal max & ave errs:\t{opts.err_arr.emax}\t{opts.err_arr.eave}")
    if opts.err_arr.emax > 0:
        parts.append(f"Find # {len(opts.err_arr.emax_locs)} max errors:")
        for i, iloc in enumerate(opts.err_arr.emax_locs):
            parts.append(f"# {i+1}\t{opts.get_err_loc(iloc)}")

    # Print and write to loge
    err_info = "\n".join(parts)
    logger.info(err_info)
    opts.loge.write("\n\n" + err_info)
    opts.loge.close()

    # Write all errors if small
    if not opts.wid:
        logger.info("%s", _get_err_info(opts))


def _get_err_info(opts: RunSmlSetting) -> str:
    """Record all errors for each simulation."""
    if opts.nerr < 1:
        return ""

    if not opts.err_arr:
        raise ValueError("Error array is not available.")

    parts = ["\n\nAll Errors for " + " ".join(opts.err_sps)]
    initids = [opts.init_sets[initid] for initid in opts.initids] if opts.init_sets else ["default"]

    for i, imech in enumerate(opts.mech_names):
        parts.append(f"\n{imech}")
        for j1, inml in enumerate(opts.labels):
            for j2, initid in enumerate(initids):
                j = j1 * opts.ninit + j2
                for k in range(opts.nerr):
                    line = [f"  {inml}\t{initid}\terr {k}"]
                    line.extend(f"{opts.err_arr.arrs[i, j, k, s]:.5f}" for s in range(opts.nsps))
                    parts.append("\t".join(line))

    return "\n".join(parts)


def write_einfo_to_loge(opts: RunSmlSetting) -> None:
    """Save detailed error information in tab-separated format."""

    if opts.nerr < 1 or not opts.loge or opts.loge.muted:
        return
    if not opts.err_arr:
        raise ValueError("Error array is not available.")

    # Get concentrations
    gconcs = get_genoa_concs(opts)

    # Write header
    opts.loge.write(f"\nRun id\tMech\tNml\tInit\tConcs of {'; '.join(opts.err_sps)}\t|\tErrors (# {opts.nerr})\n")

    if opts.init_sets:
        initids = [opts.init_sets[initid] for initid in opts.initids]
    else:
        initids = ["default"]
    for i, mech in enumerate(opts.mech_names):
        n = 0  # Run id
        pinfos = []
        for j1, inml in enumerate(opts.labels):
            pinfo0 = f"{n}\t{mech}\t{inml}"  # Run id, mechanism, nml
            for j2, initid in enumerate(initids):
                j = j1 * opts.ninit + j2
                # "\t" to separate
                parts = [pinfo0, initid, ", ".join(f"{v:6.3E}" for v in gconcs[i, j1, j2, :]), "|"]
                for k in range(opts.nerr):
                    parts.append(", ".join(f"{k:.4f}" for k in opts.err_arr.arrs[i, j, k, :]))
                pinfos.append("\t".join(parts))
                n += 1
        opts.loge.list_to_file(pinfos)
