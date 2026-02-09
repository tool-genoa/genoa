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
This module provides functions to read concentration files used in the GENOA reductions.
"""

import numpy as np

from .logger import setup_logger
from .simulation_init import RunSmlSetting
from .simulation_gecko import get_gck_concs_paths
from .simulation_ssh import get_ssh_concs_paths
from .utils import get_data_from_concs_files


# Logger
logger = setup_logger(__name__)


def get_genoa_concs(opts: RunSmlSetting) -> np.ndarray:
    """Get genoa concs with input settings."""

    path_func = {
        "GECKO": get_gck_concs_paths,
        "SSH": get_ssh_concs_paths,
    }
    if opts.box not in path_func:
        raise ValueError(f"Cannot find box model type {opts.box} for genoa concs. Check.")
    # Update SOA path if needed
    if not opts.soa_path:
        opts.soa_path = opts.paths.res
    # Get paths to genoa concs
    paths = path_func[opts.box](opts)
    return _get_gconcs_from_paths(opts, paths)


def _get_gconcs_from_paths(opts: RunSmlSetting, paths: list) -> np.ndarray:
    """Get and reshape genoa concs with input paths."""

    gconcs = get_data_from_concs_files(paths)  # In shapes of [nchem * nml * init, nsps]
    nall, nsps = gconcs.shape

    # check size
    if nall != opts.nchem * opts.nnml * opts.ninit or nsps != opts.nsps:
        raise ValueError(
            f"Invalid genoa concs shape: {gconcs.shape}."
            + f"  Should be w/ {opts.nchem} chems, {opts.nnml} nmls, {opts.ninit} inits, and {opts.nsps} species ->"
            + f"  {opts.nchem * opts.nnml * opts.ninit} * {opts.nsps}."
        )
    return gconcs.reshape(opts.nchem, opts.nnml, opts.ninit, opts.nsps)
