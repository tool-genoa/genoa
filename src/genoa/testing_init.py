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
This module contains functions to initialize testing options and settings for
running simulations with the testing process.
"""

import os
from typing import Union

from .constants import SNAME
from .folder_path import FolderPath
from .logger import setup_logger
from .record import record_initialization
from .setting_init import TestingOption, check_and_set_attr
from .setting_update import setup_4action


# Logger
logger = setup_logger(__name__)


def setup_testing(settings_map: dict, settings_in: dict) -> TestingOption:
    """Setup testing options and logfile."""

    topt = get_testing_options(settings_map, settings_in)  # Get testing options
    record_initialization(settings_map, topt, SNAME.TST, topt.loge)  # Record

    return topt


def get_testing_options(settings_map: dict, opt_in: Union[TestingOption, dict, None] = None) -> TestingOption:
    """Get testing options from settings or input."""

    # Build new testing options if needed
    if isinstance(opt_in, TestingOption):
        opt = opt_in
    else:
        opt = settings_map[SNAME.TST]()
        if isinstance(opt_in, dict):
            logger.info("Update testing options with input: %s", opt_in)
            for k, v in opt_in.items():
                check_and_set_attr(opt, k, v)

    # Update testing options
    setup_4action(SNAME.TST, opt, settings_map[SNAME.SETUP]())

    # Check if reuse namelists & reference results
    manage_reuse_ref(opt)

    return opt


def manage_reuse_ref(opt: TestingOption) -> None:
    """Manage settings to reuse ref results for testing."""

    if opt.read_ref:  # Read reference results
        if os.path.isdir(opt.read_ref):
            opt.read_ref = os.path.abspath(opt.read_ref)
        else:
            logger.warning("Reference results: %s does not exist!", opt.read_ref)
            opt.read_ref = None


def get_sml_setting_dict(opt: TestingOption, mode: str = "ref") -> dict:
    """Build settings for running simulations with the testing mechanisms."""

    # Get dictionary to initialize settings for running simulations
    sml_dict = {
        "mech_path": opt.mech_path,
        "path_cond": opt.path_cond,
        "conds": opt.conds,
        "wid": True,
    }

    if mode == "ref":  # For reference mechanism
        sml_dict["mech_names"] = [opt.ref_mech_name]
        sml_dict["ref_files_str"] = ""

    elif mode == "test":  # For testing mechanism
        sml_dict["mech_names"] = [opt.mech_name]
        # Add reference results for error calculation
        sml_dict["ref_files_str"] = opt.read_ref
        sml_dict["loge"] = opt.loge

    else:
        raise ValueError(f"Unknown mode {mode}")

    # Get paths
    sml_dict["paths"] = FolderPath(work=os.path.join(opt.path_work, mode))
    sml_dict["paths"].init_with_path_work()

    return sml_dict
