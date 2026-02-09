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
This module contains the main function to run the testing process.
"""

import os

from typing import Union

from .logger import setup_logger
from .setting_global import get_all_settings_map
from .setting_init import TestingOption
from .simulation_run import run_simulation
from .testing_init import get_sml_setting_dict, setup_testing


# Logger
logger = setup_logger(__name__)


def run_testing(settings_in: Union[TestingOption, dict, None] = None) -> TestingOption:
    """Run GENOA v3 testing or pre-testing process."""

    # Get settings
    settings_map = get_all_settings_map()

    # Testing preparation
    topt = setup_testing(settings_map, settings_in)

    # Run reference mechanism
    if not topt.read_ref:
        logger.info("Running simulations with the reference mechanism: %s ...", topt.ref_mech_name)
        rsml = run_simulation(get_sml_setting_dict(topt, "ref"))
        topt.read_ref = os.path.join(rsml.paths.res, topt.ref_mech_name)  # Ref result folder
        logger.info("Finished reference simulations and saved results in: %s", topt.read_ref)

    # Run targeted mechanism
    topt.sml = run_simulation(topt.sml if topt.sml else get_sml_setting_dict(topt, "test"))

    return topt
