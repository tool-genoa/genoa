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
GENOA v3 main driver module.
"""

import os
import time

from .build import build_new_mechanism
from .constants import SNAME
from .logger import setup_logger
from .reduction_run import run_reduction
from .setting_update import set_genoa_parameters
from .simulation_run import run_simulation
from .postprocess_run import run_postprocess


# Logger
logger = setup_logger(__name__)


def _record_genoa_info() -> str:
    """Record the program information."""

    pinfos = ["\n\n" + "=" * 75 + "\n"]
    pinfos.append(" " * 5 + "GENOA v3.0: the GENerator of reduced Organic Aerosol mechanism\n")
    pinfos.append(" " * 5 + "Copyright (C) 2025 CE-CERT (UCR) - ACOM (NCAR) - CEREA (ENPC) - INERIS")
    pinfos.append(" " * 5 + "GENOA is distributed under GPLv3.\n")
    pinfos.append("Documentation and source:\n" + " " * 5 + "https://zhizhaow.github.io/genoa-v3/\n")
    pinfos.append("=" * 75 + "\n")
    return "".join(pinfos)


def _execute_basic_action(actions: list) -> list:
    """Executes the basic actions."""

    # Action dictionary
    action_dict = {
        SNAME.NEW: build_new_mechanism,
        SNAME.SML: run_simulation,
        SNAME.PST: run_postprocess,
    }

    remain = []
    for act in actions:
        if act in action_dict:
            logger.info("Start action: %s ...", act)
            t0 = time.perf_counter()  # Start time
            action_dict[act]()
            logger.info("Action %s is finished in %.1f seconds.", act, time.perf_counter() - t0)
        else:
            remain.append(act)

    return remain


def run_genoa() -> None:
    """Main function for running the GENOA program."""

    # Record the program information
    logger.info("%sCurrent job id: %s", _record_genoa_info(), os.getpid())

    # Get starting time
    tbeg_genoa = time.perf_counter()

    # Read the configuration file
    actions = set_genoa_parameters()

    if not actions:
        logger.warning("No GENOA action is specified. Stop.")
        return

    # Check basic actions in order
    actions = _execute_basic_action(actions)

    # Perform reduction
    if actions:
        actions = run_reduction(actions)

    # Check remaining actions
    if actions:
        logger.error("Unknown actions: %s", actions)

    # Get ending time
    logger.info("GENOA is finished. Total time used: %.1f seconds.", time.perf_counter() - tbeg_genoa)
