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
This module is the main script to run the reduction process.
"""
import time

from .constants import SNAME
from .logger import setup_logger
from .threshold_reduction import run_tbr
from .training_run import run_para_training
from .testing_run import run_testing


# Logger
logger = setup_logger(__name__)


def run_reduction(actions: list) -> list:
    """
    Perform reduction based on input action dictionary
    """
    # Record total reduction time
    tbeg = time.perf_counter()
    rdc_dict = None
    logger.info("Start GENOA Reduction with actions: %s", actions)

    # Threshold-based reduction
    if SNAME.TBR in actions:
        logger.info("Start threshold-based reduction ...")
        t0 = time.perf_counter()
        rdc_dict = run_tbr()
        t1 = time.perf_counter()
        logger.info("Threshold-based reduction is finished. Time used: %.1f seconds.", t1 - t0)
        actions.remove(SNAME.TBR)
    else:
        logger.info("No threshold-based reduction.")

    # Simulation-based reduction (training)
    if SNAME.TRN in actions:
        logger.info("Start training (simulation-based reduction) ...")
        t0 = time.perf_counter()
        rdc_dict = run_para_training(rdc_dict)
        t1 = time.perf_counter()
        logger.info("Training is finished. Time used: %.1f seconds.", t1 - t0)
        actions.remove(SNAME.TRN)
    else:
        logger.info("No training (simulation-based reduction).")

    # Mechanism evaluation (testing)
    if SNAME.TST in actions:
        logger.info("Start testing (mechanism evaluation) ...")
        t0 = time.perf_counter()
        run_testing(rdc_dict)
        t1 = time.perf_counter()
        logger.info("Testing is finished. Time used: %.1f seconds.", t1 - t0)
        actions.remove(SNAME.TST)
    else:
        logger.info("No Testing.")

    tend = time.perf_counter()
    logger.info("GENOA reduction is finished. Total time used: %.1f seconds.", tend - tbeg)
    return actions
