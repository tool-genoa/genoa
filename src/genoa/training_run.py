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
This module contains the main function to run training process for reduction.
"""
from typing import Optional

from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_all_settings_map
from .training_cycle import run_reduction_cycle
from .training_init import setup_training, clean_training


# Logger
logger = setup_logger(__name__)


def run_para_training(settings_dict: Optional[dict] = None) -> dict:
    """Run GENOA v3 training reduction"""

    # Get settings
    settings_map = get_all_settings_map()
    trn_opt = settings_map[SNAME.TRN]()

    # Initialization for training
    setup_training(settings_map, settings_dict)

    icyl_opt = None  # Options for current cycle
    while icyl_opt is None or not icyl_opt.is_stop:
        icyl_opt = run_reduction_cycle(settings_map, icyl_opt)

    # Clean temporary files
    clean_training(settings_map)

    # Return final reduced mechanism
    return {"mech_name": icyl_opt.mechv.name, "mech_path": trn_opt.path_mech}
