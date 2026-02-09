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
This module contains the main function to run simulations using the given mechanism(s).
"""

from typing import Union

from .logger import setup_logger
from .simulation_cmp import set_boxmodel
from .simulation_core import write_namelists, launch_smls_in_para
from .simulation_init import RunSmlSetting, setup_run_sml
from .simulation_print import record_simulation_info
from .simulation_gecko import get_gck_sml_inputs
from .simulation_ssh import get_ssh_sml_inputs


# Logger
logger = setup_logger(__name__)


def run_simulation(run_settings: Union[RunSmlSetting, None, dict] = None) -> RunSmlSetting:
    """Run a chemical mechanism using the GENOA v3.0 algorithm."""

    sml_opts = setup_run_sml(run_settings)
    run_dict = {
        "SSH": run_gecko_or_ssh,
        "GECKO": run_gecko_or_ssh,
    }
    if sml_opts.box not in run_dict:
        raise ValueError(f"Unknown box model mode: {sml_opts.box}")

    sml_opts = run_dict[sml_opts.box](sml_opts)

    record_simulation_info(sml_opts)

    return sml_opts


def run_gecko_or_ssh(opts: RunSmlSetting) -> RunSmlSetting:
    """Run 0D simulations with box models: GECKO-A or SSH-aerosol."""

    # Compilation if need
    if not opts.box_exec:
        opts.box_exec = set_boxmodel(opts.paths.work, opts.box)

    # Update namelist files if need
    if not opts.labels:
        write_namelists(opts)

    # Get simulation inputs
    input_dict = {
        "GECKO": get_gck_sml_inputs,
        "SSH": get_ssh_sml_inputs,
    }
    sml_inputs_chem = input_dict[opts.box](opts)

    # Launch simulations
    errs = launch_smls_in_para(opts, sml_inputs_chem)

    # Record errors if need
    if opts.nerr:
        opts.update_errors(errs)

    return opts
