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
This module provides functions to read reference concentrations from simulation results.
"""

import os
from pprint import pformat
from typing import Optional

import numpy as np

from .build import build_fake_mechanism
from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_all_settings_map
from .simulation_run import run_simulation
from .postprocess_gecko import get_gecko_concs
from .postprocess_ssh import get_ssh_concs


# Logger
logger = setup_logger(__name__)


def get_refconcs_from_files(refconc_path: str, options: Optional[dict] = None) -> dict:
    """
    Get reference concentrations from path.
    If path is not given, run & analyze results for fake mechanism
    """

    # Get settings
    settings = get_all_settings_map()
    gnl = settings[SNAME.SETUP]()

    # Check path
    if not refconc_path:
        refconc_path = _get_res_path_from_fake_mech()

    # Read reference concentrations from npz file
    if os.path.isfile(refconc_path):
        logger.info("Loading reference concentrations from npz file: %s ...", refconc_path)
        return np.load(refconc_path, allow_pickle=True)

    # Read reference concentrations from result folders
    if os.path.isdir(refconc_path):
        refconc_dict = {
            "GECKO": get_gecko_concs,
            "SSH": get_ssh_concs,
        }
        if gnl.boxmodel not in refconc_dict:
            raise ValueError(f"Mode {gnl.boxmodel} is not supported.")

        # Update options with defaults
        opts = {"time": "ave", "use_fake": gnl.refconc_wfake, "winit": False, "total": "with"}
        if options:
            opts.update(options)
        logger.info("Options for reference concentrations: %s", opts)

        # Get reference concentrations
        refconc = refconc_dict[gnl.boxmodel](refconc_path, opts)

        # Save as npz file
        if gnl.refconc_save:
            if opts.get("time", None):
                filename = f"{refconc_path}/ref_conc_{opts['time']}.npz"
            else:
                filename = f"{refconc_path}/ref_conc.npz"
            np.savez(filename, **refconc)
            logger.info("Saved reference concentrations to %s.", filename)

        return refconc

    raise ValueError(f"Not able to find reference files or folder: {refconc_path}")


def _get_res_path_from_fake_mech() -> str:
    """
    If need:
    Build fake mechanism (with fake radicals recording radical production)
    Run fake mechanism and return result folder.
    Those results are used as reference concentrations.
    """

    # Initialize settings
    general = get_all_settings_map()[SNAME.SETUP]()

    # Get fake mechanism info
    # mechansim name
    fmechname = f"{general.mech_name}FA"
    # mechanism result folder
    fmech_res_path = os.path.join(general.path_sav_res, f"Results_{fmechname}")

    # Check if fake mechanism result folder exists
    if os.path.isdir(fmech_res_path):
        logger.info("Fake mechanism result folder exists: %s", fmech_res_path)
        return fmech_res_path

    # Build fake mechanism if need
    build_fake_mechanism(fmechname, True)

    # Get simulation settings
    run_setting = {
        "mech_names": [fmechname],
        "out_mode": 2,
    }
    # Run simulation with fake mechanism
    run_simulation(run_setting)

    # Return result folder
    return fmech_res_path


def get_concs_for_kinetics(concs_in: dict, species: list, basic_dict: dict) -> dict:
    """Get the reference concentration with ro2 pool concentrations and basic species."""

    concs_k = {}  # Concentrations for kinetics

    # Make sure contains concs for basic species
    empty_sp = []
    for s in basic_dict:
        if s in concs_in:
            concs_k[s] = concs_in[s]
        else:
            concs_k[s] = 0.0
            empty_sp.append(s)

    logger.info("Get reference concentrations for kinetics:\n%s", pformat(concs_k))
    if empty_sp:
        logger.info("Added reference concentrations for # %s basic species:\n%s", len(empty_sp), empty_sp)

    # Add ro2 pool concentrations
    concs_k.update(_get_ro2_concs(concs_in, species))

    return concs_k


def _get_ro2_concs(concs_in: dict, species: list) -> dict:
    """Get ro2 pool concentrations and counts"""
    nro2 = 9  # Number of RO2 groups
    concs_ro2 = {i + 1: 0.0 for i in range(nro2)}
    counts_ro2 = {i + 1: 0 for i in range(nro2)}
    for s in species:
        if s.status < 1 or s.RO2 < 1 or s.name not in concs_in:  # Invalid
            continue
        concs_ro2[s.RO2] += concs_in[s.name]
        counts_ro2[s.RO2] += 1

    tot_ro2 = sum(concs_ro2.values())  # Total ro2 concentration
    nro2 = sum(counts_ro2.values())  # Total number of ro2 species

    # Update ro2 pool names in dict
    concs_ro2 = {f"ro2s{i}": conc for i, conc in concs_ro2.items()}
    concs_ro2.update({"ro2s": tot_ro2})  # Add total

    return concs_ro2


def get_concs_for_species(concs_in: dict, species: list) -> dict:
    """Get the reference concentration for current species list."""
    if not concs_in:
        return {}

    concs_s = {}
    for s in species:
        if s.status <= 0:  # Invalid species
            continue
        sname = s.name
        concs = concs_in.get(sname, 0.0)
        for i, sps in s.reductions.items():
            if i not in ["lp", "rp"]:
                continue
            for p in sps:
                concs += concs_in.get(p, 0.0)
        concs_s[sname] = concs

    # Make sure contains ro2 concs
    concs_s.update(_get_ro2_concs(concs_in, species))

    return concs_s
