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
This module contains basic functions to run a simulation and get error from the output file.
"""

import os
import subprocess
from multiprocessing import Pool, Manager
from typing import Optional

import numpy as np

from .constants import SNAME
from .logger import setup_logger
from .simulation_basic import get_default_namelist, write_new_namelist
from .simulation_init import RunSmlSetting
from .setting_global import get_attrs_from_settings
from .simulation_gecko import prepare_gck_nmls
from .simulation_ssh import prepare_ssh_nmls


# Logger
logger = setup_logger(__name__)


def launch_smls_in_para(opts: RunSmlSetting, sml_in_chem: list) -> list:
    """Launch GECKO-A simulations with the given simulation inputs and return errors."""

    if opts.wlock:  # Run w/ lock
        return run_smls_wlock(opts, sml_in_chem)

    # Run w/o lock
    sml_inputs = []
    for sml_in in sml_in_chem:
        sml_inputs.extend(sml_in)
    with Pool(opts.npara) as pool:
        errs = pool.starmap(run_one_simulation, sml_inputs)
    return errs


def write_namelists(opts: RunSmlSetting) -> None:
    """Write GECKO-A namelists and update operation ids (labels, resids, nmlids) in the settings."""

    nml_func = {
        "GECKO": prepare_gck_nmls,
        "SSH": prepare_ssh_nmls,
    }
    if opts.box not in nml_func:
        raise ValueError(f"Cannot find box model type {opts.box} for namelist generation. Check.")

    # Get default GECKO-A namelist
    namelist = get_attrs_from_settings({SNAME.BOX[opts.box]: "namelist"})
    nml_dict = get_default_namelist(namelist)

    # Check nmlpath
    if not os.path.exists(opts.paths.nml):
        os.makedirs(opts.paths.nml)

    # Get namelist inputs for different model
    pool_inputs = nml_func[opts.box](nml_dict, opts)

    # Writting new namelists
    with Pool(opts.npara) as pool:
        pool.starmap(write_new_namelist, pool_inputs)

    logger.info("Finished updating # %d %s namelist files in %s.", len(opts.labels), opts.box, opts.paths.nml)


def run_smls_wlock(opts: RunSmlSetting, sml_in_chem: list) -> list:
    """Run GECKO-A simulations with lock and return errors."""

    if not opts.wid:
        raise ValueError("Cannot run w/ lock without id. Check.")

    errs, chemids, is_run = [], [], 1
    # Group mechanisms by scores
    for mechids in sort_mechid_by_scores(opts):
        if not mechids:
            raise ValueError("No mechanism ids found. Check.")
        if is_run > 0:  # Run simulations
            # logger.info("Running sml_wlock w/ %d mechanism(s): %s ...", len(mechids), mechids)
            # Get sml inputs
            sml_inputs = []
            for i in mechids:
                chemids.append(opts.chemids[i])
                sml_inputs.extend(sml_in_chem[i])

            # Run w/ lock
            with Manager() as manager:
                stop_dict = manager.dict()
                lock = manager.Lock()
                sml_inputs_wlock = [i + (stop_dict, lock) for i in sml_inputs]
                with Pool(opts.npara) as pool:
                    ierrs = pool.starmap(run_one_simulation, sml_inputs_wlock)
                    errs.extend(ierrs)
                    # Find mechanism & check if need more runs
                    if find_valid_mechs(mechids, ierrs, stop_dict, opts):
                        is_run -= 1

        else:  # Add errors as -1 for non-run simulations
            # logger.info("Skip mechanism(s): %s.", mechids)
            chemids.extend([opts.chemids[i] for i in mechids])
            ierrs = np.zeros((len(sml_in_chem[0]) * len(mechids),) + errs[0].shape)
            errs.extend(np.full_like(ierrs, -1))

    # Update chemids
    opts.chemids = chemids
    return errs


def run_one_simulation(
    cmd: str,  # Command to run simulation
    totof: str,  # Output file w/ print info
    outf: str,  # Output file for errors
    nerrs: Optional[list] = None,  # [nerr, nsps] or None if nerr = 0
    # For stop
    chemid: Optional[str] = None,
    err_checks: Optional[list] = None,
    stop_dict: Optional[dict] = None,
    lock: Optional[bool] = None,
) -> np.ndarray:
    """Run one simulation and return errors."""

    # Check stop signal before starting the simulation
    if nerrs and lock is not None and chemid in stop_dict:
        return np.full(nerrs, -1)  # Return a placeholder number

    # Execute linux command to run a simluation
    try:
        with open(totof, "w", encoding="utf-8") as f:
            f.write(f"Run command: {cmd}\n")
            subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("Error running simulation: %s", e)
        return np.full(nerrs, np.nan)

    if not nerrs:  # No errors expected
        return np.array([])

    # Extract errors from output file
    errs = []
    try:
        # Read errors from output error file
        with open(outf, "r", encoding="utf-8") as f:
            for line in f.readlines():
                if line.strip():
                    ierrs = [float(i) for i in line.split() if i.strip()]
                    if ierrs:
                        errs.append(ierrs)
        # Convert errors
        errs = np.array(errs, dtype=float)
    except (FileNotFoundError, ValueError) as e:
        logger.error("Error processing error file %s or errors %s: %s", outf, errs, e)
        return np.full(nerrs, np.nan)  # Return stopper

    # Update stop signal
    if lock is not None and err_checks is not None:
        with lock:
            # Check if err is valid
            if np.any(np.isnan(errs) | (errs == -1)):
                stop_dict[chemid] = 1
                # logger.info("Stop smls for %s: Found -1 or np.nan.", chemid)
                return errs
            err_check = err_checks[0]  # For Max err

            # Check max err tolerance
            err_maxs = np.max(errs, axis=1)
            exceed_max = err_maxs > err_check[0]
            if np.any(exceed_max):
                stop_dict[chemid] = 2
                # logger.info("Stop smls for %s: Max err %s > %s: %s.", chemid, err_maxs, err_check[0], exceed_max)
                return errs
            # Check max err increase
            derr = err_maxs - err_check[2]
            exceed_derr = derr > err_check[1]
            if np.any(exceed_derr):
                stop_dict[chemid] = 3
                # logger.info("Stop smls for %s: Max err up %s > %s: %s.", chemid, derr, err_check[1], exceed_derr)
                return errs

    return errs


def find_valid_mechs(mechids: list, ierrs: list, stop_dict: dict, opts: RunSmlSetting) -> bool:
    """Check mean errors to find valid mechanisms."""

    # logger.info("Checking mean errors for mechanisms w/ ids: %s and stop_dict: %s ...", mechids, stop_dict)

    # Return if all mechanisms are in stop_dict
    if len(stop_dict) == len(mechids):
        # logger.info("All mechanisms are in stop_dict. Skip.")
        return False

    # Error checks
    merr_tol, derr_tol, merr0 = opts.err_checks[1][0], opts.err_checks[1][1], opts.err_checks[1][2]
    # Reshape errors
    merrs = np.moveaxis(np.array(ierrs), 2, 1)
    n1, n2, n3 = merrs.shape
    merrs = merrs.reshape((len(mechids), n1 * n2 // len(mechids), n3))
    # Mean errors
    merrs = np.mean(merrs, axis=1)

    for i, mid in enumerate(mechids):
        chemid = opts.chemids[mid]
        if chemid in stop_dict:
            # logger.info("Skip mechanism %s in stop_dict.", chemid)
            continue
        # Check mean errors
        exceed_merr = merrs[i] > merr_tol
        if np.any(exceed_merr):
            stop_dict[chemid] = 4
            # logger.info("Stop mechanism %s: Mean error %s > %s: %s", chemid, merrs[i], merr_tol, exceed_merr)
            continue
        # Check mean error increase
        derr = merrs[i] - merr0
        exceed_derr = derr > derr_tol
        if np.any(exceed_derr):
            stop_dict[chemid] = 5
            # logger.info("Stop mechanism %s: Mean error increase %s > %s: %s", chemid, derr, derr_tol, exceed_derr)
            continue

        # logger.info("Find Mechanism %s is valid!", chemid)
        return True

    return False


def sort_mechid_by_scores(opts: RunSmlSetting) -> list:
    """Sort mechanisms by scores and return grouped mechanism ids"""

    if not (opts.scores and opts.err_checks):
        raise ValueError("No scores or error check found. Check.")

    if opts.nchem == 1:  # wid must be True
        return [[0]]

    # Sort indices of scores in descending order
    sorted_indices = sorted(range(len(opts.scores)), key=lambda i: opts.scores[i], reverse=True)

    # Group mechanisms by scores
    mechs_list, score = [], -999
    for idx in sorted_indices:
        s = opts.scores[idx]
        if s != score:  # New group
            mechs_list.append([idx])
            # if score != -999:
            # logger.info("Grouped # %d mechs: %s  with score %s.", len(mechs_list[-2]), mechs_list[-2], score)
            score = s
        else:  # Add to the same group
            mechs_list[-1].append(idx)
    # logger.info("Grouped # %d mechs: %s  with score %s.", len(mechs_list[-1]), mechs_list[-1], score)
    return mechs_list
