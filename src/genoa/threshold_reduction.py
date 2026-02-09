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
This module run threshold-based reduction depending on the given reduction
 parameters and options.
"""

import os
import json
import logging
import time

from .constants import SNAME
from .logger import isout, setup_logger
from .mechanism_in import read_and_update_genoa_mech
from .mechanism_out import mech_output
from .mechanism_update import update_mechanism, balance_carbon_in_reactions
from .rate_constant import update_rate_cst_values
from .reduction_basic import ReductionListInfo
from .reduction_setting import rdc_checks
from .reduction_strategy import reduction_on_mechanism
from .setting_init import GlobalSetting, TbrOption
from .setting_global import get_all_settings_map
from .setting_update import setup_4action
from .simulation_refconc import get_refconcs_from_files
from .simulation_run import run_simulation
from .species import remove_pvoc_partitioning
from .utils import get_size_diff_pinfo, size_difference


# Logger
logger = setup_logger(__name__)


def run_tbr() -> dict:
    """Threshold-based reduction (TBR) on the specified mechanism."""

    # Get threshold-based options
    settings_map = get_all_settings_map()
    gnl = settings_map[SNAME.SETUP]()
    tbr_opt = settings_map[SNAME.TBR]()

    # Update settings
    setup_4action(SNAME.TBR, tbr_opt, gnl)

    # Read mechanism for threshold-based reduction
    read_opt = {"clean": "wsps", "check_rate": True}
    logger.info("Reading mechanism %s for threshold-based reductions w/ %s ...", tbr_opt.ref_mech_name, read_opt)
    reactions, species, _ = read_and_update_genoa_mech(tbr_opt.ref_mech_path, tbr_opt.ref_mech_name, read_opt)

    # Conduct threshold-based reduction
    reactions, species = tbr_on_mechanism(reactions, species, tbr_opt, gnl)

    # Save threshold-reduced mechanism
    _save_tbr_mechanism(reactions, species, gnl.path_sav_mech, tbr_opt.mech_name)

    # Evaluate the threshold-reduced mechanism using simulation if needed
    if tbr_opt.tag_sim:

        t0 = time.perf_counter()  # Start time
        # Update settings for running simulation
        gnl.path_read_mech = gnl.path_sav_mech  # Read from the saved mechanism path
        settings = {"mech_names": [tbr_opt.mech_name], "loge": tbr_opt.loge}
        logger.info("Run simulation for threshold-reduced  mechanism %s w/ settings: %s", tbr_opt.mech_name, settings)

        # run_option
        run_simulation(settings)
        logger.info("Simulations completed in %.1f seconds.", time.perf_counter() - t0)

    return {"mech_name": tbr_opt.mech_name, "mech_path": gnl.path_sav_mech}


def tbr_on_mechanism(reactions: list, species: list, popt: TbrOption, gnl: GlobalSetting) -> list:
    """Threshold-based reduction on the mechanism scheme"""

    logger.info("\nConducting threshold-based reduction on the mechanism scheme ...")

    # Remove the gas-particle partitioning of primary vocs
    if popt.remove_pvoc_gp:
        remove_pvoc_partitioning(species, gnl.primary_vocs)

    # Get reference concentrations
    refconcs = get_refconcs_from_files(gnl.refconc_file)

    # Initialize kvalues
    update_rate_cst_values(reactions)

    # Update mechanism after initialization & record size infos
    sizes, prdc0 = [], "Initialization"
    update_opt = {"clean": "wgen", "merge": True, "tracer": gnl.tracers}
    reactions, species, mech_info = update_mechanism(reactions, species, update_opt)
    sizes.append(mech_info.read_size())

    # Update group list if given
    if popt.lumping:
        _update_lumping_within_tbr(species, popt, gnl)

    # Conduct threshold-based reduction  in order
    prdcs_run, prdcs_tried = popt.get_valid_in_order(), []
    ntry = {prdc[0]: popt.ntry for prdc in prdcs_run}  # No.retry times
    logger.info("Threshold-based reductions with options: %s & ntry: %s", prdcs_run, popt.ntry)

    while prdcs_run:

        # Get current strategy
        prdc = prdcs_run.pop(0)
        logger.info("Threshold reduction with %s after %s ...", prdc, prdc0)
        logger.info("Remaining reduction options: %s, tried: %s", prdcs_run, prdcs_tried)

        # Conduct reduction
        rdc_infos = ReductionListInfo(
            targeted_sps=popt.target_species,
            ref_concs=refconcs,
            given_groups=popt.given_groups,
        )
        reduction_on_mechanism(reactions, species, prdc, rdc_infos)
        prdcs_tried.append(prdc)

        # Update mechanism after reduction
        reactions, species, mech_info = update_mechanism(reactions, species, update_opt)
        sizes.append(mech_info.read_size())

        # Save mechanism if needed
        if popt.tag_save_all:
            savpath = os.path.join(gnl.path_sav_mech, f"tbr_{popt.runid}")
            savname = f"{len(prdcs_tried)}-{prdc[0]}"
            _save_tbr_mechanism(reactions, species, savpath, savname)

        # Update for next reduction
        prdc0 = prdc[0]
        ntry[prdc0] = _update_prdc_list(prdc, prdcs_run, sizes, ntry[prdc0])

    # Carbon balance
    if popt.carbon_balance:
        lost_c = popt.carbon_balance if isinstance(popt.carbon_balance) else "XCLOST"
        logger.info("Balancing carbon in reactions with lost carbon species: %s ...", lost_c)
        balance_carbon_in_reactions(reactions, species, lost_c)
        reactions, species, _ = update_mechanism(reactions, species, update_opt)

    # Record
    if isout(logger, logging.INFO):
        pinfos = [f"\nThreshold-based reduction w/ # {len(prdcs_tried)}: {prdcs_tried}\n  Mechanism size changes:"]

        for i, prdc in enumerate(prdcs_tried):
            pinfos.append(f"\n{prdc}\n{get_size_diff_pinfo(sizes[i], sizes[i + 1])}")

        # Final size
        pinfos.append(f"\nTotal size changes:\n{get_size_diff_pinfo(sizes[0], sizes[-1])}")

        logger.info("\n".join(pinfos))

    return [reactions, species]


def _update_prdc_list(prdc: list, prdc_run: list, sizes: list, ntry: int) -> int:
    """If retry, add back threshold-based reduction to the list. Return the number of retry."""

    ntry -= 1  # Decrease the number of retry
    if ntry > 0:
        sdiff = size_difference(sizes[-2], sizes[-1])
        if sum(sdiff) == 0:  # No more retry if no size changes
            logger.info("No more reduction w/ %s. Size before & now: %s & %s.", prdc, sizes[-2], sizes[-1])
            ntry = 0
        else:
            logger.info("Retry %s. ntry: %d", prdc, ntry)
            logger.info("Size before & now (diff): %s & %s (%s)", sizes[-2], sizes[-1], sdiff)
            prdc_run.append(prdc)

    return ntry


def _update_lumping_within_tbr(species: list, popt: TbrOption, gnl: GlobalSetting) -> None:
    """Update lumping settings for threshold-based reduction."""

    if not popt.lumping:  # No lumping
        return

    if not os.path.isfile(popt.lumping):
        raise FileNotFoundError(f"Lumping file {popt.lumping} not found.")

    srdc_dict = {s.name: s for s in species if s.status == 1}  # Active species only
    for s in gnl.primary_vocs:
        srdc_dict.pop(s, None)  # Remove primary VOCs from the dict

    # Read raw groups
    logger.info("Reading group info for lumping from %s", popt.lumping)
    with open(popt.lumping, "r", encoding="utf-8") as f:
        groups_in = json.load(f)

    # Process
    if isinstance(groups_in, list):
        lgroups = _update_lumping_groups_from_list(groups_in, srdc_dict)
    elif isinstance(groups_in, dict):
        lgroups = _update_lumping_groups_from_dict(groups_in, srdc_dict)
    else:
        raise ValueError(f"Invalid lumping groups format in {popt.lumping}. Expected list or dict.")

    # Update groups in lumping
    popt.target_species = lgroups["sps"]  # Target species for reduction
    popt.given_groups = lgroups  # All groups with species indices
    if not (popt.target_species and popt.given_groups):
        raise ValueError("No target species or lumping groups found for threshold-based reduction.")

    # Lumping options
    lp_opt = rdc_checks["lp"]
    if popt.nsmax and popt.nsmax > 0:  # Limit on number of species in a group
        lp_opt["nsmax"] = popt.nsmax

    if popt.ngmax and popt.ngmax > 0:  # Limit on number of groups
        lp_opt["ngmax"] = popt.ngmax
    logger.info("Lumping in threshold-based reduction with options: %s", lp_opt)


def _update_lumping_groups_from_list(groups_in: list, rdc_dict: dict) -> dict:
    """Update lumping groups from a list of species groups."""

    if not groups_in:
        return {}

    sid_dict, target_sps, target_grps = {}, [], []
    for grp_in in groups_in:
        grp = [s for s in grp_in if s in rdc_dict]
        if not grp:  # Skip empty groups
            continue
        ind = len(target_grps)
        for s in grp:
            if s in sid_dict:
                raise ValueError(f"Species {s} appears in multiple groups.")
            sid_dict[s] = ind  # Record species index in the group
        # Update
        target_grps.append(grp)
        target_sps.extend(grp)

    logger.info("Lumping groups read for %d groups & %d species from lists.", len(target_grps), len(target_sps))
    sid_dict["grps"] = target_grps  # All groups
    sid_dict["sps"] = target_sps  # All target species
    return sid_dict


def _update_lumping_groups_from_dict(groups_in: dict, rdc_dict: dict) -> dict:
    """
    Update lumping groups from a dictionary of species groups
    in the format {targeted_species: [lumpable_species]}
    """

    if not groups_in:
        return {}

    target_sps = []
    grp_dict = {}  # Species index in the group

    for k, grp in groups_in.items():
        if k not in rdc_dict:  # No target species
            continue
        grp = [s for s in grp if s in rdc_dict]
        if not grp:  # No lumpable species
            continue
        grp_dict[k] = grp
        target_sps.append(k)

    logger.info("Lumping groups read for %d targeted species from dict.", len(target_sps))
    grp_dict["sps"] = target_sps  # All target species
    return grp_dict


def _save_tbr_mechanism(reactions: list, species: list, savpath: str, savname: str) -> None:
    """Save threshold-reduced mechanism to the specified path."""
    logger.info("Saving threshold-reduced mechanism %s to %s", savname, savpath)
    mech_output(savpath, savname, reactions, species)
    logger.info("Saved threshold-reduced mechanism %s to %s successfully.", savname, savpath)
