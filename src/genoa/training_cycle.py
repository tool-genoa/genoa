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
This module contains functions to run a reduction cycle in training.
"""

import os
import json
import shutil
import time
from multiprocessing import Pool
from typing import Optional

import numpy as np

from .constants import SNAME, LINE_SEP2, TRN_MODE
from .folder_path import move_to_path
from .jumping import change_scheme_via_jumping
from .logger import setup_logger, mute_all_loggers, unmute_all_loggers
from .lumping import update_scheme_via_merging
from .mechanism_pack import Mechanism
from .reduction_print import get_reduction_pinfos, get_simple_reduction_pline
from .reduction_setting import rdc_checks
from .reduction_strategy import find_reduction_via_strategy
from .setting_init import TrainingOption
from .simulation_run import run_simulation
from .testing_run import run_testing
from .training_init import CycleOption, ReductionItem, setup_cycle
from .training_order import species_list_to_group, read_group_from_file


# Logger
logger = setup_logger(__name__)


def run_reduction_cycle(settings_map: dict, iopt: Optional[CycleOption]) -> CycleOption:
    """Run reduction cycle"""

    # Get settings
    trn_opt = settings_map[SNAME.TRN]()

    # Initialize or update cycle options
    iopt = setup_cycle(trn_opt, iopt)

    # Prepare reduction candidates in order
    load_reduction_candidates(trn_opt, iopt)

    # Loop over reduction candidate groups
    while iopt.candidates and not iopt.is_stop:
        if fetch_targets(trn_opt, iopt):  # Get targeted sps & rcn
            run_training_reduction(trn_opt, iopt)  # Run reduction on targeted sps & rcn

    # Save reduced mechanism
    save_mechanism(trn_opt, iopt)

    # Pre-testing
    run_pre_testing(trn_opt, iopt)

    # Process the end of reduction cycle
    cycle_end_processing(trn_opt, iopt)

    return iopt


def save_mechanism(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Save the final reduced mechanism if any."""
    if not (iopt.run and iopt.run.nval):
        logger.info("No valid reduction to save for # %d reduction cycle.", iopt.icycle)
        return

    # Save the final reduced mechanism
    logger.info("Save the final reduced mechanism for # %d reduction cycle as %s ...", iopt.icycle, iopt.iname)
    iopt.mechv.name, iopt.mechv.path = iopt.iname, trn_opt.path_mech
    iopt.mechv.output()


def run_pre_testing(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Run pre-testing for the current reduction cycle"""

    if iopt.is_stop or not trn_opt.pretesting:
        logger.info("Skip pre-testing due to stop or no pre-testing options.")
        return
    if iopt.run.nval == 0:
        if iopt.icycle > 0:
            logger.info("Skip pre-testing for # %d reduction cycle due to no valid reductions.", iopt.icycle)
            return
        mech_name = iopt.mechv.name
    else:
        mech_name = iopt.iname

    t0 = time.perf_counter()
    logger.info("%sPretesting on %s for cycle # %d (%d val)...", LINE_SEP2, mech_name, iopt.icycle, iopt.run.nval)
    pinfos = [f"Pre-testing for reduction cycle # {iopt.icycle} with mechanism {iopt.mechv.name}."]

    # Update pre-testing options
    trn_opt.pretesting.mech_name = mech_name
    trn_opt.pretesting.sml = None

    tst_opt = run_testing(trn_opt.pretesting)  # Run pre-testing
    pinfos.append(f"  Pre-testing results:\n  Max err: {tst_opt.sml.err_arr.emax} at {tst_opt.sml.err_arr.emax_locs}")
    pinfos.append(f"  Ave err: {tst_opt.sml.err_arr.eave}\n  Finished in {time.perf_counter() - t0:.2f} s.")
    iopt.logt.write("\n".join(pinfos))

    # Check reduction cycle termination due to pre-testing results
    if tst_opt.sml.err_arr.emax > trn_opt.stop_at_emax or tst_opt.sml.err_arr.eave > trn_opt.stop_at_eave:
        iopt.stop(
            "Stop reduction cycle due to pre-testing results exceed the limit:\n"
            + f"  Max/Ave Err (Limit): {tst_opt.sml.err_arr.emax} ({trn_opt.stop_at_emax})/"
            + f"{tst_opt.sml.err_arr.eave} ({trn_opt.stop_at_eave})."
        )


def cycle_end_processing(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Post-processing for the end of reduction cycle"""

    if not iopt.run:
        return

    # Record the end of reduction cycle
    iopt.logt.write(iopt.get_pinfo_end())
    iopt.logt.flush()
    logger.info("End of reduction cycle # %d with %d valid reductions.", iopt.icycle, iopt.run.nval)

    # Reset training options
    if iopt.run.nval > 0:
        trn_opt.group_order_file = None  # Reset group order file
    if trn_opt.restart_from:
        trn_opt.restart_from = 0  # Reset restart input

    # Save for restart & reset cycle options
    save_for_restart(trn_opt, iopt)
    iopt.candidates = None  # Reset candidates


def save_for_restart(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Save the current cycle options for restart."""
    if not iopt.check_time(trn_opt):
        return

    if trn_opt.tlim > 0:  # No need restart
        data = {}
        logger.info("Saving empty restart file to %s ...", trn_opt.restart_to)
    else:
        logger.info("Saving current cycle options for restart to %s ...", trn_opt.restart_to)
        data = iopt.save_to_dict()
        # Get previous valid mechanism
        data["mechv"] = [trn_opt.path_mech, iopt.mechv.name]
        logger.info("Save previous valid mechanism: %s", iopt.mechv.name)
        # Add trn_opt settings
        if trn_opt.group_order_file:
            data["group_order_file"] = trn_opt.group_order_file
        # Add nval count for restart
        data["nval"] = iopt.run.nval
        data["ntry"] = iopt.run.ntry

    # Save to file
    with open(trn_opt.restart_to, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
    trn_opt.restart_to = None  # Reset restart file path


def fetch_targets(trn_opt: TrainingOption, iopt: CycleOption) -> bool:
    """Get targeted species and reactions for reduction from the list."""

    if iopt.is_stop or iopt.check_time(trn_opt):
        return False

    icans, irun = iopt.candidates, iopt.run
    if not icans:
        iopt.stop("No more candidates for reduction found in the mechanism.")
        return False

    logger.info("=> %d cands 2/ %d sps (%d valid / # %d tried rdcs)", len(icans), len(icans[-1]), irun.nval, irun.ntry)

    # Initialize
    targeted_sps, targeted_rcn = [], []
    s_info, p_info, c_info = iopt.infos.sps_info, iopt.infos.parent_info, iopt.infos.children_info
    # Get targeted species set
    for s in icans[-1]:
        if s not in s_info:
            continue
        if s in irun.checked_sps:
            continue
        targeted_sps.append(s)

    if not targeted_sps:  # No targeted species
        logger.info("No targeted species found. Searched: %s", icans[-1])
        icans.pop(-1)
        return False
    targeted_sps = set(targeted_sps)

    # Get related reaction set
    rcns = set()
    for s in targeted_sps:
        if s in p_info:
            rcns.update(p_info[s]["ircns"])
        if s in c_info:
            rcns.update(c_info[s]["ircns"])
    for ircn in rcns:
        rcn = iopt.mechv.reactions[ircn]
        if not rcn.status == 1:
            continue
        if ircn in irun.checked_rcn:
            continue
        for s in trn_opt.frozen_species:  # Skip w/ frozen species
            if s in rcn.reactants or s in rcn.products:
                continue

        # Check relationship with targeted species
        if set(rcn.products) & targeted_sps:
            targeted_rcn.append(ircn)
        elif set(rcn.reactants) & targeted_sps:  # If no product, include it in current run
            if not [s for s in rcn.products if s in s_info]:
                targeted_rcn.append(ircn)
        else:
            pinfos = [f"Reaction {ircn} does not have targeted species: {rcn.to_rcn()}"]
            pinfos.append(f" Targeted species: {targeted_sps}")
            for s in targeted_sps:
                if s in p_info:
                    if ircn in p_info[s]["ircns"]:
                        pinfos.append(f"  Parent of {s}: {p_info[s]['ircns']}")
                if s in c_info:
                    if ircn in c_info[s]["ircns"]:
                        pinfos.append(f"  Children of {s}: {c_info[s]['ircns']}")
            logger.error("\n".join(pinfos))
            raise ValueError(f"Reaction {ircn} does not have targeted species: {rcn.to_rcn()}")

    if not targeted_rcn:
        logger.info("No targeted reactions found. Searched: %s", rcns)
        icans.pop(-1)
        return False

    # Update reduction options
    iopt.infos.targeted_sps = targeted_sps
    iopt.infos.targeted_rcn = targeted_rcn

    return True


def load_reduction_candidates(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Group and order items for searching for reductions"""

    if iopt.is_stop:
        return

    if iopt.candidates:
        logger.info("Reduction candidates are already loaded.")
        return

    # Get candidate list
    if trn_opt.group_order_file:  # Read from file
        iopt.candidates = read_group_from_file(trn_opt.group_order_file, iopt.mechv.species, trn_opt.kept_all)
    elif trn_opt.group_order_mode in TRN_MODE:  # Generate order based on mode
        order_file = os.path.join(trn_opt.path_mech, f"training_order_{trn_opt.group_order_mode}-{iopt.icycle}.txt")
        imode = TRN_MODE[trn_opt.group_order_mode]
        if imode == "gprop":  # Based on species properties
            iopt.candidates = species_list_to_group(iopt.mechv.species, iopt.infos, order_file)
        elif imode == "gdel":  # Based on species deletion
            iopt.candidates = rank_species_with_deletion(trn_opt, iopt, order_file)
            # Updtae nsmax in rdc_checks
            if rdc_checks["lp"]["nsmax"] != 2:
                logger.info("Update nsmax for lp from %s to 2 for gen_delete order mode.", rdc_checks["lp"]["nsmax"])
                rdc_checks["lp"]["nsmax"] = 2
        else:
            raise ValueError(f"Invalid training order mode: {imode}. Expected key from {TRN_MODE}.")
        logger.info("Generated # %d candidates => %s.", len(iopt.candidates), order_file)
        trn_opt.group_order_file = order_file  # Save for next time
    else:
        raise ValueError("No group order file or mode is set.")

    # Check
    if not iopt.candidates:
        iopt.stop("No species for reduction found in the mechanism.")
    else:
        nsps = sum(len(i) for i in iopt.candidates)
        logger.info("Reducible elements # %d are in # %d groups.", nsps, len(iopt.candidates))


def run_training_reduction(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Find and evaluate reduction candidates in training."""

    if find_reductions(iopt):

        if write_reductions(trn_opt, iopt):

            evaluate_reductions_and_find_the_best(trn_opt, iopt)


def rank_species_with_deletion(trn_opt: TrainingOption, iopt: CycleOption, savfile: Optional[str] = None) -> list:
    """
    Find, write down and evaluate deleting species one by one. Obtained reduction errors can be used to rank species.
    Based on {species: redution error obtained when deleting that species} to build species group list
    """

    # Find species to delete
    sps_to_test = [s.name for s in iopt.mechv.species if s.status == 1 and s.name not in trn_opt.kept_all]
    if not sps_to_test:  # No deletion
        logger.info("No species deletion candidates found in the mechanism.")
        return []
    logger.info("Preparing # %d species deletion out of # %d species.", len(sps_to_test), len(iopt.mechv.species))
    iopt.run.rdcs = [ReductionItem(irdc=i, stgy="rs", rdc=[s]) for i, s in enumerate(sps_to_test)]

    # Write species deletion candidates
    if not write_reductions(trn_opt, iopt):
        logger.info("No valid species deletion candidates found.")
        return []

    # Evaluate species deletion candidates
    logger.info("Evaluating species deletion candidates ...")
    sps_del_errs = evaluate_all_reductions_and_return_errors(iopt)

    # Rank species by minimum errors
    logger.info("Ranking # %d species based on deletion errors.", len(sps_del_errs))
    sps_ranked = sorted(sps_del_errs, key=sps_del_errs.get)

    # Move species ends with "O" to the front of the list if needed
    if trn_opt.group_order_mode.endswith("O"):
        sps0 = [s for s in sps_ranked if s.endswith("O")]
        if sps0:
            logger.info("Moving # %d species ending with 'O' to the front of the list.", len(sps0))
            sps_ranked = sps0 + [s for s in sps_ranked if not s.endswith("O")]
        else:
            logger.info("No species ending with 'O' found in the list. No reordering applied.")

    # Save ranked species to file if needed
    if savfile:
        with open(savfile, "w", encoding="utf-8") as f:
            for i, s in enumerate(sps_ranked):
                f.write(f"No.{i} {sps_del_errs[s]}: {s}\n")
        logger.info("Species deletion ranking is saved to %s.", savfile)

    return [[s] for s in sps_ranked]  # Return as a list of lists for compatibility


def find_reductions(iopt: CycleOption) -> bool:
    """Find reduction candidates for each applied reduction strategy."""

    rdcs_all = []
    mechv, nmcomb, nmrdc, infos = iopt.mechv, iopt.stage.nmax_eval_comb, iopt.stage.nmax_eval_rdc, iopt.infos
    for isgy in iopt.stage.strategies:  # Find per strategy
        rdcs = find_reduction_via_strategy(mechv.reactions, mechv.species, isgy, nmcomb, nmrdc, infos)
        if not rdcs:  # No reduction found
            continue
        # logger.info("Found # %d reduction candidates for strategy %s.", len(rdcs), isgy)

        # For each reduction possibility, output a reduced mechanism
        for rdc in rdcs:
            rdc_item = ReductionItem(irdc=len(rdcs_all), stgy=isgy, rdc=rdc)
            rdcs_all.append(rdc_item)

    iopt.run.rdcs = rdcs_all  # Update

    # No reduction found, go to next
    if not iopt.run.rdcs:
        iopt.go_next("No candidate is found.")
        return False

    # Find reduction posibility
    iopt.run.ntry += 1  # Count
    # logger.info("Found # %d rdcs in Cycle # %d - Try # %d", len(rdcs_all), iopt.icycle, iopt.run.ntry)
    return True


def write_reductions(trn_opt: TrainingOption, iopt: CycleOption) -> bool:
    """Write all reduction candidates to files."""

    # Update and write mechanisms
    mute_all_loggers()  # Mute all loggers
    if len(iopt.run.rdcs) > 1 and trn_opt.npara > 1:  # Parallel

        logger.info("Write # %d mechanisms in parallel with # %d processes ...", len(iopt.run.rdcs), trn_opt.npara)
        iopt.logt.mute()  # Mute log for parallel writing
        pool_inputs = [(rdc, iopt) for rdc in iopt.run.rdcs]
        with Pool(trn_opt.npara) as pool:
            rdc_mechs = pool.starmap(write_reduction_candidate, pool_inputs)
        iopt.logt.unmute()  # Unmute log after parallel

    else:  # Serial
        iopt.logt.fout.flush()
        logger.info("Write # %d mechanisms in serial ...", len(iopt.run.rdcs))
        rdc_mechs = [write_reduction_candidate(rdc, iopt) for rdc in iopt.run.rdcs]
    unmute_all_loggers()  # Unmute loggers

    # Assess mechanisms
    nval = 0
    for rdc, mech in zip(iopt.run.rdcs, rdc_mechs):
        if mech and mech.minfo.is_valid:
            # Compute reduction score
            rdc.score = compute_reduction_score(mech.get_size(), iopt.mechv.get_size())
            nval += 1
        else:  # Record to checked dict to avoid searching again
            iopt.run.to_checked(rdc)
        rdc.mech = mech

    if not nval:
        iopt.go_next("No valid reduction is written.")
        return False

    # logger.info("Wrote # %d valid reduction candidates.", nval)
    return nval


def evaluate_reductions_and_find_the_best(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """
    Evaluate reduction candidates with reduction scores and error tolerances.
    The evaluation will stop if the best reduction is found.
    """

    # Set error check thresholds
    err_max = np.maximum(iopt.stage.err_max, iopt.mechv.errs[0])
    err_ave = np.maximum(iopt.stage.err_ave, iopt.mechv.errs[1])
    iopt.sml.err_checks = [
        [err_max, iopt.stage.delta_err, iopt.mechv.errs[0]],  # Maximum error
        [err_ave, iopt.stage.delta_err, iopt.mechv.errs[1]],  # Average error
    ]
    # logger.info("Evaluate w/ error tolerances: %s", iopt.sml.err_checks)

    # Record valid mechanisms
    scores, mech_names, rdcs_run, rinfos = [], [], {}, []
    for rdc in iopt.run.rdcs:
        iname = rdc.mech.name
        if rdc.mech.minfo.is_valid:
            scores.append(rdc.score)
            mech_names.append(iname)
            rdcs_run[iname] = rdc
            if rdc.score < 0:
                logger.warning("Rdc %s w/ %s & score %s < 0!", iname, rdc.stgy, rdc.score)
        else:
            logger.info("Invalid rdc %s w/ strategy %s", iname, rdc.stgy)

    iopt.sml.scores = scores  # Update scores
    iopt.sml.update_mech_names(mech_names)  # Update mechanism names

    # Run simulations for valid reductions
    logger.info("Running simulations for # %d valid reductions ...", len(iopt.sml.chemids))
    iopt.sml = run_simulation(iopt.sml)

    # Evaluate simulation results
    nval = 0  # Count valid reductions
    for i, s in enumerate(iopt.sml.chemids):
        rdc = rdcs_run[s]
        ierrs = iopt.sml.err_arr.arrs[i]
        max_err = np.max(ierrs, axis=(0, 2))

        # Check and compute errors
        if np.any(ierrs == -1 | np.isnan(ierrs)):
            rdc.mech.minfo.is_valid = False  # Invalid mechanism
            iopt.run.to_checked(rdc)  # Record to checked dict
            if max(max_err) == -1:
                continue
            masked_errs = np.ma.masked_equal(ierrs, -1)
            ipinfo = f"{max_err}*\t{np.ma.mean(masked_errs, axis=(0, 2)).filled(-1)}*"
        else:
            nval += 1
            rdc.errs = [max_err, np.mean(ierrs, axis=(0, 2))]
            ipinfo = f"{rdc.errs[0]}\t{rdc.errs[1]}"
        rinfos.append(
            f" {rdc.mech.name}\t{rdc.stgy}\t{rdc.mech.get_size()}\t{rdc.score}\t{ipinfo}\t"
            + get_simple_reduction_pline(rdc.stgy, rdc.rdc)
        )

    # Write to file
    winfos = [
        f"{LINE_SEP2}Cycle # {iopt.icycle} - Try # {iopt.run.ntry} - Pre Valid # {iopt.run.nval}"
        + f" - Size of previous mechanism: {iopt.mechv.get_size()}\n"
        + f"No.rdcs after searching/ writing/ running: {len(iopt.run.rdcs)} / {len(iopt.sml.chemids)} / {nval}"
    ]

    if rinfos:
        winfos.append("  Rdc ID\tStgy\tSize\tScore\tMax Errs\tAve Errs (* for partial smls)\t# Rdc: Details")

    # Reset simulation options
    iopt.sml.reset_for_training()

    # Get the best reduction & update settings
    update_with_brdc(get_best_reduction(iopt), iopt, trn_opt, winfos + rinfos)


def evaluate_all_reductions_and_return_errors(iopt: CycleOption) -> dict:
    """
    Evaluate all reduction candidates without checking reduction scores or error tolerances.
    Return a dictionary in the format of {species: reduction error}
    """

    # Update valid mechanism names
    mechs = {rdc.mech.name: rdc.rdc[0] for rdc in iopt.run.rdcs if rdc.mech.minfo.is_valid}
    iopt.sml.update_mech_names(list(mechs.keys()))

    # Run all simulations w/o lock
    iopt.sml.wlock = False
    iopt.sml = run_simulation(iopt.sml)

    # Record reduction errors
    errs = {mechs[s]: np.mean(iopt.sml.err_arr.arrs[i]) for i, s in enumerate(iopt.sml.chemids)}

    # Reset simulation options
    iopt.sml.reset_for_training()

    return errs


def update_with_brdc(brdc: Optional[ReductionItem], iopt: CycleOption, trn_opt: TrainingOption, pinfos: list) -> None:
    """Update cycle options with the selected best reduction candidate."""

    # Record targeted species
    pinfos.append(f"\nTarget # {len(iopt.infos.targeted_sps)} species: {iopt.infos.targeted_sps}")

    if brdc is None:  # No find valid reduction
        iopt.go_next()
        pinfos.append("No valid reduction found.")
    else:  # Find valid reduction
        pinfos.append(f"Accept: {brdc.irdc} with score {brdc.score:.2f} via {brdc.stgy}")
        pinfos.append(after_found(trn_opt, iopt, brdc))

    # Record
    iopt.logt.write("\n".join(pinfos))

    # Reset
    iopt.run.rdcs = None


def get_best_reduction(iopt: CycleOption) -> Optional[ReductionItem]:
    """Get the best reduction candidate that meets the error tolerances and with highest score."""

    brdc, bave, bmax = None, None, None  # Best reduction

    for rdc in iopt.run.rdcs:
        # Not qualified
        if not (rdc.mech.minfo.is_valid and rdc.errs):
            continue
        if not meet_error_tolerances(rdc, iopt):
            continue
        if not meet_aerosol_treatment(rdc, iopt):
            continue
        # Get ave & max errors
        imax, iave = np.max(rdc.errs[0]), np.mean(rdc.errs[1])
        # Get the reduction w/ highest score & lower errors
        if brdc is None or rdc.score > brdc.score:  # With higher score
            brdc = rdc
            bave, bmax = iave, imax
        elif rdc.score == brdc.score:  # With the same score, choose the one with lower errors
            if iave < bave or (iave == bave and imax < bmax):  # W/ lower ave & max error
                brdc = rdc
                bave, bmax = iave, imax

    return brdc


def after_found(trn_opt: TrainingOption, iopt: CycleOption, brdc: ReductionItem) -> str:
    """
    After finding the best reduction candidate, backup the valid reduction mechanism and update settings if need.
    """

    # Update counters
    iopt.run.nval += 1
    iopt.run.npres[brdc.stgy] += 1

    # Update previous mechanism
    iopt.mechv = brdc.mech
    iopt.mechv.errs = brdc.errs

    # Backup reaction/species files
    newf = os.path.join(iopt.pth_mech, f"{iopt.run.nval}")
    for i in brdc.mech.get_files():
        shutil.move(i, f"{newf}.{i.rsplit('.', 1)[1]}")
    # logger.info("Backup mechanism %s as %s.", brdc.mech.name, newf)

    # Backup results if needed
    if iopt.err_files_str and "pre" in iopt.err_files_str:
        move_to_path(os.path.join(iopt.paths.res, f"{brdc.irdc}"), iopt.paths.ref, "pre", True)

    # Record reduction info
    pinfos = get_reduction_pinfos(brdc.stgy, brdc.rdc, iopt.mechv.reactions, iopt.mechv.species, iopt.infos)

    # Update for searching next
    iopt.run.reset_checked()
    iopt.update_infos(trn_opt)

    return "    " + "\n".join(pinfos)


def meet_error_tolerances(ritem: ReductionItem, iopt: CycleOption) -> bool:
    """Check if errors within the tolerances."""

    err_now = ritem.errs  # Current errors
    err_pre = iopt.mechv.errs  # Previous errors

    # Maximum error
    err_max = np.maximum(iopt.stage.err_max, err_pre[0])
    for err, tol in zip(err_now[0], err_max):
        if err > tol:
            return False
    # Average error
    err_ave = np.maximum(iopt.stage.err_ave, err_pre[1])
    for err, tol in zip(err_now[1], err_ave):
        if err > tol:
            return False
    # Skip error increase check if current error is small

    # Error increase
    for n in range(2):
        for j, (err, err0) in enumerate(zip(ritem.errs[n], iopt.mechv.errs[n])):
            if err - err0 > iopt.stage.delta_err[j]:
                return False

    return True


def meet_aerosol_treatment(ritem: ReductionItem, iopt: CycleOption) -> bool:
    """If aerosol-oriented treatment is applied, only accept reductions reducing aerosol size or error decrease."""

    if iopt.stage.flg_aerosol == 0:  # No aerosol treatment
        return True

    # Check aerosol size change
    if iopt.mechv.minfo.naero - ritem.mech.minfo.naero > 0:
        return True

    # Check error change
    return with_error_decrease(ritem.errs, iopt.mechv.errs)


def with_error_decrease(errs: list, errs0: list) -> bool:
    """Return True if all errors are less than or equal to the previous ones."""

    for err, err0 in zip(errs, errs0):
        if np.any(np.array(err) > np.array(err0)):
            return False

    return True


def compute_reduction_score(cur_size, org_size) -> float:
    """Compute reudction score based on the number differences in reaction/species/aerosols."""

    dsize = [s - p for s, p in zip(org_size, cur_size)]
    if min(dsize) < 0:  # Check if size reduced
        raise ValueError("Size increased after reduction. Get: ", cur_size, org_size, dsize)
    if sum(dsize) == 0:  # No change
        return 0.0
    dsize[-1] *= 10  # Weight for aerosol size

    return sum(dsize)


def write_reduction_candidate(ritem: ReductionItem, iopt: CycleOption) -> Mechanism:
    """Write reduced mechanism to file."""

    rdc_func_dict = {
        "rs": change_mech_via_removal,
        "rm": change_mech_via_removal,
        "da": change_mech_via_removal,
        "lp": change_mech_via_merging,
        "jp": change_mech_via_jumping,
        "rp": change_mech_via_merging,
    }

    # Check if reduction strategy is valid
    if ritem.stgy not in rdc_func_dict:
        raise ValueError(f"Invalid reduction strategy: {ritem.stgy}")
    mech_lists = rdc_func_dict[ritem.stgy](ritem, iopt)
    return output_training_mechanism(ritem, mech_lists, iopt)


def output_training_mechanism(ritem: ReductionItem, mech_list: list, iopt: CycleOption) -> Mechanism:
    """Save temperoary mechanism with irdc index for training."""

    # Save redcued mechanism
    imech = Mechanism(name=str(ritem.irdc), path=iopt.paths.mech, reactions=mech_list[0], species=mech_list[1])
    imech.mopt = iopt.mechv.mopt
    imech.trim()  # Trim mechanism
    imech.output()  # Output mechanism if valid
    if not imech.minfo.is_valid:  # Clean memory if not output
        imech.reactions, imech.species = None, None

    return imech


def _copy_reactions_and_species(iopt: CycleOption) -> tuple:
    """Copy reactions and species from the current mechanism."""
    return ([r.copy() for r in iopt.mechv.reactions], [s.copy() for s in iopt.mechv.species])


def change_mech_via_removal(ritem: ReductionItem, iopt: CycleOption) -> list:
    """Change mechanism via removal of species/reactions."""

    reactions, species = _copy_reactions_and_species(iopt)
    if ritem.stgy == "rm":  # Remove reactions
        for r in ritem.rdc:
            reactions[r].status = 0
    else:
        sinfo = iopt.infos.sps_info
        if ritem.stgy == "rs":  # Remove species
            for s in ritem.rdc:
                species[sinfo[s]].status = 0

        elif ritem.stgy == "da":  # Remove gas-particle partitioning
            for s in ritem.rdc:
                isp = sinfo[s]
                species[isp].condensable = False
                species[isp].update_for_reduction()
        else:
            raise ValueError(f"Invalid reduction strategy for removal: {ritem.stgy}")

    return [reactions, species]


def change_mech_via_merging(ritem: ReductionItem, iopt: CycleOption) -> list:
    """Change mechanism via lumping or replacement."""

    reactions, species = _copy_reactions_and_species(iopt)
    update_scheme_via_merging(reactions, species, ritem.rdc, iopt.infos.sps_info, ritem.stgy)

    return [reactions, species]


def change_mech_via_jumping(ritem: ReductionItem, iopt: CycleOption) -> list:
    """Change mechanism via jumping."""

    reactions, species = _copy_reactions_and_species(iopt)
    for s in ritem.rdc:
        change_scheme_via_jumping(reactions, species, s, iopt.infos.sps_info)
    return [reactions, species]
