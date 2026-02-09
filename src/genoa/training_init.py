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
This module contains functions to initialize the training options and setup the reduction cycle options.
"""

import os
import json
import shutil
import time
from typing import Optional, Dict

import numpy as np
from attrs import define, field, asdict

from .constants import SNAME, LINE_SEP, LINE_SEP2
from .folder_path import FolderPath, update_output_path, copy_to_path, link_to_path
from .gecko_cst import REF_FILES_GCK
from .logger import setup_logger
from .mechanism_pack import Mechanism
from .record import RecordOption, build_log, record_initialization
from .reduction_basic import ReductionListInfo
from .setting_init import TrainingOption, TestingOption
from .setting_update import setup_4action
from .simulation_init import RunSmlSetting
from .simulation_refconc import get_refconcs_from_files
from .simulation_run import run_simulation
from .ssh_cst import REF_FILES_SSH
from .training_stage import StageOption, setup_stage, read_training_param_dict
from .utils import restart_filename


# Logger
logger = setup_logger(__name__)


@define(slots=True)
class ReductionItem:
    """Class to store information related to a reduction."""

    irdc: int = -1  # Reduction index
    stgy: str = ""  # Reduction strategy
    rdc: Optional[list] = None  # Reduction list
    mech: Optional[Mechanism] = None  # Reduced mechanism object
    score: Optional[float] = None  # Reduction score
    errs: Optional[list] = None  # Errors


@define(slots=True)
class CycleRun:
    """Options to record current reduction in the cycle."""

    # Counter
    ntry: int = 0  # Number of reduction candidates tried
    nval: int = 0  # Number of valid reductions found
    npres: Dict[str, int] = field(factory=dict)  # Counters for strategies

    checked_sps: set = field(factory=set)  # Checked species for current candidates
    checked_rcn: set = field(factory=set)  # Checked reactions for current candidates
    checked: dict = field(factory=dict)  # Checked reduction elements for current strategies

    # Search
    t0: float = 0  # Start time
    rdcs: Optional[ReductionItem] = None  # Reduction items for current candidates

    def init_npres(self, strategies: list) -> None:
        """Initialize npres with strategies."""
        self.npres = {s: 0 for s in strategies}

    def reset_checked(self) -> None:
        """Reset checked species and reactions."""
        self.checked_sps.clear()
        self.checked_rcn.clear()
        self.checked.clear()

    def to_checked(self, rdc: ReductionItem) -> None:
        """Record to checked dictionary."""
        if rdc.stgy not in self.checked:
            self.checked[rdc.stgy] = []
        self.checked[rdc.stgy].append(rdc.rdc)


@define(slots=True)
class CycleOption:
    """Options can be updated for each reduction cycle."""

    icycle: int = 0  # Current reduction cycle index
    iname: str = ""  # Prefix for mechanisms saved in the cycle

    # Reduction stage
    istage: int = 0  # Current reduction stage index
    stage: Optional[StageOption] = None

    # Mechanisms
    mechv: Optional[Mechanism] = None  # Previous valid mechanism

    refs: Optional[list] = None  # Reference mechanisms
    err_files_str: str = ""  # Reference error files string
    sml: Optional[RunSmlSetting] = None  # Settings for running simulations

    # Paths
    paths: Optional[FolderPath] = None
    pth_mech: str = ""

    logt: Optional[RecordOption] = None  # Track file

    # Termination
    is_stop: bool = False  # Flag for termination
    stop_infos: list = field(factory=list)  # Stop information

    # Update
    update_stage: bool = True  # If need to update parameters
    update_ref: Optional[int] = None  # If need to run reference simulations - updated in setup_cycle
    update_kinetic: Optional[int] = None  # If need to update kinetics - updated in setup_cycle
    update_mech: bool = False  # If need to evaluate difference between ref and pre

    # Reduction
    infos: Optional[ReductionListInfo] = None
    run: Optional[CycleRun] = None
    candidates: Optional[list] = None

    def init_wtraining(self, topt: TrainingOption) -> str:
        """Init cycle options by training options"""
        pinfos = [f"Initialize the first reduction cycle # {self.icycle} with training options.\n"]

        # Initialize logt file
        if self.logt is None:
            self.logt = topt.logt

        # Reload cycle option if needed
        pinfo = self.reload_from_restart(topt)
        if pinfo:
            pinfos.append(pinfo)

        # Initialize starting mechanism
        if self.mechv is None:
            self.mechv = topt.mechs[1]

        # Initialize paths
        if self.paths is None:
            if topt.paths:
                self.paths = topt.paths
            else:
                if not topt.path_work:
                    raise ValueError("No valid path for training work folder.")
                self.paths = FolderPath(work=topt.path_work)
                self.paths.init_in_training()
        return "\n".join(pinfos)

    def update_for_next(self) -> str:
        """Update for a new reduction cycle."""

        self.icycle += 1  # Update cycle index
        pinfos = ["Reduction settings:"]

        # Check if go to next
        if self.run.nval == 0 or self.run.nval <= self.stage.nmin_retry:
            pinfos.append(f"Go to next due to limit no. valid reductions {self.run.nval}")
            if self.stage.flg_aerosol == 2:
                pinfos.append("  - Retry current stage without aerosol-oriented treatment.")
                self.stage.flg_aerosol = 0  # Reset
            else:  # Go to next stage
                self.istage += 1
                self.update_stage = True
                pinfos.append(f"  - Go to next stage # {self.istage}.")
        else:  # Continue current stage
            pinfos.append(f"No change. (Previous valid run # {self.run.nval})")

        return "\n".join(pinfos)

    def update_infos(self, topt: TrainingOption) -> None:
        """Update reduction info for current reduction cycle"""

        if not self.mechv:
            raise ValueError("No valid mechanism for updating reduction info.")

        # Update reduction info
        rdc_settings = {"frozen": topt.kept_all, "refconcs": topt.ref_concs}
        self.infos = ReductionListInfo()
        self.infos.update(self.mechv.reactions, self.mechv.species, rdc_settings)
        self.infos.checked = self.run.checked

    def stop(self, msg: str = "") -> None:
        """Update stop flag and info."""
        self.is_stop = True
        if msg:
            self.stop_infos.append(msg)

        pinfos = [f"{LINE_SEP2}Stop reduction cycle # {self.icycle}."]
        if self.stop_infos:
            pinfos.extend(self.stop_infos)
        self.logt.write("\n".join(pinfos))

    def get_pinfo(self) -> str:
        """Get record info."""
        pinfos = [f"\nReduction tired # {self.run.ntry} times and found # {self.run.nval} valid reductions."]
        if self.mechv:
            pinfos.append(f"Previous valid mechanism: {self.mechv.name} with errs: {self.mechv.errs}.")
        return "\n".join(pinfos)

    def get_pinfo_end(self) -> str:
        """Get record info at the end of reduction cycle."""
        pinfos = [
            f"{LINE_SEP2}END OF CYCLE # {self.icycle}. Total run # {self.run.ntry}\tValid run # {self.run.nval}\t"
            + f"Cycle run time: {time.perf_counter() - self.run.t0:.2f} s"
        ]
        if self.run.nval:
            pinfos.append(f"  Accepted reductions by strategies: {self.run.npres}")
        if self.mechv:
            pinfos.append(f"  Final mechanism: {self.mechv.get_minfo()}")
            pinfos.append(f"    Current training errors: {self.mechv.errs}")
        else:
            pinfos.append("  No valid mechanism found.")

        return "\n".join(pinfos)

    def check_time(self, topt: TrainingOption) -> bool:
        """Check if need to stop due to time limit. Return True if need to record restart file."""
        if not topt.tlim:  # No time limit
            return False

        if not topt.restart_to:  # Update restart file
            topt.restart_to = restart_filename(topt.path_mech, is_new=True, only_num=False)

        if self.is_stop:  # Already stopped
            return True

        # Check time limit
        tnow = time.perf_counter() - topt.tstart
        if tnow < topt.tlim:  # Time limit not reached
            return False

        pinfo = f"Stop training due to time limit: {tnow:.1f} > {topt.tlim} s."
        topt.tlim = -abs(topt.tlim)  # Mark as checked with negative value
        logger.info(pinfo)
        self.stop(pinfo)
        return True

    def save_to_dict(self) -> dict:
        """Save cycle option to file fo reloading."""
        # Get attributes
        tosave = ["istage", "candidates", "icycle"]
        return asdict(self, filter=lambda a, v: a.name in tosave)

    def reload_from_restart(self, topt: TrainingOption) -> str:
        """Reload cycle option from the restart_from file if given."""
        if not topt.restart_from:
            return ""

        # Read data
        logger.info("Reload cycle options from restart file: %s", topt.restart_from)
        with open(topt.restart_from, "r", encoding="utf-8") as f:
            topt.restart_from = json.load(f)
        # If empty, return
        if not topt.restart_from:
            pinfo = f"Empty restart file: {topt.restart_from}. Stop training."
            logger.info(pinfo)
            self.stop(pinfo)
            return ""

        # Load data
        pinfos = []
        # Load mechanism
        mechv = topt.restart_from.pop("mechv", None)  # in ["mech_path", "mech_name"]
        if mechv:
            self.mechv = Mechanism(name=mechv[1], path=mechv[0], runid="pre", mopt=topt.trim_opts)
            pinfos.append(f"previous mechanism: {self.mechv.name} in {self.mechv.path}")

        # Load cycle number
        self.icycle = topt.restart_from.pop("icycle", 0)
        # Load nval for checking later
        nval = topt.restart_from.pop("nval", 0)
        ntry = topt.restart_from.pop("ntry", 0)
        logger.info("Reload: %s. Previous # %d valid out of # %d try.", pinfos, nval, ntry)

        # Load other cycle options
        for k, v in topt.restart_from.items():
            if hasattr(self, k):  # Check in copt
                logger.info("Reload: %s with value: %s to coption.", k, v)
                setattr(self, k, v)
            elif hasattr(topt, k):  # Check in topt
                logger.info("Reload: %s with value: %s to training option.", k, v)
                setattr(topt, k, v)
            else:
                raise AttributeError(f"Invalid attribute: {k} for reloading cycle option.")
            pinfos.append(k)

        if pinfos:
            pinfo = f"Reload cycle options from {topt.restart_from}. Got: {', '.join(pinfos)}"
        else:
            pinfo = f"No cycle options to reload from {topt.restart_from}. Got {topt.restart_from}"
            logger.warning(pinfo)

        # Set flag for restarting cycle
        topt.restart_from = [nval, ntry] if nval + ntry > 0 else 0
        return pinfo

    def go_next(self, msg: str = "") -> None:
        """Go to next reduction candidate."""
        logger.info("%s Move to next.", msg)
        self.run.checked_sps.update(self.infos.targeted_sps)
        self.run.checked_rcn.update(self.infos.targeted_rcn)
        self.candidates.pop(-1)


def setup_training(settings_map: dict, updates: Optional[dict]) -> None:
    """Initialize the options for all reduction cycles."""

    # Get settings
    gnl = settings_map[SNAME.SETUP]()
    trn_opt = settings_map[SNAME.TRN]()
    trn_opt.tstart = time.perf_counter()  # Start time

    # Update settings
    setup_4action(SNAME.TRN, trn_opt, gnl)

    # Initialize track file logt
    logf = os.path.join(trn_opt.path_mech, f"record_{trn_opt.runid}.txt")  # Default log file
    if trn_opt.restart_from and trn_opt.tlim:  # Update log file if restart
        fname, fext = os.path.splitext(logf)
        logf = f"{fname}_r{restart_filename(trn_opt.path_mech, is_new=False, only_num=True)}{fext}"
    logger.info("Set up training log file: %s", logf)
    trn_opt.logt = build_log(logf, f"Tracking file for training with id: {trn_opt.runid} and PID: {os.getpid()}\n")

    # Update training mechanism
    if updates:
        logger.info("Use updated mechanism to set up training: %s", updates)
        if "mech_name" not in updates or "mech_path" not in updates:
            raise ValueError("Invalid update settings for training mechanism.")
        trn_opt.mech_name = updates.pop("mech_name")
        trn_opt.mech_path = updates.pop("mech_path")

    # Read training parameters
    read_training_param_dict(trn_opt)

    # Initialize training paths
    trn_opt.paths = FolderPath(work=trn_opt.path_work)
    trn_opt.paths.init_in_training()

    # Prepare other training options
    _prepare_mechanisms(trn_opt)
    trn_opt.ref_concs = get_refconcs_from_files(gnl.refconc_file)

    # Record initialization
    record_initialization(settings_map, trn_opt, SNAME.TRN, trn_opt.logt)

    # Initialize pre-testing settings
    _init_pre_testing(trn_opt)


def clean_training(settings_map: dict) -> None:
    """Clean temporary files after training if needed."""

    trn_opt = settings_map[SNAME.TRN]()
    if trn_opt.clean_tmp:
        logger.info("Clean temporary files in training work folder: %s", trn_opt.paths.work)
        shutil.rmtree(trn_opt.paths.work, ignore_errors=True)  # Remove work folder
    else:
        logger.info("Keep temporary files in training work folder: %s", trn_opt.paths.work)


def _prepare_mechanisms(trn_opt: TrainingOption) -> None:
    """Copy reference and starting mechanism to the training folder."""

    # Get mechanisms
    imechv = Mechanism(name=trn_opt.mech_name, path=trn_opt.mech_path, runid="pre", mopt=trn_opt.trim_opts)
    imech_ref = Mechanism(name=trn_opt.ref_mech_name, path=trn_opt.ref_mech_path, runid="ref", mopt=trn_opt.trim_opts)

    # Check if changed during reading
    add_mopt = {"merge": trn_opt.with_merge, "split": False}
    logger.info("Trimming mechanism %s w/ additional options: %s ...", imechv.name, add_mopt)
    imechv.trim(add_mopt)  # w/ additional options
    if imechv.minfo.is_changed:  # Save trimmed mechanism
        logger.warning("Starting mechanism changed! Rename to %s.", imechv.name)
        imechv.save_new(f"{imechv.name}n")  # n: new
    else:
        logger.info("No change in starting mechanism %s.", imechv.name)

    # Save mechanisms
    if imechv.isdiff(imech_ref) or imechv.minfo.is_changed:  # Copy pre if diff from ref or if trimmed
        if imechv.name == imech_ref.name:
            imechv.name += "N"
            trn_opt.logt.write(f"Renamed starting mechanism: {imechv.name} due to same name with reference.")
        imechv.copy_and_link(trn_opt.path_mech, trn_opt.paths.mech)
    imech_ref.copy_and_link(trn_opt.path_mech, trn_opt.paths.mech)  # Copy and link ref
    trn_opt.mechs = [imech_ref, imechv]  # Reference and starting mechanisms


def _init_pre_testing(trn_opt: TrainingOption) -> None:
    """Initialize pre-testing settings."""

    if not (trn_opt.pretesting_path and os.path.exists(trn_opt.pretesting_path)):
        trn_opt.pretesting = None
        return  # No pre-testing

    runid = "PTST"
    if trn_opt.restart_from and trn_opt.tlim:
        runid += f"r{restart_filename(trn_opt.path_mech, is_new=False, only_num=True)}"

    pinfos = [
        f"Pre-testing {runid} with reference mechanism: {trn_opt.ref_mech_name}",
        f"Stop at max_err: {trn_opt.stop_at_emax}, ave_err: {trn_opt.stop_at_eave}",
        f"Pre-testing path: {trn_opt.pretesting_path}",
    ]
    trn_opt.logt.write("\n".join(pinfos))
    logger.info("%s", "; ".join(pinfos))

    # Get pre-testing options
    trn_opt.pretesting = TestingOption(
        runid=runid,
        path_cond=trn_opt.pretesting_path,
        mech_path=trn_opt.path_mech,
        ref_mech_name=trn_opt.ref_mech_name,
        ref_mech_path=trn_opt.path_mech,
        path_work=os.path.join(trn_opt.path_work, runid),
    )

    # Get log file for pre-testing
    trn_opt.pretesting.loge = build_log(
        os.path.join(trn_opt.path_mech, f"pretest_{runid}.txt"),
        f"Pre-testing with id: {runid} and PID: {os.getpid()}\n",
    )


def setup_cycle(trn_opt: TrainingOption, iopt: Optional[CycleOption]) -> CycleOption:
    """Initialize parameters and options for the current reduction cycle"""

    t0 = time.perf_counter()  # Start time

    if iopt is None:  # Initialize options
        iopt = CycleOption()
        pinfo = iopt.init_wtraining(trn_opt)
    else:  # Update for next cycle
        pinfo = iopt.update_for_next()
    if iopt.is_stop:  # If stopped, return
        return iopt

    # Initialize run options
    nval, ntry = trn_opt.restart_from if isinstance(trn_opt.restart_from, list) else (0, 0)
    iopt.run = CycleRun(t0=t0, nval=nval, ntry=ntry)

    # Saved mechanism name prefix
    iopt.iname = f"{trn_opt.runame}-{iopt.icycle}"
    linfo = "Start" if ntry + nval == 0 else "Continue"
    logger.info("%s%s reduction cycle # %d with mechanism prefix: %s ...", LINE_SEP, linfo, iopt.icycle, iopt.iname)
    iopt.logt.write(f"{LINE_SEP}Start reduction cycle # {iopt.icycle} with mechanism prefix: {iopt.iname}\n{pinfo}")

    # Update stage parameters if update_stage is True
    setup_stage(trn_opt, iopt)
    if iopt.is_stop:
        return iopt
    iopt.logt.write(iopt.stage.get_sginfo(iopt.istage))

    # Set reduction cycle run options w/ strategies
    iopt.run.init_npres(iopt.stage.strategies)

    # Update paths
    iopt.pth_mech = update_output_path(os.path.join(trn_opt.path_mech, f"{iopt.icycle}_sav"))

    # Update reference simulations and get sml settings
    _run_reference_simulations(trn_opt, iopt)

    # Update mechanism
    _setup_mechanism_for_cycle(trn_opt, iopt)

    # Get previous mechanism errors if needed
    _check_mech_diff(trn_opt, iopt)

    # Record
    if not iopt.is_stop:
        logger.info("Finished setting up reduction cycle # %d", iopt.icycle)

    return iopt


def _run_reference_simulations(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Run simulations with reference mechanisms (ref_cases) under training conditions (conds)"""

    if not iopt.update_ref:  # No run
        logger.info("No need to run reference simulations w/ update_ref: %d", iopt.update_ref)
        return

    # Update working folder
    iopt.paths.clean_in_training()

    # Update error evaluation files & reference mechanisms
    _update_ref_mechanisms(trn_opt, iopt)
    _update_mech_diff(iopt)

    # Prepare simulation settings
    rfiles_key = {"GECKO": REF_FILES_GCK, "SSH": REF_FILES_SSH}  # ref_files_str variable names
    settings = {
        "mech_names": [imech.name for imech in iopt.refs],
        "mech_path": iopt.paths.mech,
        "path_cond": iopt.stage.path_cond,
        "conds": iopt.stage.conds,
        "ref_files_str": "",  # iopt.err_files_str,
        "wid": True,
        "paths": iopt.paths,
        # "out_mode": 0,
    }
    settings["paths"].init_with_one_path(iopt.paths.ref)

    iopt.sml = run_simulation(settings)
    iopt.sml.reset_for_training()

    # Copy reference files
    for imech in iopt.refs:
        f = os.path.join(iopt.sml.paths.res, imech.name)
        if not os.path.exists(f):
            raise FileNotFoundError(f"Results for reference mechanism {imech.name} not found: {f}")
        copy_to_path(f, iopt.paths.ref, imech.runid)

    # Update namelist for training
    for f in os.listdir(iopt.sml.paths.nml):
        if not f.startswith("nml"):
            continue
        # Update reference files
        if iopt.sml.box not in rfiles_key:
            raise ValueError(f"Invalid box name in reference files string: {iopt.sml.box} in {rfiles_key}")
        iref = rfiles_key[iopt.sml.box]  # Get variable name
        f = os.path.join(iopt.sml.paths.nml, f)
        cmd = f"sed -i '/{iref}/c\\{iref} = \"{iopt.err_files_str}\"' {f}"
        os.system(cmd)
    iopt.sml.nerr = len(iopt.refs)  # Update error files number

    iopt.update_ref = 0  # Reset


def _update_mech_diff(iopt: CycleOption) -> None:
    """Check if need to examine mechanism differences between reference and pre mechanisms."""

    iopt.update_mech = False  # Reset

    if not iopt.update_ref:  # No need to check
        logger.info("No need to check mechanism differences w/ update_ref: %d", iopt.update_ref)
        return

    if iopt.update_ref != 3:  # Have to check
        iopt.update_mech = True
        logger.info("Check mech diffs between refs and pre mechanisms w/ update_ref: %d", iopt.update_ref)
        return

    # Check mechanism differences, no check if mechs are the same
    if not iopt.refs:
        raise ValueError("No reference mechanisms to check.")
    for imech in iopt.refs:
        if iopt.mechv.isdiff(imech):
            iopt.update_mech = True
            logger.info("Ref mechanism %s is different from pre valid mechanism %s.", imech.name, iopt.mechv.name)
            return


def _update_ref_mechanisms(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Update error files for the current reduction stage."""

    # Parse reference cases
    err_files, mech_run = [], []
    for s in iopt.stage.ref_cases:
        if s.lower() in ["ref", "pre"]:
            mech_id = s.lower()
            imech = trn_opt.mechs[0] if mech_id == "ref" else iopt.mechv
        elif os.path.exists(s):  # Read mechanism
            mech_id = os.path.basename(s)  # Get mechanism name
            logger.info("Read reference mechanism from %s and save as %s.", s, mech_id)
            imech = Mechanism(name=mech_id, path=s, wfolder=False, mopt=iopt.mechv.mopt)
            imech.copy_and_link(trn_opt.path_mech)
        else:
            raise FileNotFoundError(f"Reference mechanism not found: {s}")
        err_files.append(os.path.join(iopt.paths.ref, mech_id))
        mech_run.append(imech)

    if not err_files:
        raise ValueError(f"No reference mechanisms to run. Check {iopt.stage.ref_cases}")

    iopt.err_files_str = ",".join(err_files)  # Update error files string
    iopt.refs = mech_run
    logger.info("Current reference for reduction evaluation: %s from %s", iopt.err_files_str, iopt.stage.ref_cases)


def _setup_mechanism_for_cycle(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Setup mechanism for the current reduction cycle."""

    # Elementory-like treatment - update reaction list if needed
    _elementory_treatment(trn_opt, iopt)

    # Prepare mechanism for reduction
    iopt.mechv.set_for_training(iopt.update_kinetic)
    iopt.update_kinetic = False  # Reset

    # Update reduction info
    iopt.update_infos(trn_opt)

    # Record starting mechanism
    iopt.logt.write(f"\nStarting mechanism: {iopt.mechv.get_minfo()}\n")


def _elementory_treatment(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Active elementory-like treatment if No.reactions changes"""

    size0 = iopt.mechv.get_size()  # Initial size

    # Trim mechanism w/ additional options
    if trn_opt.restart_from:  # No add opts
        add_mopt = {}
    else:  # Merge or split
        if iopt.stage.flg_split_prod:  # Elementory-like treatment
            add_mopt = {"split": True, "merge": False}
        else:
            add_mopt = {"split": False, "merge": trn_opt.with_merge}
    logger.info("Trimming mechanism %s w/ addiitonal options: %s", iopt.mechv.name, add_mopt)
    iopt.mechv.trim(add_mopt)  # Trim mechanism

    if not iopt.mechv.minfo.is_valid:
        iopt.stop("Mechanism is not valid after trimming. Stop training.")
        return

    # Update mechanism name if needed
    if iopt.mechv.minfo.is_changed:
        new_name = f"{iopt.mechv.name}m" if add_mopt["merge"] else f"{iopt.mechv.name}s"  # m: merge, s: split
        iopt.mechv.save_new(new_name, trn_opt.path_mech)  # Save new mechanism
        iopt.logt.write(f"Mechanism is changed after trimming w/ {add_mopt}. Original size: {size0}")
        iopt.logt.write(f"New mechanism saved as: {iopt.mechv.name}.\n{iopt.mechv.minfo.get_pinfo()}")
        iopt.update_mech = True
    else:  # No change
        logger.info("Mechanism size %s does not change. No new mechanism saved.", iopt.mechv.get_size())


def _check_mech_diff(trn_opt: TrainingOption, iopt: CycleOption) -> None:
    """Check how different between the previous mechanism from the reference mechanism."""

    logger.info("Current mechanism errors: %s", iopt.mechv.errs)
    if not iopt.update_mech:
        logger.info("Skip checking differences in pre/ref mechs as update_mech is False.")
        if iopt.mechv.errs is None:
            iopt.mechv.errs = [np.zeros(len(iopt.refs)), np.zeros(len(iopt.refs))]
            logger.info("Set up default errors for mechanism %s: %s", iopt.mechv.name, iopt.mechv.errs)
        return

    logger.info("Checking new differences between previous & reference mechanisms ...")

    # Link mech files if not exist
    if iopt.mechv.name not in os.listdir(iopt.paths.mech):
        logger.info("Link mechanism files for mechanism %s to %s.", iopt.mechv.name, iopt.paths.mech)
        link_to_path(os.path.join(trn_opt.path_mech, iopt.mechv.name), iopt.paths.mech)
    # Run simulations
    iopt.sml.update_mech_names([iopt.mechv.name])
    iopt.sml.wlock = False
    iopt.sml = run_simulation(iopt.sml)

    if not iopt.sml.err_arr:
        iopt.stop("No valid error files found for checking mechanism differences.")
        return

    # Update errors
    ierrs = iopt.sml.err_arr.arrs[0]
    if np.any(np.isnan(ierrs) | (ierrs == -1)):
        iopt.stop(f"Invalid errors found for checking mechanism differences. Got: {ierrs}")
        return
    iopt.mechv.errs = [np.max(ierrs, axis=(0, 2)), np.mean(ierrs, axis=(0, 2))]
    logger.info("Update errors for mechanism %s: %s", iopt.mechv.name, iopt.mechv.errs)

    # Check errors for stopping
    if trn_opt.stop_at_emax and max(iopt.mechv.errs[0]) > trn_opt.stop_at_emax:
        iopt.stop(f"Stop training due to mechv max err {iopt.mechv.errs[0]} > {trn_opt.stop_at_emax}.")
        return
    if trn_opt.stop_at_eave and max(iopt.mechv.errs[1]) > trn_opt.stop_at_eave:
        iopt.stop(f"Stop training due to mechv ave err {iopt.mechv.errs[1]} > {trn_opt.stop_at_eave}.")
        return

    # Update flag
    iopt.update_mech = False
    iopt.sml.reset_for_training()
    logger.info("Finished checking mech_diff for mechanism %s.", iopt.mechv.name)
