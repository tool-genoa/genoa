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
This module contains classes and functions for recording outputs to screen and/or file.
"""

import os
from typing import Optional, TextIO, List, Any

import attrs

from .constants import SNAME, LINE_SEP
from .logger import setup_logger


# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class RecordOption:
    """
    Options for recording the output to screen and/or file.
    """

    muted: bool = False  # Mute all output: used for multiprocessing
    fname: str = ""  # Path for the record file
    fout: Optional[TextIO] = None  # File object to write to

    def __attrs_post_init__(self) -> None:
        """Post initialization for the RecordOption class."""
        if self.fout:
            if not self.fname:
                self.fname = self.fout.name
            elif self.fname != self.fout.name:
                raise ValueError("Inconsistent record file and fout: {self.fname} vs {self.fout.name}")

    def mute(self) -> None:
        """Mute all output if needed."""
        if self.muted:
            return
        self.muted = True
        self.close()

    def unmute(self) -> None:
        """Unmute all output if needed."""
        self.muted = False

    def open(self) -> None:
        """Open the record file for writing if needed."""
        if self.fout or self.muted:
            return
        if not (self.fname and os.path.exists(self.fname)):
            raise FileNotFoundError(f"Cannot find record file: {self.fname}")
        # logger.info("Open record file: %s ...", self.fname)
        self.fout = open(self.fname, "a", encoding="utf-8")
        if not self.fout:
            raise IOError(f"Cannot open record file: {self.fname}")

    def close(self) -> None:
        """Close the record file if opened."""
        if not self.fout:
            return
        # logger.info("Close record file: %s ...", self.fout.name)
        self.fout.close()
        self.fout = None

    def write(self, msg: str) -> None:
        """Write the message to the file and/or log it."""
        if self.muted:
            return
        if not self.fout:
            self.open()
        self.fout.write(f"{msg}\n")

    def flush(self) -> None:
        """Flush the record file if opened."""
        if self.fout:
            self.fout.flush()

    def list_to_file(self, msg_list: List[str], msg: str = "", sep: str = "\n") -> None:
        """Write the message with a list of strings to the file and/or log it."""

        if self.muted:
            return
        if not self.fout:
            self.open()
        if msg:
            self.fout.write(msg + sep)
        self.fout.writelines(line + sep for line in msg_list)


def build_log(log_file: str, msg: str = "") -> RecordOption:
    """Build a new log object."""
    prepare_logfile(log_file, msg)
    log = RecordOption(fname=log_file)
    return log


def prepare_logfile(logfile: str, header: str = "") -> None:
    """Prepare log file for writing later."""

    if not logfile:
        raise ValueError("No log file is given.")

    # Check if exists
    if os.path.exists(logfile):
        os.remove(logfile)

    # Create the folder if not exists
    else:
        log_dir = os.path.dirname(logfile)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)

    # Initialize the log file with the given header
    with open(logfile, "w", encoding="utf-8") as f:
        f.write(header)


def setup_loge_4action(act: str, gnl: "GlobalSetting", opt: Any) -> None:
    """Initialize loge file for different actions if needed. Return log file name."""

    if not act or act not in [SNAME.SML, SNAME.TBR, SNAME.TST] or opt.loge:
        return

    logger.info("Setup new loge file for action %s ...", act)

    if act == SNAME.SML:
        if not opt.ref_files_str:
            return  # No file
        filename = f"error_sml_{gnl.mech_name}.dat"
        pinfo = f"# Error file for simulations w/ {gnl.mech_name} and error_ref_str: {opt.ref_files_str}"
    else:
        if act == SNAME.TBR and not opt.tag_sim:
            return  # No file
        filename = f"error_{act}_{opt.ref_mech_name}_{opt.runid}.dat"
        pinfo = f"# Error file for {act} w/ ref {opt.ref_mech_name} & test {opt.mech_name}"
        if gnl.err_ref_str:
            pinfo += f" and error_ref_str: {gnl.err_ref_str}"

    opt.loge = build_log(os.path.join(gnl.path_sav_mech, filename), f"{pinfo}. PID: {os.getpid()}\n")

    logger.info("Record error for action %s to loge file: %s", act, opt.loge.fname)


def record_initialization(settings_map: dict, opt: Any, act: str, flog: RecordOption) -> None:
    """Record the initialization of training & testing in loge files."""

    if not flog or flog.muted:
        return

    if not act:
        raise ValueError("Invalid input action.")

    # Get settings
    gnl = settings_map[SNAME.SETUP]()
    sml = settings_map[SNAME.SML]()

    act0 = act.capitalize()

    pinfos = [
        f"{LINE_SEP}{act0} with run id {opt.runid} and box model {gnl.boxmodel}.",
        f"{act0} mechanism: {opt.mech_name} in {opt.mech_path}",
        f"Reference mechanism: {opt.ref_mech_name} in {opt.ref_mech_path}",
        f"Time setting(s) # {len(gnl.time_settings)}: {gnl.time_settings}",
    ]

    add_dict = {
        SNAME.TRN: get_pinfo_for_training,
        SNAME.TST: get_pinfo_for_testing,
    }

    if act in add_dict:
        pinfos.append(add_dict[act](opt))

    # Add common information
    if sml.error_species_str:
        pinfos.append(f"Error species string: {sml.error_species_str}")

    if sml.init_set_str:
        pinfos.append(f"Initialization set string: {sml.init_set_str} w/ start id {sml.init_id0}")
    flog.write("\n".join(pinfos))


def get_pinfo_for_testing(opt: "TestingOption") -> str:
    """Get information for testing."""
    pinfos = []

    # Read reference results
    if opt.read_ref:
        pinfos.append(f"Read reference results from {opt.read_ref}")
    else:
        pinfos.append("No reading reference results. Generate during testing.")

    # Conditions
    pinfos.append(f"Testing # {len(opt.conds) - 1} condition(s) in {opt.path_cond}")

    # Log files
    if opt.loge:
        pinfos.append(f"Error log in {opt.loge.fname}")

    return "\n".join(pinfos)


def get_pinfo_for_training(opt: "TrainingOption") -> str:
    """Get information for training."""

    pinfos = [f"\nTraining stages in total: {opt.params_dict['nstage']}"]

    for k, v in opt.params_dict.items():
        if k == "nstage":
            continue
        if isinstance(v[0], list) and len(v[0]) > 1:
            v_str = "\n  ".join(f"# {i+1} {j}" for i, j in enumerate(v))
            pinfos.append(f"{k}: {v_str}")
        else:
            pinfos.append(f"{k}: {v}")

    pinfos.append(f"\nSpecies need to be kept: {opt.kept_all}")
    pinfos.append(f"Generated mechanisms are saved in {opt.path_mech}")

    return "\n".join(pinfos)
