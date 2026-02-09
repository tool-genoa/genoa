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
This module contains classes and functions for initializing and running simulations.
"""
import os
from typing import List, Optional, Union

import numpy as np
from attrs import define, field

from .constants import SNAME
from .folder_path import FolderPath
from .logger import setup_logger
from .setting_global import get_all_settings_map
from .record import RecordOption, setup_loge_4action
from .utils import get_condition_list


# Logger
logger = setup_logger(__name__)


@define(slots=True)
class ErrorArray:
    """Record and process simulation errors."""

    # Error list input from simulations in shape of [nml(nt, ncond)*chem*init, nerr, nsps]
    arr0: list = field(factory=list)
    nchem: int = 0  # Number of chemical mechanisms
    has_nan: bool = False  # If has np.nan

    # Error array in shape of [nerr, nchem, nml(nt, ncond)*init, nsps] after reshaping
    arrs: Optional[np.ndarray] = None

    # Error statistics
    emax: Optional[float] = None
    eave: Optional[float] = None
    emax_locs: Optional[List[tuple]] = None

    def __attrs_post_init__(self):
        """Check the input error list and initialize the error array."""
        if not self.arr0:
            raise ValueError("Input error list is empty.")
        self.set_errors()  # Reshape error array

    def set_errors(self) -> None:
        """Update errors"""

        arrs = np.array(self.arr0)

        nsim, nerr, nsps = arrs.shape  # Check size consistency
        if self.nchem < 1 or nsim % self.nchem != 0:
            raise ValueError(f"Cannot reshape arr0: {nsim} with {self.nchem} mechanisms.")

        # Reshape: (chem * nml * init, nerr, nsps) -> (nchem, nml * init, nerr, nsps)
        try:
            arrs = arrs.reshape(self.nchem, nsim // self.nchem, nerr, nsps)
        except ValueError as e:
            raise ValueError(f"Error reshaping array: {e}. With shape: {arrs.shape}") from e

        # Check if has np.nan
        if np.any(np.isnan(arrs)):
            self.has_nan = True
        else:
            self.has_nan = False

        self.arrs = arrs

    def update_statistics(self) -> None:
        """Update error statistics."""
        self.emax = np.max(self.arrs)  # Max error
        self.eave = np.mean(self.arrs)  # Average error
        self.emax_locs = list(zip(*np.where(self.arrs == self.emax)))  # Locations of max error


@define(slots=True)
class RunSmlSetting:
    """Settings for running simulations."""

    # Mechanism
    mech_names: List[str] = field(factory=list)
    nchem: int = 0  # Number of mechanisms

    # Paths
    mech_path: str = ""  # Chemical mechanism path
    paths: Optional[FolderPath] = None  # Paths for simulation results

    # Conditions
    path_cond: str = ""  # Path to initial conditions
    conds: List[str] = field(factory=list)  # List: path + conditions
    ncond: int = 0  # Number of conditions

    init_sets: Optional[List[str]] = None  # Initial condition sets
    ninit: int = 1  # Number of initial conditions

    # Error computation
    ref_files_str: Optional[str] = None  # Reference files/folders
    nerr: int = 0  # Number of errors

    err_sps: Optional[List[str]] = None  # Error species
    nsps: int = 1  # Number of species in error computation

    err_arr: Optional[ErrorArray] = None  # Save obtained error array

    # Output mode: 0 - no output, 1 - fast (default) output, 2 - full output
    out_mode: Optional[int] = None
    loge: Optional[RecordOption] = None  # For outputing error tsv file
    soa_path: Optional[str] = None

    # Operators
    box: Optional[str] = None  # Box model
    npara: Optional[int] = None  # Number of parallel simulations
    box_exec: Optional[str] = None  # If the box model needs to be compiled
    labels: List[str] = field(factory=list)  # If namelists need to be written
    nnml: int = 0  # Number of namelists

    wid: bool = False  # If run commands with ids
    chemids: List[str] = field(factory=list)  # Chemical mechanism ids
    initids: List[str] = field(factory=list)  # Initial condition ids
    nmlids: List[str] = field(factory=list)  # Namelist ids - updated later by labels
    resids: List[str] = field(factory=list)  # Result ids - updated later by labels

    # For parallel simulations
    err_checks: Optional[list] = None  # Max err tolerance: [max_err, delta_err, pre_max_err]
    scores: Optional[list] = None  # Scores used for locking
    wlock: bool = False  # If run parallel simulations with lock

    def update(self) -> None:
        """Update settings with default or given values."""

        # Get settings
        settings_map = get_all_settings_map()
        gnl = settings_map[SNAME.SETUP]()
        sml_opt = settings_map[SNAME.SML]()

        if not self.mech_names:  # Default mechanism name
            self.mech_names = [gnl.mech_name]
        if not self.mech_path:  # Default mechanism path
            self.mech_path = gnl.path_read_mech
        if self.out_mode is None:  # Output mode
            self.out_mode = sml_opt.result_mode
        if self.ref_files_str is None:  # Error files
            self.ref_files_str = gnl.err_ref_str
        if self.loge is None:  # Error output
            setup_loge_4action(SNAME.SML, gnl, self)

        if len(self.mech_names) > 1 and not self.wid:  # Update wid
            self.wid = True

        self.update_err_sps(sml_opt.error_species_str)  # Error species
        self.update_cond(gnl.path_init_cond)  # Update conditions
        self.update_init_sets(sml_opt.init_set_str, sml_opt.init_id0)  # Update initial conditions
        self.update_counter_and_runid()  # Update ids w/ nchem, ncond, nsps, nerr, ninit

        # Folder paths
        if not self.paths:
            self.paths = FolderPath()
            self.paths.init_with_general(gnl, self.wid, self.mech_names)
        # box model
        if self.box is None:
            self.box = gnl.boxmodel
        # Number of parallel simulations
        if self.npara is None:
            mdict = {"SSH": SNAME.SSH, "GECKO": SNAME.GCK}
            if self.box not in mdict:
                raise ValueError(f"Invalid box model: {self.box}.")
            self.npara = settings_map[mdict[self.box]]().nsim

    def update_err_sps(self, err_sps_str: str) -> None:
        """Update error species with input string."""
        if self.err_sps:  # Already set
            return

        self.err_sps = ["SOA"]
        if err_sps_str:
            for s in err_sps_str.split(";"):
                s = s.strip()
                if s:
                    self.err_sps.append(s)

    def update_cond(self, path_cond: str) -> None:
        """Update condition list with input path."""
        if self.conds:  # Already set
            return

        if not self.path_cond:
            self.path_cond = path_cond

        self.conds = get_condition_list(self.path_cond)

    def update_init_sets(self, init_set_str: str, i0: int) -> None:
        """Update initial condition sets and ids (initids) with input string."""
        if self.init_sets and self.initids:  # Already set
            return

        if init_set_str:
            self.init_sets = ["baisc"] + [s.strip() for s in init_set_str.split(";") if s.strip()]
            ninit = len(self.init_sets)
            if i0 < 0 or i0 >= ninit:
                raise ValueError(f"Invalid init_id0: {i0} with {ninit} sets.")
            self.initids = list(range(i0, ninit))
        elif self.wid:
            self.init_sets = ["basic"]
            self.initids = [0]
        else:
            self.init_sets = None
            self.initids = [""]

    def update_counter_and_runid(self) -> None:
        """Update counters and ids for running simulations."""
        # Update counters
        self.nchem = len(self.mech_names)
        self.ncond = len(self.conds) - 1
        self.nsps = len(self.err_sps)
        self.nerr = self.ref_files_str.count(",") + 1 if self.ref_files_str else 0
        self.ninit = len(self.init_sets) if self.init_sets else 1

        # chem id
        if self.wid:  # With ids
            self.chemids = self.mech_names
        else:
            self.chemids = [""] if self.initids[0] == "" else ["-"]

    def update_errors(self, arr0: list) -> None:
        """Update error array with input list."""
        self.err_arr = ErrorArray(arr0=arr0, nchem=self.nchem)
        self.check_error_shape()

    def check_error_shape(self) -> None:
        """Check error shape."""
        nchem, nsml, nerr, nsps = self.err_arr.arrs.shape
        if nchem != self.nchem or nerr != self.nerr or nsps != self.nsps or nsml != self.nnml * self.ninit:
            raise ValueError(
                f"Invalid error array shape: {self.err_arr.arrs.shape}."
                + f"  Should be w/ {self.nchem} chems, {self.nerr} errors, {self.nsps} species, "
                + f"and {self.nnml * self.ninit} nsml."
            )

    def get_err_loc(self, iloc: tuple) -> str:
        """Get error location string."""
        imech = self.mech_names[iloc[0]]
        inml = self.labels[iloc[1] % self.nnml]
        isp = self.err_sps[iloc[3]]
        return f"{imech} - err {iloc[2]} - {inml} - init {iloc[1] // self.nnml} - {isp}\n"

    def reset_for_training(self) -> None:
        """Reset settings for training."""

        self.err_arr = None  # Reset error array

        # For parallel simulations
        self.err_checks = None
        self.scores = None
        if not self.wlock:
            self.wlock = True  # Set Lock for parallel simulations

        # For mechanism
        if not self.wid:
            self.wid = True
        self.mech_names = None
        self.nchem = 0
        self.chemids = None  # == mech_names

        if self.paths.mech and self.paths.mech != self.mech_path:
            self.mech_path = self.paths.mech  # Reset mechanism path

    def update_mech_names(self, mech_names: List[str]) -> None:
        """Update mechanism names with input list."""
        if not mech_names:
            raise ValueError("Invalid mechanism names.")
        self.mech_names = mech_names
        self.nchem = len(mech_names)
        self.chemids = mech_names if self.wid else [""]


def setup_run_sml(run_settings: Union[RunSmlSetting, None, dict] = None) -> RunSmlSetting:
    """Get an instance of the settings for running simulations."""

    if isinstance(run_settings, RunSmlSetting):
        return run_settings

    # Build new settings
    if run_settings is None:  # Default settings
        settings = RunSmlSetting()
    elif isinstance(run_settings, dict):  # From dictionary
        settings = RunSmlSetting(**run_settings)
    else:  # Invalid input
        raise ValueError(f"Invalid simulation input settings: {run_settings}")

    # Update settings
    settings.update()

    return settings


def check_condition_list(cond_list: list) -> bool:
    """Check if input is a valid condition list."""

    if not (cond_list and isinstance(cond_list, list)):
        return False

    if len(cond_list) < 2:
        logger.error("Invalid condition list with len < 2: %s.", cond_list)
        return False

    if not os.path.isdir(cond_list[0]):
        logger.error("Cannot find condition folder: %s.", cond_list[0])
        return False

    for f in cond_list[1:]:
        if not os.path.exists(os.path.join(cond_list[0], f)):
            logger.error("Cannot find condition file: %s in %s.", f, cond_list[0])
            return False

    return True
