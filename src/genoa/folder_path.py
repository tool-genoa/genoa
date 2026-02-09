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
This module contains classes and functions for managing folder paths used in running simulations.
"""

import os
import shutil
from typing import Optional

import attrs

from .logger import setup_logger


# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class FolderPath:
    """Folder paths for saving simulation results."""

    work: str = ""  # Root folder for running testing
    res: str = ""  # Root folder for saving results
    rec: str = ""  # Root folder for saving log files
    nml: str = ""  # Root folder for saving namelist files

    # Used only for training
    mech: str = ""  # Root folder for saving mechanism files (training)
    ref: str = ""  # Root folder for saving reference files (training)

    def init_with_general(self, gnl, wid, mech_names) -> None:
        """Initialize paths with global settings."""
        self.work = gnl.path_workplace
        if wid:
            self.res = gnl.path_sav_res
        else:  # Result folder
            self.res = update_output_path(f"{gnl.path_sav_res}/Results_{mech_names[0]}", True)
        self.nml = self.res
        self.rec = self.res

    def init_with_path_work(self) -> None:
        """Initialize paths for using in pre-testing or testing process."""

        if not self.work:
            raise FileNotFoundError("No working path for initialization.")

        # Update paths with path_work
        self.work = update_output_path(self.work)
        self.res = update_output_path(os.path.join(self.work, "RES"))
        self.rec = update_output_path(os.path.join(self.work, "LOG"))
        self.nml = update_output_path(os.path.join(self.work, "NML"))

    def init_in_training(self) -> None:
        """Initialize paths for using in training process."""

        if not self.work:
            raise FileNotFoundError("No working path for training.")

        # Update paths with path_work
        self.work = update_output_path(self.work)
        self.mech = update_output_path(os.path.join(self.work, "CHEM"))
        self.nml = update_output_path(os.path.join(self.work, "NML"))
        self.res = update_output_path(os.path.join(self.work, "RES"))
        self.rec = update_output_path(os.path.join(self.work, "LOG"))
        self.ref = update_output_path(os.path.join(self.work, "REF"))

    def init_with_one_path(self, path: str) -> None:
        """Initialize paths for using in testing cycle."""

        if not path:
            raise FileNotFoundError("No path for initialization.")

        # Update empty paths with the given path
        inew = update_output_path(path)
        paths = [self.work, self.res, self.rec, self.nml]
        for i, path0 in enumerate(paths):
            if not path0:
                paths[i] = inew

    def update(self):
        """Update paths."""
        for path in [self.work, self.res, self.rec, self.nml, self.mech, self.ref]:
            if path and not os.path.exists(path):
                update_output_path(path)

    def update_path(self, new_path: str, arr: str) -> None:
        """Update new_path for given attribute."""
        if not new_path:
            raise ValueError("Empty new path for updating.")
        if not hasattr(self, arr):
            raise AttributeError(f"Attribute {arr} not found in {self}")
        # Update path
        setattr(self, arr, new_path)
        update_output_path(new_path)

    def clean_in_training(self) -> None:
        """Clean paths and prepare for next reduction cycle."""
        for path in [self.res, self.rec, self.nml, self.ref]:
            update_output_path(path, del_exist=True)


def copy_to_path(source: str, targeted_path: str, new_name: str = "", remove_exist: bool = False) -> str:
    """Copy source file or folder to the targeted path."""
    if not (source and os.path.exists(source)) or not targeted_path:
        raise ValueError("Empty source or target path for copying.")

    # Get file/folder name
    fname = new_name if new_name else os.path.basename(source)
    target = os.path.join(targeted_path, fname)

    # Remove existing file/folder
    if os.path.exists(target):
        if not remove_exist:
            raise FileExistsError(f"File/folder already exists: {target}")
        shutil.rmtree(target, ignore_errors=True)

    # Copy file or directory
    if os.path.isdir(source):
        shutil.copytree(source, target)
    else:
        shutil.copyfile(source, target)
    logger.info("Copied %s to %s", source, target)
    return target


def move_to_path(source: str, targeted_path: str, new_name: str = "", remove_exist: bool = False) -> str:
    """Move source file or folder to the targeted path."""
    if not (source and os.path.exists(source)) or not targeted_path:
        raise ValueError("Empty source or target path for moving.")

    # Get the new file/folder name
    fname = new_name if new_name else os.path.basename(source)
    target = os.path.join(targeted_path, fname)

    # Remove existing file/folder
    if os.path.exists(target):
        if not remove_exist:
            raise FileExistsError(f"File/folder already exists: {target}")
        shutil.rmtree(target, ignore_errors=True)

    # Move file or directory
    shutil.move(source, target)

    return target


def link_to_path(source: str, targeted_path: str, new_name: str = "", remove_exist: bool = False) -> str:
    """Create a symbolic link to the source file or folder in the targeted path."""
    if not (source and os.path.exists(source)) or not targeted_path:
        raise ValueError("Empty source or target path for linking.")

    # Get the new file/folder name
    fname = new_name if new_name else os.path.basename(source)
    target = os.path.join(targeted_path, fname)

    # Remove existing file/folder
    if os.path.exists(target):
        if not remove_exist:
            raise FileExistsError(f"File/folder already exists: {target}")
        shutil.rmtree(target, ignore_errors=True)

    # Create a symbolic link
    os.symlink(source, target)

    return target


def update_input_path(path: str, check_exist: bool = False) -> str:
    """Update input path to an absolute path and optionally check for existence."""

    if check_exist:
        if not path:
            raise ValueError("Empty path.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input path does not exist: {path}")

    return os.path.abspath(path)


def update_output_path(path: str, del_exist: bool = False) -> str:
    """Update output path to an absolute path and optionally delete existing content."""

    if not path:
        raise ValueError("Output path is empty.")

    if del_exist and os.path.exists(path):
        shutil.rmtree(path, ignore_errors=True)

    os.makedirs(path, exist_ok=True)

    return os.path.abspath(path)


def get_common_files(fpath: str, fnames: Optional[list], suffix: str = ".png") -> list:
    """Output a list of file names if they exist in the folder and input fnames list"""

    if not suffix:
        return []

    if not (fpath and os.path.exists(fpath) and os.path.isdir(fpath)):
        logger.error("Invalid path for updating file names with %s : %s", suffix, fpath)
        return []

    pics = [f for f in os.listdir(fpath) if f.endswith(suffix)]
    if not pics:
        logger.warning("No figures found in %s with suffix %s", fpath, suffix)
        return []

    if fnames is None:
        return pics
    return [f for f in pics if f in fnames]
