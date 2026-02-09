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
This module contains the class Mechanism that stores the information of a mechanism.
"""


import os
from typing import List, Optional

import attrs

from .folder_path import link_to_path  # , copy_to_path
from .logger import setup_logger
from .mechanism_in import read_and_update_genoa_mech
from .mechanism_out import mech_output
from .mechanism_update import MechInfo, update_mechanism
from .rate_constant import update_rate_cst_values
from .species import get_species_category
from .utils import get_mech_files


# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class Mechanism:
    """Mechanism information"""

    # Default
    runid: str = ""  # Run id
    name: str = ""  # Mechanism name
    path: str = ""  # Mechanism path
    wfolder: bool = True  # If with or without folder

    # Read
    files: Optional[list] = None  # List of mechanism files
    reactions: Optional[list] = None  # List of reactions
    species: Optional[list] = None  # List of species

    # Update or output mechanism
    minfo: Optional[MechInfo] = None  # Mechanism information
    sinfo: Optional[dict] = None  # For postprocessing
    mopt: dict = attrs.field(factory=dict)  # Mechanism update options

    # Simulation
    errs: Optional[list] = None  # List of errors

    def __attrs_post_init__(self) -> None:
        """
        Post initialization.
          Default read options: {with_folder": True, "check_rate": False, "clean": "wgen",
          "merge": False, "split": False, "check_kept": None}
          Default output options:  {"out_modes": [gnl.boxmodel], "add_fake": False, "with_folder": True}
        """
        if not (self.name and self.path):
            raise ValueError("Invalid mechanism name or path")
        if not self.reactions:
            self.read_mech()

    def read_mech(self, nopt: Optional[dict] = None) -> None:
        """Read mechanism species and reactions from the mechanism file."""
        if not (self.name and self.path and os.path.exists(self.path)):
            raise ValueError(f"Invalid mechanism name or path: {self.name}, {self.path}")

        # Read from files
        mopt = _update_options(self.mopt, nopt)
        if not self.wfolder:
            mopt["with_folder"] = self.wfolder
        self.reactions, self.species, self.minfo = read_and_update_genoa_mech(self.path, self.name, mopt)

    def save_new(self, name: str = "", path: str = "") -> None:
        """Save the mechanism with a new name and path."""

        if name:
            self.name = name
        if path:
            self.path = path
        self.errs = None  # Reset errors
        self.output()  # Save mechanism

    def set_for_training(self, update_kinetic: bool) -> None:
        """Update mechanism for training."""

        logger.info("Setting mechanism %s for training ...", self.name)
        if not self.reactions:
            self.read_mech()
        for s in self.species:
            s.update_for_reduction()
        if update_kinetic:
            update_rate_cst_values(self.reactions)

    def isdiff(self, other: "Mechanism") -> bool:
        """Check if the current mechanism is different from the other mechanism."""
        try:
            self.get_files()
            other.get_files()
            for x, y in zip(self.files, other.files):
                if os.path.abspath(x) != os.path.abspath(y):
                    return True
            return False
        except ValueError as e:
            logger.warning("Cannot get files for both mechs. Got: %s", e)
            return True

    def copy_and_link(self, cpath: str, lpath: Optional[str] = None, new_name: Optional[str] = None) -> None:
        """Copy the mechanism to the given cpath and link to lpath if given."""

        if not cpath or not os.path.exists(cpath):
            logger.error("\nInvalid path in copying mechanism. Got: %s", cpath)
            return
        if os.path.abspath(self.path) == os.path.abspath(cpath):
            logger.error("\nCannot copy mechanism %s to the same path %s.", self.name, cpath)
            return

        # Copy mechanism
        # if self.wfolder:  # Copy mechanism folder
        #     copy_to_path(os.path.join(self.path, self.name), cpath, new_name, True)
        # else:  # Copy mechanism file
        #     for ifile in self.get_files():
        #         copy_to_path(ifile, cpath, new_name, True)
        #     self.wfolder = True  # Reset after copy

        self.path = cpath
        self.name = new_name if new_name else self.name
        self.output()
        self.wfolder = True  # Reset after copy

        # Link mechanism
        if lpath:
            if not os.path.exists(lpath):
                logger.error("\nInvalid link path in linking mechanism. Got: %s", lpath)
            if self.wfolder:
                link_to_path(os.path.join(self.path, self.name), lpath, new_name, True)
            else:
                for ifile in self.get_files():
                    link_to_path(ifile, lpath, new_name, True)
            self.path = lpath

    def output(self) -> None:
        """Save mechanism."""

        if not self.reactions:  # Read mechanism
            logger.info("No reactions found for saving mechanism %s. Reading mechanism ...", self.name)
            self.read_mech()

        if not self.minfo.is_valid:  # No output
            logger.info("No output for mechanism %s (is_valid = 0). Got: %s", self.name, self.minfo.get_pinfo())
            return

        # Save mechanism
        self.minfo.is_valid = mech_output(self.path, self.name, self.reactions, self.species)

    def trim(self, nopt: Optional[dict] = None) -> None:
        """Clean and merge the mechanism if needed."""

        # Update trimming options
        mopt = _update_options(self.mopt, nopt)
        if not mopt:
            raise ValueError(f"Invalid trimming options for mechanism {self.name}")
        logger.info("Trimming mechanism %s with options: %s", self.name, mopt)
        if self.reactions:  # Update mechanism
            self.reactions, self.species, self.minfo = update_mechanism(self.reactions, self.species, mopt)
        else:  # Read & update
            self.read_mech(mopt)

    def get_runid(self) -> str:
        """Get runid or mech name"""
        return self.runid if self.runid else self.name

    def get_files(self) -> List[str]:
        """Get mechanism files."""
        if not self.files:
            self.files = get_mech_files(self.path, self.name, self.wfolder)
        return self.files

    def get_size(self) -> List[int]:
        """Get size of the mechanism."""
        return self.minfo.read_size()

    def get_minfo(self) -> str:
        """Get size infos."""

        # Name
        minfos = [f"{self.name}"]
        if self.runid:
            minfos.append(f"Runid: {self.runid}")

        # Size
        minfo = ",".join(str(s) for s in self.get_size())
        minfos.append(f"Mechanism size (# reducible reactions, species, condensables: {minfo})")

        # Errors
        if self.errs:
            minfos.append(f"Training Errors: {self.errs}")
        return "\n".join(minfos)


def _update_options(idefault: dict, inew: Optional[dict] = None) -> dict:
    """Return a new dictionary with updated options."""
    if not inew:
        return idefault
    updated = idefault.copy()
    updated.update(inew)
    return updated


def load_mechanisms(mech_names: list, mech_paths: list, update_kinetic: bool, wpvocs: Optional[set]) -> list:
    """Get a list of Mechanism objects for the given mechanism names and paths."""

    logger.info("Loading # %d mechanisms ...", len(mech_names))
    mechs = []  # List of mechanisms
    for name, path in zip(mech_names, mech_paths):
        mechs.append(Mechanism(name=name, path=path))
        mechs[-1].read_mech()

    # Update kinetic rate values
    if update_kinetic:
        logger.info("Updating kinetic rate values ...")
        for mech in mechs:
            update_rate_cst_values(mech.reactions)

    # Add species category
    for mech in mechs:
        mech.sinfo = get_species_category(mech.species, wpvocs)

    logger.info("Finished loading mechanisms: %s.", mech_names)
    return mechs
