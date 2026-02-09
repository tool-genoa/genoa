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
This module contains the settings for box models used for simulations.
"""

import os

from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_attrs_from_settings
from .utils import is_same_link, build_link


# Logger
logger = setup_logger(__name__)


# Settings for box model compilation
COMPILE_GCK = {
    # "cmd": "./build.sh > tmptoto 2>&1",  # Command to compile - read from settings
    "source": "PROG/boxmod",  # Source file to check if compilation is needed
    "core": "PROG/boxmod.f90",  # Core file to check if can be compiled
    "toto": "tmptoto",  # Temporary file to check
    "mesgs": ["NO ERROR"],  # Messages to check if compilation is successful
    "links_in": ["PROG/boxmod", "INTERP30/gecko2box"],  # Links needed to build for running the model
    "links_out": ["boxmod", "gecko2box"],  # Name in the output directory
}

COMPILE_SSH = {
    # "cmd": "./clean > tmptoto 2>&1 && ./compile > tmptoto 2>&1",
    "source": "src/ssh-aerosol",
    "core": "src/ssh-aerosol.f90",
    "toto": "tmptoto",
    "mesgs": ["scons: done building targets."],
    "links_in": ["src/ssh-aerosol"],
    "links_out": ["ssh-aerosol"],
}


def _get_cmp_settings(box: str) -> dict:
    """Get the compilation settings for a given model type."""

    cmp_dict = {
        "GECKO": COMPILE_GCK,
        "SSH": COMPILE_SSH,
    }

    if box not in cmp_dict:
        raise ValueError(f"Box model type {box} not found in compile settings.")
    cmp_d = cmp_dict[box]

    # Add section name
    cmp_d["sec"] = SNAME.BOX[box]

    return cmp_d


def set_boxmodel(dir_out: str, box: str, out_exec: str = "") -> str:
    """Compile box model"""

    # Get compile settings
    cmp_d = _get_cmp_settings(box)
    # Get path to box and command to build
    dir_in, cmd_b = get_attrs_from_settings({cmp_d["sec"]: ["path", "cmd_build"]})

    icompile = False  # If force to recompile

    if not os.path.exists(dir_in):
        raise FileNotFoundError(f"Cannot find box model path: {dir_in}. Check.")

    source = os.path.join(dir_in, cmp_d["source"])
    logger.debug("Setting up box model from %s to  %s ...", dir_in, dir_out)

    # Check if needs to recompile
    if not icompile and not (os.path.exists(source) and os.path.isfile(source)):
        icompile = True

    # Compile box model if need
    if icompile:

        logger.debug("Compiling box model %s ...", box)

        # Check core file exists
        if not os.path.exists(cmp_d["core"]):
            raise FileNotFoundError(f"Cannot find core file {cmp_d['core']}. Check.")

        # Go to dir_in
        wd = os.getcwd()
        os.chdir(dir_in)

        # Compile
        os.system(cmd_b)

        # Check if message in toto file after compilation
        with open(cmp_d["toto"], "r", encoding="utf-8") as f:
            info = f.read()
        for merr in cmp_d["mesgs"]:
            if merr not in info:
                raise ValueError(f"Cannot compile. Not found {merr}. Check file: {dir_in}/{cmp_d['toto']}.")

        # Return from dir_in
        os.chdir(wd)

    # Build links
    os.makedirs(dir_out, exist_ok=True)
    for s_in, s_out in zip(cmp_d["links_in"], cmp_d["links_out"]):
        src = os.path.join(dir_in, s_in)
        tgt = os.path.join(dir_out, s_out)
        if not is_same_link(src, tgt):
            build_link(src, tgt)
        logger.debug("Finished linking %s %s.", src, tgt)

    # Output executable path
    if not out_exec:
        return os.path.join(dir_out, cmp_d["links_out"][0])

    if out_exec not in cmp_d["links_out"]:
        raise ValueError(f"Output executable {out_exec} not found in {cmp_d['links_out']}.")

    return os.path.join(dir_out, out_exec)
