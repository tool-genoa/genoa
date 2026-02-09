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
This module provides functions to input and output of chemical mechanisms.

"""

import os
from typing import Optional

from .constants import SNAME
from .gecko_out import mech_output_to_gecko
from .logger import setup_logger
from .setting_global import get_attrs_from_settings
from .ssh_aerosol import mech_output_to_ssh
from .utils import get_mechanism_size_info


# Logger
logger = setup_logger(__name__)


def _get_fake_species(species: list) -> dict:
    """Get fake species for output"""
    sps_set, fake_species = set(), {}
    for isp in species:
        if isp.status <= 0:
            continue
        sps_set.add(isp.name)
        if isp.status == 1 and isp.radical:
            fake_species[isp.name] = isp.mass

    # Check if fake species exists
    for prefix in ["FA", "f"]:
        for s, isp in fake_species.items():
            if f"{prefix}{s}" in sps_set:
                raise ValueError(f"Fake species {prefix}{s} already exists.")
    logger.info("Find %d fake species.", len(fake_species))

    return fake_species


def _get_output_modes() -> list:
    """Get output modes from settings"""
    box = get_attrs_from_settings({SNAME.SETUP: "boxmodel"})
    return [box]


def mech_output(path: str, chem: str, reactions: list, species: list, options: Optional[dict] = None) -> bool:
    """
    Generate output files. Return True if successful.
    Output options contains:
    - out_modes: [GECKO, viz, SSH, SSH_all]
    - add_fake: True, False (add fake radical species)
    - with_folder: True, False (create folder {path}/{chem} for output)
    """

    # Get options
    # Default: {"out_modes": [gnl.boxmodel], "add_fake": False, "with_folder": True}
    if options is None:
        options = {}

    # Check output folder
    if options.get("with_folder", True):
        path = os.path.join(path, chem)

    # Tag to check if default reaction & species list is generated
    has_default = False

    # Check if add fake species
    if options.get("add_fake", False):
        fakes = _get_fake_species(species)
    else:
        fakes = None

    # Get output modes
    out_modes = options.get("out_modes", _get_output_modes())
    for mode in out_modes:
        is_out = True  # Set output flag
        if mode == "SSH":
            is_out = mech_output_to_ssh(path, chem, reactions, species, fakes)
            # Default species list is generated in SSH mode
            has_default = True
        elif mode == "GECKO":
            is_out = mech_output_to_gecko(path, chem, reactions, species, fakes)
        else:
            logger.error("Unknown output mode %s", mode)
        if not is_out:
            return False

    # Generated default reaction & species list
    if not has_default:
        mech_output_to_default(path, chem, reactions, species)

    # Print info
    logger.info(get_mechanism_size_info(reactions, species))

    return True


def mech_output_to_default(path: str, chem: str, reactions: list, species: list) -> None:
    """Generate default output files"""

    # Default / SSH-aerosol reaction list
    with open(os.path.join(path, f"{chem}.reactions"), "w", encoding="utf-8") as f:
        f.writelines(rcn.to_rcn(i) for i, rcn in enumerate(reactions) if rcn.status > 0)
        f.write("END\n")

    # Default species list
    with open(os.path.join(path, f"{chem}.mol"), "w", encoding="utf-8") as f:
        for isp in species:
            if isp.status > 0:
                f.writelines(isp.output_to_lines())
