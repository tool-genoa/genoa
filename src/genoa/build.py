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
This module contains the main function to build a chemical mechanism using the GENOA algorithm.
"""

from .constants import SNAME
from .gecko_in import rename_species
from .logger import setup_logger
from .mechanism_in import read_and_update_mechanism, read_and_update_genoa_mech
from .mechanism_out import mech_output
from .setting_global import get_all_settings_map
from .species import update_species_condensibility
from .ssh_aerosol import update_soap_strucs_from_file
from .utils import is_mech_exist


logger = setup_logger(__name__)


def build_new_mechanism() -> None:
    """Build a chemical mechanism using the GENOA v3.0 algorithm."""

    # Get settings
    settings_map = get_all_settings_map()
    gnl = settings_map[SNAME.SETUP]()
    nopt = settings_map[SNAME.NEW]()

    # Read reactions and species lists
    ropt = {"check_rate": True, "clean": nopt.update_mode, "merge": nopt.if_merge, "tracer": gnl.tracers}
    logger.info("Generating mechanism from %s with options: %s", nopt.reactionfile, ropt)
    reactions, species, _ = read_and_update_mechanism(
        nopt.reactionfile,
        nopt.speciesfile,
        nopt.reactiontype,
        nopt.speciestype,
        ropt,
    )

    # Update soap structures
    if nopt.soapfile:
        update_soap_strucs_from_file(species, nopt.soapfile)

    # Update species condensibility
    update_species_condensibility(species, nopt.psat_svoc, nopt.psat_nvoc)

    # Rename species if needed for GECKO-A
    if gnl.boxmodel == "GECKO":
        logger.info("Checking species names for using in GECKO-A...")
        rename_species(reactions, species)

    # Save new mechanism
    logger.info("Saving new mechanism %s to %s", gnl.mech_name, gnl.path_sav_mech)
    out_options = {"add_fake": nopt.add_fake_radical}
    if nopt.output_modes:  # Add output modes if specified
        out_options["out_modes"] = nopt.output_modes
    mech_output(gnl.path_sav_mech, gnl.mech_name, reactions, species, out_options)


def build_fake_mechanism(mech_name: str, check_exist: bool) -> None:
    """Build a mechanism with fake radical species."""

    # Get settings
    gnl = get_all_settings_map()[SNAME.SETUP]()

    # Return if fake mechanism already exists
    if check_exist and is_mech_exist(gnl.path_read_mech, mech_name, True):
        logger.info("Fake mechanism %s already exists.", mech_name)
        return
    logger.info("Generating fake mechanism %s ...", mech_name)

    # Read original mechanism
    ropt = {"with_folder": True, "check_rate": True, "clean": "wsps", "tracer": gnl.tracers}
    reactions, species, _ = read_and_update_genoa_mech(gnl.path_read_mech, gnl.mech_name, ropt)

    # Output fake mechanism
    out_options = {"add_fake": True}
    mech_output(gnl.path_sav_mech, mech_name, reactions, species, out_options)
    logger.info("Fake mechanism is generated.")
    return
