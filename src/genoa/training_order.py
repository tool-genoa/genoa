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
This module is used to reorganize species list to groups for reduction.
"""

import os
import math
from typing import Optional, TextIO

import numpy as np

from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_attrs_from_settings
from .species import Species
from .reduction_basic import ReductionListInfo
from .reduction_setting import group_order_build


# Logger
logger = setup_logger(__name__)


def read_group_from_file(sfile: str, species: list, frozen: set) -> list:
    """Read the ordered species list from a file."""

    # Check if file exists
    logger.info("Reading species groups for training order from file: %s", sfile)
    if not (sfile and os.path.exists(sfile)):
        raise ValueError(f"Species group file not found: {sfile}")

    # Read species and ranks
    sps_in, irank = {}, 0
    with open(sfile, "r", encoding="utf-8") as f:
        for line in f.readlines():

            # Examine line
            line = line.strip()
            if not line or line.startswith("#"):  # Skip empty or comments
                continue
            if ":" in line:  # Remove index part
                line = line.split(":", 1)[1].strip()

            # Get species from line
            isps = [s.strip() for s in line.split(" ") if s.strip()]
            if not isps:
                continue

            # Record
            sps_in.update({s: irank for s in isps})
            irank += 1
    logger.info("Read # %d species from # %d lines.", len(sps_in), irank + 1)

    # Check species and update rank info for current run
    sps_out, unsps = {}, []
    for isp in species:
        if isp.status != 1 or isp.name in frozen:
            continue
        s = isp.name
        # Get highest rank for reducible species
        irank = sps_in.get(s, -1)
        for sp in isp.reductions:
            for p in sp:
                irank = max(irank, sps_in.get(p, -1))
        if irank == -1:  # Not in any group
            logger.warning("Species %s not found in any group in file %s.", s, sfile)
            unsps.append(s)
        else:  # Record species current rank
            sps_out[s] = irank
    if not sps_out:
        raise ValueError(f"No valid species found in file {sfile}.")

    logger.info("Reorganized # %d / # %d current valid / invalid species.", len(sps_out), len(unsps))

    # Get new group list
    sorted_sps = sorted(sps_out, key=sps_out.get)
    groups, irank = [], -1
    for s in sorted_sps:
        i = sps_out[s]  # Get rank of species
        if i != irank:  # New group
            groups.append([])
            irank = i
        groups[-1].append(s)
    groups.append(unsps)  # Add unsorted species

    logger.info("Updated # %d valid species groups for # %d species.", len(groups), len(species))

    return groups


def species_list_to_group(species: list, infos: ReductionListInfo, savfile: str) -> list:
    """
    Reorginize species list to species groups for reduction
    based on input group tyeps
    """

    # Get settings
    settings = get_attrs_from_settings({SNAME.SETUP: ["primary_vocs"], SNAME.TRN: ["nsps_group"]})
    primary_vocs = settings[SNAME.SETUP]
    nsmax = settings[SNAME.TRN]
    gtype, rtype = group_order_build["group_types"], group_order_build["order_types"]

    logger.info("Generating training order with group types: %s and order types: %s.", gtype, rtype)
    if nsmax > 0:
        logger.info("Max number of species in a group: %s", nsmax)
    else:
        logger.info("No max number of species limit for grouping.")

    # Generate a group list for reduction
    if rtype:

        # Get pvoc mass & psat
        pvocs = [species[infos.sps_info[s]] for s in primary_vocs]
        mass_pvoc, psat_pvoc = _get_pvoc_properties(pvocs)

        group_info, rank_info = {}, {}  # Record info used for group & rank species
        # Record species by group and rank infos
        for s in infos.targeted_sps:
            if s not in infos.sps_info:
                continue  # Not valid species
            isp = species[infos.sps_info[s]]
            igp = _get_group_name(isp, gtype if gtype else ["name"])
            if igp not in group_info:  # Check group id
                group_info[igp] = []  # Get new group
                # Get group rank: the larger the LESS important
                rank_info[igp] = _get_group_rank(isp, rtype, mass_pvoc, psat_pvoc)
            # Record species to its group
            group_info[igp].append(s)

        # Order species group by their mass
        sorted_keys = sorted(group_info.keys(), key=lambda k: rank_info[k])
        sorted_group = [group_info[s] for s in sorted_keys]

    elif gtype:
        raise ValueError(f"rtype can not be empty if gtype has been defined: {gtype}.\n")

    else:  # No gtype or rtype
        sorted_group = [[s] for s in infos.targeted_sps]
    logger.info("# %d species are grouped for reduction.", len(infos.targeted_sps))

    if nsmax > 0:  # Break large groups
        logger.info("Break large groups with more than %d species.", nsmax)
        new_group = []
        for s in sorted_group:
            ns = len(s)
            if ns > nsmax:
                new_group.extend([s[i : i + nsmax] for i in range(0, ns, nsmax)])
            else:
                new_group.append(s)
        sorted_group = new_group

    if savfile:  # Save to file
        logger.info("Save species groups to file: %s", savfile)
        _save_group_to_file(sorted_group, savfile, rtype, gtype)

    return sorted_group


def _get_pvoc_properties(primary_vocs: Optional[list]) -> list:
    """Calculate mass and saturation pressure values for primary VOCs."""

    if primary_vocs:
        mass_pvoc = np.average([isp.mass for isp in primary_vocs])
        psat_pvoc = math.exp(np.average([math.log(isp.psat_atm) for isp in primary_vocs]))
        logger.info("PVOC mass: %.2f, Geometric mean psat: %.2e", mass_pvoc, psat_pvoc)
    else:  # No primary_vocs, use default values
        mass_pvoc, psat_pvoc = 150, 1e-3
        logger.info("No primary VOCs found. Use default values: mass_pvoc=150, psat_pvoc=1e-3.")
    return [mass_pvoc, psat_pvoc]


def _save_group_to_file(group: list, savfile: str, rtype: list, gtype: list) -> None:
    """Save the ordered species list to a file."""
    with open(savfile, "w", encoding="utf-8") as f:
        if rtype and gtype:
            _write_grouped_species_to_file(f, group)
        elif rtype:
            for i, s in enumerate(group):
                f.write(f"{i}: {' '.join(s)}\n")
        else:
            for s in group:
                f.write(f"{' '.join(s)}\n")


def _write_grouped_species_to_file(f: TextIO, sorted_species: list) -> None:
    """Write grouped species with rank info to file."""
    nmax, ndis, imax = 1, {}, None
    for i, s in enumerate(sorted_species):
        n = len(s)
        nmax = max(n, nmax)
        if n == nmax:
            imax = i
        if n in ndis:
            ndis[n] += 1
        else:
            ndis[n] = 1
        f.write(f"No.{i} contains # {n}: " + " ".join(s) + "\n")
    if imax is None:
        raise ValueError("No group with more than 1 species.")
    logger.info("Max group No.%d with %d species. Distribution: %s.", imax, nmax, ndis)


def _get_group_name(isp: Species, gtype: list) -> str:
    """Get group infomation to rank and distribute species"""

    gid = []
    # Get group id & rk - the larger the LESS important
    for s in gtype:
        if s == "formula":
            gid.append(f"{isp.formula}")
        elif s == "gen":
            gid.append(f"{isp.generation}")
        elif s == "note":
            gid.append(isp.note)
        elif s == "stype":
            gid.append(isp.rdc_id)
        elif s == "basic":
            gid.append(isp.rdc_id[0])
        elif s == "name":
            gid.append(isp.name)
        else:
            raise ValueError(f"Unknown mode: {s} for group determination.")

    return "-".join(gid)


def _get_group_rank(isp: Species, rtype: list, mass_pvoc: float, psat_pvoc: float) -> float:
    """Get group infomation to rank and distribute species"""

    grk = 1.0
    for s in rtype:
        if s == "mass":
            grk *= abs(isp.mass - mass_pvoc) / (isp.mass + mass_pvoc)
        elif s == "gen":
            grk *= abs(isp.generation)
        elif s == "psat":
            grk /= abs(math.log(isp.psat_atm) - psat_pvoc)
        elif s in ["note", "type"]:
            if isp.radical:
                grk *= 10
            elif isp.condensable:
                grk *= 0.1
        else:
            raise ValueError(f"Unknown mode: {s} for order determination.")

    return grk
