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
This module contains functions to perform reduction based on the reaction tree.
"""

import os
import logging
import numpy as np

from .jumping import find_reduction_via_jumping
from .kinetics import compute_bratio, compute_product_yield, compute_yield_in_tree
from .logger import setup_logger, list_to_log, isout
from .mechanism_relation import get_downstream_species, update_sps_list_with_gen
from .reduction_basic import ReductionListInfo
from .reduction_print import get_jp_pinfos
from .unit_conversion import ppb_to_molec


# Logger
logger = setup_logger(__name__)


def find_reduction_via_reaction_tree(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Reduction based on the reaction tree in species and reactions."""

    # Perform the reduction based on reduction type [function, code for infos.update()]
    rdc_fnc_dict = {
        "nvoc": [nvoc_reduction, "s"],
        "nvoc_tree": [nvoc_reduction, "sc"],
        "svoc": [svoc_reduction, "s"],
        "generation": [generation_reduction, "s"],
        "svoc_tree": [svoc_tree_reduction, "sc"],
        "conc_tree": [conc_tree_reduction, "scr"],
        "pyield": [pyield_reduction, "s"],
        "pyield_tree": [pyield_reduction_in_tree, "srck"],
        "bratio": [bratio_reduction, "srk"],
        "lifetime": [lifetime_reduction, "spckrl"],
        "extra_sps_rdc": [extra_species_reduction, "s"],  # Extra removal for species
    }

    if rdc_set[0] not in rdc_fnc_dict:
        raise ValueError(f"Unknown reduction type: {rdc_set[0]} in reduction set: {rdc_set}")

    # Update reduction info
    infos.update(reactions, species, {"renew": rdc_fnc_dict[rdc_set[0]][1]})
    # Conduct reduction
    return rdc_fnc_dict[rdc_set[0]][0](reactions, species, rdc_set, infos)


def nvoc_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Set non-volatile species based on Psat threshold."""
    if rdc_set[1] < 0 or rdc_set[1] > 1e-3:
        raise ValueError(f"Invalid Psat threshold for NVOCs: {rdc_set[1]}")

    valids, values = [], []
    for s in infos.targeted_sps:
        isp = species[infos.sps_info[s]]
        if isp.condensable and isp.psat_atm < rdc_set[1]:
            isp.non_volatile = True
            valids.append(s)
            values.append(isp.psat_atm)

    _log_reduction(reactions, species, rdc_set, valids, values, infos)

    # Remove reactions derived from NVOCs
    if rdc_set[0] == "nvoc_tree":
        # Remove reactions
        for s in valids:
            if s not in infos.children_info:
                continue
            for ircn in infos.children_info[s]["ircns"]:
                reactions[ircn].status = 0  # Inactive


def svoc_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Set condensable species to volatile based on Psat threshold."""

    # Check Psat threshold
    threshold = rdc_set[1]
    if threshold <= 0.0 or threshold > 1.0:
        raise ValueError(f"Invalid Psat threshold for SVOCs: {threshold}")

    valids, values = [], []
    for s in infos.targeted_sps:
        isp = species[infos.sps_info[s]]
        if isp.condensable and isp.psat_atm > threshold:
            isp.condensable = False  # Set to volatile
            valids.append(s)
            values.append(isp.psat_atm)

    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def generation_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Remove species with generation > ngen."""

    # Check generation threshold
    ngen = rdc_set[1]
    if ngen <= 0:
        raise ValueError(f"Invalid generation threshold: {ngen}")

    # Update generation number
    update_sps_list_with_gen(reactions, species)

    valids, values = [], []
    nmax = 0  # Max generation
    for s in infos.targeted_sps:
        isp = species[infos.sps_info[s]]
        nmax = max(nmax, isp.generation)  # Update max generation
        if isp.generation > ngen:
            valids.append(s)
            values.append(isp.generation)
            isp.status = 0  # Inactive

    _log_reduction(reactions, species, rdc_set, valids, values, infos, f"Max gen of targeted species: {nmax}")


def svoc_tree_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Remove species with min(psat, childen_psat) > threshold."""

    # Check Psat threshold
    if rdc_set[1] <= 0 or rdc_set[1] > 1:
        raise ValueError(f"Invalid Psat threshold for SVOCs: {rdc_set[1]}")

    valids, values = [], []
    for s in infos.targeted_sps:
        isp = species[infos.sps_info[s]]
        children = get_downstream_species(s, infos.children_info)
        ipsat = isp.psat_atm if isp.condensable else 1e9
        # Get min Psat of children
        for child in children:
            if child in infos.sps_info:
                jsp = species[infos.sps_info[child]]
                if jsp.condensable:
                    ipsat = min(ipsat, jsp.psat_atm)
            if ipsat <= rdc_set[1]:  # Already below threshold
                break
        # May need to add treatment for species produces many non-condensable species

        # Find species to remove
        if ipsat > rdc_set[1]:
            valids.append(s)
            values.append(ipsat)
            isp.status = 0  # Inactive

    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def conc_tree_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Perform reduciton based on the reference concentrations of species and its children."""

    # Check reference concentrations
    if not infos.ref_concs:
        raise ValueError("Reference concentrations not provided.")
    tree_concs = infos.ref_concs["total"]

    # Check concentration threshold
    threshold = rdc_set[1]
    if threshold <= 0:
        raise ValueError(f"Invalid concentration threshold: {threshold}")
    # Convert ppb to molec/cm3
    threshold = ppb_to_molec(threshold, 298, 1)
    logger.info("Convert threshold from %s ppb to %s molec/cm3", rdc_set[1], threshold)

    valids, values, terminal_count = [], [], 0
    for s in infos.targeted_sps:
        isp = species[infos.sps_info[s]]
        children = get_downstream_species(s, infos.children_info)
        if not children:  # Terminal species, do not remove
            terminal_count += 1
            continue
        # Get reference concentrations
        if s not in tree_concs:
            logger.warning("Reference concentration not found for species: %s", s)
        iconcs = tree_concs[s]
        for p in children:
            if p in infos.sps_info:
                if p not in tree_concs:
                    logger.warning("Reference concentration not found for species: %s", p)
                iconcs += tree_concs[p]
            if iconcs >= threshold:  # Already above threshold
                break
        # Find species to remove
        if iconcs < threshold:
            valids.append(s)
            values.append(iconcs)
            isp.status = 0  # Inactive

    _log_reduction(reactions, species, rdc_set, valids, values, infos, f"Terminal species count: {terminal_count}")


def pyield_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Perform reduction based on the product yield and generation number."""

    # Check threshold
    threshold = rdc_set[1]
    if threshold < 0 or threshold > 0.5:
        raise ValueError(f"Invalid product yield threshold: {threshold}")

    # Get product yields: {rdc_infos: [species, yield ratios, ircns]}
    yields = compute_product_yield(reactions)
    # Get species with small yields {species: [rcns]}
    smalls = {}
    for rinfo, [sps, ylds, ircns] in yields.items():
        for i, s in enumerate(sps):
            if s not in infos.sps_info:
                continue
            igen = species[infos.sps_info[s]].generation
            if igen < 1:
                continue
            y = max(ylds[i]) / igen  # Max yield/ No.gen
            if y < threshold:
                if s not in smalls:
                    smalls[s] = []
                smalls[s].append([y, ylds[i], igen, ircns, rinfo])

    # Find products to remove
    valids, values = [], []
    for s in infos.targeted_sps:
        if s not in smalls:
            continue
        valids.append(s)
        pinfos = []
        for [y, y0, igen, ircns, rinfo] in smalls[s]:
            pinfos.append(f" pyield: {y0}, gen: {igen}, tyield: {y} in rcns: {ircns} w/ info: {rinfo}")
            for ircn in ircns:  # Remove products
                del reactions[ircn].products[s]
        values.append(";".join(pinfos))

    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def pyield_reduction_in_tree(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Perform reduction based on the product yield and reaction time."""

    if not infos.ref_concs:
        raise ValueError("Reference concentrations not provided.")

    # Check threshold
    threshold = rdc_set[1]
    if threshold < 0 or threshold > 0.5:
        raise ValueError(f"Invalid product yield threshold: {threshold}")

    # max_yld: {species: max yield}
    max_yld = compute_yield_in_tree(reactions, infos.ref_concs["gas"], infos.children_info, infos.kunis)

    # Find products to remove
    valids, values = [], []
    for s in infos.targeted_sps:
        if s not in max_yld:
            continue
        if max_yld[s] < threshold:
            valids.append(s)
            values.append(max_yld[s])
            species[infos.sps_info[s]].status = 0  # Inactive species

    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def bratio_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Perform reduction based on the branching ratios of species in reactions."""

    # Check reference concentrations
    if not infos.ref_concs:
        raise ValueError("Reference concentrations not provided.")

    # Check threshold
    threshold = rdc_set[1]
    if threshold < 0 or threshold > 0.5:
        raise ValueError(f"Invalid branching ratio threshold: {threshold}")

    # Get branching ratios {species: [bratios, ircns]}
    bratios = compute_bratio(reactions, infos.ref_concs["gas"], infos.kunis)

    valids, values = [], []
    for s in infos.targeted_sps:
        if s not in bratios:  # No branching ratio
            continue
        ratios, ircns = bratios[s]
        for i, ircn in enumerate(ircns):
            max_rt = max(ratios[i])
            if max_rt < threshold:  # All bratio below threshold
                valids.append(ircn)
                values.append(max_rt)
                reactions[ircn].status = 0  # Inactive

    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def lifetime_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Perform reduction if the lifetime of species is faster than the threshold."""

    # Check threshold
    threshold = rdc_set[1]
    if threshold < 0 or threshold > 1e6:
        raise ValueError(f"Invalid lifetime threshold: {threshold} s")

    valids, values = [], []
    for s in infos.targeted_sps:
        if s not in infos.ltimes:  # Lifetime is infinite
            continue
        ltime = np.min(infos.ltimes[s])  # Get lifetime at different temperatures
        if ltime < threshold:  # Can be jumped
            valids.append(s)
            values.append(ltime)

    _log_reduction(reactions, species, rdc_set, valids, values, infos)

    if valids:  # Apply jumping
        infos.targeted_sps = valids
        reduction_info = find_reduction_via_jumping(reactions, species, infos, True)
        list_to_log(logger, get_jp_pinfos(reduction_info), f"\n-- Jumping with # {len(reduction_info)} reduction(s). ")


def extra_species_reduction(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Remove species that meets input criteria."""

    einfo_func_dict = {
        "maxyield": _read_max_yield,
        "gridvalues": _read_grid_values,
    }

    # Get species info from files
    info_files = rdc_set[1]
    if not (info_files and isinstance(info_files, list)):
        raise ValueError("Invalid extra species reduction info: empty or not a list.")

    # Read extra species reduction info
    extra_sps, extra_infos = set(), {}
    for i, info in enumerate(info_files):
        # Check input
        if not isinstance(info, list) or len(info) != 3:
            raise ValueError(f"Invalid extra species reduction info: {info}")
        ifile, itype, lim = info

        # Check file type and read
        if itype not in einfo_func_dict:
            raise ValueError(f"Unknown extra species reduction type: {itype} in info: {info}")
        einfos = einfo_func_dict[itype](ifile, lim)
        extra_infos[itype] = einfos
        extra_sps.update(einfos.keys())  # Add species to remove

    # Remove species based on extra info
    valids, values = [], []
    for s in infos.targeted_sps:
        if s not in extra_sps:
            continue
        valids.append(s)
        values.append("; ".join([f"{i}: {v[s]}" for i, v in enumerate(extra_infos) if s in v]))
        # Remove species
        species[infos.sps_info[s]].status = 0

    logger.info("Extra species reduction with # %d species.", len(valids))
    _log_reduction(reactions, species, rdc_set, valids, values, infos)


def _read_max_yield(ifile: str, lim: float) -> dict:
    """Read species can be reduced based on maxyield.dat file by gecko-a."""

    if not (ifile and os.path.isfile(ifile)):
        raise ValueError(f"Invalid file: {ifile} for reading max yield.")

    esps, nsps = {}, 0
    with open(ifile, "r", encoding="utf-8") as f:
        lines = f.readlines()
        for line in lines:
            parts = [i.strip() for i in line.split() if i.strip()]
            val = float(parts[0])
            if val < lim:  # Below threshold
                esps[parts[1]] = val
            nsps += 1

    # Print
    list_to_log(
        logger,
        [f"  # {i}: {s} with value: {v}" for i, (s, v) in enumerate(esps.items())],
        f"\n-- Read max yield from file: {ifile} of {nsps} species < {lim}",
    )
    logger.info("Read max yield from file: %s of %s species < %s", ifile, nsps, lim)
    return esps


def _read_grid_values(ifile: str, lim: float) -> set:
    """Read species can be reduced based on gridvalues.csv file from Dan."""

    if not (ifile and os.path.isfile(ifile)):
        raise ValueError(f"Invalid file: {ifile} for reading grid values.")

    esps, nsps = {}, 0
    with open(ifile, "r", encoding="utf-8") as f:
        lines = f.readlines()[1:]  # Skip header
        for line in lines:
            parts = [i.strip() for i in line.split(",") if i.strip()]
            if len(parts) < 4:
                continue
            val = float(parts[2]) * float(parts[3])  # Psat * gamma
            if val > lim:
                esps[parts[0]] = val
            nsps += 1

    # Print
    list_to_log(
        logger,
        [f"  # {i}: {s} with value: {v}" for i, (s, v) in enumerate(esps.items())],
        f"\n-- Read grid values from file: {ifile} of {nsps} species < {lim}",
    )
    logger.info("Read grid values from file: %s of %s species < %s", ifile, nsps, lim)
    return esps


def _log_reduction(
    reactions: list, species: list, rdc_set: list, valids: list, values: list, infos: ReductionListInfo, add: str = ""
) -> None:
    """Record the reduction details to the file or console."""
    if not isout(logger, logging.INFO):
        return

    pinfos = [f"\n\n-- Searched reduction with type {rdc_set[0]} and value {rdc_set[1]}:"]
    pinfos.append(f"  Target # {len(infos.targeted_sps)} from {len(species)} species & {len(reactions)} reactions.")
    pinfos.append(f"  *-- Find {len(valids)} valids.")

    if add:  # Add extra info if needed
        pinfos.append(add)

    list_to_log(logger, [f"  # {i}: {s} with value: {values[i]}" for i, s in enumerate(valids)], "\n".join(pinfos))
