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
This module analyzes the relationships among reactions and species within the same reaction tree.
These relationships helps to define candidate reductions.
"""


from collections import deque
from typing import Optional

from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_attrs_from_settings


# Logger
logger = setup_logger(__name__)


def update_sps_list_with_gen(reactions: list, species: list, savfile: Optional[str] = None) -> list:
    """Update species status, generation, precursors by primary VOCs."""

    # Get dict of generation and precursors
    gens = get_source_and_generation(reactions, savfile)

    # Update species list with generation and precursors
    logger.info("Updating species list with generation and precursors from primary VOCs.")
    # Loop over species list
    sps_no_gen, nval = [], 0
    for s in species:
        if s.name in gens:
            s.generation, precursors = gens[s.name]
            # Convert to list for json output
            s.precursors = sorted(precursors)
            nval += 1
        elif s.status == 1:
            # Reducible species not from primary VOCs
            sps_no_gen.append(s)
    logger.info("Updated # %d species. # %d species w/o generation info.", nval, len(sps_no_gen))

    # Inactivate species not from primary VOCs
    for s in sps_no_gen:
        s.status = 0

    return sps_no_gen


def get_downstream_infos(sname: str, infos: dict) -> dict:
    """Breadth-first search to find all downstream products of a species."""

    if sname not in infos:
        return {"species": set(), "ircns": set()}
    # Initialize queue with the starting species: sname
    queue = deque([sname])
    visited = set()
    all_products = set()
    all_reactions = set()

    while queue:
        sname = queue.popleft()
        if sname in visited:
            continue

        # Get products and reaction info
        pdts = infos.get(sname, {}).get("species", set())
        all_products.update(pdts)
        all_reactions.update(infos.get(sname, {}).get("ircns", set()))

        # Add products to queue
        for s in pdts:
            if s not in visited:
                queue.append(s)
        visited.add(sname)

    return {"species": all_products, "ircns": all_reactions}


def get_downstream_species(sname: str, infos: dict) -> set:
    """Breadth-first search to find all downstream products of a species."""

    if sname not in infos:
        return set()

    queue = deque([sname])
    visited = set()
    all_products = set()
    while queue:
        sname = queue.popleft()
        if sname in visited:
            continue
        pdts = infos.get(sname, {}).get("species", set())
        all_products.update(pdts)
        for s in pdts:
            if s not in visited:
                queue.append(s)
        visited.add(sname)

    return all_products


def get_species_relation(reactions: list, is_reversed: bool = False) -> dict:
    """
    Generates a dictionary (or a tree) of information related to
      the destruction or production (reverse) of species from the reaction list.
    """
    # Get settings
    basic_dict = get_attrs_from_settings({SNAME.SETUP: ["basicspecies_dict"]})

    rcn_info = {}
    for i, rcn in enumerate(reactions):
        if rcn.status <= 0:
            continue
        # Get sets
        if is_reversed:
            species_set = [s for s in rcn.products if s not in basic_dict]
            related_set = [s for s in rcn.reactants if s not in basic_dict]
        else:
            species_set = [s for s in rcn.reactants if s not in basic_dict]
            related_set = [s for s in rcn.products if s not in basic_dict]
        # Record sets
        for s in species_set:
            if s not in rcn_info:
                rcn_info[s] = {"species": set(), "ircns": set()}
            rcn_info[s]["species"].update(related_set)
            rcn_info[s]["ircns"].add(i)
    return rcn_info


def get_source_and_generation(reactions: list, savefile: Optional[str] = None) -> dict:
    """Compute generations for all nodes in the tree starting from the start species."""

    # Get settings
    primary_vocs = get_attrs_from_settings({SNAME.SETUP: ["primary_vocs"]})

    # Get species relationships
    chd_info = get_species_relation(reactions)

    # Get paths
    paths, sources = {}, {}
    for start in primary_vocs:
        paths[start] = [start]  # Record paths from precursors
        sources[start] = {start}  # Record precursors
    queue = deque([(start, [start]) for start in primary_vocs])
    while queue:
        s0, path = queue.popleft()
        if s0 not in chd_info:
            continue  # End species
        npath = len(path) + 1
        for s1 in chd_info[s0]["species"]:
            # Check and record sources
            if s1 not in sources:
                sources[s1] = set()
            if path[0] not in sources[s1]:
                sources[s1].add(path[0])
            # Check and record the shortest path
            # Add to generation list
            if s1 not in paths or npath < len(paths[s1]):
                paths[s1] = path + [s1]
                queue.append((s1, path + [s1]))

    # Obtain the generation and sources for each species
    outs, maxgen = {}, 0
    for p, vals in paths.items():
        igen = len(vals)
        maxgen = max(maxgen, igen)
        outs[p] = [igen, sources[p]]

    # Save to file
    if savefile:
        pinfos = ["species,generation,sources,path"]
        for p, vals in paths.items():
            pinfos.append(f"{p},{len(vals)},{'|'.join(sources[p])},{' -> '.join(vals)}")
        with open(savefile, "w", encoding="utf-8") as f:
            f.write("\n".join(pinfos))
        logger.info("Saved generations and sources to file: %s", savefile)

    logger.info("Find generations for # %d species with max No.gen at %d.", len(outs), maxgen)

    return outs


def get_source(reactions: list, chd_info: Optional[dict] = None) -> dict:
    """Compute sources for all nodes in the tree starting from the start species."""

    # Get settings
    primary_vocs = get_attrs_from_settings({SNAME.SETUP: ["primary_vocs"]})

    if chd_info is None:
        chd_info = get_species_relation(reactions)

    # Get paths
    sources = {}
    for start in primary_vocs:
        sources[start] = {start}  # Record precursors
    queue = deque([(start, [start]) for start in primary_vocs])
    while queue:
        s0, path = queue.popleft()
        if s0 not in chd_info:
            continue  # End species
        for s1 in chd_info[s0]["species"]:
            # Check and record sources
            if s1 not in sources:
                sources[s1] = set()
            if path[0] not in sources[s1]:
                sources[s1].add(path[0])
            # Check and record the shortest path
            queue.append((s1, path + [s1]))

    return sources
