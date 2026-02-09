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
This module contains functions to find reduction candidates based on jumping strategy.
"""

from typing import Optional

import numpy as np

from .logger import setup_logger
from .reduction_basic import ReductionListInfo
from .reduction_setting import rdc_checks
from .utils import calculate_normalized_ratios


# Logger
logger = setup_logger(__name__)


def find_reduction_via_jumping(
    reactions: list, species: list, infos: ReductionListInfo, change_scheme: bool = False
) -> list:
    """Find reduction candidates based on jumping strategy."""

    # Get checking for jumping
    jp_checks = rdc_checks["jp"]
    ngmax, nsmax = jp_checks["ngmax"], jp_checks["nsmax"]
    rtlim = jp_checks.get("lrtmin", None)  # Ratio threshold for negligible species

    reduction_info, jumped = [], set()  # Record reduction info

    # Find reduction candidates over targeted species
    for s in infos.targeted_sps:
        if len(reduction_info) >= ngmax:
            break  # Break if enough groups

        if s in jumped:
            continue  # Skip if already jumped

        # Get children for s to jump
        children = set()  # Record children species names
        if s in infos.children_info:
            for ircn in infos.children_info[s]["ircns"]:
                children.update(reactions[ircn].products)  # Get all products

        if not children:
            # logger.info("No children species for jumping %s", s)
            continue  # No children to jump

        children = list(children)  # Convert to list
        crcns = {ircn: reactions[ircn] for ircn in infos.children_info[s]["ircns"]}

        # Remove negligible children species and reactions
        remove_negligibles(children, crcns, rtlim)

        if len(children) > nsmax:
            continue  # Too many children species

        children_ratios = _get_jumping_ratios(children, crcns)

        if children_ratios is None:  # No ratios
            # logger.info("No ratios for jumping %s to %s", s, children)
            continue

        # Get reduction info
        rdc_info = {
            "jumped": s,
            "species": children,
            "ratios": children_ratios,
            "crcns": list(crcns.keys()),
            "prcns": infos.parent_info[s]["ircns"] if s in infos.parent_info else [],
        }

        # Record reduction info
        reduction_info.append(rdc_info)
        jumped.add(s)

        # Update reduction info after jumping
        parents = infos.parent_info.get(s, {"species": []})["species"]
        update_infos_after_jumping(s, children, parents, infos)

        # Change scheme via jumping
        if change_scheme:
            change_scheme_via_jumping(reactions, species, rdc_info, infos.sps_info)

    # Resume infos if not change_scheme:
    if not change_scheme:
        infos.update(reactions, species, {"renew": "pc"})

    return reduction_info


def _get_jumping_ratios(children: list, crcns: list) -> Optional[list]:
    """Get jumping ratios for children species."""

    # Get reduction ids
    crcns_ids = get_reaction_ids(crcns)

    if len(crcns_ids) == 1:  # Only one type of reactions
        return get_children_ratio_wsame_type(children, crcns)

    # Multiple types of reactions
    return compute_complex_jumping_ratios(children, crcns, crcns_ids)


def get_reaction_ids(rcns: dict) -> dict:
    """Return reaction information from input reaction list."""

    rinfos = {}  # Record rcn ids
    for rcn in rcns.values():
        rcn_type = rcn.rcn_id

        if rcn_type not in rinfos:
            rinfos[rcn_type] = []

        rinfos[rcn_type].append(rcn)

    return rinfos


def get_children_ratio_wsame_type(spnames: list, rcns: dict) -> Optional[list]:
    """Return ratios for species in given reactions with the same type."""

    if not rcns:
        return None
    if len(rcns) == 1:
        return next(iter(rcns.values())).get_ratios(slist=spnames)

    # Multiple reactions
    rcn_rts = {i: rcn.rate.get_ratio() for i, rcn in rcns.items()}

    # Normalize ratios
    rcn_rts = dict(zip(rcn_rts.keys(), calculate_normalized_ratios(rcn_rts.values())))

    sps_rts = np.zeros(len(spnames))  # Get species ratios
    for i, rcn in rcns.items():
        sps_rts += rcn.get_ratios(slist=spnames) * rcn_rts[i]

    return sps_rts


def get_max_children_ratios(spnames: list, rcns: dict) -> list:
    """Return ratios for species in given reactions."""

    rcn_list = list(rcns.keys())
    rcn_rts = np.array([rcns[i].kuni[1] for i in rcn_list]).T
    rcn_nrts = np.array([calculate_normalized_ratios(rt) for rt in rcn_rts]).T

    # Get max ratios for children species
    sps_rts = np.zeros(len(spnames))
    for i, s in enumerate(rcn_list):
        sps_rts += rcns[s].get_ratios(slist=spnames) * np.max(rcn_nrts[i])

    return sps_rts


def remove_negligibles(spnames: list, rcns: dict, rtlim: Optional[float]) -> list:
    """
    Remove species from jumping list if all ratios are below threshold.
    They are considered to be removed directly.
    """

    # Check ratio threshold
    if not rtlim or rtlim <= 0.0:
        return [spnames, rcns]

    if rtlim >= 0.5:
        raise ValueError(f"Invalid ratio threshold: {rtlim}")

    # Get max ratios for children species in shape of (nsp, nrcn)
    ratios = get_max_children_ratios(spnames, rcns)

    # Remove species if all ratios are below threshold
    sps_rm = [spnames[i] for i, ratio in enumerate(ratios) if ratio < rtlim]

    if sps_rm:
        logger.info("Find # %d negligible species with threshold %s: %s in %s", len(sps_rm), rtlim, sps_rm, spnames)
        # Remove species from list
        for sp in sps_rm:
            spnames.remove(sp)

    return sps_rm


def compute_complex_jumping_ratios(children: list, crcns: dict, crcns_ids: dict) -> Optional[list]:
    """Compute jumping ratios for children species with multiple types of reactions."""

    # When one temp & rh, should be the same
    # !!! return get_max_children_ratios(children, crcns)  # Shape: (nsp, nrcn)

    return None  # Not checked


def update_infos_after_jumping(spname: str, children: list, parents: list, infos: ReductionListInfo) -> None:
    """Update reduction information for next jumping."""

    # Update parent info
    if spname in infos.parent_info:
        info = infos.parent_info.pop(spname)
        for p in children:
            if p in infos.parent_info:
                infos.parent_info[p]["species"].update(info["species"])
                infos.parent_info[p]["ircns"].update(info["ircns"])

    # Update children info
    if spname in infos.children_info:
        info = infos.children_info.pop(spname)
        for p in parents:
            if p in infos.children_info:
                infos.children_info[p]["species"].update(info["species"])
                infos.children_info[p]["ircns"].update(info["ircns"])


def change_scheme_via_jumping(reactions: list, species: list, rdc_info: dict, sps_info: dict) -> None:
    """Change mechanism based on jumping strategy."""

    # Update species list
    jumped = rdc_info["jumped"]
    species[sps_info[jumped]].status = 0

    # Update reaction list
    children = rdc_info["species"]
    ratios = rdc_info["ratios"]

    # Change parents -> targeted species * ratios
    for ircn in rdc_info["prcns"]:
        rcn = reactions[ircn]
        # Update products
        if jumped in rcn.products:
            irt = rcn.products[jumped]
            for ip, p in enumerate(children):
                rcn.products[p] += irt * ratios[ip]
            del rcn.products[jumped]
            rcn.clean_products()

    # Change targeted -> children
    for ircn in rdc_info["crcns"]:
        rcn = reactions[ircn]
        if jumped in rcn.reactants:
            rcn.status = 0
