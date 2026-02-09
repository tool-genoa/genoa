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
This module contains functions to compute kinatic parameters for reactions.
"""


from collections import deque
from typing import Optional

import numpy as np

from .constants import SNAME
from .mechanism_relation import get_species_relation
from .reaction import get_reduction_info_from_list
from .setting_global import get_attrs_from_settings
from .utils import calculate_normalized_ratios


def compute_kuni(reactions: list, ref_concs: dict, update_rcn: bool = False) -> dict:
    """Compute pseudo-unimolecular rate constant (kuni) for reactions based on reference concentrations."""
    # Get settings
    basic_dict = get_attrs_from_settings({SNAME.SETUP: ["basicspecies_dict"]})

    kuni = {}  # reaction index: [reactant, kunis]
    for ircn, rcn in enumerate(reactions):
        if rcn.status != 1:  # Find reducible reactions
            continue
        # Get species and concentrations to convert k to kuni
        sps = rcn.get_sps_for_rate(basic_dict)
        iconc = 1.0
        for s in sps[1:]:
            iconc *= ref_concs[s]
        kuni[ircn] = [sps[0], rcn.rate.kvalues * iconc]

    if update_rcn:
        for ircn, ks in kuni.items():
            reactions[ircn].kuni = ks

    return kuni


def reshape_kuni_wsps(kuni: dict) -> dict:
    """Reorganize kuni by reactants"""
    if not kuni:
        return {}

    kuni_sps = {}
    for ircn, [rct, ks] in kuni.items():
        if rct not in kuni_sps:
            kuni_sps[rct] = [[], []]  # rcns, kunis
        kuni_sps[rct][0].append(ircn)
        kuni_sps[rct][1].append(ks)

    return kuni_sps


def compute_product_yield(reactions: list) -> dict:
    """Compute yields for products in the reactions with the same reduction info."""

    # Get settings
    basic_dict = get_attrs_from_settings({SNAME.SETUP: ["basicspecies_dict"]})

    yields = {}  # Dict: info: [species, yield ratios, ircns]
    rdc_info = get_reduction_info_from_list(reactions)

    for info, ircns in rdc_info.items():
        ylds = {}  # Get yields for non-basic products
        for ircn in ircns:
            rcn = reactions[ircn]
            if rcn.status != 1:
                continue
            for p, rt in rcn.products.items():
                if p in basic_dict:  # Skip basic species
                    continue
                if p not in ylds:
                    ylds[p] = [np.zeros(rcn.rate.kvalues.shape), []]
                ylds[p][0] += rt * rcn.rate.kvalues
                ylds[p][1].append(ircn)

        if ylds:
            sps, vals = zip(*ylds.items())
            rts = np.array([v[0] for v in vals]).T
            ircns = set()
            for v in vals:
                ircns.update(v[1])
            norm_ylds = np.array([calculate_normalized_ratios(v) for v in rts]).T
            yields[info] = [sps, norm_ylds, ircns]

    return yields


def compute_yield_in_tree(reactions: list, ref_concs: dict, chd_info: dict = None, kuni: dict = None) -> dict:
    """Compute the estimated maximum product yields in the reaction tree."""
    # Get settings
    primary_vocs = get_attrs_from_settings({SNAME.SETUP: ["primary_vocs"]})

    kuni = compute_kuni(reactions, ref_concs)
    if chd_info is None:
        chd_info = get_species_relation(reactions, is_reversed=False)

    max_yld = {s: 1.0 for s in primary_vocs}
    rct_queue = deque(primary_vocs)  # Queue to check products

    while rct_queue:
        s = rct_queue.popleft()
        if s not in max_yld:
            raise ValueError(f"Species {s} not in max_yld")
        if not (s in chd_info and chd_info[s]["species"]):  # No products
            continue

        # Get max_yld in 1-step reactions
        ylds = {p: 0.0 for p in chd_info[s]["species"]}
        for ircn in chd_info[s]["ircns"]:
            if ircn not in kuni:  # Skip non-reducible reactions
                continue
            if s != kuni[ircn][0]:
                raise ValueError(f"Reactant {s} mismatch in kuni: {kuni[ircn]}")
            for p, rt in reactions[ircn].products.items():
                if p in ylds:
                    ylds[p] += kuni[ircn][1] * rt

        # Normalize
        pds, ylds = zip(*ylds.items())
        nylds = [calculate_normalized_ratios(v) for v in np.array(ylds).T]
        nyld_max = np.max(nylds, axis=0)  # Get max yield
        for i, p in enumerate(pds):
            yld = max_yld[s] * nyld_max[i]
            if p not in max_yld or yld > max_yld[p]:
                max_yld[p] = yld
                rct_queue.append(p)

    return max_yld


def compute_yield_with_rxntime(reactions: list, ref_concs: dict, rxntime: float = 21600) -> dict:
    """
    Compute the estimated product yields with reaction time in the environment:
    Reaction correction = 1 - exp(-kuni * rxntime)
    Yield of product = yield of reactant forming it * reaction correction * yield of product from reactant
    """
    # Get settings
    primary_vocs = get_attrs_from_settings({SNAME.SETUP: ["primary_vocs"]})

    kuni = compute_kuni(reactions, ref_concs)
    chd_info = get_species_relation(reactions, is_reversed=False)

    yields = {s: 1.0 for s in primary_vocs}
    rct_queue = deque(primary_vocs)  # Queue to check products

    while rct_queue:
        s = rct_queue.popleft()
        if not (s in chd_info and chd_info[s]["species"]):  # No products
            continue
        if s not in yields:
            raise ValueError(f"Species {s} not in yields")

        ylds = {p: -1 for p in chd_info[s]["species"]}

        for ircn in chd_info[s]["ircns"]:  # Update yields per reaction
            if ircn not in kuni:  # Skip non-reducible reactions
                continue
            if s != kuni[ircn][0]:
                raise ValueError(f"Reactant {s} mismatch in kuni: {kuni[ircn]}")
            yld = yields[s] * (1 - np.exp(-kuni[ircn][1] * rxntime))
            for p, rt in reactions[ircn].products.items():  # Update yields for products
                if p in ylds:
                    ylds[p] = max(ylds[p], yld * rt)

        for p, yld in ylds.items():
            if yld > 0:  # Update yield
                if p not in yields or yld > yields[p]:
                    yields[p] = yld
                    rct_queue.append(p)

    return yields


def compute_lifetime(reactions: list, ref_concs: dict, kuni: Optional[dict] = None) -> dict:
    """Compute lifetimes for reactants in the reaction tree."""

    # Get kuni order by reactants
    if not kuni:
        kuni = compute_kuni(reactions, ref_concs)
    kuni_sps = reshape_kuni_wsps(kuni)

    lifetime = {}
    for s, [_, ks] in kuni_sps.items():
        totals = np.sum(ks, axis=0)  # Sum over reactions for each kvalue
        lifetime[s] = np.where(totals == 0, np.inf, 1.0 / totals)
    return lifetime


def compute_bratio(reactions: list, ref_concs: dict, kuni: Optional[dict] = None) -> dict:
    """Compute branching ratios for 1-step reactions."""
    # Get kuni order by reactants
    if not kuni:
        kuni = compute_kuni(reactions, ref_concs)
    kuni_sps = reshape_kuni_wsps(kuni)

    bratio = {}
    for s, [ircns, ks] in kuni_sps.items():
        rts = np.array([calculate_normalized_ratios(v) for v in np.array(ks).T]).T
        bratio[s] = [rts, ircns]
    return bratio
