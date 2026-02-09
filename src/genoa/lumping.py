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
This module contains functions for finding reduction candidates based on lumping or replacement strategy.
"""

import math

from .mechanism_relation import get_downstream_species
from .logger import setup_logger
from .reduction_basic import is_two_similar, is_within_scale, ReductionListInfo
from .reduction_setting import rdc_checks, N_LUMP_SIMPLE
from .utils import calculate_normalized_ratios


# Logger
logger = setup_logger(__name__)


def find_reduction_via_replacement(reactions: list, species: list, infos: ReductionListInfo) -> list:
    """Find reduction candidates on species based on replacement strategy."""

    # Get checking for replacement
    rp_checks = rdc_checks["rp"]

    # Get reduction candidates
    reduction_info = []  # Record found reduction candidates
    changed_species_all = set()  # Record all species that have been changed

    for s in infos.targeted_sps:
        if len(reduction_info) >= rp_checks["ngmax"]:
            break  # Break if enough groups
        if s in changed_species_all:
            continue  # Skip if already included

        # Get parent reactions
        rcns = infos.parent_info.get(s, {}).get("ircns", [])
        isp = species[infos.sps_info[s]]

        # Find species for replacement
        changed_species = {s}
        for ircn in rcns:  # Find from products
            rcn = reactions[ircn]
            # Find speices can replace targeted
            for p in rcn.products:
                # Find species p that can replace s
                if p in infos.sps_info:
                    if p in changed_species_all or p in changed_species:
                        continue
                    jsp = species[infos.sps_info[p]]
                    if is_two_similar(isp, jsp, rp_checks):
                        changed_species.add(p)  # Add to replace set
                        if len(changed_species) >= rp_checks["nsmax"]:
                            break  # Break if enough species

        # Record replacement info
        if len(changed_species) > 1:
            rdc_info = {
                "species": [s] + [p for p in changed_species if p != s],
                "ratios": [1.0] + [0.0] * (len(changed_species) - 1),
                "rcns": get_related_reactions(changed_species, infos),
            }

            reduction_info.append(rdc_info)
            changed_species_all.update(changed_species)

    return reduction_info


def get_lumpable_list(spname: str, lumped: set, species: list, infos: ReductionListInfo, checks: dict) -> list:
    """Get and rank species to be checked for lumping."""

    # Get nearby species
    pinfo = get_downstream_species(spname, infos.parent_info)  # Get parent species
    cinfo = get_downstream_species(spname, infos.children_info)  # Get children species

    # Get targting species
    targeting = _get_sorted_targeting_sps(spname, infos)
    # logger.info("Found # %d species for %s.", len(targeting), spname)

    # Check for lumping
    lumpables = [spname]
    isp = species[infos.sps_info[spname]]

    for s in targeting:

        if s in lumped or s in lumpables or s not in infos.sps_info:
            continue
        if s in pinfo or s in cinfo:
            continue

        jsp = species[infos.sps_info[s]]
        if not is_two_similar(isp, jsp, checks):
            continue  # Skip if not similar

        # Check lifetime if needed
        if "tau" in checks:
            tau_out = False
            if spname in infos.ltimes and s in infos.ltimes:
                for t0, t1 in zip(infos.ltimes[spname], infos.ltimes[s]):
                    if not is_within_scale(t0, t1, checks["tau"]):
                        tau_out = True
                        break
            elif spname in infos.ltimes or s in infos.ltimes:
                continue  # Skip if only one species has lifetime
            if tau_out:  # Skip if lifetime is out of scale
                continue

        # Record a valid lumpable species
        lumpables.append(s)

        # Break if enough lumpable species
        if len(lumpables) >= checks["nsmax"]:
            break

        # Update parent and children info
        pinfo.update(get_downstream_species(s, infos.parent_info))
        cinfo.update(get_downstream_species(s, infos.children_info))

    return lumpables


def get_related_reactions(targeted: list, infos: ReductionListInfo) -> set:
    """Get related reactions for targeted species."""

    rcns = set()
    for s in targeted:
        if s in infos.children_info:
            rcns.update(infos.children_info[s]["ircns"])
        if s in infos.parent_info:
            rcns.update(infos.parent_info[s]["ircns"])
    return rcns


def get_lumping_group_and_ratio(lumpables: list, concs: dict, lrmin: float) -> list:
    """Get lumping group and ratio based on the lumpable species."""

    # Compute lumping ratios based on reference concentrations
    ratios = [concs[s] for s in lumpables]

    # Add weight for lifetime if needed
    # ratios = [ratios[i]/lifetime[s] for i, s in enumerate(lumpables)]

    # Normalize ratios
    ratios = calculate_normalized_ratios(ratios, iround=4, limit=lrmin)

    # Sort species based on ratios
    lumpables, ratios = zip(*sorted(zip(lumpables, ratios), key=lambda x: x[1], reverse=True))

    return [lumpables, ratios]


def _get_sorted_targeted_sps(infos: ReductionListInfo) -> list:
    """Output the targeted species for lumping."""

    # Get sorted species based on reference concentrations
    rconcs = infos.ref_concs["total"]
    return list(sorted(infos.targeted_sps, key=lambda s: rconcs[s], reverse=True))


def _get_sorted_targeting_sps(sname: str, infos: ReductionListInfo) -> list:
    """Output the targeting species for lumping."""

    # Settings
    rconcs = infos.ref_concs["total"]

    # Use given groups if available
    if infos.given_groups:

        if sname not in infos.given_groups:
            return []

        if isinstance(infos.given_groups, dict):

            # Get lumpable species from given groups
            igx = infos.given_groups[sname]
            ikey = "grps"

            # igx is group index w/ targeted species
            if ikey in infos.given_groups:
                igx_sps = infos.given_groups[ikey][igx]
                return [s for s in igx_sps if s in infos.sps_info]
            # igx is lumpable species list
            return [s for s in igx if s in infos.sps_info]

        sgroups = [s for s in infos.given_groups if s in rconcs]

    # Use targeted species as lumpable species
    elif len(infos.sps_info) > N_LUMP_SIMPLE:
        sgroups = [s for s in infos.targeted_sps if s in rconcs]

    # Use all species in reference concentrations
    else:
        sgroups = rconcs.keys()

    # Use reference concentrations to find targeting species
    rconcs = infos.ref_concs["total"]
    iconc = rconcs.get(sname, 0.0)
    return sorted(sgroups, key=lambda x: abs(rconcs[x] - iconc))


def find_reduction_via_lumping(reactions: list, species: list, infos: ReductionListInfo) -> list:
    """Find reduction candidates on species based on lumping strategy."""

    # Get checking for lumping
    lp_checks = rdc_checks["lp"]

    # Re-rank targeted and tartgeting species
    targeted = _get_sorted_targeted_sps(infos)

    # Start to search for lumping in order
    lumped_all = set()  # Store lumped species
    reduction_info = []  # Store lumping informations

    for s in targeted:

        if len(reduction_info) >= lp_checks["ngmax"]:
            break  # Break if enough groups

        if s in lumped_all:
            continue  # Skip if already lumped

        # Get all possible lumpable species
        lumped = get_lumpable_list(s, lumped_all, species, infos, lp_checks)
        if len(lumped) < 2:
            continue  # Break if no lumpable is found

        # Get lumping species in order and ratio
        lumped, ratios = get_lumping_group_and_ratio(lumped, infos.ref_concs["total"], lp_checks.get("lrtmin", 0.0))
        lumped_all.update(lumped)

        # Record lumping info
        rdc_info = {
            "species": lumped,  # Lumpable species names
            "ratios": ratios,
            "rcns": get_related_reactions(lumped, infos),  # Related reactions
        }

        reduction_info.append(rdc_info)

    return reduction_info


def update_scheme_via_merging(reactions: list, species: list, rdcs: list, sps_info: dict, rdc_type: str) -> None:
    """Update the mechanism for lumping or replacement based on the reduction sets."""

    for rdc in rdcs:
        # Prepare for merging
        merged_species = [species[sps_info[s]] for s in rdc["species"]]
        merged_ratios = rdc["ratios"]
        merged_rcns = [reactions[s] for s in rdc["rcns"]]  # List of reactions involved

        # Update species list
        update_species_via_merging(merged_species, merged_ratios, rdc_type)
        # Update reactions list
        update_reactions_via_merging(merged_rcns, merged_species, merged_ratios)


def update_reactions_via_merging(rcns: list, merged_species: list, merged_ratios: list) -> None:
    """Modify a reaction based on the specified merged parameters."""

    # Get species names
    snames = [s.name for s in merged_species]
    newname = snames[0]

    for rcn in rcns:
        if rcn.status != 1:
            continue  # Skip if inactive

        # Check reactants
        lump_rc = [s for s in rcn.reactants if s in snames]
        if lump_rc:

            # Get ratio for new kinetic rate
            rt = sum(merged_ratios[snames.index(s)] for s in lump_rc)

            # Update kinetic rate
            if not rcn.rate.multiply_rate_by_ratio(rt):
                rcn.status = 0
                continue

            # Change reactant names
            for s in lump_rc:
                if s == newname:
                    continue
                irt = rcn.reactants.pop(s)
                rcn.reactants[newname] += irt

            # Update reaction id
            rcn.update_rcn_id()

        # Check products
        lump_pd = [s for s in rcn.products if s in snames]
        if lump_pd:
            if lump_rc:  # Should not lumping reactants with products
                # logger.warning("Merged in reactant & product: %s & %s (rcn: %s)", lump_rc, lump_pd, rcn.to_rcn())
                rcn.clean_duplicates()
                # if set(lump_rc) != set(lump_pd):
                #    raise ValueError("Stop cuz lump_rc != lump_pd")

            # Change product names
            for s in lump_pd:
                if s == newname:
                    continue
                irt = rcn.products.pop(s)  # Find old ratio
                rcn.products[newname] += irt


def update_species_via_merging(merged_species: list, merged_ratios: list, mode: str) -> None:
    """Update species properties after lumping or replacement."""

    # Check mode
    if mode not in ["lp", "rp"]:
        raise ValueError(f"Unknown mode: {mode}")
    # Check species
    if len(merged_species) < 2:
        raise ValueError(f"Not enough species for merging. Got {merged_species}")

    # Get new species
    newsp = merged_species[0]

    # Update properties if needed
    if sum(merged_ratios[1:]) > 0.0:
        newsp.update_with_dict(get_merged_property(merged_species, merged_ratios))

    # Discard merged species & update reductions info
    if mode not in newsp.reductions:  # Initialize
        newsp.reductions[mode] = []
    for s in merged_species[1:]:
        s.status = 0  # Discard lumped species
        newsp.reductions[mode].append(s.name)


def get_merged_property(merged_species: list, merged_ratios: list) -> dict:
    """Get new species properties after merging."""

    # Initialize property dictionary for new species
    new = {"mass": 0.0, "fgroups": {}, "reductions": {}, "precursors": []}
    ncond = 0  # Count condensables

    # Calculate new properties
    for i, s in enumerate(merged_species):
        rt = merged_ratios[i]  # Get ratio
        new["mass"] += s.mass * rt  # Update mass

        for j in s.fgroups:  # Update functional groups
            if j not in new["fgroups"]:
                new["fgroups"][j] = 0.0
            new["fgroups"][j] += s.fgroups[j] * rt

        for k, v in s.reductions.items():  # Update reductions
            if k not in new["reductions"]:
                new["reductions"][k] = []
            for j in v:
                new["reductions"][k].append(j)

        for j in s.precursors:  # Update precursors
            if j not in new["precursors"]:
                new["precursors"].append(j)

        if s.condensable:  # Count condensables
            ncond += 1

    # Update properties for condensables
    if ncond:
        # Check if all species are condensables
        if ncond != len(merged_species):
            raise ValueError(f"Not all condensables: find # {ncond} out of # {len(merged_species)}")
        # Update condensable properties
        new.update(get_merged_condensed_property(merged_species, merged_ratios))

    return new


def get_merged_condensed_property(merged_species: list, merged_ratios: list) -> dict:
    """Get new species properties for condensables after merging."""

    # Initialize property dictionary for new condensed species
    log_psat, ihvap, istruc = 0.0, 0.0, {}

    # Calculate new properties
    for i, s in enumerate(merged_species):
        if not s.condensable:
            raise ValueError(f"Species {s.name} is not condensable")

        rt = merged_ratios[i]  # Get ratio
        log_psat += math.log10(s.psat_atm) * rt
        ihvap += s.dhvap_kj * rt

        for k, v in s.soap_strucs.items():  # Update soap structures
            if k not in istruc:
                istruc[k] = 0.0
            istruc[k] += v * rt

    return {"psat_atm": 10**log_psat, "dhvap_kj": ihvap, "soap_strucs": istruc}
