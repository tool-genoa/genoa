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
This module contains functions for reduction strategies based on the reduction type and value.
"""

from copy import deepcopy
from itertools import combinations, product

from .jumping import find_reduction_via_jumping
from .logger import setup_logger, list_to_log
from .lumping import find_reduction_via_lumping, find_reduction_via_replacement, update_scheme_via_merging
from .reduction_basic import get_size_for_combination, ReductionListInfo
from .reduction_setting import N_LUMP_SIMPLE, NSP_EXPAND_LIM
from .reduction_tree import find_reduction_via_reaction_tree
from .reduction_print import get_reduction_pinfos
from .utils import calculate_normalized_ratios


# Logger
logger = setup_logger(__name__)


def reduction_on_mechanism(reactions: list, species: list, rdc_set: list, infos: ReductionListInfo) -> None:
    """Reduce the mechanism based on the specified reduction sets (reduction type and value)."""

    # Get settings
    rdc_type = rdc_set[0]
    logger.info("Conducting reduction based on <%s> with value <%s> ...", rdc_type, rdc_set[1])

    # Function list for complex reductions & code for updating infos
    complex_rdc_dict = {
        "lp": [find_reduction_via_lumping, "spckrl"],
        "rp": [find_reduction_via_replacement, "spc"],
        "jp": [find_reduction_via_jumping, "spck"],
    }

    # Simple reductions
    if rdc_type not in complex_rdc_dict:
        return find_reduction_via_reaction_tree(reactions, species, rdc_set, infos)

    # Update for reduction
    infos.update(reactions, species, {"renew": complex_rdc_dict[rdc_type][1]})
    logger.info("Updating SIMPLE species types for reduction ...")
    for s in species:
        if s.status > 0:
            s.update_for_reduction(is_simple=True)  # CAUTION: simple for threshold-based reduction

    # Complex reductions
    if rdc_type == "jp":
        # Search for jumping candidates & update the mechanism
        rdcs = complex_rdc_dict[rdc_type][0](reactions, species, infos, change_scheme=False)

    elif rdc_type in ["lp", "rp"]:
        # Search for reduction candidates
        rdcs = complex_rdc_dict[rdc_type][0](reactions, species, infos)
        # Update the mechanism based on the reduction sets
        update_scheme_via_merging(reactions, species, rdcs, infos.sps_info, rdc_type)

    else:
        raise ValueError(f"Unknown complex reduction type: {rdc_type}")

    # Record reduction information
    list_to_log(
        logger,
        get_reduction_pinfos(rdc_type, rdcs, reactions, species, infos),
        f"\nFind # {len(rdcs)} {rdc_type} candidates.",
    )

    return None


def find_reduction_via_strategy(
    reactions: list, species: list, strategy: str, nmcomb: int, nmrdc: int, infos: ReductionListInfo
) -> list:
    """Find reduction candidates on species or reactions based on specified reduction strategy and criteria."""

    rdc_func_dict = {
        "rs": find_reduction_via_removal,
        "rm": find_reduction_via_removal,
        "da": find_reduction_via_removal,
        "lp": find_reduction_via_lumping,
        "jp": find_reduction_via_jumping,
        "rp": find_reduction_via_replacement,
    }

    # Check if the strategy is valid
    if strategy not in rdc_func_dict:
        raise ValueError(f"Reduction strategy {strategy} not recognized.")

    if strategy in ["rs", "rm", "da"]:  # Removal
        targets = rdc_func_dict[strategy](strategy, species, infos)

    else:  # Lumping, jumping, replacement
        targets = rdc_func_dict[strategy](reactions, species, infos)

    # No reduction found
    if not targets:
        return None

    # Generate reduction sets based on found targets
    return generate_reduction_sets(strategy, targets, nmcomb, nmrdc, infos)


def find_reduction_via_removal(strategy: str, species: list, infos: ReductionListInfo) -> list:
    """Find reduction candidates via removal"""

    # Get checked items
    checked = infos.checked.get(strategy, [])

    # Remove species
    if strategy == "rs":
        return [s for s in infos.targeted_sps if s not in checked]

    # Remove reactions
    if strategy == "rm":
        return [s for s in infos.targeted_rcn if s not in checked]

    # Remove gas-particle partitioning
    if strategy == "da":
        return [s for s in infos.targeted_sps if (s not in checked and species[infos.sps_info[s]].condensable)]

    raise ValueError(f"removal strategy {strategy} not recognized.")


def generate_diverse_combinations(targets: list, nmax_comb: int, nmax_rdc: int) -> list:
    """Generate combinations of target objects."""

    n = len(targets)
    if nmax_comb <= 0 or n <= 1:  # No combinations
        return [targets]

    selected = []  # List of selected combinations

    # Return all combinations if n is small
    if n <= nmax_comb:
        for i in range(1, n + 1):
            selected.extend([list(c) for c in combinations(targets, i)])
        return selected

    # Default n+1 list if n is large
    selected = [targets] + [[c] for c in targets]  # Add single sets of targets

    # Return default
    if len(selected) >= nmax_rdc:
        return selected

    # Get possible sizes for combinations
    poss_size = get_size_for_combination(n)

    while len(selected) < nmax_rdc and poss_size:
        size = poss_size.pop(0)  # Pop the most diverse
        for c in combinations(targets, size):
            selected.append(list(c))
            if len(selected) >= nmax_rdc:
                break

    return selected


def _get_new_ratio_sets(species: list, infos: ReductionListInfo) -> list:
    """Get new ratio sets based on the species and ratios."""

    # Get new ratio based oo gas concentrations
    rts_gas = [infos.ref_concs["gas"].get(s, 0.0) for s in species]

    # Get new ratios based on reactivity: gas_concs/lifetime
    rts_flux = [infos.ref_concs["gas"].get(s, 0.0) / infos.ltimes.get(s, [1.0])[0] for s in species]

    # Get new ratios based on aerosol concentrations
    rts_aero = [infos.ref_concs["aero"].get(s, 0.0) for s in species]

    # Normalize ratios
    return [calculate_normalized_ratios(rts, 4) for rts in (rts_gas, rts_flux, rts_aero) if sum(rts)]


def generate_reduction_sets(strategy: str, targets: list, nmcomb: int, nmrdc: int, infos: ReductionListInfo) -> list:
    """Generate reductions based on the strategy and targets."""

    rdc_all = generate_diverse_combinations(targets, nmcomb, nmrdc)

    if len(rdc_all) >= nmrdc or strategy not in ("lp", "rp") or len(infos.sps_info) > N_LUMP_SIMPLE:
        return rdc_all

    # Check if w/ further expansion
    expand_sps = [s["species"] for tars in rdc_all for s in tars]
    nsps = max(len(s) for s in expand_sps)
    # No expansion if too many species
    if nsps > NSP_EXPAND_LIM:
        logger.info("Skip expansion w/ %s: too many species (%d) from %s", strategy, nsps, expand_sps)
        return rdc_all

    # Expand reduction sets
    logger.info("Expand reduction sets w/ %s: %s", strategy, expand_sps)
    if strategy == "rp":
        return expand_replacement_targets(rdc_all, nmrdc)

    if strategy == "lp":
        return expand_lumping_targets_wratios(rdc_all, nmrdc, infos)

    raise ValueError(f"Unknown reduction strategy: {strategy}")


def expand_replacement_targets(rdc_sets: list, nmax: int) -> list:
    """Expand replacement targets by replacing one species with another."""

    new_sets, nrdc = [], len(rdc_sets)
    for tars in rdc_sets:

        ntars = [len(s["species"]) - 1 for s in tars]
        if max(ntars) <= 1:  # No need to regroup
            continue

        # Get all species combinations for replacement
        sps_combs = _get_species_combs([s["species"][1:] for s in tars], 1)
        if not sps_combs:  # No combinations available
            continue

        for comb in product(*sps_combs):

            # Skip original
            if all(len(comb[i]) == ntar for i, ntar in enumerate(ntars)):
                continue

            new_group = []
            for i, tar in enumerate(tars):
                new_tar = deepcopy(tar)
                new_tar["species"] = [tar["species"][0]] + list(comb[i])
                new_tar["ratios"] = [1.0] + [0.0] * len(comb[i])
                new_group.append(new_tar)

            if new_group:
                new_sets.append(new_group)
                logger.info("Expand rp %s => %s", [s["species"] for s in tars], [s["species"] for s in new_group])
                nadd = len(new_sets)
                if nadd + nrdc >= nmax:
                    return rdc_sets + new_sets

    return rdc_sets + new_sets


def expand_lumping_targets_wratios(rdc_sets: list, nmax: int, infos: ReductionListInfo) -> list:
    """Expand lumping targets by applying ratio separation to each species group."""

    new_sets, nrdc = [], len(rdc_sets)
    for tars in rdc_sets:

        # Get original ratios and new ratios for each target
        org_rt = [tar["ratios"] for tar in tars]
        new_rts_all = [_get_new_ratio_sets(tar["species"], infos) for tar in tars]
        min_len = min(len(nrts) for nrts in new_rts_all)
        ntars = len(tars)

        for i in range(min_len):
            new_group = [new_rts_all[j][i] for j in range(ntars)]

            # Similar ratios
            if all(_is_two_ratios_similar(new_group[j], org_rt[j]) for j in range(ntars)):
                continue

            # Create new set w/ updated ratios
            new_set = [deepcopy(tar) for tar in tars]
            for j, rt in enumerate(new_group):
                new_set[j]["ratios"] = rt
            new_sets.append(new_set)
            logger.info("Expand lumping ratios for %s: %s (pre: %s)", new_set[0]["species"], new_group, org_rt)

            nadd = len(new_sets)
            if nadd + nrdc >= nmax:
                return rdc_sets + new_sets

    return rdc_sets + new_sets


def expand_lumping_targets_wspecies(rdc_sets: list, nmax: int) -> list:
    """Expand lumping targets by applying species number changes from each group"""

    new_sets, nrdc = [], len(rdc_sets)
    for tars in rdc_sets:

        # Check if there are enough species to relump
        ntars = [len(s["species"]) for s in tars]
        if max(ntars) <= 2:
            continue
        # Get all species combinations
        sps_combs = _get_species_combs([s["species"] for s in tars], 2)
        if not sps_combs:
            continue

        # Get species and ratios
        sps_rts = [dict(zip(s["species"], s["ratios"])) for s in tars]

        # Cross-product of combinations
        for comb in product(*sps_combs):

            # Skip original species by checking len
            if all(len(comb[i]) == ntar for i, ntar in enumerate(ntars)):
                continue

            new_group = []
            for i, tar in enumerate(tars):
                nrts, nsps = zip(*sorted(zip([sps_rts[i][s] for s in comb[i]], comb[i]), reverse=True))
                new_tar = deepcopy(tar)
                new_tar["species"] = list(nsps)
                new_tar["ratios"] = calculate_normalized_ratios(list(nrts), 4)
                new_group.append(new_tar)

            if new_group:
                new_sets.append(new_group)
                logger.info("Expand lumping species for %s: %s", new_group[0]["species"], new_group[0]["ratios"])

                nadd = len(new_sets)
                if nadd + nrdc >= nmax:
                    return rdc_sets + new_sets

    return rdc_sets + new_sets


def _get_species_combs(groups: list, nmin: int) -> list:
    """Get all combinations of species (n >2) from the groups."""

    species_combs = []
    with_comb = False
    for sps in groups:
        ns = len(sps)
        if ns <= nmin:
            species_combs.append([sps])
        else:
            with_comb = True
            species_combs.append([list(c) for i in range(nmin, ns + 1) for c in combinations(sps, i)])

    if not with_comb:
        return []

    return species_combs


def _is_two_ratios_similar(ratio1: list, ratio2: list, threshold: float = 0.01) -> bool:
    """Check if two normalized ratios are similar within a given threshold."""

    if len(ratio1) != len(ratio2):
        return False
    return all(abs(r1 - r2) < threshold for r1, r2 in zip(ratio1, ratio2))
