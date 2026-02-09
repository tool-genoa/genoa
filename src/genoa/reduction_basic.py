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
This module contains basic functions for species and reaction reduction.
"""

import math
import random
from typing import Optional

from attrs import define, field

from .kinetics import compute_kuni, compute_lifetime
from .logger import setup_logger, list_to_log
from .mechanism_relation import get_species_relation
from .simulation_refconc import get_concs_for_species


# Logger
logger = setup_logger(__name__)


@define(slots=True)
class ReductionListInfo:
    """Information required for reduction."""

    # Targeted species
    targeted_sps: list = field(factory=list)
    # Targeted reactions
    targeted_rcn: list = field(factory=list)

    # Frozen species
    frozen: set = field(factory=set)
    # Checked reduction elements per strategy
    checked: dict = field(factory=dict)

    # Information used for reduction
    # Species information
    sps_info: dict = field(factory=dict)
    # Reaction information: parent
    parent_info: dict = field(factory=dict)
    # Reaction information: children
    children_info: dict = field(factory=dict)

    # Reference concentrations
    ref_concs: dict = field(factory=dict)
    # Pseudo-unimolecular rate constant
    kunis: dict = field(factory=dict)
    # Lifetime
    ltimes: dict = field(factory=dict)

    # Targeting groups of species for reduction
    given_groups: Optional[dict] = None

    def update(self, reactions: list, species: list, inputs: dict) -> None:
        """
        Update the reduction information with the species and reactions, following the inputs dictionary:
        {"frozen": set/list, "refconcs": dict, "renew": str, "given_groups": list/dict}
        """
        # Get flag for renewing the information
        renew = inputs.get("renew", "spckrl")  # Deault renew all

        # Update frozen species
        if inputs.get("frozen", None):
            self.frozen.update(inputs["frozen"])

        # Update given groups
        if inputs.get("given_groups", None):
            self.given_groups = inputs["given_groups"]

        # Update species information
        if "s" in renew:
            self.sps_info = {s.name: i for i, s in enumerate(species) if s.status == 1}

        # Update targeted species
        if not self.targeted_sps:
            self.init_targeted_species(species)
        else:
            self.update_targeted_species()

        # Update reaction information
        if "p" in renew:
            self.parent_info = get_species_relation(reactions, is_reversed=True)

        # Not need for replacement
        if "c" in renew:
            self.children_info = get_species_relation(reactions, is_reversed=False)

        # Update reference concentrations
        if inputs.get("refconcs", None) or "r" in renew:
            ref_concs = inputs.get("refconcs", self.ref_concs)
            for k in ["gas", "aero", "total"]:
                ref_concs[k] = get_concs_for_species(ref_concs[k], species)
            self.ref_concs = ref_concs

        # Update pseudo-unimolecular rate constants
        if "k" in renew:
            self.kunis = compute_kuni(reactions, self.ref_concs["gas"], update_rcn=True)

        # Update lifetimes
        if "l" in renew:
            self.ltimes = compute_lifetime(reactions, self.ref_concs["gas"], self.kunis)

    def init_targeted_species(self, species: list, add_rule: bool = False) -> None:
        """Get targeted species for reduction."""
        pinfos = ["Initializing targeted species for reduction ..."]

        # Get all reducible species in reverse order of the species list
        targeted = [s.name for s in reversed(species) if s.status == 1 and s.name not in self.frozen]

        if not add_rule:
            # pinfos.append(f"INFOS: {len(targeted)} species are targeted for reduction.")
            self.targeted_sps = targeted
            return

        # Add rules for targeted species
        gen_max = 10  # Maximum generation number
        pinfos.append(f"Add rules for targeted species with generation number <= {gen_max}.")

        # Get new targeted species
        new_targeted = []
        for s in targeted:
            isp = species[self.sps_info[s]]
            if isp.generation > gen_max:
                continue
            new_targeted.append(s)

        # Update targeted species
        pinfos.append(f"Update {len(new_targeted)} targeted species from {len(targeted)} species")
        self.targeted_sps = new_targeted

        list_to_log(logger, self.targeted_sps, "\n".join(pinfos))

    def update_targeted_species(self) -> None:
        """Remove invalid species from the targeted species."""

        self.targeted_sps = [s for s in self.targeted_sps if s in self.sps_info and s not in self.frozen]


def is_within_scale(v1: float, v2: float, factor: float = 2.0) -> bool:
    """Check if the ratio of two values (v1 and v2) is within a specified scale factor."""

    # Handle zero values
    if v1 == 0.0 or v2 == 0.0 or factor == 0.0 or factor == 1.0:
        return v1 == v2
    # Check if the ratio is within the scale
    s1, s2 = (factor, 1.0 / factor) if factor < 1.0 else (1.0 / factor, factor)
    return s1 <= v1 / v2 <= s2


def is_two_similar(x, y, lims: dict) -> bool:
    """Determine if a target species is similar to a single peer species based on certain criteria."""

    # Comparing species type
    if "stype" in lims and x.rdc_id != y.rdc_id:
        return False

    # Comparing formula
    if "formula" in lims and x.formula != y.formula:
        return False

    # Comparing GECKO type
    if "note" in lims and x.note != y.note:
        return False

    # Comparing generation number - limit is a float
    if "gen" in lims and abs(x.generation - y.generation) > lims["gen"]:
        return False

    # Comparing fractional mass differences - limit is a float
    if "frc_mass" in lims and abs(x.mass - y.mass) / (x.mass + y.mass) * 2 > lims["frc_mass"]:
        return False

    # Comparing absolute mass differences - limit is a float
    if "abs_mass" in lims and abs(x.mass - y.mass) > lims["abs_mass"]:
        return False

    # Check functional groups - limits is a dictionary with keys as functional groups and values as limits
    if "fgroup" in lims and _is_dict_diff(x.fgroups, y.fgroups, lims["fgroup"]):
        return False

    # Check ratios - limits is a dictionary
    if "ratios" in lims and _is_dict_diff(x.ratios, y.ratios, lims["ratios"]):
        return False

    # Check condensable properties
    if x.condensable and y.condensable:
        # Saturation pressure - limits is a float
        if "psat" in lims and abs(math.log10(x.psat_atm) - math.log10(y.psat_atm)) > lims["psat"]:
            return False
        # SOAP structure - limits is a dictionary
        if "soap_strucs" in lims and x.soap_strucs and y.soap_strucs:
            if _is_dict_diff(x.soap_strucs, y.soap_strucs, lims["ratios"]):
                return False

    # Return True if all criteria are met
    return True


def _is_dict_diff(dict1: dict, dict2: dict, limits: dict) -> bool:
    """Check if the differences between two dictionaries are within specified limits."""

    for k, lim in limits.items():
        v1 = dict1.get(k, None)
        v2 = dict2.get(k, None)
        if v1 is None and v2 is None:  # No key in both
            continue
        if v1 is None or v2 is None:  # Key is missing in one
            return True
        if abs(v1 - v2) > lim:  # Limit exceeded
            return True

    return False  # No differences found within limits


def is_all_similar(target, peers: list, checks: dict) -> bool:
    """Determine if a target species is similar to a group of peer species"""

    for peer in peers:
        if not is_two_similar(target, peer, checks):
            return False
    return True


def order_size_with_log(n: int) -> list:
    """Generate combination sizes using logarithmic spacing (in the range 2 to n-1)."""
    return sorted({max(2, math.ceil(n / (2**i))) for i in range(1, int(math.log(n, 2)) + 1)})


def order_size_via_half(n: int) -> list:
    """Generate combination sizes by halving the size progressively (in the range 2 to n-1)."""
    sizes = []
    current_size = n // 2
    while current_size >= 2:
        sizes.append(current_size)
        current_size = current_size // 2
    return sorted(sizes)


def order_size_random(n: int) -> list:
    """Generate combination sizes randomly (in the range 2 to n-1)."""
    sizes = list(range(2, n))
    random.shuffle(sizes)  # Randomize the order of sizes
    return sizes


def order_size_weighted(n: int) -> list:
    """Generate combination sizes with a bias towards middle sizes (in the range 2 to n-1)."""
    sizes = list(range(2, n))
    sizes.sort(key=lambda x: abs(x - n / 2))  # Prioritize sizes closer to the middle
    return sizes


def get_size_for_combination(n: int, mode: str = "log"):
    """Generate size for combination in rnage (2, n-1)"""

    psize_dict = {
        "log": order_size_with_log,
        "half": order_size_via_half,
        "random": order_size_random,
        "weighted": order_size_weighted,
    }
    if mode in psize_dict:
        return psize_dict[mode](n)

    # Reverse order as default
    return list(range(2, n + 1))[::-1]
