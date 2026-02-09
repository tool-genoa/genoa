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
This module contains the settings for the reduction strategies.
The user can modify here the settings to customize the reduction strategies.
"""

from typing import Any, Dict, Optional

import attrs

from .logger import setup_logger


# Logger
logger = setup_logger(__name__)


# Number of species over which simple lumping is applied
N_LUMP_SIMPLE = 4000

# Number of species over which reduction expansion is NOT applied
NSP_EXPAND_LIM = 5

# Threshold for building reduction possibilities during training
rdc_variant_build = {
    "ncomb": 0,  # Maximum number of reduction candidate sets: 3 with 7 combinations
    "nrdc": 0,  # Maximum number of reduction candidate combinations
}

# Order to run threshold-based strategies
tbr_order = [
    "lifetime",
    "generation",
    "nvoc",
    "nvoc_tree",
    "svoc",
    "svoc_tree",
    "pyield",
    "jumping",
    "conc_tree",
    "pyield_tree",
    "bratio",
    "extra_sps_rdc",
    "lumping",
    "replacement",
]

# Config for training order if generated from species properties
group_order_build = {
    "nsmax": 0,  # Max number of species in a group if > 0
    "group_types": ["gen", "note"],  # Keywrods for grouping similar species
    "order_types": ["gen", "mass", "note"],  # Keywrods for ordering groups
}

# Settings for lumping to find lumpable pairs. See class FindPairOption for details and default values.
lumping_settings = {
    # "nsmax": 2,
    # "lrtmin": 0.0,
    "tau": 10,  # 5,
    # "gen": 1,
    # "formula": 1,
    "stype": 1,
    "note": 1,
    "fgroup": {"C": 2, "O": 3},  # Functional groups
    # "C": 2,
    # "N": 1,
    # "O": 3,  # 1
    # "H": None,
    "frc_mass": 0.1,
    "abs_mass": 20.0,
    # "ratios": {"O/C": 0.2, "N/C": 0.2, "OM/OC": 0.2},  # Ratios
    # "O/C": 0.1,
    # "N/C": 0.1,
    # "OM/OC": 0.1,
    "psat": 2.0,
    # "sgroup": 0.0,
}

# Settings for replacement
replacement_settings = {
    # "formula": 1,
    # "stype": 1,
    "frc_mass": 0.5,
    # "abs_mass": 100.0,
    "fgroup": {"C": 3},
    # "C": 3,  # 2,
    # "O/C": 0.2,
    # "N/C": 0.2,
}

# Settings for jumping
# Available keys: "ngmax", "nsmax", "lrtmin" (ratio threshold for negligible children species)
jumping_settings = {}


@attrs.define(slots=True)
class FindPairOption:
    """Settings for finding reductions via lumping & replacement."""

    # Max number to build reduction groups (and species in a group)
    ngmax: Optional[int] = 1e9  # Maximum number of merging groups, None for no limit
    nsmax: Optional[int] = 1e9  # Maximum number of species in a merging group, None for no limit

    # Check strings - in is_two_similar()
    stype: Optional[int] = None  # species type1 == species type2
    formula: Optional[int] = None  # formula1 == formula2 in CHNO
    note: Optional[int] = None  # GECKO type1 == GECKO type2

    # Check differences w/ input float - in is_two_similar()
    gen: Optional[int] = None  # Generation/Node depth: abs(gen1-gen2) <= iGen, None for no check
    frc_mass: Optional[float] = None  # abs(mass1-mass2)/(mass1+mass2) <= frc_mass
    abs_mass: Optional[float] = None  # abs(mass1-mass2) <= abs_mass
    psat: Optional[float] = None  # abs(log10(psat1)-log10(psat2)) <= psat - for condensable species ONLY

    # Check differences w/ input dictionary - in is_two_similar()
    # Available keys: "C", "N", "O", "H"
    fgroup: Optional[dict] = None  # Functional groups: abs(fgroup1-fgroup2) <= fgroup
    # Available keys: "OM/OC", "O/C", "N/C", "H/C"
    ratios: Optional[dict] = None  # Ratios: abs(ratio1-ratio2) <= ratios
    soap_strucs: Optional[dict] = None  # SOAP structure: abs(soap1-soap2) <= soap_strucs

    # Check in lumping
    tau: Optional[float] = None  # Lifetime threshold

    # Check in lumping & jumping
    lrtmin: Optional[float] = None  # Minimum lumping ratio | jumping threshold for negligible children species


def _get_reduction_settings(settings: Dict[str, Any]) -> Dict[str, Any]:
    """Return the settings for the reduction strategies."""
    # Get default options
    options = attrs.asdict(FindPairOption())

    # Update the default options with the given settings
    for k, v in settings.items():
        if k not in options:
            raise ValueError(f"Invalid option '{k}' in the settings.")
        options[k] = v

    # Remove none values
    options = {k: v for k, v in options.items() if v is not None}

    # Check types
    for k, v in options.items():
        if k in ["fgroup", "ratios", "soap_strucs"]:
            if not isinstance(options[k], dict):
                raise TypeError(f"Option '{k}' must be a dictionary.")
            vals = options[k].values()
            for v in vals:
                if not isinstance(v, (int, float)):
                    raise TypeError(f"Values in option '{k}' must be numeric. Got {vals}.")
        elif not isinstance(v, (int, float)):
            raise TypeError(f"Option '{k}' must be numeric. Got {v}.")

    return options


def log_reduction_settings(checks: dict) -> str:
    """Record the reduction settings."""

    pinfos = ["Reduction settings:"]
    pinfos.append(f"  Lumping: {checks['lp']}")
    pinfos.append(f"  Replacement: {checks['rp']}")
    pinfos.append(f"  Jumping: {checks['jp']}")
    return "\n".join(pinfos)


# Get settings for using in the reduction strategies
rdc_checks = {
    "lp": _get_reduction_settings(lumping_settings),
    "rp": _get_reduction_settings(replacement_settings),
    "jp": _get_reduction_settings(jumping_settings),
}

logger.debug(log_reduction_settings(rdc_checks))
