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
This module contains the function to convert mechanism from MechGen to SMILES format.
"""

from .logger import setup_logger
from .smiles import string_to_smiles
from .utils import isfloat


# Logger
logger = setup_logger(__name__)

# Convert functional groups (this might need to be expanded for all functional groups)
MECHGEN_TO_SMILES = {
    # C
    "aC": "c",  # aromatic
    "pC": "C",  # allylic
    # H (ignored)
    "H4": "",
    "H3": "",
    "H2": "",
    "H": "",
    # Single bonds
    "-": "",
    "-^": "/",  # Cis and trans isomerizatio
    "-v": "\\",
    # O
    "OH": "O",
    # 'CO' may overlap with COO!!!
    "CO": "C(=O)",
    "C1O": "C1(=O)",
    "C2O": "C2(=O)",
    # C=O
    "CHO": "C(=O)",
    # NO2
    "NO2": "N(=O)=O",
    # radicals
    "C[.]": "[C]",
    "CH[.]": "[CH]",
    "CO[.]": "C[O]",
    "[O.]": "[O]",
    "[OO.]": "O[O]",
    "[OO]": "O[O]",
    # Syn/Anti isomerization not used
}

# Add user-defined SMILES for specific MechGen groups HERE
MECHGEN_SMILES = {}

# Rank mechgen groups from the longest to shortest
_MECHGEN_GRP = sorted(MECHGEN_TO_SMILES, key=len, reverse=True)


def mechgen_to_smiles(mechgen: str, tag_canonical: bool = True) -> str:
    """Convert input mechgen format to its canonical (optional) smiles structure"""

    mechgen = mechgen.strip()

    if not mechgen:
        logger.WARNING("Empty input string.")
        return ""

    # Check if mechgen with specific SMILES
    if mechgen in MECHGEN_SMILES:
        return MECHGEN_SMILES[mechgen]

    # Check if contain non-smiles structure
    for i in ["syn", "anti", "{"]:  # Some products may contain "{...}"
        if i in mechgen:
            raise NameError("Contain non-smiles structure: ", mechgen)

    # Replace mechgen groups step by step
    for i in _MECHGEN_GRP:
        if i in mechgen:
            mechgen = mechgen.replace(i, "+" + MECHGEN_TO_SMILES[i] + "+")
    mechgen = mechgen.replace("+", "")

    # Add number 1 after C for monocyclic ("*" MUST followed by a number)
    mechgen = _add_monocyclic_number(mechgen)

    # Remove * if contains any
    mechgen = mechgen.replace("*", "")

    # Check string and convert to canonical smiles if need
    return string_to_smiles(mechgen, tag_canonical)


def _add_monocyclic_number(mg_in: str) -> str:
    """
    Add ring closure number '1' for monocyclic structures indicated by '*' without a number.
    The number is inserted after the nearest preceding carbon atom ('C').
    """

    # Find all "*" w/o number - monocyclic
    star_inds = [i for i in range(len(mg_in) - 1) if (mg_in[i] == "*" and not isfloat(mg_in[i + 1]))]

    if not star_inds:  # No modification needed
        return mg_in

    updated, i_cut = [], 0  # Record new string

    for ind in star_inds:

        # Get substring
        subline = mg_in[i_cut:ind]

        # Find the last 'C'
        cind = subline.rfind("C")
        if cind == -1:
            raise ValueError(f"No C found before * in subline: '{subline}' from input: {mg_in}")

        # Insert '1' after the 'C'
        updated.append(subline[: cind + 1] + "1" + subline[cind + 1 :])

        i_cut = ind + 1  # Update index for next subline

    # Append the remaining part of the string after the last '*'
    updated.append(mg_in[i_cut:])

    # New string
    new_mg = "".join(updated).strip()
    # logger.info("Process new smiles for monocyclics: %s", new_mg")
    return new_mg
