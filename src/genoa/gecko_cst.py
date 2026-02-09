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
This module contains constants related to GECKO-A generator and box model.
"""

# === Reaction & Kinetic related constants ===

# Max number of lines for a reaction
MAX_LINES = 8
# Max length of a line for a reaction
MAX_LINE_LEN = 120

# Default arrhenius parameters
GARR3 = ["1.0E+00", "0.0", "0"]
# EXTRA keyword for extra reactions
EXTRA_LABEL = {100, 200, 500, 501, 502, 550}
EXTRA_KEYS = {f"EXTRA {i}": i for i in EXTRA_LABEL}

# Peroxy radical group keyword
GRO2_GRP = ["PERO1", "PERO2", "PERO3", "PERO4", "PERO5", "MEPERO", "PERO7", "PERO8", "PERO9"]
# Number of peroxy groups
NPERO = len(GRO2_GRP)

# Keyword to add extra reactant name
GRCT_DICT = {"ARR": None, "FALLOFF": "FALLOFF", "PHOT": "HV", "EXTRA 200": "ISOM"}

# Keyword conversion GECKO -> SSH
GCK_TO_SSH_KEYW = {
    "HV": "PHOT",
    "TBODY": "TB M",
    "PERO": "RO2",  # MEPERO? 01-09
    "EXTRA": "EXTRA",
    "OXYGEN": "TB O2",
    "FALLOFF": "FALLOFF",
    "ISOM": "EXTRA 200",
}

# All kinetic types and codes
GTYPES = list(GCK_TO_SSH_KEYW.keys())
GTYPES.extend(GRO2_GRP)
# notes = ['HABS', 'MJ19', 'RA07', 'OLDH', 'IUPA', 'LV09', 'RXEX']

# === Species related constants ===

# Length of species name
LENSP = 8
# Default number of speces to separate kinetic info
NSPACE = " " * 20

# Convert functional groups (this might need to be expanded for all functional groups)
GCK_TO_SMILE = {
    # C
    "Cd": "C",
    "c": "c",
    # H
    "H4": "",
    "H3": "",
    "H2": "",
    "H": "",
    # C=O
    "CHO": "C(=O)",
    # Warning: "CO" may overlap with "COO"
    "CO": "C(=O)",
    "C1O": "C1(=O)",
    "C2O": "C2(=O)",
    # NO2
    # "NO2": "[N+](=O)[O-]",
    "NO2": "N(=O)=O",
    # O
    "OH": "O",
    # Crieege radicals
    "EOO.": "OO.",
    "ZOO.": "OO.",
    # Handwritten
    "#": "",
    "#mm": "",
    "-": "",
}

# GECKO-A notes to functional groups
GCODE_DICT = {
    "A": "-CO(OH)",
    "B": "-Br",
    "D": "-CHO",
    "E": "R-O-R",
    "F": "-F",
    "G": "-CO(OOH)",
    "H": "-OOH",
    "K": ">CO",
    "L": "-Cl",
    "N": "-ONO2",
    "O": "-OH",
    "P": "-CO(OONO2)",
    "R": "aromatic rings",
    "T": "aliphatic rings",
    "U": ">C=C<",
    "V": "-NO2",
    "X": "-C=C=O",
    "1": "->C(O.)",
    "2": "->C(OO.)",
    "3": "-CO(OO.)",
    "4": ">C.(OO.)",
    "S": "-OSO3",
}

# === Namelist related constants ===

# Namelist variables and default values - not used, for information only
GCK_NML_ITEMS = {
    "fg_out_col": 1,  # If > 0, output "output.genoa" file
    "fg_output": 0,  # 0: no print and no normal output; 1: print + bof; 2: all outputs
    "fg_interp": 0,  # If = 0, no product ratio info
    "err_sps_list_in": "",  # List of species for error calculation
    "ref_conc_list_in": "",  # List of files for reference concentrations
    "init_conc_list_in": "",  # List of species with specified initial concentrations
    "input_dir_chem": "",  # Input directory to read chem files
    "input_cond": "",  # Input directory to read initial condition file: [].key
    "output_dir": "",  # Output directory to write results
    "phot_file": "",  # Photolysis file
}

# Variable name for the reference concentration list - update in constants.py
REF_FILES_GCK = "ref_conc_list_in"
