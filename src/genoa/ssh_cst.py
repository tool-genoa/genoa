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
This module contains constants related to SSH-aerosol box model.
"""

# === Namelist related constants ===

# Namelist variables and default values - not used, for information only
SSH_NML_ITEMS = {
    "initial_time": 0,  # Initial time in seconds
    "final_time": 0,  # Final time in seconds
    "delta_t": 0,  # Time step in seconds
    "latitude": 0,  # Latitude
    "longitude": 0,  # Longitude
    "meteo_file": "",  # Meteorological file
    "photolysis_file": "",  # Photolysis file
    "reaction_list_file": "",  # Reaction list file
    "species_list_file": "",  # Species list file
    "aerosol_species_list_file": "",  # Aerosol species list file
    "nout_soa": 0,  # Number of output SOA
    "RO2_list_file": "",  # RO2 list file
    "init_gas_conc_file": "",  # Initial gas phase concentrations
    "init_aero_conc_mass_file": "",  # Initial aerosol phase concentrations
    "cst_gas_file": "",  # Constant gas file
    "cst_aero_file": "",  # Constant aerosol file
    "init_species_file": "",  # Initial species file
    "err_species_list": "",  # Error species list
    "ref_conc_files_in": "",  # Reference concentrations files
    "output_directory": "",  # Output directory
    "output_type": 0,  # Output type
    "output_aero_list": "",  # Output aerosol list
    "output_gas_list": "",  # Output gas list
    # "particles_composition_file": "",  # Particles composition file
}

# Variables can be read from condition files
SCOND_ITEMS = {
    "latitude": "latitude",
    "longitude": "longitude",
    "meteo_file": "meteo_file",
    "init_gas": "init_gas_conc_file",
    "init_aero": "init_aero_conc_mass_file",
    "cst_gas": "cst_gas_file",
    "cst_aero": "cst_aero_file",
}

SCOND_HOURLY = {"init_gas_conc_file", "init_aero_conc_mass_file", "cst_gas_file", "cst_aero_file"}

# Variables & file suffixes can be read from mech folder
SMECH_SUFFIX = {
    "reaction_list_file": ".reactions",
    "species_list_file": ".species",
    "aerosol_species_list_file": ".aer.vec",
    "RO2_list_file": ".RO2",
}

# Variable name for the reference concentration list - update in constants.py
REF_FILES_SSH = "ref_conc_files_in"


# === Species related constants ===

# Total number of functional groups in SOAP
NGP_SOAP: int = 60

# SSH-aerosol functional group names
SSH_FGP_NAME = (
    ["RC"] * 12
    + ["RC"] * 4
    + ["C=C"] * 5
    + ["RC"] * 5
    + ["OH"] * 3
    + ["RCO"] * 3
    + ["RCOO"] * 2
    + ["COC"] * 3
    + ["COOH"]
    + ["NO3"] * 4
    + ["CO-OH"] * 3
    + ["CO-OC"] * 9
    + ["PAN"]
    + ["COOOH"]
    + ["O=COC=O"]
    + ["CHNO2"] * 3
)

# SSH-aerosol functional group index
SSH_FGP_INDEX = {
    "RC": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22, 23, 24, 25],
    "C=C": [16, 17, 18, 19, 20],
    "OH": [26, 27, 28],
    "RCO": [29, 30, 31],
    "RCOO": [32, 33],
    "COC": [34, 35, 36],
    "COOH": [37],
    "NO3": [38, 39, 40, 41],
    "CO-OH": [42, 43, 44],
    "CO-OC": [45, 46, 47, 48, 49, 50, 51, 52, 53],
    "PAN": [54],
    "COOOH": [55],
    "O=COC=O": [56],
    "CHNO2": [57, 58, 59],
}

# SSH-aerosol functional group mass
SSH_FGP_MASS = [
    15.0,
    14.0,
    13.0,
    12.0,
    15.0,
    14.0,
    13.0,
    12.0,
    15.0,
    14.0,
    13.0,
    12.0,
    15.0,
    14.0,
    13.0,
    12.0,
    27.0,
    26.0,
    26.0,
    25.0,
    24.0,
    13.0,
    12.0,
    27.0,
    26.0,
    25.0,
    17.0,
    18.0,
    29.0,
    43.0,
    42.0,
    29.0,
    59.0,
    58.0,
    31.0,
    30.0,
    29.0,
    45.0,
    58.0,
    76.0,
    75.0,
    74.0,
    47.0,
    46.0,
    45.0,
    61.0,
    60.0,
    59.0,
    60.0,
    59.0,
    58.0,
    58.0,
    57.0,
    56.0,
    106.0,
    61.0,
    72.0,
    61.0,
    60.0,
    59.0,
]
