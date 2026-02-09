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
This module contains the constants used in multiuse modules.
"""

import logging


# === Writing & Printing related constants === #

NOUT_LIM = 20  # Maximum number to output from a list

# Line separators for printing
LINE_SEP = "#" + "=" * 40 + "#\n"
LINE_SEP2 = "%" + "-" * 40 + "%\n"

# Precision for formatting stoichiometric coefficients
PRC_ST: str = ".5f"
# Precision for formatting rate constants
PRC_FL: str = "6.3E"


# === Reduction related constants ===

# All reduction strategies can be used in training
RDC_SGY = ["rm", "rm1", "jp", "lp", "rp", "rs", "da"]

# Mode for training order generation
# gprop: based on species properties; gdel: based on species deletion
TRN_MODE = {"gen_prop": "gprop", "gen_delete": "gdel", "gen_delete_O": "gdel"}


# === Species & Reaction Processing related constants ===

# Solar zenith angles
SZAS = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 78.0, 86.0, 90.0]
# Number of solar zenith angles
NSZA: int = len(SZAS)

# Reactant dictionary - species name: molar mass
RCT_DICT = {"HO": 17, "OH": 17, "O3": 48, "NO3": 62, "NO": 30, "NO2": 46, "HO2": 33, "SO2": 64, "CO": 28}

# No rate string
NORATE_STR = "!!!!"


# === Post-processing related constants ===

# Concentration phase
PHASES = ["gas", "aero", "total"]

# Available units
UNITS_DICT = {"ug": "ug/m3", "ppb": "ppb", "molec": "molec/cm3", "ppbC": "ppbC"}

# Ratios
RATIOS = ["OM/OC", "O/C", "H/C", "N/C"]

# Psat_atm upper limits for SVOCs, LVOCs, ELVOCs
PSAT_BIN_LIMS = {"SVOC": 1e-6, "LVOC": 1e-9, "ELVOC": 1e-13, "NVOC": 0}

NTOP_MAX: int = 10  # Maximum number of top species to display

EPS: float = 1e-30  # Epsilon to avoid division by zero


# === Configuration file section names === #
class SNAME:
    """Section name in the configuration file."""

    # Read activated actions
    ACN = "action"

    # General settings
    SETUP = "general"

    # Other settings
    NEW = "build_mechanism"
    SML = "run_simulation"
    PST = "post-processing"
    TBR = "threshold_reduction"
    TRN = "training"
    TST = "testing"

    GCK = "GECKO-A"
    SSH = "SSH-aerosol"
    ENV = "environment"

    # Available action with settings
    ACNS = [NEW, SML, TRN, TST, TBR, PST]
    # Read from configuration file by default
    READ = [SML]
    # Box models
    BOX = {"SSH": SSH, "GECKO": GCK}


# === Logging settings === #
class LOG:
    """Settings for logging."""

    FLAG: int = 3  # 0: no record; 1: record to console; 2: record to file; 3: record to both
    FILE: str = "./runtime/logs/log"
    LV_LOG: int = logging.INFO  # Default logging level for the logger (INFO, DEBUG, WARNING, ERROR, CRITICAL)
    LV_CON: int = logging.INFO  # Default logging level for the console
    LV_FIL: int = logging.WARNING  # Default logging level for the file
