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
 This module contains functions to convert units between ug/m3, ppb, molec/cm3.

For Unit conversion: PV = nRT
   [ug/m3] = [molec/cm3] * MW * MULTIMASS
   [ppb] = [ug/m3] * R * T / (P * mw)
   [molec/cm3] = [ppb] * P * Na / (R * T)
"""

import numpy as np

# Avogadro's number in mol-1
N_AVO = 6.02214076e23
# Conversion factors
# Ideal gas constant in L* atm/(mol * K)
R_ATM = 0.082057
# 1E12/avogadro (6.022e23): factor involved for molec/cm3 to microg/m3
MULTIMASS = 1e12 / N_AVO
# MULTIMASS in GECKO
# molecules per cubic centimeter (molec/cm3) to micrograms per cubic meter (ug/m3)
STD_PRESS_PA = 101325.0  # Standard pressure in Pa


def ug_to_ppb(val: float, mw: float, temp: float, fpres: float) -> float:
    """
    Convert unit ug/m3 to unit ppb.
    Args:
     - val: value in ug/m3
     - mw: molecular weight in g/mol
     - temp: temperature in K
     - fpres: factor for pressure to standard pressure
    """
    return val * R_ATM * temp / (fpres * mw)


def ppb_to_ug(val: float, mw: float, temp: float, fpres: float) -> float:
    """Convert unit ppb to ug/m3"""
    return val * mw * fpres / (temp * R_ATM)


def ppb_to_molec(val: float, temp: float, fpres: float) -> float:
    """Convert unit ppb to molecule/cm3"""
    return val * fpres / (R_ATM * temp * MULTIMASS)


def molec_to_ppb(val: float, temp: float, fpres: float) -> float:
    """Convert unit molecule/cm3 to ppb"""
    return val * R_ATM * temp * MULTIMASS / fpres


def get_carbon_number(species: list, only_hc: bool = True) -> dict:
    """Get the number of carbon atoms in each species."""
    nc_dict = {s.name: s.fgroups.get("C", s.name.count("C")) for s in species}
    # Remove inorganic and methane hydrocarbons
    if only_hc:
        for s in ["CO", "CO2", "CH4", "XCLOST"]:
            if s in nc_dict:
                nc_dict[s] = 0
    return nc_dict


def get_unit_conv_factors(species, unit_in="ug", unit_out="molec", temp=298.0, nc_dict=None):
    """Calculate unit conversion factors for species. unit_in * fac = unit_out"""

    # No conversion needed
    if unit_in == unit_out:
        return {s.name: 1 for s in species}

    # Factor to pressure
    fpres = 1.0
    # Get basic factor based on the type of temperature: either np list or a float
    f0 = 1 if isinstance(temp, (int, float)) else np.ones(len(temp))

    # Get number of carbon atoms if needed
    if (unit_in == "ppbC" or unit_out == "ppbC") and nc_dict is None:
        nc_dict = get_carbon_number(species)

    # Define the conversion functions
    conversion_funcs = {
        "ug": {
            "ppb": lambda s: ug_to_ppb(f0, s.mass, temp, fpres),
            "ppbC": lambda s: ug_to_ppb(f0, s.mass, temp, fpres) * nc_dict[s.name],
            "molec": lambda s: f0 / (s.mass * MULTIMASS),
        },
        "ppb": {
            "ug": lambda s: ppb_to_ug(f0, s.mass, temp, fpres),
            "ppbC": lambda s: f0 * nc_dict[s.name],
            "molec": lambda s: ppb_to_molec(f0, temp, fpres),
        },
        "ppbC": {
            "ug": lambda s: (ppb_to_ug(f0, s.mass, temp, fpres) / nc_dict[s.name] if nc_dict[s.name] != 0 else 0),
            "ppb": lambda s: f0 / nc_dict[s.name] if nc_dict[s.name] != 0 else 0,
            "molec": lambda s: (ppb_to_molec(f0, temp, fpres) / nc_dict[s.name] if nc_dict[s.name] != 0 else 0),
        },
        "molec": {
            "ug": lambda s: s.mass * MULTIMASS,
            "ppb": lambda s: molec_to_ppb(f0, temp, fpres),
            "ppbC": lambda s: molec_to_ppb(f0, temp, fpres) * nc_dict[s.name],
        },
    }

    # Check if the conversion is supported
    if unit_in not in conversion_funcs or unit_out not in conversion_funcs[unit_in]:
        raise ValueError(f"Unsupported conversion from {unit_in} to {unit_out}.")

    # Perform the conversion
    unit_conv = {s.name: conversion_funcs[unit_in][unit_out](s) for s in species}

    return unit_conv
