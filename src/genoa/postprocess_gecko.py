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
This module contains functions to read GECKO-A results from binary files.
"""

import os
import struct
from typing import BinaryIO

import numpy as np

from .logger import setup_logger
from .utils import get_average_concs, get_data_based_on_mode, get_files_from_folder, replace_species_in_conc
from .unit_conversion import MULTIMASS, molec_to_ppb


# Logger
logger = setup_logger(__name__)


def get_gecko_concs(folderpath: str, conc_opts: dict) -> dict:
    """
    Get reference concentrations from multiple GECKO-A simulation results.
    Output a dict with two keys: "gas" and "aero" for {species: concentrations}
    """

    # Get file paths
    filepaths = get_files_from_folder(folderpath, ".bof")
    if not filepaths:
        raise ValueError(f"No .bof files in {folderpath}. Got: {folderpath}")

    # Read concs from files
    logger.info("Read concentrations in GECKO-A results from %d files ...", len(filepaths))
    concs_all = []
    for filepath in filepaths:
        # Extract raw results
        raw = _extract_gecko_bof(filepath)
        # Obtain processed results
        concs = _process_gecko_bof(raw, conc_opts)
        concs_all.append(concs)

    return get_average_concs(concs_all)


def _get_species_names(bin_obj: BinaryIO, numsp: int, lensp: int) -> list:
    """Read species names from the GECKO-A binary file."""
    species_names = []
    bin_obj.read(lensp)  # Skip first lensp bytes
    for _ in range(numsp):
        species_name = bin_obj.read(lensp).decode("latin-1").strip()
        species_names.append(species_name)
    return species_names


def _get_molar_mass(bin_obj: BinaryIO, numsp: int, lensp: int) -> np.ndarray:
    """Read molar mass from the GECKO-A binary file."""
    bin_obj.read(lensp)  # Skip first lensp bytes
    return np.frombuffer(bin_obj.read(numsp * 8), dtype=np.float64)


def _get_time_temp_conc(bin_obj: BinaryIO, numsp: int, ndat: int) -> list:
    """Read time, temperature, and concentration data from the GECKO-A binary file."""

    # Initialization
    time = np.zeros(ndat)
    temp = np.zeros(ndat)
    conc = np.zeros((ndat, numsp))

    # Read data
    for i in range(ndat):
        _, time[i], temp[i] = struct.unpack("3d", bin_obj.read(24))  # Double precision
        conc[i, :] = np.frombuffer(bin_obj.read(numsp * 8), dtype=np.float64)
    return [time, temp, conc]


def _get_gecko_unit_conv(raw_result: dict, unit: str) -> float:
    """Get unit conversion factor for molecule/cm3 to the desired unit."""

    if unit == "molec":  # No conversion
        return np.ones(len(raw_result["name"]))

    if unit == "ug":
        return raw_result["mw"] * MULTIMASS

    if unit == "ppb":
        fpres = 1.0  # Factor to standard pressure
        if np.all(raw_result["temp"] == raw_result["temp"][0]):  # Single temperature
            return molec_to_ppb(np.ones(len(raw_result["name"])), raw_result["temp"][0], fpres)

        # Multiple temperatures
        fac = np.array([molec_to_ppb(1, t, fpres) for t in raw_result["temp"]])
        return np.tile(fac, (len(raw_result["name"]), 1))

    raise ValueError(f"Unknown unit: {unit}")


def _process_gecko_bof(raw_result: dict, options: dict) -> dict:
    """
    Process the gecko-a binary results and return concentration dictionary.
    Available options:
    - winit: True/False, whether with initial concentrations (at step 0)
    - unit: molec, ug, ppb. Not support ppbC.
    - time: all/ave/max/last, time points to be included
    - use_fake: True/False, whether to use fake concentrations
    - total: "only", "with", None. Total concentrations output mode.
    """

    # Settings
    prefix_dict = {"G": "gas", "A": "aero"}  # Perfix for gas, aerosol species name
    unit_in = "molec"  # Input unit is molec/cm3
    fakep = "f"  # Fake species prefix

    # Default options - to get average concentrations per species
    # {"time": "ave", "unit": unit_in, "use_fake": False, "winit": False, "total": None}

    # Unit conversion factor
    unit_convs = _get_gecko_unit_conv(raw_result, options.get("unit", unit_in))

    # Start time index
    t0 = 0 if options.get("winit", False) else 1

    # Output concentration mode
    cmode = options.get("time", "ave")

    # Extract concentrations by species
    aconcs = raw_result["conc"]
    conc_dict = {s: {} for s in prefix_dict.values()}
    for i, s in enumerate(raw_result["name"]):
        for k, v in prefix_dict.items():
            if s.startswith(k):
                # Save concentrations for gas/aerosol species found
                p = s[len(k) :]
                conc_dict[v][p] = get_data_based_on_mode(aconcs[t0:, i] * unit_convs[i], cmode)
                break

    # Process fake species concentration if needed
    if options.get("use_fake", False):
        replace_species_in_conc(conc_dict["gas"], fakep)

    # Output total concentrations if needed
    total_mode = options.get("total", None)
    if total_mode in ["only", "with"]:
        tconcs = conc_dict["gas"].copy()
        for s, v in conc_dict["aero"].items():
            if s in tconcs:
                tconcs[s] += v
            else:
                tconcs[s] = v

        if total_mode == "only":  # Otput: {species: total concs}
            return tconcs
        if total_mode == "with":  # Add to the output dict
            conc_dict["total"] = tconcs

    elif total_mode is not None:
        raise ValueError(f"Unknown total mode: {total_mode}")

    # Output: {"gas": {species: concs}, "aero": {species: concs}}
    return conc_dict


def _extract_gecko_bof(bof_filepath: str) -> dict:
    """Read gecko-a box simulation results in the binary file format."""

    # Check file
    if not (bof_filepath and os.path.exists(bof_filepath)):
        raise FileNotFoundError(f"File {bof_filepath} not found for reading binary file.")

    # Read file
    with open(bof_filepath, "rb") as f:

        # Read integers numsp, lensp, ndat
        _, numsp, lensp, ndat = struct.unpack("4i", f.read(16))

        # Check lensp
        if lensp != 8:
            raise ValueError(f"Length of species does not match: {lensp} != 8")

        logger.debug("Read # %d species with %d max species length & %d time points.", numsp, lensp, ndat)

        # Read species names and molar mass
        snames = _get_species_names(f, numsp, lensp)

        # Read molar mass - ZW added for genoa
        mws = _get_molar_mass(f, numsp, lensp)

        # Read time, temperature, and concentrations
        times, temps, concs = _get_time_temp_conc(f, numsp, ndat)

    # Build result dictionary
    raw_result = {"name": snames, "mw": mws, "time": times, "temp": temps, "conc": concs}

    return raw_result
