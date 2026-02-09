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
This module contains functions to post-process the SSH-aerosol simulation results.
"""

import os

import numpy as np
from netCDF4 import Dataset

from .logger import setup_logger
from .utils import get_average_concs, get_data_based_on_mode, get_files_from_folder, replace_species_in_conc
from .unit_conversion import MULTIMASS, STD_PRESS_PA, ug_to_ppb


# Logger
logger = setup_logger(__name__)


# Global variables
NC_FILES = {"gas": "gas.nc", "aero": "aero.nc"}


def get_ssh_concs(folderpath: str, conc_opts: dict) -> dict:
    """
    Get reference concentrations from multiple SSH-aerosol simulation results.
    output a dict with two keys: "gas" and "aero" for {species: concentrations}
    """

    # Get file paths
    gas_paths = get_files_from_folder(folderpath, NC_FILES["gas"])
    gas_dir = [os.path.dirname(f) for f in gas_paths]
    if not gas_dir:
        raise FileNotFoundError(f"No SSH-aerosol results found in {folderpath}")

    # Read concentrations
    logger.info("Read SSH-aerosol concentrations from # %d files w/ options: %s", len(gas_dir), conc_opts)
    concs_all = []
    for f in gas_dir:
        concs = _get_conc_from_ssh_nc(f, conc_opts)
        concs_all.append(concs)

    return get_average_concs(concs_all)


def _get_conc_from_ssh_nc(result_path: str, options: dict) -> dict:
    """
    Extract SSH-aerosol simulation results from the result netCDF files:
    output_directory/ gas.nc and aero.nc
    """

    # Settings
    unit_in = "ug"  # Input unit is ug/m3
    aerp = "P"  # Aerosol species prefix
    fakep = "FA"  # Fake species prefix

    # Default options -> average concentrations per species
    # {"time": "ave", "unit": unit_in, "use_fake": False, "winit": False, "total": None}

    # Check result folder
    if not (result_path and os.path.exists(result_path)):
        raise FileNotFoundError(f"Result folder not found: {result_path}")
    for f in NC_FILES.values():
        nc_file = os.path.join(result_path, f)
        if not os.path.exists(nc_file):
            raise FileNotFoundError(f"NetCDF file not found: {nc_file}")

    # Start time index
    t0 = 0 if options.get("winit", False) else 1

    # Output unit
    unit = options.get("unit", unit_in)

    # Output concentration mode
    cmode = options.get("time", "ave")

    # Extract gas & time data
    nc_file = os.path.join(result_path, NC_FILES["gas"])
    with Dataset(nc_file, "r") as ncf:
        # a_time = ncf.variables["Time"][:]  # Time series, nt: number of time step
        # Species names
        gnames = [b"".join(s).decode().strip() for s in ncf.variables["species_name"][:]]
        # Concentrations
        if unit == unit_in:  # No conversion
            gconcs = {s: get_data_based_on_mode(ncf.variables[s][t0:], cmode) for s in gnames}
        else:  # Convert to other units
            # Molecular weights in g/mol
            gmws = ncf.variables["molecular_weight"][:]  # ns: number of species
            # Temperature in K
            a_temp = np.array(ncf.variables["temperature"][t0:])  # nt
            # Pressure in Pa
            a_press = np.array(ncf.variables["pressure"][t0:])  # nt
            f_press = a_press / STD_PRESS_PA
            # Get unit conversion factor
            u_convs = _get_ssh_unit_conv(gnames, gmws, options.get("unit", unit_in), a_temp, f_press)
            gconcs = {s: get_data_based_on_mode(ncf.variables[s][t0:] * u_convs[s], cmode) for s in gnames}

    # Extract aerosol data
    nc_file = os.path.join(result_path, NC_FILES["aero"])
    with Dataset(nc_file, "r") as ncf:
        # Aerosol species names
        anames = [b"".join(s).decode().strip() for s in ncf.variables["species_name"][:]]
        # Concentrations
        if unit == unit_in:  # No conversion
            aconcs = {
                s[1:]: get_data_based_on_mode(np.sum(ncf.variables[s][:, t0:], axis=0), cmode)
                for s in anames
                if (s[0] == aerp and s[1:] in gnames)
            }
        else:  # Convert to other units
            aconcs = {
                s[1:]: get_data_based_on_mode(np.sum(ncf.variables[s][:, t0:], axis=0) * u_convs[s[1:]], cmode)
                for s in anames
                if (s[0] == aerp and s[1:] in gnames)
            }

    # Process fake species concentration if needed
    if options.get("use_fake", False):
        replace_species_in_conc(gconcs, fakep)

    # Output total concentrations if needed
    total_mode = options.get("total", None)
    if total_mode in ["only", "with"]:
        tconcs = gconcs.copy()
        for s, v in aconcs.items():
            if s in tconcs:
                tconcs[s] += v
            else:
                tconcs[s] = v
        if total_mode == "only":  # Output: {species: total concs}
            return tconcs
        if total_mode == "with":
            return {"gas": gconcs, "aero": aconcs, "total": tconcs}
    elif total_mode is not None:
        raise ValueError(f"Unknown total mode: {total_mode}")

    # Output: {"gas": {species: concs}, "aero": {species: concs}}
    return {"gas": gconcs, "aero": aconcs}


def _get_ssh_unit_conv(names: list, mws: list, unit: str, temps: list, press: list) -> dict:
    """
    Get unit conversion factor for SSH-aerosol simulation results.
    Input unit is ug/m3, unit out is unit
    """

    if unit == "ug":  # No conversion
        return {s: 1.0 for s in names}

    if unit == "molec":  # ug/m3 to molec/cm3
        return {s: 1 / (mws[i] * MULTIMASS) for i, s in enumerate(names)}

    if unit == "ppb":  # ug/m3 to ppb
        return {s: ug_to_ppb(1, mws[i], temps, press) for i, s in enumerate(names)}

    raise ValueError(f"Unsupport or unknown output unit: {unit}. Use: ug, molec, or ppb")
