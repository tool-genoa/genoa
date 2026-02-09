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
This module computes the environment parameters for estimating kinetic rates.
"""

import os
import json
import math

import attrs
import numpy as np

from .constants import SZAS
from .logger import setup_logger
from .setting_global import get_env_option
from .utils import get_data_based_on_mode


# Logger
logger = setup_logger(__name__)

# Global variable
_ENV_PARAMS_LIST = None


@attrs.define(slots=True)
class EnvSettingGlobal:
    """
    Environment parameters for GENOA.
    """

    # Temperature in K
    # Temperature sets [253, 273, 293, 313]
    temp: float = 298.15
    # Pressure
    press: float = 2.55e19
    # Air mass in molec/cm3
    airm: float = 0.0
    # Oxygen mass: airm * 0.2 (or 0.2095)
    # Nitrogen mass: airm * 0.8 (or 0.7809)
    # Hydrogen mass: airm * 5.8e-7

    # Water mass
    h2o: float = 0.0
    # Relative humidity
    rh: float = 0.0
    # Absolute humidity
    sh: float = 0.0

    # Solar zenith angle
    azi: float = 0.0
    # sza zone in szas
    izone: int = -1

    def update_air_mass(self) -> bool:
        """Update the air mass based on the reactants."""
        if self.press == 0.0 or self.temp == 0.0:
            logger.error("Can't update air mass %f with P %f and T %f.", self.airm, self.press, self.temp)
            return False
        self.airm = self.press * 7.242e16 / self.temp
        return True

    def update_humidity(self) -> bool:
        """Update the humidity based on the absolute humidity."""
        if self.temp == 0.0 or self.rh == 0.0 or self.press == 0.0:
            logger.error("Can't update SH %f with T %f, RH %f, P %f.", self.sh, self.temp, self.rh, self.press)
            return False
        self.sh = compute_sh(self.rh, self.temp, self.press)
        self.h2o = compute_gas_phase_water(self.temp, self.rh)
        return True

    def update_basic_parameters(self, temp: float = 298.0, press: float = 2.55e19, rh: float = 0.70) -> bool:
        """Update basic the environment parameters."""
        self.temp = temp
        self.press = press
        self.rh = rh
        self.update_air_mass()
        self.update_humidity()

        return True

    def update_azi(self, azi: float, szas: list) -> bool:
        """Update the solar zenith angle and izone in szas."""

        if azi < 0.0 or azi > 180.0:
            logger.error("\nInvalid solar zenith angle %f.", azi)
            return False

        # Get index in szas
        if azi < szas[0]:
            izone = 0
        elif azi < szas[-1]:
            for i in range(len(szas) - 1):
                if szas[i] <= azi < szas[i + 1]:
                    izone = i
                    break
        else:
            izone = -1

        self.azi = azi
        self.izone = izone

        return True


def get_env_params_list() -> list:
    """
    Returns a list of environment parameters to compute rate coefficients under different conditions.
    """
    global _ENV_PARAMS_LIST
    if _ENV_PARAMS_LIST is None:
        env_option = get_env_option()
        env_list = []
        for temp in env_option.temps:
            for press in env_option.press:
                for rh in env_option.rhs:
                    for azi in env_option.azis:
                        env = EnvSettingGlobal()
                        env.update_basic_parameters(temp, press, rh)
                        env.update_azi(azi, SZAS)
                        env_list.append(env)

        _ENV_PARAMS_LIST = env_list
        logger.info(
            "Generated # %d environmental sets from T in K: %s, Press in atm: %s, RH in %%: %s, SZA in degree: %s",
            len(env_list),
            env_option.temps,
            env_option.press,
            env_option.rhs,
            env_option.azis,
        )

    return _ENV_PARAMS_LIST


def compute_sh(rh: float, temp: float, pres: float) -> float:
    """
    Computes specific humidity.
    Converted from SSH-aerosol: Meteorology.f

    Args:
    rh - Relative humidity
    temp - Temperature (K)
    pres - Pressure (Pa)

    Returns:
    sh - Specific humidity (kg/kg).
    """
    # Calculate saturation vapor pressure (psat)
    psat = 611.2 * math.exp(17.67 * (temp - 273.15) / (temp - 29.65))

    # Calculate specific humidity (sh)
    sh = 1.0 / (1.0 + pres / (psat * rh * 0.62197) - 1.0 / 0.62197)

    return sh


def compute_gas_phase_water(temp0: float, rh0: float) -> float:
    """
    Computes the concentration of water in the gas phase.
    Used in both GECKO and SSH-aerosol

    Args:
    temp0 - Temperature (K)
    rh0 - Relative Humidity

    Returns:
    water - Water concentration in molecules
            per cubic centimeter (molec/cm3)

    """
    # H2O = 6.1078*math.exp(-1.0E0*(597.3-0.57*(TEMP-273.16))*
    #       18./1.986*(1./TEMP-1./273.16))*10./(1.38E-16*TEMP)*
    #       meteo_ave['RelativeHumidity']
    avogadro = 6.02214e23  # avogadro number
    rgas = 8.3144621  # gas constant (J.K-1.mol-1)
    c1 = 610.94
    c2 = 17.625
    c3 = 243.04  # Magnus formula parameters

    # Convert temperature to Celsius
    t_c = temp0 - 273.16
    # Calculate saturation water pressure
    psat_h2o = c1 * math.exp(c2 * t_c / (c3 + t_c))
    # Calculate water pressure
    p_h2o = psat_h2o * rh0 / 100
    # Calculate water concentration (molec/cm3)
    water = p_h2o / (rgas * temp0 * 1e6 / avogadro)

    return water


def get_list_from_file(inval, num=0):
    """
    Use to read condition identifiers from a file or a list
    if num != 0: get number of num from inval
    """
    with open(inval, "r", encoding="utf-8") as f:
        locs = json.loads(f.read())
    nlocs = len(locs)
    if 0 < num <= nlocs:
        locs = locs[:num]
        logger.info("Read # %d out of # %d conditions from file %s.", num, nlocs, inval)
    else:
        logger.info("Read all # %d conditions from file %s.", nlocs, inval)

    return locs


def _get_locations(locs, npz_file):
    """Get a list of [latitude, longitude] pairs from locs."""
    data = np.load(npz_file)
    lats = data["latitudes"]
    lons = data["longitudes"]
    return [[round(lats[i[0]], 2), round(lons[i[1]], 2)] for i in locs]


def get_list_conditions(inval, npz_file):
    """Read input and output a list of condition indentifiers"""

    if not inval:
        raise ValueError(f"Check input {inval} for locations.")

    if isinstance(inval, list):
        locs = inval
    elif "," in inval:  # separate file name and number
        nloc = int(inval.split(",", 1)[1])
        inval1 = inval.split(",", 1)[0]
        locs = get_list_from_file(inval1, nloc)
    elif os.path.exists(inval):
        locs = get_list_from_file(inval)
    else:
        raise ValueError(f"Not able to read locations from {inval}.")
    return locs, _get_locations(locs, npz_file)


def _get_data_from_meteo_file(paths, items, mode="ave"):
    """
    Return the average value of given items from the given paths.

    Args:
        paths (list of str): List of file paths to read data from.
        items (list of str): List of items (columns) to calculate average for.

    Returns:
        dict: A dictionary mapping each item to its average value.
    """

    data = {i: [] for i in items}
    for path in paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Error: File not found - {path}")
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()
            titles = lines[0].strip().split(",")

            # Find index of required items in the file
            item_inds = {i: titles.index(i) for i in items if i in titles}

            for line in lines[1:]:  # Process lines except the title
                values = line.strip().split(",")
                for key, val in item_inds.items():
                    try:  # Read data
                        data[key].append(float(values[val]))
                    except ValueError:
                        logger.warning("Warning: Non-numeric value encountered in %s for %s", path, items)

    # Output based on mode
    for s in items:
        if data[s] == []:
            data[s] = None
        else:
            data[s] = get_data_based_on_mode(np.array(data[s]), mode=mode)
    return data


def update_env_params_from_files(locs_info, path_locations, path_initfile):
    """Compute environmental paramters from given files"""

    env_params = get_env_params_list()
    locs, locations = locs_info

    if locs[0]:
        new_env = EnvSettingGlobal()
        meteo_paths = []
        logger.info("Get environmental parameters for locations: %s", locs)
        # Get locations if not provided
        if locations is None:
            get_list_conditions(locs, path_locations)
        for il in locs:  # y,x,m = iloc
            ipath = f"{path_initfile}/m{il[2]}/y{il[0]}/x{il[1]}/meteo.dat"
            meteo_paths.append(ipath)

        # Get data from meteo files
        # Time,Temperature,Pressure,RelativeHumidity,SpecificHumidity
        meteo_ave = _get_data_from_meteo_file(meteo_paths, ["Temperature", "Pressure", "RelativeHumidity"])
        # Update parameters from meteo data
        new_env.update_basic_parameters(
            temp=meteo_ave["Temperature"],
            rh=meteo_ave["RelativeHumidity"],
            press=meteo_ave["Pressure"],
        )

        # Get solar zenith angle
        # new_env.azi = sza_max(locs_info)

        env_params.append(new_env)

    # Print to check
    logger.info("Updated environmental parameters: %s", env_params)

    return env_params
