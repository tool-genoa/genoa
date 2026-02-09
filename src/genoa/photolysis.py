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
This module contains functions to compute photolysis rates in GENOA.
"""

import os
from datetime import date, datetime

import numpy as np

from .constants import SNAME, SZAS, NSZA, NORATE_STR
from .utils import isfloat
from .setting_global import get_attrs_from_settings
from .logger import setup_logger


# Logger
logger = setup_logger(__name__)

# Photolysis coefficients from the photolysis file
_PHOTO_COEFS_DICT = None

# Photolysis attenuation factor
ATTENUATION = 1.0

# Need to keep update with MCM website
MCM_PHOT_COEF_DICT = {
    1: [6.073e-05, 1.743, 0.474],
    2: [4.775e-04, 0.298, 0.080],
    3: [1.041e-05, 0.723, 0.279],
    4: [1.165e-02, 0.244, 0.267],
    5: [2.485e-02, 0.168, 0.108],
    6: [1.747e-01, 0.155, 0.125],
    7: [2.644e-03, 0.261, 0.288],
    8: [9.312e-07, 1.230, 0.307],
    11: [4.642e-05, 0.762, 0.353],
    12: [6.853e-05, 0.477, 0.323],
    13: [7.344e-06, 1.202, 0.417],
    14: [2.879e-05, 1.067, 0.358],
    15: [2.792e-05, 0.805, 0.338],
    16: [1.675e-05, 0.805, 0.338],
    17: [7.914e-05, 0.764, 0.364],
    18: [1.140e-05, 0.396, 0.298],
    19: [1.140e-05, 0.396, 0.298],
    20: [7.600e-04, 0.396, 0.298],
    21: [7.992e-07, 1.578, 0.271],
    22: [5.804e-06, 1.092, 0.377],
    23: [1.836e-05, 0.395, 0.296],
    24: [1.836e-05, 0.395, 0.296],
    31: [6.845e-05, 0.130, 0.201],
    32: [1.032e-05, 0.130, 0.201],
    33: [3.802e-05, 0.644, 0.312],
    34: [1.537e-04, 0.170, 0.208],
    35: [3.326e-04, 0.148, 0.215],
    41: [7.649e-06, 0.682, 0.279],
    51: [1.588e-06, 1.154, 0.318],
    52: [1.907e-06, 1.244, 0.335],
    53: [2.485e-06, 1.196, 0.328],
    54: [4.095e-06, 1.111, 0.316],
    55: [1.135e-05, 0.974, 0.309],
    56: [7.549e-06, 1.015, 0.324],
    57: [3.363e-06, 1.296, 0.322],
    61: [7.537e-04, 0.499, 0.266],
}


def days_in_month(year: int, month: int) -> int:
    """Returns the number of days in the given month of the given year."""
    if month == 12:
        next_month = 1
        year += 1
    else:
        next_month = month + 1

    return (date(year, next_month, 1) - date(year, month, 1)).days


def sza_max(locs_info: list, year: int = 2015) -> float:
    """Computes the average of the maximum Solar Zenith Angle (SZA) at given locations"""

    max_sza = 0.0  # Initialize max_sza
    # Get locs
    locs, locations = locs_info
    for i, iloc in enumerate(locs):
        latitude, longitude = locations[i]
        month = iloc[2] + 1  # Get current month
        num_days = days_in_month(year, month)
        max_at_loc = 0.0
        # Compute max sza during this month
        for day in range(num_days):
            for hour in range(24):
                # Calculate seconds since the start of the year
                tu = (datetime(year, month, day + 1, hour) - datetime(year, 1, 1)).total_seconds()
                max_at_loc = max(max_at_loc, ssh_muzero(tu, latitude, longitude))

        max_sza = max(max_sza, max_at_loc)

    if max_sza == 0.0:
        raise ValueError("Maximum sza is zero. Check the input data and calculations.")

    logger.info("Got max SZA at %s from %s locations.", max_sza, len(locs))

    return max_sza


def mcm_photolysis_coefs(n: int) -> str:
    """
    Convert MCM photolysis index to photolysis coefficients in the format [l, m, n].

    These coefficients (read from function mcm_manual_rates()) can be used in the formula:
    J = l * np.cos(X)**m * EXP(-n * (1 / np.cos(X)))
    """

    # Get coefficients
    if n not in MCM_PHOT_COEF_DICT:
        logger.error("Non-recognized MCM photolysis index: %s", n)
        return None

    coefs = MCM_PHOT_COEF_DICT[n]
    return f"{coefs[0]} {coefs[1]:5.3f} {coefs[2]:5.3f}"


def mcm_to_ssh_photolysis(jin: str) -> str:
    """Convert photolysis kinetic rate from MCM format to a format readable by SSH-aerosol."""

    # Check format
    if not jin.startswith("J<"):
        return None

    # Extract photolysis index from jin in the format "J<1>*0.03*2*..."
    values = jin.split("*")

    # Photolysis index
    ind = values[0].split("J<")[1].split(">")[0]

    # Factor to photolysis rate
    factor = 1.0
    if len(values) > 1:
        for i in values[1:]:
            factor *= float(i)

    if not ind.isdigit():
        return f"KINETIC photolysis {jin} with non-integer index: {ind} {NORATE_STR}"
    # Get photolysis coefficients
    phot_str = mcm_photolysis_coefs(int(ind))
    if not phot_str:
        return f"KINETIC photolysis {jin} with invalid index: {ind} {NORATE_STR}"
    return f"KINETIC EXTRA 91 {phot_str} {factor:6.3E}"


def read_photo_coefs_from_file(pfile) -> dict:
    """Read photolysis tabulation from a file in the GECKO format"""

    if not pfile or not os.path.exists(pfile):
        return {}

    photo_coefs_all = {}
    with open(pfile, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("/") or line.startswith("END"):
                continue
            if "PHOT" in line:
                parts = [i for i in line.split(" ") if i != ""]
                if parts[-1] != str(NSZA):
                    raise ValueError(f"Not the same NSZA: {NSZA}", line, parts[-1])
                ikey = int(parts[2])
                if ikey in photo_coefs_all:
                    raise ValueError(f"Find duplicated key: {ikey}")
                photo_coefs_all[ikey] = []

                iphoto = photo_coefs_all[ikey]
            else:  # Read data
                parts = [float(i) for i in line.split(" ") if isfloat(i)]
                if len(parts) == 2:
                    iphoto.append(parts[1])  # sza, coef

    return photo_coefs_all


def get_photo_coefs_dict() -> dict:
    """Get photolysis coefficients from the photolysis file"""

    global _PHOTO_COEFS_DICT
    if _PHOTO_COEFS_DICT is None:
        photo_file = get_attrs_from_settings({SNAME.ENV: "photolysis_file"})
        _PHOTO_COEFS_DICT = read_photo_coefs_from_file(photo_file)
    return _PHOTO_COEFS_DICT


def get_photolysis_rate_spack(azi: float, izone: int, photo_in: list) -> float:
    """Compute photolysis rate using spack mode."""
    # if azi >= 90:
    #    return 0.0  # No photolysis

    photo_out = ssh_spl3(NSZA, SZAS, photo_in)
    photo = photo_out[3, izone]
    photo = photo_out[2, izone] + (azi - SZAS[izone]) * photo
    photo = photo_out[1, izone] + (azi - SZAS[izone]) * photo
    photo = photo_out[0, izone] + (azi - SZAS[izone]) * photo

    return max(photo * ATTENUATION, 0.0)


def get_photolysis_rate_gecko(azi: float, izone: int, photo_in: list) -> float:
    """Compute photolysis rate using spack mode."""
    if azi >= 90:
        return 0.0  # No photolysis

    photo_out = ssh_gck_cspline(NSZA, SZAS, photo_in)
    zp = (azi - SZAS[izone]) / (SZAS[izone + 1] - SZAS[izone])
    photo = ((zp * photo_out[0, izone] + photo_out[1, izone]) * zp + photo_out[2, izone]) * zp + photo_out[3, izone]
    photo = abs(photo)

    return max(photo * ATTENUATION, 0.0)


def ssh_gck_cspline(ndat, xp, yp):
    """
    Convert from cubicsplinf.f90 in ssh-aerosol
    Calculate cubic spline coefficients.

    Args:
    xp : numpy array of x-coordinates.
    yp : numpy array of y-coordinates.

    Returns:
    yout : numpy array of cubic spline coefficients.
    """
    # initialize
    yout = np.zeros((4, ndat - 1))
    w = np.zeros(ndat * 3)
    coe = np.zeros(ndat * 3)
    deriv = np.array([0, 0])

    j = 1  # derivative at end point provided
    xdel = xp[1] - xp[0]
    ydel = yp[1] - yp[0]

    if j == 2:
        coe[1] = 0.0
        w[1] = 0.5 * xdel * xdel * deriv[0]
    else:
        coe[0] = xdel * deriv[0]
        coe[1] = 1.0
        w[1] = ydel - coe[0]

    m = ndat - 2
    if m > 0:
        for i in range(1, m + 1):
            xdel1 = xdel
            xdel = xp[i + 1] - xp[i]
            ratio = xdel1 / xdel
            idx = 3 * i
            coe[idx] = -ratio / (2.0 - coe[idx - 1])
            w[idx] = (-ydel - w[idx - 1]) / (2.0 - coe[idx - 1])
            coe[idx + 1] = -ratio * ratio / (ratio - coe[idx])
            w[idx + 1] = (ydel - w[idx]) / (ratio - coe[idx])
            ydel = yp[i + 1] - yp[i]
            coe[idx + 2] = 1.0 / (1.0 - coe[idx + 1])
            w[idx + 2] = (ydel - w[idx + 1]) / (1.0 - coe[idx + 1])

    if j == 1:
        coe[3 * ndat - 3] = (xdel * deriv[1] - ydel - w[3 * ndat - 4]) / (2.0 - coe[3 * ndat - 4])
    else:
        coe[3 * ndat - 3] = (xdel * xdel * deriv[1] / 2.0 - w[3 * ndat - 4]) / (3.0 - coe[3 * ndat - 4])

    m = 3 * ndat - 6
    if m > 0:
        for ii in range(1, m + 1):
            i = m - ii + 3
            coe[i] = w[i] - coe[i] * coe[i + 1]

    if j == 1:
        coe[1] = w[1] - coe[2]
    else:
        coe[0] = yp[1] - yp[0] - w[1] - coe[2]
        coe[1] = w[1]

    for i in range(1, ndat):
        xdel = xp[i] - xp[i - 1]
        for j in range(1, 4):
            yout[j - 1, i - 1] = coe[3 * i - 2 - j]
        yout[3, i - 1] = yp[i - 1]

    return yout


def ssh_spl3(n1, a, b):
    """
    Calculate spline coefficients.

    Args:
    a : numpy array of x-coordinates.
    b : numpy array of y-coordinates.

    Returns:
    c : numpy array
        Array of calculated spline coefficients.
    """

    # Initialize array sizes based on input
    n0 = n1 - 1
    c = np.zeros((4, n0))
    d = np.zeros(n1)
    c20 = np.zeros(n1)
    c21 = np.zeros(n1)
    c30 = np.zeros(n1)
    c31 = np.zeros(n1)
    c40 = np.zeros(n1)
    c41 = np.zeros(n1)

    for i in range(n0):
        c[0, i] = b[i]
        d[i] = a[i + 1] - a[i]

    # Initial coefficients setup
    c31[0] = 1.0
    c40[0] = (b[1] - c[0, 0] - c20[0] * d[0] - c30[0] * d[0] ** 2) / d[0] ** 3
    c41[0] = (-c21[0] * d[0] - c31[0] * d[0] ** 2) / d[0] ** 3

    # Calculate spline coefficients
    for i in range(1, n0):
        c20[i] = c20[i - 1] + 2.0 * c30[i - 1] * d[i - 1] + 3.0 * c40[i - 1] * d[i - 1] ** 2
        c21[i] = c21[i - 1] + 2.0 * c31[i - 1] * d[i - 1] + 3.0 * c41[i - 1] * d[i - 1] ** 2
        c30[i] = c30[i - 1] + 3.0 * c40[i - 1] * d[i - 1]
        c31[i] = c31[i - 1] + 3.0 * c41[i - 1] * d[i - 1]
        c40[i] = (b[i + 1] - c[0, i] - c20[i] * d[i] - c30[i] * d[i] ** 2) / d[i] ** 3
        c41[i] = (-c21[i] * d[i] - c31[i] * d[i] ** 2) / d[i] ** 3

    cc31 = c20[n0 - 1] + 2.0 * c30[n0 - 1] * d[n0 - 1] + 3.0 * c40[n0 - 1] * d[n0 - 1] ** 2
    cc31 = -cc31 / (c21[n0 - 1] + 2.0 * c31[n0 - 1] * d[n0 - 1] + 3.0 * c41[n0 - 1] * d[n0 - 1] ** 2)
    # Compute c matrix
    for i in range(n0):
        c[1, i] = c20[i] + cc31 * c21[i]
        c[2, i] = c30[i] + cc31 * c31[i]
        c[3, i] = c40[i] + cc31 * c41[i]

    return c


# Convert from angzen.edg.f in ssh-aerosol
def ssh_tsolv(tu, long):
    """
    calcul de l'heure solaire "vraie"
    fonction de l'heure tu (gmt)
    exprimee en secondes dans l'annee
    et de la longitude en degres

    retour = valeur en heures
    """

    t = tu / 86400.0
    dlaent = ssh_aent(t)
    t = t - dlaent
    t = t * 24.0
    dlcorheu = ssh_corheu(tu)
    tsolv = t + long / 15.0 + dlcorheu
    tsolv = ssh_amodp(tsolv, 24.0)
    tsolv = tsolv - 12.0
    tsolv = tsolv * np.pi / 12.0
    return tsolv


def ssh_declin(tsec):
    """
    calcul de la declinaison du soleil fonction
    du temps tu en secondes dans l'annee
    """
    t = 1.0 + (tsec / 86400.0)
    t = t + 0.1
    t = 2 * np.pi * t / 365.0
    declin = 0.006918
    declin = declin - 0.399912 * np.cos(t) + 0.070257 * np.sin(t)
    declin = declin - 0.006758 * np.cos(2 * t) + 0.000907 * np.sin(2 * t)
    declin = declin - 0.002697 * np.cos(3 * t) + 0.001480 * np.sin(3 * t)
    return declin


def ssh_muzero(tu, long, lat):
    """
    calcul du cosinus de l'angle zenithal
    fonction de l'heure tu (gmt)
    exprimee en secondes dans l'annee
    et de la longitude,latitude en degres
    """
    tul = tu % 31536000 if tu > 31536000 else tu
    hr = ssh_tsolv(tul, long)
    decl = ssh_declin(tul)
    flat = lat * np.pi / 180.0
    muzero = np.sin(decl) * np.sin(flat) + np.cos(decl) * np.cos(flat) * np.cos(hr)
    azi = np.degrees(np.arccos(np.clip(muzero, -1.0, 1.0)))
    return azi


def ssh_corheu(tsec):
    """
    calcul de la difference temps "vrai" - temps "moyen"
    fonction du temps tu en secondes dans l'annee

    retour = valeur en heures
    """
    t = 1.0 + (tsec / 86400.0)
    t = t + 0.1
    t = 2.0 * np.pi * t / 365.0
    corheu = 0.000075
    corheu = corheu + 0.001868 * np.cos(t) - 0.032077 * np.sin(t)
    corheu = corheu - 0.014615 * np.cos(2.0 * t) - 0.040849 * np.sin(2.0 * t)
    corheu = corheu * 12.0 / np.pi
    return corheu


def ssh_aent(x):
    """
    PARTIE ENTIERE "USUELLE" : UNIQUE ENTIER K TEL QUE K =<X < K+1
    """
    if x >= 0:
        return np.floor(x)

    return np.floor(x) - 1


def ssh_amodp(x, y):
    """
    L'UNIQUE REEL F TEL QUE X= N*|Y| +F, 0<= F < |Y|
    """
    z = np.max([y, -y])
    if z > 0:
        dlaent = ssh_aent(x / z)
        return x - z * dlaent

    raise ValueError("Invalid input for ssh_amodp")
