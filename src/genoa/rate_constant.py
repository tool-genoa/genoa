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
This module computes knietic rate constants at different environmental condition.
"""
import math
from typing import List, Any

import numpy as np

from .gecko_cst import EXTRA_LABEL
from .env_parameters import get_env_params_list
from .logger import setup_logger
from .photolysis import get_photo_coefs_dict, get_photolysis_rate_spack, get_photolysis_rate_gecko


# Logger
logger = setup_logger(__name__)


def update_rate_cst_values(reactions: list, force_update: bool = True) -> None:
    """Update rate constants (w/o set to unimolecular reactions)"""

    # Get environmental parameters
    env_params = get_env_params_list()
    nenv = len(env_params)

    # Get rate lines and bimolacular species if any
    to_process = []
    for rcn in reactions:

        if rcn.status < 1:  # Only update valid rcn w/o kvalues
            continue

        # Already updated
        if rcn.rate.kvalues is not None and not force_update:
            continue

        to_process.append(rcn.rate)
        rcn.rate.kvalues = np.zeros(nenv)

    if not to_process:
        logger.warning("Not updating kinetic values: no reactions meet the update criteria.")
        return

    logger.info("Updating kvalues with # %d environments for %d reactions ...", nenv, len(to_process))

    # Get photolysis coefs dict
    photo_in_all = get_photo_coefs_dict()

    # Function mappings for EXTRA label cases
    extra_funcs = {
        91: ssh_mcm_rate_91,
        92: ssh_mcm_rate_92,
        99: ssh_genoa_spec,  # Add NEW LABEL HERE
        10: ssh_spec_cb05_10,
        # 20: ssh_spec_racm2_20,
        # 30: ssh_spec_melchior2_30
    }

    # Get values for different environmental condition
    for ienv, env in enumerate(env_params):

        # Init dicts for mcm93
        init_mcm93_list = ssh_mcm_rate_93_init(env)

        for rte in to_process:

            # Init with ratio
            qfor = rte.ratio

            # Check the main keyword
            mkey, coefs = rte.mainkey, rte.coefs
            if mkey == "ARR":
                qfor *= arrhenius_law(coefs, env.temp)
            elif mkey == "FALLOFF":
                qfor *= ssh_gck_forate(coefs, env)
            elif mkey == "PHOT":
                qfor *= get_photolysis_rate_gecko(env.azi, env.izone, photo_in_all[int(coefs[0])])
            elif mkey == "PHOTOLYSIS":
                qfor *= get_photolysis_rate_spack(env.azi, env.izone, coefs)
            elif "EXTRA" in mkey:
                mlabel = int(mkey.rsplit(" ", 1)[-1])
                if mlabel in EXTRA_LABEL:  # GECKO-A rate constant
                    qfor *= ssh_gck_extrarate(mlabel, coefs, env)
                elif mlabel == 93:
                    qfor *= ssh_mcm_rate_93(coefs, init_mcm93_list)
                else:
                    qfor *= extra_funcs[mlabel](coefs, env)
            else:
                raise ValueError(f"Can't process KINETIC keyword {mkey}")

            # Check 2nd keywords
            if rte.tb:  # 3rd body reactions
                rlabel = rte.tb
                if rlabel == "O2":
                    qfor *= env.airm * 0.2
                elif rlabel == "H2O":
                    qfor *= env.h2o
                elif rlabel == "M":
                    qfor *= env.airm
                elif rlabel == "N2":
                    qfor *= env.airm * 0.8
                elif rlabel == "H2":
                    qfor *= env.airm * 5.8e-7
                else:
                    raise ValueError(f"3rd body rcn type unknown: {rlabel}")

            # if rte.ro2:
            #    qfor *= ref_concs.get(rte.ro2, 0.0)  # RO2 concentration

            # Update kvalues
            rte.kvalues[ienv] = qfor


def arrhenius_law(coefs: list, temp: float) -> float:
    """Get the rate constant using the Arrhenius law"""
    return coefs[0] * (temp ** coefs[1]) * math.exp(-coefs[2] / temp)


def ssh_gck_forate(coefs: list, envp: Any) -> float:
    """
    Converted from GECKO-A to compute kinetic rates for GECKO-A FALLOFF reactions.
    Extracted from forate(). Late update 2025-05-22
    """

    # Get coefs
    arrs3 = coefs[0:3]
    fall_coefs = coefs[3:]

    xki = arrs3[0] * ((envp.temp / 300) ** arrs3[1]) * math.exp(-arrs3[2] / envp.temp)
    xkom = fall_coefs[0] * ((envp.temp / 300) ** fall_coefs[1]) * math.exp(-fall_coefs[2] / envp.temp) * envp.airm
    kratio = xkom / xki
    factor = 1 / (1 + (math.log10(kratio)) ** 2)

    # Usual Troe equation
    if fall_coefs[3] != 0 and fall_coefs[4] == 0:
        qfor = (xkom / (1 + kratio)) * (fall_coefs[3] ** factor)
    # Equation for MCM rates constant
    elif fall_coefs[3] != 0 and fall_coefs[4] == 1:
        factor2 = 1 / (1 + (math.log10(kratio) / (0.75 - 1.27 * math.log10(fall_coefs[3]))) ** 2)
        qfor = (xkom / (1 + kratio)) * (fall_coefs[3] ** factor2)
    elif fall_coefs[3] == 0 and fall_coefs[4] == 204 and fall_coefs[5] == 0.17 and fall_coefs[6] == 51:
        fcent = math.exp(-envp.temp / fall_coefs[4]) + fall_coefs[5] + math.exp(-fall_coefs[6] / envp.temp)
        factor2 = 1 / (1 + (math.log10(kratio) / (0.75 - 1.27 * math.log10(fcent))) ** 2)
        qfor = (xkom / (1 + kratio)) * (fcent**factor2)
    else:
        logger.error("MCk: unexpected setup in falloff rxn. Got fall_coefs: %s", fall_coefs)
        raise ValueError("MCk: unexpected setup in falloff rxn")

    return qfor


def ssh_gck_extrarate(label: int, coefs: list, envp: Any) -> float:
    """
    Converted from GECKO-A to compute kinetic rates for GECKO-A EXTRA reactions.
    Extracted from extrarate(). Late update 2025-05-22
    """

    # Get coefs
    arrs3 = coefs[0:3]
    extra_coefs = coefs[3:]

    qfor = arrhenius_law(arrs3, envp.temp)

    # O+O2+M=>O3 // 3rd order reaction
    if label == 100:
        qfor = qfor * envp.airm * envp.airm * 0.2
    # reaction with water vapor
    elif label == 500:
        qfor = qfor * envp.h2o
    # reaction with H20 (xk1) and H2O+M (xk2) // specific to HO2+HO2 => H2O2
    elif label == 501:
        xk1 = qfor * envp.h2o
        xk2 = extra_coefs[0] * math.exp(-extra_coefs[2] / envp.temp) * envp.h2o * envp.airm
        qfor = xk1 + xk2
    # reaction with water dimer (Kd according to Scribano et al., 2006)
    elif label == 502:
        # dimerization constant [atm-1]
        kd = 4.7856e-4 * math.exp(1851 / envp.temp - 5.10485e-3 * envp.temp)
        # convert atm-1 to molec-1 cm3
        kd = kd * 8.314 * envp.temp * 1e6 / (1.01325e5 * 6.02e23)
        water_dimer = kd * (envp.h2o**2)
        qfor = qfor * water_dimer
    elif label == 550:  # OH+HNO3 reaction: k=k0 + k3*M/(1+K3*M/K2)
        xk2 = extra_coefs[0] * (envp.temp ** extra_coefs[1]) * math.exp(-extra_coefs[2] / envp.temp)
        xk3 = extra_coefs[3] * (envp.temp ** extra_coefs[4]) * math.exp(-extra_coefs[5] / envp.temp)
        qfor = qfor + xk3 * envp.airm / (1 + ((xk3 * envp.airm) / xk2))
    elif label == 200:  # Add: ISOM for isomerisation reaction
        qfor = qfor * (
            extra_coefs[0] * (envp.temp**4)
            + extra_coefs[1] * (envp.temp**3)
            + extra_coefs[2] * (envp.temp**2)
            + extra_coefs[3] * envp.temp
            + extra_coefs[4]
        )
    else:
        raise ValueError("EXTRA label unknown for GECKO-A", label)

    return qfor


def ssh_mcm_rate_91(coefs: List[float], envp: Any) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for MCM reactions.
    Label EXTRA 91: Photolysis Rates
    """

    qfor = 0.0  # Default no light
    if envp.azi < 90.0:  # with light
        cosx = max(0.0, math.cos(envp.azi / 180.0 * math.pi))  # cosx
        if cosx > 1.0e-10:
            secx = 1.0 / cosx  # secx
            # J = l * np.cos(X)**m * EXP(-n * (1 / np.cos(X)))
            qfor = coefs[0] * (cosx ** coefs[1]) * math.exp(-coefs[2] * secx)

    return qfor


def ssh_mcm_rate_92(coefs: list, envp: Any) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for MCM reactions.
    Label EXTRA 92 for Generic Rate Coefficients. Last update: 2025-05-22
    """

    # Get sublabel
    sublabel = int(coefs[0])

    mcm92 = {
        1: [2.7e-12, 0, -360],  # KRO2NO
        2: [2.91e-13, 0, -1300],  # KRO2HO2
        3: [5.2e-13, 0, -980],  # KAPHO2
        4: [7.5e-12, 0, -290],  # KAPNO
        5: [2.3e-12, 0, 0],  # KRO2NO3
        6: [1.44e-12, 0, 1862],  # KNO3AL
        7: [1.0e6, 0, 0],  # KDEC
        8: [2.5e-14, 0, 300],  # KROPRIM
        9: [2.5e-14, 0, 300],  # KROSEC
        10: [1.03e-13, 0, -365],  # KCH3O2
        11: [3.5e-13, 0, 0],  # K298CH3O2
        12: [3.0e7, 0, 5300],  # K14ISOM1
    }

    if sublabel not in mcm92:
        raise ValueError(f"Error: MCM Generic Rate coeff id not known: {sublabel}")

    return arrhenius_law(mcm92[sublabel], envp.temp)


def _calculate_kd(ka, kb, log_val):
    """Compute factors for complex MCM rate coefficients"""
    return 10 ** (math.log10(log_val) / (1 + (math.log10(ka / kb) / (0.75 - 1.27 * math.log10(log_val))) ** 2))


def ssh_mcm_rate_93_init(envp: Any) -> list:
    """
    Converted from SSH-aerosol to compute kinetic rates for MCM reactions.
    Label EXTRA 93 for Complex Rate Coefficients. Last update: 2025-05-22
    """

    ka_kb_val_values = {
        # KMT01 to KMT18
        1: (
            1.0e-31 * envp.airm * (envp.temp / 300) ** (-1.6),
            5.0e-11 * (envp.temp / 300) ** (-0.3),
            0.85,
        ),
        2: (
            1.3e-31 * envp.airm * (envp.temp / 300) ** (-1.5),
            2.3e-11 * (envp.temp / 300) ** 0.24,
            0.6,
        ),
        3: (
            3.6e-30 * envp.airm * (envp.temp / 300) ** (-4.1),
            1.9e-12 * (envp.temp / 300) ** 0.2,
            0.35,
        ),
        4: (
            1.3e-3 * envp.airm * (envp.temp / 300) ** (-3.5) * math.exp(-11000 / envp.temp),
            9.7e14 * (envp.temp / 300) ** 0.1 * math.exp(-11080 / envp.temp),
            0.35,
        ),
        7: (
            7.4e-31 * envp.airm * (envp.temp / 300) ** (-2.4),
            3.3e-11 * (envp.temp / 300) ** (-0.3),
            0.81,
        ),
        8: (3.2e-30 * envp.airm * (envp.temp / 300) ** (-4.5), 3.0e-11, 0.41),
        9: (1.4e-31 * envp.airm * (envp.temp / 300) ** (-3.1), 4.0e-12, 0.4),
        10: (
            4.10e-5 * envp.airm * math.exp(-10650 / envp.temp),
            6.0e15 * math.exp(-11170 / envp.temp),
            0.4,
        ),
        12: (2.5e-31 * envp.airm * (envp.temp / 300) ** (-2.6), 2.0e-12, 0.53),
        13: (2.5e-30 * envp.airm * (envp.temp / 300) ** (-5.5), 1.8e-11, 0.36),
        14: (
            9.0e-5 * math.exp(-9690 / envp.temp) * envp.airm,
            1.1e16 * math.exp(-10560 / envp.temp),
            0.36,
        ),
        15: (
            8.6e-29 * envp.airm * (envp.temp / 300) ** (-3.1),
            9.0e-12 * (envp.temp / 300) ** (-0.85),
            0.48,
        ),
        16: (
            8e-27 * envp.airm * (envp.temp / 300) ** (-3.5),
            3.0e-11 * (envp.temp / 300) ** (-1),
            0.5,
        ),
        17: (
            5.0e-30 * envp.airm * (envp.temp / 300) ** (-1.5),
            1.0e-12,
            0.17 * math.exp(-51 / envp.temp) + math.exp(-envp.temp / 204),
        ),
        # KFPAN
        21: (
            3.28e-28 * envp.airm * (envp.temp / 300) ** (-6.87),
            1.125e-11 * (envp.temp / 300) ** (-1.105),
            0.30,
        ),
        # KBPAN
        22: (
            1.10e-5 * envp.airm * math.exp(-10100 / envp.temp),
            1.90e17 * math.exp(-14100 / envp.temp),
            0.30,
        ),
        # KBPPN
        23: (
            1.7e-3 * math.exp(-11280 / envp.temp) * envp.airm,
            8.3e16 * math.exp(-13940 / envp.temp),
            0.36,
        ),
    }

    qfor_special_cases = {
        5: 1.44e-13 * (1 + envp.airm / 4.2e19),
        6: 1 + (1.40e-21 * math.exp(2200 / envp.temp) * envp.h2o),
        11: 2.40e-14 * math.exp(460 / envp.temp)
        + (6.50e-34 * math.exp(1335 / envp.temp) * envp.airm)
        / (1 + (6.50e-34 * math.exp(1335 / envp.temp) * envp.airm) / (2.70e-17 * math.exp(2199 / envp.temp))),
        18: 9.5e-39
        * envp.airm
        * 2e-1
        * math.exp(5270 / envp.temp)
        / (1 + 7.5e-29 * envp.airm * 2e-1 * math.exp(5610 / envp.temp)),
    }

    return [ka_kb_val_values, qfor_special_cases]


def ssh_mcm_rate_93(coefs: list, mcm93_list) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for MCM reactions.
    Label EXTRA 93 for Complex Rate Coefficients
    """
    # Get sublabel
    sublabel = int(coefs[0])
    # Get computed dict
    ka_kb_val_values, qfor_special_cases = mcm93_list

    if sublabel in ka_kb_val_values:
        ka, kb, val = ka_kb_val_values[sublabel]
        kd = _calculate_kd(ka, kb, val)
        qfor = (ka * kb) * kd / (ka + kb)
    elif sublabel in qfor_special_cases:
        qfor = qfor_special_cases[sublabel]
    else:
        raise ValueError(f"MCM Complex Rate coef id not known: {sublabel}")

    return qfor


def ssh_spec_cb05_10(coefs: list, envp: Any) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for SPEC reactions.
    Label EXTRA 10. Last update: 2025-05-22
    """

    # Get sublabel
    sublabel = int(coefs[0])

    if sublabel == 1:
        qfor = envp.airm * 6.0e-34 * (envp.temp / 3e2) ** (-2.4) * envp.airm * 0.2
    elif sublabel == 2:
        qfor = 2.3e-13 * math.exp(600 / envp.temp) + 1.7e-33 * envp.airm * math.exp(1000 / envp.temp)
    elif sublabel == 3:
        qfor = (3.22e-34 * math.exp(2800 / envp.temp) + 2.38e-54 * envp.airm * math.exp(3200 / envp.temp)) * envp.h2o
    elif sublabel == 4:
        ka = 2.4e-14 * math.exp(460 / envp.temp)
        kb = 2.7e-17 * math.exp(2199 / envp.temp)
        kd = 6.5e-34 * math.exp(1335 / envp.temp) * envp.airm
        qfor = ka + kd / (1 + kd / kb)
    elif sublabel == 5:
        qfor = 1.44e-13 * (1 + 2.381e-20 * 0.8 * envp.airm)
    elif sublabel == 6:
        ka = 2e-30 * (envp.temp / 300) ** -3
        kb = 2.5e-11
        kd = ka * envp.airm / (1.0 + ka * envp.airm / kb)
        qfor = kd * 0.6 ** (1.0 / (1.0 + math.log10(kd / kb) ** 2)) * 0.885
    elif sublabel == 7:
        qfor = 1.8e-39 * envp.h2o * envp.h2o
    else:
        raise ValueError("Error: CB05 sublabel non known: ", sublabel)
    return qfor


def ssh_spec_racm2_20(coefs: list, envp: Any) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for SPEC reactions.
    Label EXTRA 20. Need update!
    """

    # Get sublabel
    sublabel = int(coefs[0])

    logger.warning("ssh_spec_racm2_20 not update!")

    if sublabel == 1:
        qfor = envp.airm * 6.0e-34 * (envp.temp / 3e2) ** (-2.3)
    elif sublabel == 2:
        qfor = 2.3e-13 * math.exp(600 / envp.temp) + 1.73e-33 * envp.airm * math.exp(1000 / envp.temp)
    elif sublabel == 3:
        qfor = 3.22e-34 * math.exp(2800 / envp.temp) + 2.38e-54 * envp.airm * math.exp(3200 / envp.temp)
    elif sublabel == 4:
        ka = 7.2e-15 * math.exp(785 / envp.temp)
        kb = 4.1e-16 * math.exp(1440 / envp.temp)
        kd = 1.9e-33 * math.exp(725 / envp.temp) * envp.airm
        qfor = ka + kd / (1 + kd / kb)
    elif sublabel == 5:
        qfor = 1.5e-13 * (1 + 2.439e-20 * envp.airm)
    elif sublabel == 6:
        ka = 3.4e-30 * (3e2 / envp.temp) ** 3.2 * envp.airm
        kb = ka / (4.77e-11 * (3e2 / envp.temp) ** 1.4)
        qfor = (ka / (1 + kb)) * 3e-1 ** (1 / (1 + ((math.log10(kb) - 0.12) / 1.2) ** 2))
    elif sublabel == 7:
        qfor = 2.0e-39 * envp.h2o * envp.h2o
    else:
        raise ValueError("Error: RACM2 sublabel non known: ", sublabel)
    return qfor


def ssh_genoa_spec(coefs: list, envp: Any) -> float:
    """
    Converted from SSH-aerosol to compute kinetic rates for additional reactions.
    Label EXTRA 99. Last update: 2025-05-22
    """

    # Get sublabel
    sublabel = int(coefs[0])

    if sublabel == 1:
        qfor = 3.8e-13 * math.exp(780 / envp.temp) * (1 - 1 / (1 + 498 * math.exp(-1160 / envp.temp)))
    elif sublabel == 2:
        qfor = 3.8e-13 * math.exp(780 / envp.temp) * (1 / (1 + 498 * math.exp(-1160 / envp.temp)))
    elif sublabel == 3:
        qfor = 1.03e-13 * math.exp(365 / envp.temp) * (1.0 - 7.18 * math.exp(-885 / envp.temp))
    elif sublabel == 4:
        qfor = 8.8e-12 * math.exp(-1320 / envp.temp) + 1.7e-14 * math.exp(423 / envp.temp)
    elif sublabel == 5:
        qfor = 5.00e-12 * 3.2 * (1.0 - math.exp(-550 / envp.temp))
    elif sublabel == 6:
        qfor = 1.03e-13 * 0.5 * math.exp(365 / envp.temp) * (1 - 7.18 * math.exp(-885 / envp.temp))
    elif sublabel == 7:
        qfor = 2.20e10 * math.exp(-8174 / envp.temp) * math.exp(1e8 / envp.temp**3)
    elif sublabel == 8:
        qfor = 8.14e9 * math.exp(-8591 / envp.temp) * math.exp(1e8 / envp.temp**3)
    elif sublabel == 9:  # ssh_IRDICARB
        ka = 611.2 * math.exp(17.67 * (envp.temp - 273.15) / (envp.temp - 29.65))
        kb = envp.rh * envp.pres / ((0.62197 * (1 - envp.rh) + envp.rh) * ka) * 100
        qfor = 1e-9 * kb**3 - 1e-7 * kb**2 + 3e-7 * kb + 0.0003
    # Add more label HERE
    else:
        raise ValueError(f"Error: SSH gnoea_spec sublabel not known: {sublabel}")

    return qfor
