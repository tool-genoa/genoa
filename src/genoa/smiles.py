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
This module contains functions to handle SMILES strings and extract properties.
"""

import sys
import math
from typing import Any

import numpy as np
from openbabel import openbabel as ob, pybel

# UManSysProp
TOOLS_PATH = "./vendor/UManSysProp_public-1.04-genoa"
if TOOLS_PATH not in sys.path:
    sys.path.append(TOOLS_PATH)
from umansysprop import boiling_points, vapour_pressures

from .logger import setup_logger


# Logger
logger = setup_logger(__name__)

# Basic SMART structures
key_pybel = {
    "C": pybel.Smarts("[#6]"),  # carbon
    "H": pybel.Smarts("[H]"),  # hydrogen
    "O": pybel.Smarts("[O]"),  # oxygen 6
    "N": pybel.Smarts("[N]"),  # nitrogen 7
    # "S":("S")
}


def _check_radical(smiles: str) -> bool:
    """Check if the molecule is a radical"""
    return "[O+]" in smiles or "[O]" in smiles or "O." in smiles


def _check_ro2(smiles: str) -> bool:
    """Check if the molecule is a RO2"""
    return any(i in smiles for i in ("O[O]", "[O]O", "OO."))


def _get_pybelmol(mol):
    """Get a molecule in pybel format"""
    if isinstance(mol, pybel.Molecule):
        return mol
    # if input is a smiles string
    if isinstance(mol, str):
        try:
            mol = pybel.readstring("smi", mol)
            return mol
        except IOError as exc:
            raise ValueError("Not a smi string: ", mol) from exc
    else:
        logger.error("Not recognize pymol data type: %s", mol)
        raise ValueError("Not recognize pymol data type, not mol and not string.", mol)


def _get_functional_groups(smiles: Any) -> dict:
    """return a dictionary of the key functional groups read by pybel
    the key functional group can be modified in Parameter.py"""

    pymol = _get_pybelmol(smiles)

    fgroups = {}

    # Count atoms
    for k, v in key_pybel.items():
        tmp = v.findall(pymol)
        if tmp:
            fgroups[k] = len(tmp) * 1.0
        else:
            fgroups[k] = 0.0

    # Count rings
    fgroups.update(_get_ring_in_dict(pymol))

    return fgroups


def _get_ring_in_dict(mol) -> dict:
    """Analyze and return ring structure info from a pybel.Molecule."""

    ring_info = {
        "R": "The number of total rings",
        "Rc": "The number of aromatic rings",
        "RC": "The number of carbocyclic rings",
        "RN": "The number of heterocyclic rings",
    }
    rdict = {k: 0 for k in ring_info}

    obmol = mol.OBMol
    rings = obmol.GetSSSR()

    for ring in rings:
        atom_ids = list(ring._path)
        atoms = [obmol.GetAtom(i + 1) for i in atom_ids]  # OBAtom is 1-based
        atomic_nums = [atom.GetAtomicNum() for atom in atoms if atom]

        is_aromatic = all(atom.IsAromatic() for atom in atoms if atom)
        is_carbocyclic = all(Z == 6 for Z in atomic_nums)
        is_heterocyclic = any(Z not in (6, 1) for Z in atomic_nums)

        rdict["R"] += 1
        if is_aromatic:
            rdict["Rc"] += 1
        if is_carbocyclic:
            rdict["RC"] += 1
        if is_heterocyclic:
            rdict["RN"] += 1

    return {k: v for k, v in rdict.items() if v > 0}


def _get_saturation_vapor_pressure(mol, mode="evap", temperature=298):
    """return (Psat,T) Psat at the certain temperature T in unit: atm, K"""

    # Boiling points [(K)] bp
    bp_methods = {
        "0": boiling_points.nannoolal,
        "1": boiling_points.stein_and_brown,
        "2": boiling_points.joback_and_reid,
    }
    # Vapour pressures [log10 (atm) at a specific temperature] vp
    vp_methods = {
        "evap": vapour_pressures.evaporation,
        "evap2": vapour_pressures.evaporation2,
        "simpol": vapour_pressures.simpol,
        "0": vapour_pressures.nannoolal,
        "1": vapour_pressures.myrdal_and_yalkowsky,
    }

    pymol = _get_pybelmol(mol)

    if mode in vp_methods:
        vp = vp_methods[mode](pymol, temperature)
    elif "VP" in mode and "BP" in mode:
        bp = bp_methods[mode[5]](pymol)
        vp = vp_methods[mode[2]](pymol, temperature, bp)
    else:
        raise ValueError("Not recognize Psat mode: ", mode)

    return 10**vp


def _get_estimated_psat(mol, mode="ave", temperature=298.0):
    """return estimated Psat at the certain temperature T in unit: atm, K"""

    basic_methods = ["evap", "evap2", "simpol", "VP0BP0", "VP0BP1", "VP0BP2", "VP1BP0", "VP1BP1", "VP1BP2"]
    advanced_methods = ["ave", "ave_log10", "ave_log10_3", "ave_log10_3_2", "min", "max"]

    if mode in basic_methods:
        return _get_saturation_vapor_pressure(mol, mode=mode, temperature=temperature)

    if mode not in advanced_methods:
        logger.error("Not recognize mode: %s. Available methods: %s", mode, basic_methods + advanced_methods)
        raise ValueError("Not recognize mode: ", mode)

    # Get all the Psat values
    vp_list = np.zeros(len(basic_methods))
    for i, vp in enumerate(basic_methods):
        vp_list[i] = _get_saturation_vapor_pressure(mol, mode=vp, temperature=temperature)

    # Define a mapping of modes to computations
    mode_computations = {
        "ave": np.average(vp_list),
        "ave_log10": 10 ** np.average(np.log10(vp_list)),
        "ave_log10_3": 10 ** np.average(np.sort(np.log10(vp_list))[:3]),
        "ave_log10_3_2": 10 ** (np.average(np.sort(np.log10(vp_list))[:3]) / 2.0),
        "min": np.min(vp_list),
        "max": np.max(vp_list),
    }

    return mode_computations.get(mode, 0.0)


def _get_enthalpy_vaporization(mol: Any, mode: str = "evap") -> float:
    """return Hvap in unit KJ.mol-1"""

    temp1, temp2 = 298.0, 308.0
    psat1 = _get_estimated_psat(mol, mode=mode, temperature=temp1)
    psat2 = _get_estimated_psat(mol, mode=mode, temperature=temp2)

    #  Clausius-Clapeyron equation
    return abs(math.log(psat2 / psat1) * 8.31446 / (1.0 / temp1 - 1.0 / temp2)) / 1000.0


def string_to_smiles(string: str, to_canonical: bool = False) -> str:
    """Convert input smiles to its canonical smiles format"""

    # initialize Open Babel objects
    obj = ob.OBConversion()
    obj.SetInAndOutFormats("smi", "can")

    mol = ob.OBMol()
    is_valid = obj.ReadString(mol, string)

    # check if input smiles is valid
    if not is_valid:
        logger.error("Invalid smiles string: %s", string)
        return False

    if to_canonical:  # convert to canonical smiles

        # Not convert radicals as it may change the structure
        if _check_radical(string):
            logger.info("Radical smiles not converted to canonical: %s", string)
            return string

        return obj.WriteString(mol).strip()

    return string


def get_properties_from_smiles(smiles: str, vp_type: str) -> dict:
    """Get properties from smiles string"""

    # get molecule
    pymol = _get_pybelmol(smiles)

    properties = {}

    # Check for Radicals and RO2
    if _check_radical(smiles):  # Radical
        properties["radical"] = True
        properties["condensable"] = False
        if _check_ro2(smiles):  # RO2
            properties["RO2"] = 1  # Default Pool 1

    else:  # Condensable
        properties["radical"] = False
        properties["condensable"] = True

    properties["formula"] = pymol.formula
    properties["mass"] = pymol.molwt
    properties["fgroups"] = _get_functional_groups(pymol)
    properties["psat_atm"] = _get_estimated_psat(pymol, vp_type)
    properties["dhvap_kj"] = _get_enthalpy_vaporization(pymol)
    return properties
