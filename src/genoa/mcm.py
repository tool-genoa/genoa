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
This module contains functions to process mechanism files from
the Master Chemical Mechanism (MCM) in https://mcm.york.ac.uk/MCM/
"""

import re
from collections import OrderedDict

from .constants import SNAME, NORATE_STR
from .logger import setup_logger
from .photolysis import mcm_to_ssh_photolysis
from .reaction import Reaction
from .setting_global import get_attrs_from_settings
from .species import Species, add_basic_to_species_list
from .utils import isfloat


# Logger
logger = setup_logger(__name__)

# MCM kinetic rate
MCM_SIM_COEF_DICT = {
    "KRO2NO": 1,
    "KRO2HO2": 2,
    "KAPHO2": 3,
    "KAPNO": 4,
    "KRO2NO3": 5,
    "KNO3AL": 6,
    "KDEC": 7,
    "KROPRIM": 8,
    "KROSEC": 9,
    "KCH3O2": 10,
    "K298CH3O2": 11,
    "K14ISOM1": 12,
}

# KMT01 to KMT18
MCM_CMX_COEF_DICT = {f"KMT{i:02}": i for i in range(1, 19)}
MCM_CMX_COEF_DICT.update({"KFPAN": 21, "KBPAN": 22, "KBPPN": 23})

# Manually written rate
MCM_NML_DICT = {
    # ssh_genoa_spec
    # CH3O2 + HO2 -> CH3OOH
    "3.8E-13*EXP(780/TEMP)*(1-1/(1+498*EXP(-1160/TEMP)))": "EXTRA 99 1",
    # CH3O2 + HO2 -> HCHO
    "3.8E-13*EXP(780/TEMP)*(1/(1+498*EXP(-1160/TEMP)))": "EXTRA 99 2",
    # CH3O2 + RO2 -> CH3OH + HCHO
    # 1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
    "2*KCH3O2*RO2*(1-7.18*EXP(-885/TEMP))": "RO2 1 EXTRA 99 3",
    # In case that lump CH3COCH3
    "8.8E-12*EXP(-1320/TEMP)+1.7E-14*EXP(423/TEMP)": "EXTRA 99 4",
    # HCOCO -> 1.0000 CO + 1.0000 OH
    # 5.00E-12*O2*3.2*(1-EXP(-550/TEMP))
    "5.00E-12*O2*3.2*(1-EXP(-550/TEMP))": "TB O2 EXTRA 99 5",
    # CH3O2 + RO2 -> CH3OH + HCHO
    # 1.03E-13*math.exp(365/TEMP)*(1-7.18*EXP(-885/TEMP))
    "2*KCH3O2*RO2*0.5*(1-7.18*EXP(-885/TEMP))": "RO2 1 EXTRA 99 6",
    # Add for LIMONENE species: INDO
    "2.20E+10*EXP(-8174/TEMP)*EXP(1.00E+8/TEMP@(3))": "EXTRA 99 7",
    "8.14E+9*EXP(-8591/TEMP)*EXP(1.00E+8/TEMP@(3))": "EXTRA 99 8",
    # Non-EXTRA ratio
    "2*(K298CH3O2*8.0E-12)@(0.5)*RO2": "RO2 1 ARR 3.34664e-12 0 0",
    "2*(K298CH3O2*2.9E-12*EXP(500/TEMP))@(0.5)*RO2": "RO2 1 ARR 2.01494e-12 0 -250",
    "2*(KCH3O2*7.8E-14*EXP(1000/TEMP))@(0.5)*RO2": "RO2 1 ARR 1.7926e-13 0 -682.5",
    "(KCH3O2*7.8E-14*EXP(1000/TEMP))@(0.5)": "ARR 8.963e-14 0 -682.5",
    "2*KCH3O2*RO2*7.18*EXP(-885/TEMP)": "RO2 1 ARR 1.4791e-12 0 520",
    # TOL MCM
    "2*(KCH3O2*2.4E-14*EXP(1620/TEMP))@(0.5)*RO2": "RO2 1 ARR 9.9438e-14 0 -992.5",
    # Nc12H26
    # 2*(1.03E-13*6.4E-14*math.exp(365/TEMP))**0.5 => 1.6238e-13 182.5
    "2*(KCH3O2*6.4E-14)@(0.5)*RO2": "RO2 1 ARR 1.6238e-13 0 -182.5",
    # 2*(3.5E-13*3E-13)**0.5
    "2*(K298CH3O2*3E-13)@(0.5)*RO2": "RO2 1 ARR 6.4807e-13 0 0",
    # 2 * (1.03E-13*1.6E-12)**0.5 1.6238e-12 (365-2200)*0.5 = -917.5
    "2*(KCH3O2*1.6E-12*EXP(-2200/TEMP))@(0.5)*RO2": "RO2 1 ARR 1.6238e-12 0 917.5",
    # Add for LIMONENE species: INDO
    "1.80E+13*(TEMP/298)@(1.7)*EXP(-4079/TEMP)": "ARR 1.119708E+09 1.7 4079.",
    "1.80E+13*(TEMP/298)@(1.7)*EXP(-4733/TEMP)": "ARR 1.119708E+09 1.7 4733.",
}


def mcm_manual_rates(kin):
    """
    Manually update the kinetic rate string from MCM format to SSH format.
    """
    if kin in MCM_NML_DICT:  # Find pair
        val = MCM_NML_DICT[kin]
        if "EXTRA 99" in val:
            return f"KINETIC {val} 1.0"  # Add a ratio
        return f"KINETIC {val}"

    # Check if there is a factor
    for key, val in MCM_NML_DICT.items():
        if key not in kin:
            continue
        res = kin.replace(key, "").strip()
        try:  # Update factor - 1st case
            factor = float(res[1:])
        except ValueError:
            try:  # Update factor - 2nd case
                factor = float(res[:-1])
            except ValueError:
                logger.error("MCM maunal can not process rate sting: %s\n with key: %s and res: %s", kin, key, res)
                return None
        # Process operator
        operator = res[0]
        if operator == "*":
            ratio = factor
        elif operator == "/":
            ratio = 1 / factor
        else:
            raise ValueError(f"Find unrecognized operator {operator} in {kin}")

        if "EXTRA 99" in val:  # Add ratio
            return f"KINETIC {val} {ratio:6.3E}"

        if "ARR" in val:  # Update c1
            val = val.split(" ")
            ind = val.index("ARR")
            val[ind + 1] = f"{float(val[ind+1])*ratio:6.3E}"
            return f'KINETIC {" ".join(val)}'

        raise ValueError(f"Need to manually update the rate {val} with factor {factor}\n" f"raw line: {kin}")

    return None


def parse_mcm_kinetic_rate(kin):
    """
    Parse a kinetic rate string from MCM format to a dictionary of components.
    """

    # C, return C
    if isfloat(kin):
        return {"c1": float(kin)}
    # Generic rate coefficients
    if kin in MCM_SIM_COEF_DICT:
        return {"GRC": kin}
    # Complex rate coefficients
    if kin in MCM_CMX_COEF_DICT:
        return {"CRC": kin}
    # Third body or RO2-RO2
    if kin in ["H2O", "M", "O2", "N2", "H2", "RO2"]:
        return {"TB": kin}
    # EXP(-Ea/TEMP), return Ea
    if kin.startswith("EXP(") and kin.endswith("/TEMP)"):
        val = kin.split("EXP(")[1].split("/TEMP")[0]
        if not isfloat(val):
            return {"UK": kin}
        return {"c3": float(val) * -1}

    # (TEMP/300)@(-n), return c2 = n, c1 = 1/(300 ** -n)
    if kin.startswith("(TEMP/300)@(") and kin.endswith(")"):
        val = kin.split("TEMP/300)@(")[1].split(")")[0]
        if not isfloat(val):
            return {"UK": kin}
        return {"c2": float(val), "c1": 1 / (300.0 ** float(val))}

    # TEMP@(n), return n
    if kin.startswith("TEMP@(") and kin.endswith(")"):
        val = kin.split("TEMP@(")[1].split(")")[0]
        if not isfloat(val):
            return {"UK": kin}
        return {"c2": float(val)}

    # Unknown
    return {"UK": kin}


def mcm_to_ssh_kinetic_rate(kin: str) -> str:
    """Convert a kinetic rate string from MCM format to SSH format."""

    # Remove space
    kin = kin.replace(" ", "")

    # Photolysis
    rate = mcm_to_ssh_photolysis(kin)
    if rate:
        return rate

    # Manual sets
    rate = mcm_manual_rates(kin)
    if rate:
        return rate

    # Break down to pieces separated by "*"
    df = {}
    for val in kin.split("*"):
        val = val.strip()
        if not val:
            continue

        rate_dict = parse_mcm_kinetic_rate(val)
        for key, value in rate_dict.items():
            # Record new element
            if key not in df:
                df[key] = value
            # Update c1
            elif key == "c1":
                df[key] *= value
            # Update c2 in T**c2 or c3 in exp(-c3/T)
            elif key in ["c2", "c3"]:
                df[key] += value
            # Update unknown
            elif key == "UK":
                df[key] = "\t".join([df[key], value])
            else:
                raise ValueError(
                    f"Find values for the same key {key} that can not be merged.\n"
                    f"Get {df}, {value}, from raw line: {kin}"
                )

    # Output rate
    rate_ssh = "KINETIC "

    # Check unknown elements
    if "UK" in df:
        logger.error("Find unknown elements in %s: %s", kin, df["UK"])
        return f"{rate_ssh}{NORATE_STR}"

    # Ensure consistency of components
    valid_keys = {"c1", "TB"}
    additional_keys = set(df.keys()) - valid_keys
    if len(additional_keys) > 1 and additional_keys != {"c2", "c3"}:
        logger.error("Find multiple keys in %s: %s, %s", kin, df, additional_keys)
        return f"{rate_ssh}{NORATE_STR}"

    # Build output in the SSH format
    c1 = f"{df['c1']:6.3E}" if "c1" in df else "1"

    # Third body or RO2-RO2
    if "TB" in df:
        if df["TB"] == "RO2":
            rate_ssh += "RO2 1 "
        else:
            rate_ssh += f"TB {df['TB']} "

    # MCM specific rates
    if "GRC" in df:
        return f"{rate_ssh}EXTRA 92 {MCM_SIM_COEF_DICT[df['GRC']]} {c1}"
    if "CRC" in df:
        return f"{rate_ssh}EXTRA 93 {MCM_CMX_COEF_DICT[df['CRC']]} {c1}"

    # Arrhenius rate
    c2 = f"{df['c2']:6.3E}" if "c2" in df else "0"
    c3 = f"{df['c3']:6.3E}" if "c3" in df else "0"

    return f"{rate_ssh}ARR {c1} {c2} {c3}"


def process_mcm_facsmile_reactions(reaction_file):
    """Reads MCM reactions from a .fac file and returns a list of Reaction instances."""

    with open(reaction_file, "r", encoding="utf-8") as f:  # mcm_export.fac
        lines = f.read().splitlines()

    reactions, norate = [], []
    for line in lines:
        line = line.strip()

        # Find reaction e.g., % 3.31D-11 : C722OOH + OH = C722O2 ;
        if line.startswith("%") and line.endswith(";") and ":" in line:

            logger.info("Processing line: %s", line)

            # Read kinetic and reactants/products
            parts = line[1:-1].split(":")
            if len(parts) != 2:
                raise ValueError(f"Find reaction that can not be processed: {line}")

            # Build new reaction
            rcn = Reaction()

            # Reactants and products
            rcn.from_line_wo_kinetic(parts[1], separator=" = ")

            # Update notation in the rate
            if "D" in parts[0]:
                # Replace D+ with E+, D- with E-
                parts[0] = parts[0].replace("D+", "E+").replace("D-", "E-")
                # Replace 'D' with 'E', but not if it's part of 'KDEC'
                if "D" in parts[0]:  # Still find D, change to E+ if not followed by EC
                    parts[0] = re.sub(r"D(?!EC)", "E+", parts[0])
            # Record rate
            rcn.rate.string = parts[0].strip()
            rcn.rate.ssh = mcm_to_ssh_kinetic_rate(rcn.rate.string)
            # Check ssh string
            if NORATE_STR in rcn.rate.ssh:
                norate.append(rcn.rate.string)
            elif not rcn.rate.init_rate_with_string():
                raise ValueError(f"Can not initialize rate in {rcn.rate.ssh} from {line}.")

            # Add new reaction
            reactions.append(rcn)

    logger.info("Read # %d reactions from %s.", len(reactions), reaction_file)
    if norate:
        logger.error("Find %d reactions with no rate: %s", len(norate), "\n".join(norate))

    return reactions


def process_mcm_species(species_file: str) -> list:
    """
    Processes MCM species from a species .tsv file and returns a list of Species instances.
    Current format: \t is the separator, with 9 columns
    Name\tSmiles\tInchi\tInchiKey\tFormula\tMass\tExcited\tPeroxy radical\tSynonyms

    """

    logger.info("Reading species list from %s ...", species_file)

    # Get settings
    settings = get_attrs_from_settings({SNAME.SETUP: "basicspecies_dict", SNAME.NEW: "vptype"})
    basic_dict = settings[SNAME.SETUP]
    vptype = settings[SNAME.NEW]

    # Read species info
    splist = OrderedDict()
    ncol = None  # No.columns
    with open(species_file, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    for line in lines:
        line = line.strip()
        # Skip empty lines or comments
        if line == "" or line[0] == "*":
            continue
        # Find header
        if line.startswith("Name"):
            if ncol is not None:
                raise ValueError(f"Find multiple headers in {species_file}.")
            # Check number of columns
            ncol = len(line.split("\t"))
            if ncol != 9:
                raise ValueError(f"ncol {ncol} != 9. Check {line} in {species_file}.")
        # Find species info
        elif ncol:
            # Add new species, read smiles, InChI, InChIKey (no), ...
            # ... formula, mass, excited, perox.radical, synonyms (no)
            parts = [i.strip() for i in line.split("\t") if i != ""]
            if len(parts) < 6:
                logger.warning("Not able to read line with %d columns < 6: %s", len(parts), line)
                continue
            # Check and clean smiles
            if "/" in parts[1] or "\\" in parts[1]:
                logger.info("Find / or \\ in smiles: %s", parts[1])
                parts[1] = parts[1].replace("/", "").replace("\\", "")
            # name: [smiles, formula, mass]
            splist[parts[0]] = [parts[1], parts[4], parts[5]]

    # Check if read
    if ncol is None:
        raise ValueError(f"Find no header line starting with 'Name' in {species_file}. Check the file format.")
    if not splist:
        raise ValueError(f"Find no species in {species_file}. Check the file format.")

    # Record in species list
    species, sps_dict = [], {}
    for sp, vals in splist.items():
        # Buil new species
        isp = Species(name=sp, string="", smiles=vals[0])
        # Check formula
        if "C" not in vals[1]:
            if sp not in basic_dict:
                logger.info("Find new species %s not in basic species set.", sp)
                basic_dict[sp] = float(vals[2])
            continue

        # Update species properties based on smiles
        isp.update_based_on_smiles(vptype)

        # Check radical
        if sp.endswith("O2") and not isp.RO2:
            logger.warning("Find species %s is not RO2: %s.", sp, isp.smiles)

        # Add new species
        sps_dict[sp] = isp
        species.append(isp)

    # Check functional groups
    check_mcm_fgroups(species)

    logger.info("Read # %d out of # %d species from %s.", len(species), len(splist), species_file)

    # Add basic species
    add_basic_to_species_list(species, sps_dict, basic_dict)

    return species


def check_mcm_fgroups(species: list) -> None:
    """Update functional groups in MCM species based on their name and SMILES."""

    voc_from_name = ["PAN", "OOH", "CO3H", "CO2H"]  # Read voc species types from name
    radical_from_name = ["OO", "O2", "O3"]  # Read radical species types from name
    struc_from_smiles = {"=C": "w=", "N": "wN"}  # Read structural species types from SMILES
    counts = {p: 0 for p in voc_from_name + radical_from_name + list(struc_from_smiles.values())}

    for s in species:
        if s.status != 1:
            continue
        sname = s.name
        fgps = []
        if s.radical:  # Radical species
            for p in radical_from_name:
                if sname.endswith(p):
                    fgps.append(p)
        else:  # Not radical
            for p in voc_from_name:
                if sname.endswith(p):
                    fgps.append(p)
        # Structural species
        for k, v in struc_from_smiles.items():
            if k in s.smiles:
                fgps.append(v)
        # Update functional groups
        if fgps:
            logger.info("Update fgroups %s -> %s for MCM species %s w/ %s", s.fgroups, fgps, s.name, s.smiles)
            for p in fgps:
                s.fgroups[p] = 1
                counts[p] += 1

    if sum(counts.values()) > 0:
        logger.info("Read total functional groups in MCM species: %s", counts)
