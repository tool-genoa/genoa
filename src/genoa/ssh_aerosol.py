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
This module provides functions to write species and reaction lists in the SSH-aerosol format.
"""

import os

from .constants import SNAME
from .ssh_cst import NGP_SOAP
from .logger import setup_logger
from .reaction import get_reactions_with_fake_species
from .setting_global import get_all_settings_map, get_attrs_from_settings
from .species import Species
from .utils import isfloat


# Logger
logger = setup_logger(__name__)

# Global variables
_SSH_AERO_LINES = None


def get_default_ssh_aero_lines():
    """Get the default aerosol species list for SSH-aerosol"""
    global _SSH_AERO_LINES
    if _SSH_AERO_LINES is None:
        # Get default aerosol species list
        aerosol_list = get_attrs_from_settings({SNAME.SSH: "aero_file"})
        with open(aerosol_list, "r", encoding="utf-8") as f:
            faero = f.readlines()
            # Get the header and tail
            _SSH_AERO_LINES = ["".join(faero[0:10]), faero[-1]]
    return _SSH_AERO_LINES


def to_ssh_soap_str(isps: Species) -> str:
    """
    Converts the soap_strucs attribute to a formatted string
    used for SSH-aerosol aerosol species list.
    """
    if not isps.soap_strucs:  # Empty
        return "-"

    vals = []  # Init
    for i in range(NGP_SOAP):  # Process each group
        j = str(i)
        if j in isps.soap_strucs:
            num = isps.soap_strucs[j]
        else:
            num = 0.0
        vals.append(f"&{num:5.2E}")

    return "".join(vals)  # Reformat into a single string


def to_ssh_aero_list(isps: Species, mode: str = "smiles", first_run: bool = False) -> str:
    """Formats species data into a line for the SSH-aerosol aerosol species list."""

    # Check if it is an aerosol
    if not isps.condensable:
        logger.warning("\nSpecies %s is not an aerosol for writing to SSH-aerosol aerosol list.", isps.name)
        return None

    # Structure
    if mode == "vectors":
        struc = to_ssh_soap_str(isps)
    elif mode == "smiles":
        struc = isps.smiles
    else:
        raise NameError(f"MD: Unknown mode: {mode}")

    # Define the content and formatting style
    content = [
        # Aerosol name
        f"P{isps.name}",
        # Organic in SSH-aerosol - aerosol_type
        "4",
        # group id: 2 for organic - Index_groups
        "2" if isps.grp_id is None else f"{isps.grp_id:d}",
        # MWs
        f"{isps.mass:.2f}",
        # Remove precursor if it is the 1st run
        "--" if first_run else isps.name,
        "687d0",  # coll_fac
        "8.39d0",  # mole_diam
        "30.D-03",  # surf_tens
        "1.0",  # accomod
        "1.30D-06",  # mass_dens
        "0",  # '1' if isps.non_volatile else '0',  # non_volatile
        "BOTH",  # partitioning
        struc,  # smiles/vector
        f"{isps.psat_atm * 760.:5.2E}",  # psat in torr
        f"{isps.dhvap_kj:5.2E}",  # dHvap
        "0.0",  # Henry - compute in ssh-aerosol
        "0.0",  # Tref - compute in ssh-aerosol
    ]

    return "\t".join(content) + "\n"


def update_soap_strucs_from_file(species: list, filename: str, mode: str = "SOAP") -> None:
    """Read SOAP structures from a file and update species with the structures."""

    if not (filename and os.path.exists(filename)):
        raise ValueError(f"Cannot read SOAP strctures from {filename}")

    # Get species infos
    sps_dict = {s.name: s for s in species}

    logger.info("Reading soap_strucs from %s with mode %s ...", filename, mode)

    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    mode_dict = {
        "SOAP": read_soap_strucs_from_soap,
        "SSH": read_soap_strucs_from_ssh,
        "TABLE": read_soap_strucs_from_table,
    }

    if mode not in mode_dict:
        raise ValueError(f"Unrecognized mode: {mode}")
    # Return [sps_updated, nsp_read]
    infos = mode_dict[mode](lines, sps_dict)

    logger.info("%s", _get_pinfo_soap_strucs(filename, sps_dict, infos))


def _get_pinfo_soap_strucs(filename: str, sps_dict: dict, infos: list) -> str:
    """Record the check result of the SOAP structures."""

    if not (infos[0] and infos[1]):
        raise ValueError(f"No valid SOAP structure found in the file {filename}.")

    parts = [f"\nSOAP structures are read from {filename}:"]
    if infos[0]:
        parts.append(f"  Updated # {len(infos[0])} species with SOAP structures.")
    else:
        parts.append("  No species updated with SOAP structures.")

    if infos[1]:
        parts.append(f"  Find # {infos[1]} species from the soapfile but not in the species list.")

    sps = [s.name for s in sps_dict.values() if s.condensable and not s.soap_strucs]
    if sps:
        parts.append(f"  Find # {len(sps)} condensable species still without SOAP structures.")
        parts.extend(["Condensables w/o SOAP:"] + sps)
    else:
        parts.append("  All condensable species now have SOAP structures.")

    return "\n".join(parts)


def read_soap_strucs_from_soap(lines: list, sps_dict: dict) -> list:
    """Read SOAP structures from a file in the SOAP format."""

    iden = "is constructed from smiles:"  # Used to find new species
    sep = "==="  # Separator for the end of a block
    l_checked = -1  # Last checked line
    nsp_read, sps_updated = 0, []  # Counters
    for i, line in enumerate(lines):
        if i <= l_checked:  # Skip checked lines
            continue
        if iden in line:  # Find new species name
            sname = line.split(iden)[0].strip()
            if not sname:
                raise ValueError(f"Empty species name found in line {i}: {line}")
            l_checked = _find_soap_block_ends(lines, i, sep)  # Update checked line
            if sname not in sps_dict:  # Not recorded
                nsp_read += 1
                continue
            # Read SOAP structure
            sps_dict[sname].soap_strucs = _read_next_soap_struc(lines[i + 1 : l_checked + 1], sep)
            sps_updated.append(sname)
    return [sps_updated, nsp_read]


def _find_soap_block_ends(lines: list, start: int, sep: str) -> int:
    """Find the end of a SOAP block: lines[start: end+1]"""
    i = -1
    for line in lines[start:]:
        i += 1  # Count lines
        if line.startswith(sep):  # End of the block
            return i + start
    return start


def _read_next_soap_struc(sublines: list, sep: str) -> dict:
    """Read the given SOAP structure block"""
    struc = {}
    for line in sublines:
        line = line.strip()
        if not line or line.startswith(sep):
            continue
        # Example of a line: 7 C linked to alcohol: 1 1 0
        parts = [i for i in line.split(" ") if i != ""]
        if len(parts) >= 4 and parts[0].isdigit():  # Contains group index and numbers x 3
            ind = parts[0]  # Group index
            if ind in struc:
                raise ValueError(f"Repeat group {ind} found in the SOAP structure.")
            struc[ind] = int(parts[-3])  # Group index and number
        elif not (line.startswith("W") or line.startswith("G")):  # Warning/ WARNING/ Group
            # Group O[N+](=O)[O-]) used instead
            raise ValueError(f"Unrecognized line\n    {line}")
    return struc


def read_soap_strucs_from_ssh(lines: list, sps_dict: dict) -> list:
    """Read SOAP structures from a aerosol species list file in the SSH format."""
    iden = "&"  # Used to find new species
    nsp_read, sps_updated = 0, []  # Counters
    for line in lines:
        if iden in line:
            # Find species name
            sname = line.split("\t")[0][1:]  # Remove "P"
            if not sname:
                raise ValueError(f"Empty species name found in line: {line}")
            # Find vector groups
            parts = [float(i) for i in line.split(iden) if isfloat(i)]
            if len(parts) != NGP_SOAP:
                raise ValueError(f"Unrecognized line: {line} with {len(parts)} soap groups.")
            # Record if the species is in the list
            if sname in sps_dict:
                sps_dict[sname].soap_strucs = {str(i): val for i, val in enumerate(parts) if val > 0.0}
                sps_updated.append(sname)
            else:
                nsp_read += 1
    return [sps_updated, nsp_read]


def read_soap_strucs_from_table(lines: list, sps_dict: dict) -> list:
    """Read SOAP structures from a table file."""
    sep = ","  # Separator
    nadd = 0  # Number of additional columns with numbers
    sps_updated, nsp_read = [], 0  # Counters
    for line in lines:
        if sep in line:
            # Find species name in the first column
            sname = line.split(sep, 1)[0]
            if not sname:
                raise ValueError(f"Empty species name found in line: {line}")
            # Find vector groups
            parts = [float(i) for i in line.split(sep) if isfloat(i)]
            if len(parts) != NGP_SOAP + nadd:  # Name + groups
                raise ValueError(f"Unrecognized line: {line} with {len(parts)} parts != {NGP_SOAP + nadd}.")
            # Record if the species is in the list
            if sname in sps_dict:
                sps_dict[sname].soap_strucs = {str(i): val for i, val in enumerate(parts[nadd:]) if val > 0.0}
                sps_updated.append(sname)
            else:
                nsp_read += 1
    return [sps_updated, nsp_read]


def mech_output_to_ssh(path: str, chem: str, reactions: list, species: list, fake_sps: dict) -> bool:
    """Return mechanism files in the SSH format"""

    # Settings
    ssh_opt = get_all_settings_map()[SNAME.SSH]()

    if not os.path.exists(path):
        os.makedirs(path)

    # Add fake species if needed
    if fake_sps:
        current_reactions = get_reactions_with_fake_species(reactions, fake_sps, "FA")
    else:
        current_reactions = reactions

    # Write reactions
    ro2s_rcn = set()  # Record RO2 species with RO2-RO2 reactions
    with open(os.path.join(path, f"{chem}.reactions"), "w", encoding="utf-8") as f:

        def gen_reaction_lines():
            for i, rcn in enumerate(current_reactions):
                if rcn.status <= 0:
                    continue
                yield rcn.to_rcn(i)
                if rcn.rate.ro2:  # Update used RO2 pool
                    ro2s_rcn.update(rcn.reactants)

        f.writelines(gen_reaction_lines())
        f.write("END\n")

    # Update species group id if needed
    if ssh_opt.soa_grps:
        _update_ssh_group_ids(species, ssh_opt.soa_grps)

    # Write species files
    aero_lines = get_default_ssh_aero_lines()
    valid_species = [isp for isp in species if isp.status > 0]  # Valid species

    # .mol file
    with open(os.path.join(path, f"{chem}.mol"), "w", encoding="utf-8") as f:
        for isp in valid_species:
            f.writelines(isp.output_to_lines())

    # .species file
    with open(os.path.join(path, f"{chem}.species"), "w", encoding="utf-8") as f:
        # Write header
        f.write(f"# Species list for {chem}\n")
        f.writelines(f"{isp.name}    {isp.mass:.2f}\n" for isp in valid_species)

        # Add fake species
        if fake_sps:
            f.writelines(f"FA{s}    {imw:.2f}\n" for s, imw in fake_sps.items())

    # .RO2 species
    ro2s = {isp.name: isp.RO2 for isp in valid_species if isp.RO2}
    ro2s_rcn = sorted(set(ro2s.keys()) & set(ro2s_rcn))
    with open(os.path.join(path, f"{chem}.RO2"), "w", encoding="utf-8") as f:
        f.write(f"# {len(ro2s_rcn)} out of # {len(ro2s)} RO2 species with RO2 reactions.\n")
        f.writelines(f"{s}    {ro2s[s]:d}\n" for s in ro2s_rcn)

    # .aer.vec file
    with open(os.path.join(path, f"{chem}.aer.vec"), "w", encoding="utf-8") as f:
        # Aerosol list
        f.write(aero_lines[0])
        f.writelines(to_ssh_aero_list(isp, "vectors") for isp in valid_species if isp.condensable)
        # Add water
        f.write(aero_lines[1])

    if ssh_opt.out_all_aeros:
        # .aer.1st and .aer file
        with open(os.path.join(path, f"{chem}.aer.1st"), "w", encoding="utf-8") as f0, open(
            os.path.join(path, f"{chem}.aer"), "w", encoding="utf-8"
        ) as f1:

            f0.write(aero_lines[0])
            f1.write(aero_lines[0])

            for isp in valid_species:
                if isp.condensable:
                    f0.write(to_ssh_aero_list(isp, "smiles", True))
                    f1.write(to_ssh_aero_list(isp, "smiles"))

            # Add water
            f0.write(aero_lines[1])
            f1.write(aero_lines[1])

    return True


def _update_ssh_group_ids(species: list, group_dict: dict) -> None:
    """Update grp_ids for aerosol species based on the precursor group ids"""

    if not group_dict:
        logger.warning("No SSH aerosol groups defined. Skipping group id update.")
        return

    counts = {k: 0 for k in group_dict}  # Initialize counts for each group
    outs = set()  # Species not in the group_dict

    for s in species:
        if s.status <= 0 or not s.condensable:
            continue
        igrp = ",".join(sorted(s.precursors))
        if igrp in group_dict:
            s.grp_id = group_dict[igrp]
            counts[s.grp_id] += 1
        else:
            s.grp_id = None
            outs.add(s.name)

    logger.info("Updated group_id for %s & ungroupd # %d condensables: %s", counts, len(outs), outs)
