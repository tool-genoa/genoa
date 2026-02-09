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
This module contains functions to read GECKO-A mechanisms.
"""

import os

import attrs

from .constants import SNAME, NORATE_STR
from .gecko_cst import NPERO, LENSP, GTYPES, GRO2_GRP, GCK_TO_SSH_KEYW, EXTRA_LABEL, GCK_TO_SMILE
from .logger import setup_logger
from .reaction import Reaction
from .setting_global import get_all_settings_map, get_attrs_from_settings
from .smiles import string_to_smiles, get_properties_from_smiles
from .species import Species, add_basic_to_species_list
from .utils import isfloat

# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class GCKFileNames:
    """GECKO-A file names."""

    path: str
    rcn: str = "reactions.dum"
    sps: str = "dictionary.out"
    psat: str = "pvap.nan.dat"  # VP0BP0
    primary: str = "listprimary.dat"
    ro2: str = "pero.ro2"
    tg: str = "Tg.dat"
    mdv: str = "vdiffusion.dat"
    ro2s: list = attrs.field(factory=list)

    def __attrs_post_init__(self):
        """Check and initialize GECKO-A file names from a path."""
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"Path {self.path} does not exist.")
        if not os.path.isdir(self.path):
            raise NotADirectoryError(f"Path {self.path} is not a directory.")
        files_in = os.listdir(self.path)

        # Map the attribute names to their file names
        attr_to_file = {
            "rcn": self.rcn,
            "sps": self.sps,
            "psat": self.psat,
            "primary": self.primary,
            "tg": self.tg,
            "mdv": self.mdv,
        }

        # Check if all files exist
        for attr, f in attr_to_file.items():
            if f not in files_in:
                raise FileNotFoundError(f"File {f} does not exist.")
            # Update the attribute with the full path to the file
            setattr(self, attr, os.path.join(self.path, f))

        # Check ro2 files
        ro2_files = [f"pero{i}.dat" for i in range(1, NPERO + 1)]
        if self.ro2 in files_in:
            self.ro2s = [os.path.join(self.path, self.ro2)]
        elif all(f in files_in for f in ro2_files):
            self.ro2s = [os.path.join(self.path, f) for f in ro2_files]
        else:
            raise FileNotFoundError("RO2 files do not exist.")


def process_gck_reactions_and_species(geckopath: str) -> list:
    """Processes reaction and species lists from a GECKO-A path."""

    # Get settings
    remove_empty_ro2 = False  # Remove reactions with empty ro2 groups
    settings_map = get_all_settings_map()
    gnl = settings_map[SNAME.SETUP]()
    nopt = settings_map[SNAME.NEW]()

    # Get gck file names
    gckfile = GCKFileNames(path=geckopath)
    logger.info("Reading GECKO-A mechanism from %s with files: %s", geckopath, gckfile)

    _read_and_check_primary_vocs(gckfile.primary, gnl.primary_vocs)

    # Read species list
    species, sps_dict = _read_and_check_species_dictionary(gckfile.sps)

    # Read RO2 info
    if len(gckfile.ro2s) == 1:
        ro2_counts = _update_ro2_info_from_one_file(gckfile.ro2s[0], sps_dict)
    elif len(gckfile.ro2s) == NPERO:
        ro2_counts = _update_ro2_info_from_files(gckfile.ro2s, sps_dict)
    else:
        raise ValueError(f"no. ro2 files should be 1 or {NPERO}. Read {len(gckfile.ro2s)}.")

    logger.info("%s", "\n".join([f"Read # {j} RO2s for pero{i+1}." for i, j in enumerate(ro2_counts)]))

    _check_ch3o2(sps_dict, ro2_counts)  # Check CH3O2

    # Read aerosol files
    naero = _read_aerosol_files(gckfile, sps_dict)
    logger.info("Read # %d aerosols from %s, %s, %s.", naero, gckfile.psat, gckfile.tg, gckfile.mdv)

    # Update and check species list
    if nopt.update_from_smiles:
        logger.info("Checking and updating species properties from SMILES ...")
        _check_and_update_species_list(gckfile, sps_dict)

    # Get empty ro2 pool
    if remove_empty_ro2:
        rm_ro2 = [f"{GRO2_GRP[i] }" for i, j in enumerate(ro2_counts) if j == 0]
    else:
        rm_ro2 = None

    # Read reaction list
    reactions, sps_dict_order = _read_reaction_file(gckfile.rcn, sps_dict, rm_ro2)

    # Reorder species based on their appearance in reactions
    species.sort(key=lambda s: sps_dict_order[s.name])

    return [reactions, species]


def gecko_to_smiles(gecko: str, tag_canonical: bool = True) -> str:
    """Convert input GECKO-A format to its canonical (optional) smiles structure"""

    gecko_smiles = {}  # specific smiles

    # check if gecko in certain format
    if gecko in gecko_smiles:
        return gecko_smiles[gecko]

    # rank gecko groups from the longest to shortest
    gecko_group = sorted(GCK_TO_SMILE, key=len, reverse=True)

    # replace gecko groups step by step
    for i in gecko_group:
        if i in gecko:
            gecko = gecko.replace(i, GCK_TO_SMILE[i])

    # check string and convert to canonical smiles if need
    return string_to_smiles(gecko, tag_canonical)


def gecko_to_ssh_kinetic_rate(rate_str: str) -> str:
    """Convert a kinetic rate string from GECKO format to SSH format."""

    # Update rate_type and rate_info
    if "T:" in rate_str:
        tmp = rate_str.split("T:", 1)[1].split(" ")
        # Get kinetic type
        rate_type = tmp[0].strip()
        # Get kinetic index and ratio
        rate_info = [i for i in tmp[1:] if i != ""]
    else:
        rate_type, rate_info = None, []

    # Get arrhenius's law coefficients
    arrs3 = []
    for i in rate_str.split(" "):
        if not isfloat(i):
            continue
        arrs3.append(i.strip())
        if len(arrs3) == 3:
            break
    # Build ssh kinetic
    parts = ["KINETIC"]

    # Check type
    if rate_type:
        if rate_type in GCK_TO_SSH_KEYW:
            parts.append(GCK_TO_SSH_KEYW[rate_type])
        elif rate_type in GRO2_GRP:  # Check RO2 pool index
            parts.append(f"{GCK_TO_SSH_KEYW['PERO']} {GRO2_GRP.index(rate_type)+1}")
        else:
            raise NameError(f"Not found kinetic rate type: {rate_type}")

    # Photolysis
    if rate_type == "HV":  # no arrhenius(3) is used
        # Check read len
        if len(rate_info) != 2:
            raise ValueError(f"Invalid parameters for 'HV' type: {rate_info}", rate_str)

        parts.extend(map(str, rate_info))  # index and ratio

    # FALL OFF with 11 paramters
    elif rate_type == "FALLOFF":  # arrhenius(3) + FALLOFF(7) + ratio = 11
        # Check read len
        if len(rate_info) != 7:
            raise ValueError(f"Invalid parameters for 'FALLOFF' type: {rate_info}")

        # Write 11 params
        parts.extend(map(str, arrs3 + rate_info + [1.0]))

    elif rate_type == "EXTRA":  # arrhenius(3) + EXTRA(?) + ratio = ?

        # Check type
        if int(rate_info[0]) not in EXTRA_LABEL:
            raise NameError(f"READ EXTRA, type not found: {rate_info[0]}")
        # write type
        parts.append(str(rate_info[0]))

        # Write the rest of params
        parts.extend(map(str, arrs3 + rate_info[1:] + [1.0]))

    elif rate_type == "ISOM":  # as EXTRA 200

        # Check read len isom(5)
        if len(rate_info) != 5:
            raise ValueError("READ ISOM: ", rate_str, len(rate_info), rate_info)

        # specific case: if rate info is [0,0,0,0,1], then it is a simple rate
        rate_val = [float(i) for i in rate_info]
        if rate_val == [0.0, 0.0, 0.0, 0.0, 1.0]:
            return "KINETIC ARR " + " ".join(map(str, arrs3))

        # Build rate to 11 params: '200' + arrhenius(3) + isom(5) + 1.0 = 10
        parts.extend(map(str, arrs3 + rate_info + [1.0]))

    else:
        # Write 3 params
        parts.extend(map(str, ["ARR"] + arrs3))

    return " ".join(parts)


def _read_aerosol_files(gckfile: str, sps_dict: dict) -> int:
    """Read psay, Tg, mdv files. Return number of aerosols."""

    # Read Psat
    with open(gckfile.psat, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    naero = 0  # count for SVOCs
    for line in lines:
        info = [i for i in line.split(" ") if i != ""]
        # In format species name, Past, evap
        if len(info) == 3:
            # Remove "G" from species name
            # e.g., "GGO03001" -> "GO03001"
            sname = info[0][1:]
            if sname not in sps_dict:
                raise NameError(f"Condensable {sname} from {gckfile.psat} not in {gckfile.sps}.")
            # Get speices
            isp = sps_dict[sname]
            # Update aero info
            isp.condensable = True
            isp.psat_atm = float(info[1])
            isp.dhvap_kj = float(info[2])
            naero += 1

    # Read Tg from Tg.dat
    icount = 0  # Count for Tg
    with open(gckfile.tg, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    for line in lines:
        info = [i for i in line.split(" ") if i != ""]
        if len(info) == 2:  # In format species name, Tg
            sname = info[0][1:]
            if sname not in sps_dict:
                raise NameError(f"Aerosol {sname} from {gckfile.tg} is not in {gckfile.sps}.")
            sps_dict[sname].tg = float(info[1])
            icount += 1
    # Check if all aerosols have Tg
    if icount != naero:
        raise ValueError(f"Find # {icount} Tg for aerosols but # {naero} aerosols.")

    # Read mdv from mdv.dat -> vdiffusion.dat
    icount = 0  # Count for mdv
    with open(gckfile.mdv, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    for line in lines:
        info = [i for i in line.split(" ") if i != ""]
        if len(info) == 2:
            sname = info[0][1:]
            if sname not in sps_dict:
                raise NameError(f"Aerosol {sname} from {gckfile.mdv} is not in {gckfile.sps}.")
            sps_dict[sname].dvol = float(info[1])
            icount += 1
    # Check if all aerosols have mdv
    if icount != naero:
        raise ValueError(f"Find # {icount} mdv for aerosols but # {naero} aerosols.")
    return naero


def _read_and_check_primary_vocs(filename: str, primary_set: list) -> None:
    """Read primary VOCs from a file and check if it is the same as the input set."""

    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    new_vocs = []
    for line in lines:
        info = [i for i in line.strip().split(" ") if i != ""]
        if len(info) == 2:
            pvoc = info[0].strip()
            if pvoc not in primary_set:
                raise ValueError(f"{pvoc} in {filename} but not in {primary_set}.")
            new_vocs.append(pvoc)

    for s in primary_set:
        if s not in new_vocs:
            raise ValueError(f"{s} in {primary_set} but not in {filename}.")


def _get_new_species(info: list, basic_sps: dict) -> Species():
    """Return new species generated based on input gecko-a info line"""
    n = len(info)

    # Check if valid line
    # 12 basic columns + note (3rd col.)  + no.generation (NCAR model, 4th col.)
    if n < 12 or n > 14:
        raise ValueError(f"Not 12 <= n <= 14. Read {n} for {info} in species file.")

    # Build new species with GECKO structure (2nd col.) and status
    isp = Species(name=info[0], string=info[1], status=2 if n == 12 else 1)

    # GECKO-A simplified structure (3rd col.)
    if n == 14:
        # With no.generation (4th col., NCAR version)
        isp.note = info[2] = f"{info[2]}_{info[3]}"
        ind = 3  # Index for info reading
    # Not with GECKO-A simplified structure
    elif n == 12:
        ind = 1
    else:
        isp.note = info[2]
        ind = 2

    # Molar mass (4th col.)
    isp.mass = float(info[ind + 1])

    if isp.mass <= 0.0:
        isp.status = 2
        # Update mass
        if isp.name in basic_sps:
            isp.mass = basic_sps[isp.name]
        else:
            raise ValueError(f"Read basic species {isp.mass} w/o mass.")
        return isp

    # Radical (RO2 is assigned later)
    if "O." in isp.string:
        isp.radical = True

    # Radical or not 1/0 (5th col.)
    if isp.radical and int(info[ind + 2]) == 0:
        raise ValueError(f"Radical {isp.name} with invalid string: {isp.string}. status: {isp.radical}. Read: {info}")

    # GECKO string to smiles
    isp.smiles = gecko_to_smiles(isp.string, tag_canonical=not isp.radical)

    # Decomposition: number of atoms in the order of
    # C, H, N, O, S, F, Cl, Br, and . (6th to end col.)
    fgroups = [int(i) for i in info[ind + 3 :]]
    # Not support for S, F, Cl, Br yet
    if sum(fgroups[4:]) > 0:
        raise ValueError("Decomposition not support for f{isp.name}. Read: {fgroups} from:\n{info}")

    # Get fgroups, formula, ratios
    isp.fgroups = dict(zip("CHNO", fgroups[:4]))
    isp.update_by_functional_group()

    return isp


def _read_and_check_species_dictionary(filename: str) -> list:

    # Save species objects and name for index
    species, sps_dict = [], {}
    logger.info("Reading species list from %s ...", filename)

    # Get settings
    basic_dict = get_attrs_from_settings({SNAME.SETUP: "basicspecies_dict"})

    # Example:
    # 1N3009   CH2(OH)CH(O.)CH2(ONO2) 1.NO 3  136.1  1  3  6  1  5  0  0  0  0
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    # Remove initial and END lines
    for line in lines[1:-1]:
        # Read line and split to columns
        info = [i for i in line.split(" ") if i != ""]

        # Check if string starts with "#"
        if info[1].startswith("#"):
            logger.info("Find compound %s with special formula: %s", info[0], info[1])
            raise ValueError("Species with special formula not supported yet.")
            # info[1] = info[1][1:]  # Remove "#"

        # Get new species
        isp = _get_new_species(info, basic_dict)
        if isp.name in sps_dict:  # If already record
            raise NameError(f"Repeat {isp.name} in {filename}")

        # Add to species list
        species.append(isp)
        sps_dict[isp.name] = isp
        # Add basic species
        if isp.status == 2:
            if isp.name not in basic_dict:
                logger.info("Adding basic species %s from predefined list ...", isp.name)
            basic_dict[isp.name] = isp.mass

    # Update species list with basic species
    add_basic_to_species_list(species, sps_dict, basic_dict)
    logger.info("Read # %d species from # %d lines.", len(species), len(lines))

    return [species, sps_dict]


def _update_ro2_info_from_one_file(filename: str, sps_dict: dict) -> list:
    """Read ro2 name and group id from a file (.ro2) and return number of ro2 groups."""

    num_counts = [0] * NPERO
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    # Read counts from header
    nums_read = [int(i) for i in lines[0].split("!")[0].split(" ") if isfloat(i)]
    if len(nums_read) != NPERO:
        raise ValueError(f"Find invalid header {lines[0]} in ro2s file: {filename}.")
    # Read species
    for line in lines[1:-1]:
        info = [i for i in line.split(" ") if i != ""]
        if len(info) != 2 or not info[0].startswith("G"):
            raise ValueError(f"Find invalid species line {line} in file: {filename}.")
        # Get species name without "G"
        sname, ro2_id = info[0][1:], int(info[1])
        if sname not in sps_dict:
            raise NameError(f"RO2 {sname} from {filename} not in species list.")
        if not sps_dict[sname].radical:
            raise ValueError(f"Species {sname} is not a radical in species list.")
        sps_dict[sname].RO2 = ro2_id
        num_counts[ro2_id - 1] += 1

    # Check
    for i, icount in enumerate(num_counts):
        if icount != nums_read[i]:
            raise ValueError(f"Find # {icount} RO2s for pero{i+1}. Should be: {nums_read[i]}.")

    return num_counts


def _update_ro2_info_from_files(filenames: list, sps_dict: dict) -> list:
    """read ro2 infor from multiple files (pero[i].dat)"""
    if len(filenames) != NPERO:
        raise ValueError(f"Find {len(filenames)} files but {NPERO} files are needed.")

    nums_counts = []
    for ind, filename in enumerate(filenames):
        with open(filename, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        # Read counts from header
        num_read = int(lines[0].split("!")[0])
        # Read species
        count, ro2_id = 0, ind + 1
        for line in lines[1:-1]:
            info = line.strip()
            if not info.startswith("G"):
                raise ValueError(f"Find invalid species line {line} in file: {filename}.")
            # Get species name without "G"
            sname = info[1:]
            if sname not in sps_dict:
                raise NameError(f"RO2 {sname} from {filename} not in species list.")
            if not sps_dict[sname].radical:
                raise ValueError(f"Species {sname} is not a radical in species list.")
            sps_dict[sname].RO2 = ro2_id
            count += 1
        # Check
        if count != num_read:
            raise ValueError(f"Count # {count} RO2s for {ro2_id}. Read {num_read}.")

        nums_counts.append(count)

    return nums_counts


def _check_ch3o2(sps_dict: dict, ro2_counts: list) -> None:
    """Check ch3o2 and MEPERO"""
    # Index in gecko ro2 group
    n = GRO2_GRP.index("MEPERO")
    if "CH3O2" in sps_dict:
        sps_dict["CH3O2"].RO2 = n + 1
    else:
        logger.warning("CH3O2 not found in species list.")
        if ro2_counts[n] > 0:
            raise ValueError("CH3O2 not found in species list but MEPERO is not empty.")


def _read_reaction_file(filename: str, sps_dict: dict, rm_ro2: list) -> list:
    """
    Read reaction file and update species and reactions.
    """

    # Get settings
    basic_dict = get_attrs_from_settings({SNAME.SETUP: "basicspecies_dict"})

    # Counts for empty pero, wall loss, index of checked line, gas-particle reactions
    n_empty_ro2, nwall, ichecked = 0, 0, -1

    # Basic species dict that can be updated
    basic_all = set(basic_dict.keys())

    # Order of species in reactions - default large number
    sps_dict_order = {s.name: 1e9 for s in sps_dict.values()}

    # Read reactions from reaction file
    logger.info("Reading reactions from %s ...", filename)

    reactions, pinfos = [], {"rcn": [], "sps": []}
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    for iline, line in enumerate(lines):
        # Skip checked lines
        if iline <= ichecked:
            continue

        # Get block of reactions in case the line is not finished (end with +)
        block = line.strip()
        ichecked = iline

        # Skip comments
        if not block or any(block.startswith(i) for i in ("!", "END", "REACTIONS")):
            continue

        # Check if need to read more lines for a whole block
        if block.endswith("+"):
            blocks, current_block = [block], block
            while current_block.endswith("+"):
                ichecked += 1
                current_block = lines[ichecked].strip()
                blocks.append(current_block)
            block = "".join(blocks)

        # Check if block is valid
        if "=>" not in block:
            raise ValueError(f"Find invalid block with lines # {iline} - {ichecked}: {block}")

        # Check if need to skip reaction
        # Gas-to-particle reactions
        if any(i in block for i in ["AIN ", "AOU "]):
            # sname = block.split("+")[0][1:]
            continue
        # Wall loss reactions
        if any(i in block for i in ["WIN ", "WOU "]):
            # Only update with win
            if "WOU " in block:
                continue
            sname = block.split("+", 1)[0][1:]
            if sname not in sps_dict:
                raise ValueError(f"Find {sname} in {block} but not in species list.")
            sps_dict[sname].wall_loss = True
            nwall += 1
            continue
        # RO2-RO2 reactions with empty pero groups
        if rm_ro2 and any(i in block for i in rm_ro2):
            n_empty_ro2 += 1
            continue

        # Get treatable elements
        parts = []
        for i in block.replace("=>", " => ").split(" "):
            if i == "":
                continue
            # Find 'GNO2+GNO3', not '4.670E+14' !
            if "+" in i and not isfloat(i):
                # Add species
                for j in i.split("+"):
                    if j != "":
                        parts.append(j)
            else:
                parts.append(i)

        # Index to separate reactants and products
        ind = parts.index("=>")

        # Build new reaction
        rcn, rtype, tag_basic = Reaction(), None, True

        # Reactants
        for i in parts[:ind]:
            # Find species name
            if i[0] == "G":
                s = i[1:]  # Species name
                if s not in sps_dict:
                    raise ValueError(f"Find {s} in {block} but not in species list.")
                rcn.reactants[s] += 1  # Add reactants
                # Check if all basic species
                if tag_basic and s not in basic_all:
                    tag_basic = False

                # Update order of species in reactions
                if sps_dict_order[s] > len(sps_dict_order):
                    sps_dict_order[s] = len(sps_dict_order)

            # Get reaction type
            else:
                # Unrecognized type
                if i not in GTYPES:
                    raise NameError(f"Find new type {i} from reaction: {block}")
                # Multiple reaction type
                if rtype:
                    raise NameError(f"Find multiple rtype: {rtype} for reaction {block}")
                rtype = i

        # Update reaction status
        if tag_basic:
            pinfos["rcn"].append(block)
            rcn.status = 2

        # Products
        rt = 1.0  # default ratio
        for i in parts[ind + 1 : -3]:  # Remove arrhenius coefs
            if i == "NOTHING":  # No products
                continue
            # Find number: product ratio
            if isfloat(i):
                rt = float(i)
            # Find species name for products (gas)
            elif i[0] == "G":
                s = i[1:]
                if s not in sps_dict:
                    raise ValueError(f"Find {s} in {block} but not in species list.")
                rcn.products[s] += rt
                # Update basic species
                if tag_basic and s not in basic_all:
                    pinfos["sps"].append(s)
                    # Update status
                    sps_dict[s].status = 2
                    logger.info("Basic species %s from reaction not in the list: %s", s, block)
                    basic_all.add(s)
                rt = 1.0  # Reset ratio

                # Update order of species in reactions
                if sps_dict_order[s] > len(sps_dict_order):
                    sps_dict_order[s] = len(sps_dict_order)
            else:
                raise ValueError(f"Find unknown product {i} w/o 'G' in {block}.")

        # Kinetic rate
        rate_parts = parts[-3:]
        # Read more paramters to the string if need
        if rtype:
            rate_parts.append(f"T:{rtype}")
            # Get more parameters from the next line
            if rtype in ["HV", "FALLOFF", "EXTRA", "ISOM"]:
                ichecked += 1  # Read next line
                rate_parts.append(lines[ichecked].replace("/", "").replace(rtype, "").strip())
        rate_str = " ".join(rate_parts)

        # Process string
        rcn.rate.string = rate_str
        rcn.rate.ssh = gecko_to_ssh_kinetic_rate(rate_str)
        # Check ssh string
        if NORATE_STR in rcn.rate.ssh:
            raise ValueError(f"Find {NORATE_STR} in ssh string: {rcn.rate.ssh}")
        if not rcn.rate.init_rate_with_string():
            raise ValueError(f"Cannot initialize rate with string: {rcn.rate.string} from {block}")

        # Update reaction list
        reactions.append(rcn)

    logger.info(
        "Read # %d reaction(s) from # %d lines (%d basic reaction(s) & %d basic species).",
        len(reactions),
        len(lines),
        len(pinfos["rcn"]),
        len(basic_all),
    )
    # logger.info(_get_reaction_pinfo(n_empty_ro2, rm_ro2, nwall, pinfos, basic_all))

    return [reactions, sps_dict_order]


def _get_reaction_pinfo(n_empty_ro2: int, rm_ro2: list, nwall: int, pinfos: dict, basics: set) -> str:
    """Get reaction information for recording"""
    parts = []

    if n_empty_ro2 > 0:
        parts.append(f" Skipped # {n_empty_ro2} reaction(s) with empty ro2 groups: {rm_ro2}.")

    if nwall > 0:
        parts.append(f" Updated # {nwall} wall loss reaction(s).")

    if pinfos["rcn"]:
        parts.append(f"\nRead # {len(pinfos['rcn'])} basic reaction(s):")
        parts.extend(pinfos["rcn"])

    if pinfos["sps"]:
        parts.append(f"\nAdded # {len(pinfos['sps'])} basic specie(s) from reaction list:")
        parts.extend(pinfos["sps"])

    parts.append(f"\nTotal # {len(basics)} basic species:")
    parts.extend(f"  {sp}" for sp in sorted(basics))

    return "\n".join(parts)


def rename_species(reactions: list, species: list, lsps: int = LENSP - 1) -> None:
    """Rename species in both reactions and species list to a shorter name if needed."""

    sps_dict, sps_new, n = {}, {}, 0
    for s in species:
        if not s.status:
            continue
        if len(s.name) > lsps:  # Rename species
            find_new = False
            for _ in range(10):  # Try 10 times
                sname = f"R{n}"
                sname += s.name[-(lsps - len(sname)) :]  # Add last part
                if sname not in sps_dict:
                    find_new = True
                    logger.info("Rename %s to %s due to length limit.", s.name, sname)
                    break
                n += 1  # Update index
            if not find_new:
                raise ValueError(f"Cannot find new name for {s.name} in species list.")
            if len(sname) > lsps:
                raise ValueError(f"Species name {sname} for {s.name} exceeds limit {lsps}.")
            # Update new name
            sps_dict[sname] = s
            sps_new[s.name] = sname
            s.name = sname  # Update name in species list

    if not sps_new:
        logger.info("No species need to be renamed.")
        return

    logger.info("Renamed in total # %d species.", len(sps_new))

    # Update species in reactions
    nrcn = 0
    for rcn in reactions:
        csps = [s for s in rcn.reactants + rcn.products if s in sps_new]
        if csps:
            nrcn += 1
            for s in csps:
                if s in rcn.reactants:
                    rcn.reactants[sps_new[s]] = rcn.reactants.pop(s)
                if s in rcn.products:
                    rcn.products[sps_new[s]] = rcn.products.pop(s)
            rcn.update_rcn_id()
    logger.info("Renamed species in # %d reactions.", nrcn)


def _check_and_update_species_list(gckfile: GCKFileNames, sps_dict: dict) -> None:
    """Check and update species list with properties computed by its SMILES string."""

    # Compute properties from SMILES
    props = list(get_properties_from_smiles("CC", "VP0BP0").keys())
    props += ["tg", "dvol"]  # Add Tg and dvol
    if "fgroups" not in props:
        raise ValueError(f"fgroups not found in properties. Got: {props}")
    ifg = props.index("fgroups")  # Index of functional groups

    def _update_species_and_out_diff():
        """Update species and output the differences from read values"""
        for s, isp in sps_dict.items():
            if not isp.smiles or isp.smiles == "-":
                continue

            props_sml = get_properties_from_smiles(isp.smiles, "VP0BP0")
            isml = [props_sml.get(k, None) for k in props[:-2]]  # exclude Tg, dvol

            # Update functional groups
            fg_inferred = isml[ifg]
            if fg_inferred:
                isp.fgroups.update({k: v for k, v in fg_inferred.items() if k not in isp.fgroups})

            iread = [getattr(isp, prop, None) for prop in props]

            # Compute Tg and dvol if needed
            if isp.tg is None:
                isml.extend([None, None])
            else:
                isp.update_prop_wgck()
                isml.extend([isp.tg, isp.dvol])
                if not iread[-2:] == isml[-2:]:
                    logger.warning("Tg or dvol for %s w/ %s not match: %s, %s", s, isp.smiles, iread[-2:], isml[-2:])

            yield f"{s}\t{isp.smiles}\t" + "\t".join(f"{i},{j}" for i, j in zip(isml, iread)) + "\n"

    file_smiles = f"{gckfile.sps}.genoa"
    logger.info("Saving species properties from SMILES to %s ...", file_smiles)
    with open(file_smiles, "w", encoding="utf-8") as f:
        # Write header
        f.write("Species\tSMILES\t" + "\t".join(f"{s}_read,{s}_computed" for s in props) + "\n")
        f.writelines(_update_species_and_out_diff())

    logger.info("Update species properties from SMILES => %s.genoa", gckfile.sps)
