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
This module contains utility functions for handling mechanisms, conditions, and data processing.
"""

import os
import re
from typing import Any, Optional, Union

import numpy as np
from .logger import setup_logger


# Logger
logger = setup_logger(__name__)


def isfloat(value: Any) -> bool:
    """Check if the input is a float"""
    try:
        float(value)
        return True
    except ValueError:
        return False


def get_condition_list(path_cond: str) -> list:
    """output path and a list of condition filenames used for GECKO-A simulations."""

    if not os.path.isdir(path_cond):  # Check if is folder
        raise FileNotFoundError(f"Cannot find condition folder: {path_cond}.")

    # Get condition files
    conds = [i for i in os.listdir(path_cond) if os.path.isfile(os.path.join(path_cond, i))]
    if len(conds) == 0:
        raise FileNotFoundError(f"Cannot find any file in {path_cond}.")
    conds.sort()  # Sort by name
    # list_to_log(logger, conds, f"Find {len(conds)} condition(s) in {path_cond}:", "\t")

    return [path_cond] + conds


def get_files_from_folder(folderpath: str, filename: str) -> list:
    """Return a list of file paths with the given filename  obtained from the given folder or subfolders."""

    # Check folder path: should be Results folder
    if not (filename and folderpath and os.path.isdir(folderpath)):
        raise FileNotFoundError(f"{folderpath} is not a directory or filename {filename} is empty.")

    folderin = os.listdir(folderpath)
    fileout = [os.path.join(folderpath, f) for f in folderin if filename in f]

    # File in the folder - no search in subfolders
    if fileout:
        return fileout

    # Check if in subfolders
    for subfolder in folderin:
        ipath = os.path.join(folderpath, subfolder)
        # Result folder only
        if not os.path.isdir(ipath):
            continue
        # Check if filename in subfolder
        ifind = [f for f in os.listdir(ipath) if filename in f]
        if ifind:
            fileout.extend([os.path.join(ipath, f) for f in ifind])
    return fileout


def load_common_conditions(paths: list) -> list:
    """Get common conditions from given result paths."""

    logger.info("Loading common conditions from # %d result paths ...", len(paths))
    conds = None  # Common conditions
    for path in paths:
        if not (path and os.path.isdir(path)):
            raise FileNotFoundError(f"Path not found or not a folder: {path}")
        icond = [i.strip() for i in os.listdir(path) if ".T" in i]
        if not icond:
            raise FileNotFoundError(f"No conditions containing '.T' found in {path}")
        if conds is None:
            conds = icond
        else:
            conds = [i for i in conds if i in icond]

    if not conds:
        logger.warning("No common conditions found from paths: %s", paths)
        return []

    logger.info("Found # %d common conditions: %s.", len(conds), conds)
    return conds


def get_data_from_concs_files(filenames: list, mode: str = "ave") -> np.ndarray:
    """Get concentration from simulation results"""
    concs, x0, y0 = [], None, None
    for filename in filenames:
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"File not found: {filename}")
        with open(filename, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        x = len([i for i in lines[0].split(" ") if i.strip()])
        y = len(lines) - 1
        if x0 is None:
            x0, y0 = x, y
        elif x0 != x or y0 != y:
            raise ValueError(f"File {filename} has different dimensions: {x}x{y} vs {x0}x{y0}")
        # Get concentrations
        iconc = np.zeros((x, y))
        for i, line in enumerate(lines[1:]):
            iconc[:, i] = [float(j) for j in line.split(" ") if j.strip()]
        concs.append([get_data_based_on_mode(i, mode) for i in iconc])

    return np.array(concs)  # (n, x) or (n, x, y)


def is_mech_exist(path: str, chem: str, with_folder: bool = False, mode="GENOA") -> bool:
    """Check mechanism file exists or not."""

    if with_folder:  # Update path with folder
        path = os.path.join(path, chem)

    key_filename_dict = {
        "SSH": f"{chem}.species",
        "GCK": f"{chem}.mech",
        "GENOA": f"{chem}.mol",
    }

    if mode not in key_filename_dict:
        raise ValueError(f"Unknown mode: {mode}")
    filename = os.path.join(path, key_filename_dict[mode])
    if not os.path.exists(filename):
        return False

    return True


def get_mech_files(path: str, chem: str, with_folder: bool = False) -> list:
    """Return reaction and species files of a mechanism."""

    if with_folder:  # Update path with folder
        path = os.path.join(path, chem)

    # Read reaction and species files
    mech_files = [f"{chem}.reactions", f"{chem}.mol"]
    for i, filename in enumerate(mech_files):
        if not os.path.exists(os.path.join(path, filename)):
            raise FileNotFoundError(f"Can not find {filename} in {path}")
        # Update file name
        mech_files[i] = os.path.join(path, filename)

    return mech_files


def get_info_from_mech(reactions: list, species: list) -> dict:
    """Output a size dictionary from input lists."""

    rsize, ssize = count_reactions(reactions), count_species(species)
    return {
        "nrcn": rsize[0],
        "ntcn": rsize[0] + rsize[1],
        "nsps": ssize[0] + ssize[1],
        "nrsp": ssize[0],
        "naero": ssize[2],
    }


def get_info3_from_mech(reactions: list, species: list) -> list:
    """Output No. reducible reactions, species, condensables from input lists."""

    size_dict = get_info_from_mech(reactions, species)
    return [size_dict["nrcn"], size_dict["nrsp"], size_dict["naero"]]


def count_reactions(reactions: list) -> list:
    """Count the number of reactions: [reducible, basic, invalid]."""

    reducible, basic, invalid = 0, 0, 0
    for r in reactions:
        if r.status == 1:
            reducible += 1
        elif r.status == 2:
            basic += 1
        elif r.status == 0:
            invalid += 1
        else:
            raise ValueError(f"Find rcn {r.to_rcn()} with invalid status: {r.status}")

    if len(reactions) != reducible + basic + invalid:
        sum_rcn = reducible + basic + invalid
        raise ValueError(
            f"Counting error in reactions: sum: {sum_rcn} != total: {len(reactions)}",
            f"Got: {reducible} reducibles, {basic} basics, {invalid} invalids",
        )

    return [reducible, basic, invalid]


def count_species(species: list) -> list:
    """Count the number of species: [reducible, basic, aerosol, invalid]."""

    reducible, basic, aerosol, invalid = 0, 0, 0, 0
    for s in species:
        if s.status == 1:
            reducible += 1
            if s.condensable:  # Only reducible aerosol
                aerosol += 1
        elif s.status == 2:
            basic += 1
        elif s.status == 0:
            invalid += 1
        else:
            raise ValueError(f"Find sp {s.to_ssh_spc()} with invalid status: {s.status}")

    if len(species) != reducible + basic + invalid:
        sum_sps = reducible + basic + invalid
        raise ValueError(
            f"Counting error in species: sum: {sum_sps} != total: {len(species)}",
            f"Got: {reducible} reducibles, {basic} basics, {aerosol} aerosols, {invalid} invalids",
        )
    return [reducible, basic, aerosol, invalid]


def get_mechanism_size_info(reactions: list, species: list) -> str:
    """Obtain mechanism size information"""

    # Get size info
    rcn_size, sps_size = count_reactions(reactions), count_species(species)

    # Get record info
    pinfos = [f"Reactions: {len(reactions)}, Species: {len(species)}\nSize Counts:"]

    # Reaction size
    i, j, k = rcn_size
    pinfos.append(f"No. reactions: {len(reactions)} containing:  {i} reducibles; {j} basics; {k} invalids")

    # Species size
    i, j, k, l = sps_size
    pinfos.append(f"No. species: {len(species)} containing:  {i} reducibles; {j} basics; {k} aerosols; {l} invalids")

    return "\n".join(pinfos)


def size_difference(list1: list, list2: list, mode: str = "abs") -> list:
    """compute the size difference between two lists"""

    if len(list1) != len(list2):
        raise ValueError("Two lists must have the same length.")
    if mode == "abs":
        return [abs(list1[i] - list2[i]) for i in range(len(list1))]
    if mode == "rel":
        return [abs(list1[i] - list2[i]) / list1[i] for i in range(len(list1))]
    if mode == "rel100":
        return [abs(list1[i] - list2[i]) / list1[i] * 100 for i in range(len(list1))]
    if mode == "fb":
        return [abs(list1[i] - list2[i]) / (list1[i] + list2[i]) for i in range(len(list1))]

    raise ValueError(f"Not recognize mode {mode}.")


def get_size_diff_pinfo(before: list, after: list) -> str:
    """Record the size difference between two lists."""
    sdiff = size_difference(before, after, "rel100")
    sum_sdiff = sum(sdiff)
    if sum_sdiff > 0.0:
        info = [f"Before: {', '.join([str(i) for i in before])}\nAfter: {', '.join([str(i) for i in after])}"]
        info.append(f"Diff: {', '.join([str(i) for i in size_difference(before, after)])}")
        info.append("Rel Diff in %: " + "%, ".join([str(round(i, 3)) for i in sdiff]) + "%")
        return "\n".join(info)
    if sum_sdiff == 0.0:
        return "No size difference."
    raise ValueError(f"Negative sum of size difference: {sum_sdiff}")


def calculate_normalized_ratios(f_list: list, iround: Optional[int] = None, limit: float = 0.0) -> list:
    """Calculate normalized ratios of each element in a float list relative to the sum of all elements."""

    # No value or only one value
    n = len(f_list)
    if n == 0:
        return []
    if n == 1:
        return [1.0]

    # Compute ratios
    total = sum(f_list)
    if total <= 0.0:
        return [1.0 / n] * n

    ratios = [i / total for i in f_list]

    # Remove small values if limit is specified
    if limit > 0.0:
        ratios = [i if i >= limit else 0.0 for i in ratios]
        total = sum(ratios)
        if total <= 0.0:
            return [1.0 / n] * n
        ratios = [i / total for i in ratios]

    # Optional rounding and renormalization
    if iround is not None:
        ratios = [round(i, iround) for i in ratios]
        ratios = [max(i, 1e-18) for i in ratios]  # Avoid division by zero
        total = sum(ratios)
        if total <= 0.0:
            return [1.0 / n] * n
        ratios = [i / total for i in ratios]

    return ratios


def is_same_link(link1: str, link2: str) -> bool:
    """Check if two paths (assumed to be links) point to the same file, with error handling."""

    try:
        return os.path.realpath(link1) == os.path.realpath(link2)
    except OSError:
        return False


def build_link(source: str, target: str) -> bool:
    """Build a symbolic link from the source to the target."""

    if not os.path.exists(source):
        return False
    if os.path.islink(target) or os.path.exists(target):
        os.remove(target)
    os.symlink(source, target)
    return True


def get_data_based_on_mode(data: np.ndarray, mode: str) -> Any:
    """Processes the given data based on the specified mode."""

    if len(data) == 0:
        raise ValueError("Data is empty.")

    if mode == "all":
        return data
    if mode == "ave":
        return np.mean(data)
    if mode == "max":
        return np.max(data)
    if mode == "min":
        return np.min(data)
    if mode == "last":
        return data[-1]
    if mode.isdigit():
        ind = int(mode)
        if ind >= len(data) or ind < 0:
            raise ValueError(f"Index out of range: {ind}")
        return data[ind]
    raise ValueError(f"SSU: Unknown mode: {mode}")


def add_val_to_str_list(instring: str, val: str, delimiter: str = ",") -> str:
    """
    Split the input list by the delimiter, strip spaces,
    Then add the value to each non-empty element and return a string with delimiter
    """
    out_list = [f"{item.strip()}{val}" for item in instring.split(delimiter) if item.strip()]

    return delimiter.join(out_list)


def split_and_get_sublist(inlist: list, num_splits: int, sublist_index: int) -> list:
    """Split the input list into sublists and return the specified sublist."""

    # Check if the sublist index is valid
    if not 0 <= sublist_index < num_splits:
        raise ValueError("Invalid sublist index")

    # Calculate the size of each split and the number of larger splits
    split_size, num_larger_splits = divmod(len(inlist), num_splits)

    # Calculate the start and end indices for the sublist
    start_index = (
        sublist_index * (split_size + 1)
        if sublist_index < num_larger_splits
        else sublist_index * split_size + num_larger_splits
    )
    end_index = start_index + split_size + 1 if sublist_index < num_larger_splits else start_index + split_size

    # Return the specified sublist
    return inlist[start_index:end_index]


def convert_number_format_in_string(string: str) -> str:
    """Convert scientific notation in a string from aEb to a$\\times 10^b$ format for LaTeX."""
    pattern = r"(\d+(?:\.\d+)?)E([+-]?\d+)"
    return re.sub(pattern, lambda m: f"{float(m.group(1)):.1f}$\\times 10^{{{int(m.group(2))}}}$", string)


def reformat_number_in_string(number_str):
    """Reformat number in a string separated by &. aEb to a$^{b}$ used in LaTeX"""

    numbers = number_str.split("&")

    for i, num in enumerate(numbers):
        try:
            num = float(num)
            if 0.01 < num < 1000:
                decimal_places = len(str(num).split(".")[1].rstrip("0"))
                a_str = f"{num:.{decimal_places}f}"
                numbers[i] = a_str.rstrip("0") if "." in a_str else a_str
            else:
                numbers[i] = f"{num:.2e}".replace("e", "$\\times 10^{{").replace("+", "") + "}$"
        except ValueError:
            pass

    return "&".join(numbers)


def get_average_concs(concs_all: list) -> dict:
    """Get average of a list of concentrations (from multiple simulations)."""

    if len(concs_all) == 1:  # no average needed
        return concs_all[0]

    # Find shared keys among all dicts
    shared_phases = set(concs_all[0].keys())
    for c in concs_all[1:]:
        shared_phases &= set(c.keys())

    ave_concs = {}
    for phase in shared_phases:
        shared_species = set(concs_all[0][phase].keys())
        nsps = len(shared_species)
        for conc in concs_all[1:]:
            shared_species &= set(conc[phase].keys())
            if nsps != len(shared_species):
                logger.warning("Species mismatch in %s: expected %d, found %d.", phase, nsps, len(shared_species))
        sps_ave = {}
        for s in shared_species:
            # Compute average concentration for each species
            sps_ave[s] = np.mean([c[phase][s] for c in concs_all], axis=0)
        ave_concs[phase] = sps_ave

    return ave_concs


def replace_species_in_conc(concs_in: dict, fakep: str) -> None:
    """Update the concentration dictionary with fake species concentrations."""

    sps_wfake = [s for s in concs_in if f"{fakep}{s}" in concs_in]

    # Check if any fake species found
    if not sps_wfake:
        logger.warning("No fake species found in concentrations with prefix '%s'.", fakep)
        return

    # Replace concentrations
    for s in sps_wfake:
        concs_in[s] = concs_in.pop(f"{fakep}{s}")
    logger.debug("Replaced concentrations of # %d real species with fake species", len(sps_wfake))


def restart_filename(file_dir: str, is_new: bool, only_num: bool) -> Union[str, int, None]:
    """Return the restart filename."""
    if not (file_dir and os.path.exists(file_dir)):
        return None
    pattern = re.compile(r"^restart_(\d+)\.json")
    max_n = 0
    for f in os.listdir(file_dir):
        match = pattern.match(f)
        if match:
            n = int(match.group(1))
            max_n = max(max_n, n)

    if is_new:  # Return new filename
        max_n += 1

    if max_n == 0:  # No restart file found
        return None

    return max_n if only_num else os.path.join(file_dir, f"restart_{max_n}.json")
