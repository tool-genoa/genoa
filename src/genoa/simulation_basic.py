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
This module contains basic functions to run a simulation and get error from the output file.
"""

from .logger import setup_logger


# Logger
logger = setup_logger(__name__)


def get_default_namelist(namelist: str) -> dict:
    """return namelist content and a dict containing default keys and values."""

    # Get default namelist
    with open(namelist, "r", encoding="utf-8") as f:
        info = f.readlines()

    # Find changable content
    nml_content = {"info": info}  # Record all info to a list

    for i, line in enumerate(info):
        # Find key & values separated by "="
        if line and "=" in line:
            key = line.split("=", 1)[0].strip()
            # Record line number in content
            nml_content[key] = i

    # Check if no changable values
    if len(nml_content) <= 1:
        raise ValueError(f"Not found changable values in {namelist}.")

    return nml_content


def write_new_namelist(new_file: str, nml_info: list, new_lines: dict) -> None:
    """Write down a new namelist with items updated from the default gecko namelist."""

    with open(new_file, "w", encoding="utf-8") as f:
        for i, line in enumerate(nml_info):
            if i in new_lines:  # Use new line
                iline = new_lines[i]
            else:  # Use from content
                iline = line
            f.write(iline)


def update_nml_line(nml_content: dict, items: dict) -> dict:
    """Update default namelist information"""

    for key, val in items.items():
        if key not in nml_content:
            raise ValueError(f"Cannot find {key} in namelist. Check.")

        # Update content in the list for output
        nml_content["info"][nml_content[key]] = val

    return nml_content


def get_nline_dict(items: dict) -> dict:
    """Output a new line for the namelist based on the input value type."""

    new = {}

    for k, v in items.items():
        if isinstance(v, str):  # String
            new[k] = f'{k} = "{v}" \t!Updated\n'
        else:  # Number
            new[k] = f"{k} = {v} \t!Updated\n"

    return new
