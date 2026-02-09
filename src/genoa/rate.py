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
This module contains the class to store and manage kinetic information.
"""

from typing import List, Optional

import numpy as np
from attrs import define, evolve

from .constants import PRC_FL, NSZA
from .gecko_cst import EXTRA_KEYS
from .logger import setup_logger


# Logger
logger = setup_logger(__name__)

# Symbol to record reduction ratios in reactions
SYM_RT = " ;rt:"

# Length of coefficient for different type of kinetic (w/o ratio)
SSH_RTE_LEN = {
    0: ["EXTRA 94"],
    1: ["PHOT", "EXTRA 10", "EXTRA 99", "EXTRA 92", "EXTRA 93"],
    3: ["ARR", "EXTRA 91", "EXTRA 100", "EXTRA 500", "EXTRA 502"],
    6: ["EXTRA 501"],  # arrs(3) + coefs(3)
    8: ["EXTRA 200"],  # arrs(3) + coefs(5)
    9: ["EXTRA 550"],  # arrs(3) + coefs(6)
    10: ["FALLOFF"],  # arrs(3) + coefs(7)
}
# Add photolysis
if NSZA in SSH_RTE_LEN:
    SSH_RTE_LEN[NSZA].append("PHOTOLYSIS")
else:
    SSH_RTE_LEN[NSZA] = ["PHOTOLYSIS"]  # tabs(nsza)

# Get dictionary of key to len
SSH_KYW_LEN = {}
for ncoef, keys in SSH_RTE_LEN.items():
    for k in keys:
        SSH_KYW_LEN[k] = ncoef


@define(slots=True)
class Kinetic:
    """
    Class to store and manage kinetic information.
    SSH-aerosol format as default format.
    """

    # Original input string or comments with reduction info
    string: str = ""
    # SSH-aerosol format
    ssh: str = ""
    # Kinetic rate contant values for different environmental set
    kvalues: Optional[np.ndarray] = None

    # Decomposition of the ssh string
    mainkey: str = ""  # Main keyword for the kinetic
    key_str: str = ""  # String of keys for the kinetic
    coefs: Optional[List[float]] = None  # Coefficients for the kinetic
    coefs_str: Optional[List[str]] = None  # Coefficients as strings
    w_c1: bool = False  # Whether the kinetic has C1 value
    ratio: float = 1.0  # Ratio for the kinetic
    ro2: Optional[str] = None  # RO2 index for the kinetic, default is 0
    tb: Optional[str] = None  # TB index for the kinetic, default is empty string

    def copy(self) -> "Kinetic":
        """Return a deep copy of the kinetic instance."""
        new_k = evolve(self)
        new_k.kvalues = np.copy(self.kvalues) if self.kvalues is not None else None
        new_k.coefs = self.coefs.copy() if self.coefs is not None else None
        new_k.coefs_str = self.coefs_str.copy() if self.coefs_str is not None else None
        return new_k

    def init_rate_with_string(self) -> bool:
        """Decompose the ssh string and update the kinetic attributes."""

        if not self.ssh:
            logger.error("SSH string is empty. Cannot decomposition.")
            return False

        # Get parts from ssh string
        parts = [i for i in self.ssh.split(" ") if i != ""]
        nparts = len(parts)
        if nparts < 2:
            logger.error("SSH string is too short: %s", self.ssh)
            return False

        # Get 1st keyword if any
        n, key, self.ro2, self.tb = 1, parts[1], None, None
        if key in ("TB", "RO2"):
            if key == "TB":
                self.tb = parts[2]
            elif key == "RO2":
                self.ro2 = parts[2]
            else:
                logger.error("Unknown keyword: %s from %s", key, parts)
                return False
            n += 2
            key = parts[n]  # Update to 2nd keyword

        # Get mainkey
        mainkey = ""
        if key in ("ARR", "PHOT", "FALLOFF", "PHOTOLYSIS"):
            mainkey = key  # No sublabel
            n += 1
        elif key in ("EXTRA",):
            skey = parts[n + 1]  # W/ sublabel
            mainkey = f"{key} {skey}"
            n += 2
        # Check mainkey
        if mainkey not in SSH_KYW_LEN:
            logger.error("Unknown keyword: %s from %s", mainkey, parts)
            return False

        # Get number of coefficients for the mainkey
        ncoefs = SSH_KYW_LEN[mainkey]
        if nparts < n + ncoefs:
            logger.error("No enough coefficents. Got # parts: %s, n: %s, # coefs: %s", nparts, n, ncoefs)
            return False

        # Load mainkey & coefficients
        self.mainkey = mainkey
        self.coefs = [float(i) for i in parts[n : n + ncoefs]]
        self.key_str = " ".join(parts[1:n])
        self.coefs_str = parts[n : n + ncoefs]
        n += ncoefs

        # Get ratio if any
        if nparts > n + 1:
            logger.error("Too many parts for %s from %s. Got %s, need %s or %s", mainkey, parts, nparts, n, n + 1)
            return False
        self.ratio = float(parts[n]) if nparts == n + 1 else 1.0

        # Update  C1
        self.w_c1 = mainkey in ("ARR", "FALLOFF") or mainkey in EXTRA_KEYS

        return True

    def update_ssh_with_ratio(self) -> None:
        """Update the ssh string with the current kinetic information."""

        if not self.mainkey:
            logger.error("Main key is not set. Cannot update SSH.")
            return

        if self.ratio == 1.0:
            return  # No need to update if ratio is 1.0

        if self.w_c1:
            self.coefs[0] *= self.ratio  # Update C1 with ratio
            self.coefs_str[0] = f"{self.coefs[0]:{PRC_FL}}"
            self.ssh = f"KINETIC {self.key_str} " + " ".join(self.coefs_str)
            self.ratio = 1.0  # Reset ratio for C1 kinetics
        else:
            self.ssh = f"KINETIC {self.key_str} " + " ".join(self.coefs_str) + f" {self.ratio:{PRC_FL}}"

    def multiply_rate_by_ratio(self, ratio: float) -> bool:
        """Multiply kinetic rate by a ratio. Return rcn status."""

        # Check if ratio is valid
        if ratio <= 0.0:  # Inactive
            return False

        if ratio == 1.0:  # No change
            return True

        # Update strings
        self.string += f"*{ratio:{PRC_FL}}" if SYM_RT in self.string else f"{SYM_RT}{ratio:{PRC_FL}}"
        self.ratio *= ratio
        self.update_ssh_with_ratio()

        # Update kinetic values if needed
        if self.kvalues is not None:
            self.kvalues *= ratio

        return True

    def get_ratio(self) -> float:
        """Get the ratio of the kinetic."""
        if self.w_c1:
            return self.coefs[0] * self.ratio
        return self.ratio
