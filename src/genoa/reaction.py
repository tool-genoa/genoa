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
This module contains the class Reaction to store and manage information
related to chemical reactions.
"""

from collections import defaultdict
from typing import List, Dict, Optional

import numpy as np
from attrs import define, field

from .constants import PRC_ST, PRC_FL, RCT_DICT, NORATE_STR
from .logger import setup_logger
from .rate import Kinetic


# Logger
logger = setup_logger(__name__)


@define(slots=True)
class Reaction:
    """Class to store and manage information related to chemical reactions"""

    # Reduction status
    # 0: inactive; 1: active and reducible; 2: active and non-reducible
    status: int = 1

    # List of reactants: ratio
    reactants: Dict[str, float] = field(factory=lambda: defaultdict(float))

    # List of products: ratio
    products: Dict[str, float] = field(factory=lambda: defaultdict(float))

    # Kinetic rate
    rate: Kinetic = field(factory=Kinetic)
    kuni: Optional[list] = field(factory=list)  # [reactant, [kuni_s]]

    # Product function ratio
    # prt_fnc: ProductFunction = ProductionFunction()

    # Information for merging/splitting reactions
    rcn_id: Optional[str] = None  # Reaction ID

    def copy(self) -> "Reaction":
        """Return a deep copy of the reaction instance."""
        new_r = Reaction(
            status=self.status,
            reactants=self.reactants.copy(),
            products=self.products.copy(),
            rate=self.rate.copy(),
            rcn_id=self.get_rcn_id(),
        )
        new_r.kuni = [self.kuni[0], self.kuni[1].copy()] if self.kuni else []
        return new_r

    def update_status(self, basic_dict: dict) -> None:
        """Update status of reducible reaction."""
        if self.status != 1:
            return
        if not self.reactants or self.rate.ratio <= 0.0:  # Inactive
            self.status = 0
            return
        for s in self.reactants:
            if s not in basic_dict:
                return  # Active reducible reaction
        self.status = 2  # Active non-reducible reaction

    def reset_species(self) -> None:
        """Empty reactants and products and their ratio."""
        self.reactants = defaultdict(float)
        self.products = defaultdict(float)

    def reset_products(self) -> None:
        """Empty all product and its ratio."""
        self.products = defaultdict(float)

    def clean_duplicates(self) -> bool:
        """Clean duplicates in reactants and products."""

        # Get list of duplicates in reactants and products
        dsps = set(self.reactants) & set(self.products)

        # Allow duplicates for default reactants
        if dsps:
            dsps = [s for s in dsps if s not in RCT_DICT]

        # No duplicate
        if not dsps:
            return False

        # logger.warning("Find repeated species in reaction: %s", self.to_rcn())
        # Remove duplicates
        for s in dsps:
            r_rt, p_rt = self.reactants[s], self.products[s]
            # Mass is not balanced
            if p_rt >= r_rt:
                self.status = 0
                # logger.warning("A reaction is removed cuz %s w/ %s >= %s", s, p_rt, r_rt)
                return True

            # Remove product
            if p_rt <= 0.0:
                self.products.pop(s)
                # logger.warning("A product is removed cuz %s w/ %s <= 0.0", s, p_rt)
                return True

            # Update product & kinectic rate
            new_rt = r_rt - self.products.pop(s)

            # Update ratios for other products
            for p in self.products:
                if p not in self.reactants:
                    self.products[p] /= new_rt
            # Update kinetic rate
            if not self.rate.multiply_rate_by_ratio(new_rt):
                # logger.warning("A reaction is removed cuz rate is not multiplied by ratio %s", new_rt)
                self.status = 0

            # logger.warning("A reaction has been updated with ratio %s. New rcn: %s", new_rt, self.to_rcn())

        return True

    def clean_products(self) -> None:
        """Clean products with zero ratio."""

        for s in list(self.products):
            if self.products[s] > 0.0:
                continue
            if self.products[s] < 0.0:
                raise ValueError(f"Negative ratio for product: {s} in {self.to_rcn()}")
            self.products.pop(s)

    def get_sps_for_rate(self, basic_dict: dict, only_main: bool = False) -> Optional[List[str]]:
        """Return a list containing main reactant, the rest reactants, and ro2 pool id for rate computation."""
        if self.status != 1:
            return None
        sps = [s for s, rt in self.reactants.items() if s not in basic_dict and rt > 0.0]
        if not sps:
            raise ValueError(f"No main reactant in reducible reaction: {self.to_rcn()}")
        if len(sps) > 1:
            raise ValueError(f"Multiple main reactants: {sps} in {self.to_rcn()}")

        if only_main:
            return sps

        sps = sps + [s for s in self.reactants if s not in sps]  # Add the rest reactants

        # Check if contains ro2 pool
        if self.rate.ro2 is not None:
            sps.append(f"ro2s{self.rate.ro2}")
        return sps

    def get_ratios(self, slist: list) -> list:
        """Get product ratios."""

        # Get basic ratios
        ratios = np.zeros(len(slist))
        for i, s in enumerate(slist):
            if s in self.products:
                ratios[i] = self.products[s]

        return ratios

    def update_rcn_id(self) -> None:
        """Get reduction identifier from this reaction to check."""

        info = " ".join(sorted(self.reactants)) + " " + self.rate.key_str + " "
        n = 1 if self.rate.w_c1 else 0
        self.rcn_id = info + " ".join(f"{i:{PRC_FL}}" for i in self.rate.coefs[n:])

    def get_rcn_id(self) -> str:
        """Return the reduction identifier for this reaction."""
        if self.rcn_id is None:
            self.update_rcn_id()
        return self.rcn_id

    def to_line_wo_kinetic(self, separator="->"):
        """Output reaction line w/o kinetics"""

        # Reactants
        parts = []
        for s, rt in self.reactants.items():
            for _ in range(int(rt)):
                parts.append(s)
        outline = " + ".join(sorted(parts)) + f" {separator} "

        # Products
        parts = []
        for s in sorted(self.products.keys()):
            if self.products[s] != 1.0:
                parts.append(f"{self.products[s]:{PRC_ST}} {s}")
            else:
                parts.append(s)

        return outline + " + ".join(parts)

    def from_line_wo_kinetic(self, line, separator="->"):
        """Reconstructs a reaction from line w/o kinetics."""

        # Read lines for reactants and products
        rcline, pdline = line.split(separator, 1)

        # Reset reactants and products
        self.reset_species()

        # Record reactants
        parts = []
        for s in rcline.strip().split(" "):
            if s in ("", "+"):
                continue
            parts.append(s)

        for s in parts:
            self.reactants[s] += 1

        # No read reactants
        if len(self.reactants) == 0:
            raise ValueError(f"No reactant can be read from line: {line}")

        # Products
        parts = [i for i in pdline.strip().split(" ") if i != ""]

        # No product
        if not parts:
            return

        n = len(parts)

        # Only one product
        if n == 1:
            self.products[parts[0]] += 1.0
            return

        # Multiple products
        sps = []
        # Add + to find end
        for i, s in enumerate(parts + ["+"]):
            if s == "+":  # Find end
                n1 = len(sps)
                if n1 not in (1, 2, 4):
                    raise ValueError(f"# {n1} in a product: {sps} in {line}. n1 should be 1,2,4. parts: {parts}", i, s)

                isp = sps[-1]  # Species name

                # Get ratios
                if n1 == 1:
                    rt = 1.0
                else:
                    rt = float(sps[0])

                if n1 == 4:  # Find function ratio
                    raise ValueError("Not set yet.")

                # Add ratio
                self.products[isp] += rt
                sps = []  # Reset

            else:  # Find one element: ratio or species
                sps.append(s)

    def to_rcn(self, str_in=None) -> str:
        """write down a reaction."""

        rlines = [] if str_in is None else [f"%=========={str_in}=========="]

        # Reactants and products
        rlines.append(self.to_line_wo_kinetic())
        # Kinetic comment
        rlines.append(f"%{self.rate.string}")
        # Kinetic SSH-aerosol excutable format
        rlines.append(self.rate.ssh)

        return "\n".join(rlines) + "\n"

    def from_rcn(self, strin: str) -> None:
        """Read a reaction."""

        # Get string for reaction and kinetic
        rcn_str, rate_str, comments, n = "", "", [], 0
        # Remove product line break if needed
        # line_break_str = " //\n"
        # strin = strin.replace(line_break_str, "")

        # Read reaction string
        for line in strin.split("\n"):
            line = line.strip()
            # Skip comment or empty line
            if not line:
                continue
            # Find comments
            if line.startswith("%"):
                # Record comments for kinetic
                if "===" not in line:
                    # Remove %
                    comments.append(line[1:])
            # Find reaction line
            elif "->" in line:
                rcn_str = line
                self.from_line_wo_kinetic(line)
            # Read kinetic line
            elif line.startswith("KINETIC"):
                n += 1  # Count
                if n > 1:  # Find for ratio
                    raise ValueError(f"Multiple kinetics but the no sign for function ratio in {strin}")

                rate_str = line
                self.rate.ssh = line
                if not self.rate.init_rate_with_string():
                    raise ValueError(f"Can not initialize rate with string: {line}")

            elif line not in ("END"):
                raise ValueError(f"Can not understand line for load reactions: {line}")

        # Check and update updates
        if rcn_str == "" or rate_str == "":
            raise ValueError("Can not update reaction with a string: {strin}")

        # Update rate string to preserve comments
        self.rate.string = "\n%".join(comments)


def get_reduction_info_from_list(reactions: list) -> dict:
    """Get a dictionary containing reactions with the same reduction identifier"""

    rdc_info = {}

    for i, rcn in enumerate(reactions):

        # Get reducible reaction
        if rcn.status != 1:
            continue

        # Store info
        info = rcn.get_rcn_id()
        if info in rdc_info:
            rdc_info[info].append(i)
        else:
            rdc_info[info] = [i]

    return rdc_info


def get_reactions_with_fake_species(reactions: list, fakes: dict, prefix: str) -> list:
    """Return a list of reactions with updated fake species"""
    # Check fake species
    if not fakes:
        logger.warning("No fake species found.")
        return reactions

    new_reactions = []
    for rcn in reactions:
        if rcn.status <= 0:
            continue
        # Check products in fakes
        fake_pdts = [s for s in rcn.products if s in fakes]
        if fake_pdts:
            rcn = rcn.copy()
            # Update products
            for s in fake_pdts:
                rcn.products[f"{prefix}{s}"] += rcn.products[s]
        new_reactions.append(rcn)
    return new_reactions


def with_norate_reactions(reactions: list) -> int:
    """Check if there is any reaction with no rate."""

    norate = 0  # Number of reactions with no rate

    logger.info("Checking reactions with no rate...")
    for rcn in reactions:
        if rcn.status <= 0:  # Only check active reactions
            continue
        if not rcn.rate or not rcn.rate.ssh:
            raise ValueError(f"Reaction {rcn.to_rcn()} has no rate string in the ssh format.")
        if NORATE_STR in rcn.rate.string:
            norate += 1
            logger.warning("Reaction %s has no rate (# %d).", rcn.to_rcn(), norate)

    logger.info("Total # %d reactions with no rate.", norate)

    return norate


def check_conservation(reactions: list) -> None:
    """Checks and ensures conservation of RCT_DICT species in the reaction."""

    for rcn in reactions:
        if rcn.status != 1:  # Skip non-reducible reactions
            continue

        # Filter out constant species from reactants
        rc_to_keep = {s: rcn.reactants[s] for s in rcn.reactants if s in RCT_DICT}

        # Initialize list to track matched species
        pd_keep = set()

        # Check and update products
        for s in list(rcn.products.keys()):
            if s in RCT_DICT:
                if s in rc_to_keep:
                    # Set ratio to the same in the reactants
                    rcn.products[s] = rc_to_keep[s]
                    pd_keep.add(s)
                else:
                    logger.info("Reaction conservation: delete %s from %s", s, rcn.to_rcn())
                    rcn.products.pop(s, None)

        # Add any unmatched reactant species to products
        for s, ratio in rc_to_keep.items():
            if s not in pd_keep:
                rcn.products[s] = ratio
