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
This module contains the class to store and manage information related to chemical species.
"""

import math
import json
from typing import List, Dict, Any, Optional, ClassVar

from attrs import define, field, evolve

from .logger import setup_logger
from .smiles import get_properties_from_smiles


# Logger
logger = setup_logger(__name__)


@define(slots=True)
class Species:
    """Class to store and manage information related to chemical species"""

    # Basic information
    name: str = ""  # Species name
    string: str = ""  # Input format

    # Status
    """
    Status codes:
    0: Not active
    1: Active (organics can be reduced)
    2: Basic species (should not be reduced)
    """
    status: int = 1

    # Molecular mass in g/mol
    mass: float = 0.0

    # Flags
    radical: bool = False  # Radical flag
    RO2: int = 0  # RO2 flag or group id
    condensable: bool = False  # Condensable flag
    non_volatile: bool = False  # Non-volatile flag

    # Chemical information
    formula: str = "-"  # Formula as C[i]H[i]N[i]O[i]
    smiles: str = "-"  # SMILES structure
    fgroups: Dict[str, float] = field(factory=dict)  # Functional group info
    ratios: Dict[str, float] = field(factory=dict)  # Atomic ratios: OM/OC, N/C, H/C, O/C

    # Aerosol properties (Tref: 298 K)
    psat_atm: float = 0.0  # Saturation vapor pressure at Tref (atm)
    dhvap_kj: float = 0.0  # Enthalpy of vaporization (kJ/mol)

    # SSH-aerosol related
    grp_id: Optional[int] = None  # Group id
    soap_strucs: Dict[str, float] = field(factory=dict)  # Non-zero functional groups

    # GECKO-A related
    note: str = ""  # Simplified structure code
    wall_loss: bool = False  # Wall loss flag
    tg: Optional[float] = None  # Tg
    dvol: Optional[float] = None  # Diffusion coefficient (cm^2/s)

    # Mechanism related
    generation: int = -1  # Generation number
    precursors: List[str] = field(factory=list)  # Precursors

    # Reduction related
    rdc_id: Optional[str] = None  # Species type used for reduction
    reductions: Dict[str, List] = field(factory=dict)  # Reduction records

    _default: ClassVar["Species"] = None

    @classmethod
    def _get_default(cls) -> "Species":
        """Get the default species instance."""
        if cls._default is None:
            cls._default = cls(status=0)
        return cls._default

    def check_and_set_attr(self, key: str, value: Any) -> bool:
        """Check and set an attribute if it exists."""
        if hasattr(self, key):
            setattr(self, key, value)
            return True
        return False

    def copy(self) -> "Species":
        """Create a deep copy of the species instance."""
        new_sps = evolve(self)
        new_sps.fgroups = self.fgroups.copy()
        new_sps.ratios = self.ratios.copy()
        new_sps.precursors = self.precursors.copy()
        new_sps.soap_strucs = self.soap_strucs.copy()
        new_sps.reductions = {k: v.copy() for k, v in self.reductions.items()}

        return new_sps

    def output_to_lines(self) -> list:
        """Output all attributes to a text string."""
        lines = ["\n"]
        for f in self.__attrs_attrs__:
            name = f.name
            current = getattr(self, name)
            if current != getattr(self._get_default(), name):
                lines.append(f"{name}\t{json.dumps(current)}\n")
        lines.append("\n")  # Add an empty line at the end
        return lines

    def read_from_lines(self, lines: List[str]) -> None:
        """Updates object attributes from a list of text strings."""
        for line in lines:
            if "\t" in line:
                key, val_json = line.split("\t", 1)
                self.check_and_set_attr(key, json.loads(val_json))
        # Update soap_strucs from list to dict
        # if isinstance(self.soap_strucs, list):
        #    self.soap_strucs = {v: self.soap_strucs[1][i] for i, v in enumerate(self.soap_strucs[0])}

    def update_by_functional_group(self) -> None:
        """
        Get the number of key chemical elements in the given molecule,
        return OM/OC mass ratio, H/C, O/C, N/C atomic ratios
        and the degree of unsaturation

        For a compound with formula CaHbNcOdXe where X is F, Cl, Br or I
        the degree of unsaturation is given by:
        degree of unsaturation = 1/2 (2 + 2C + N - H - X)
        int(0.5*(2+2*tmp['C']-tmp['H']+tmp['N']))
        """
        self.formula = f"C{self.fgroups['C']}H{self.fgroups['H']}" + f"N{self.fgroups['N']}O{self.fgroups['O']}"
        self.ratios["OM/OC"] = round(self.mass / (self.fgroups["C"] * 12.0), 3)
        self.ratios["H/C"] = round(self.fgroups["H"] / self.fgroups["C"], 3)
        self.ratios["O/C"] = round(self.fgroups["O"] / self.fgroups["C"], 3)
        self.ratios["N/C"] = round(self.fgroups["N"] / self.fgroups["C"], 3)

    def update_based_on_smiles(self, vptype: str) -> None:
        """Updates various properties based on its SMILES string."""

        # Update molecular properties from smiles
        # radical, RO2, mass, fgroups, psat_atm, dhvap_kj
        self.update_with_dict(get_properties_from_smiles(self.smiles, vptype))

    def update_for_reduction(self, is_simple: bool = False) -> None:
        """Update the type of species based on its characteristics."""

        itype = "R" if self.radical else "A" if self.condensable else "G"

        if is_simple:  # Only basic type
            self.rdc_id = itype
            return

        if self.radical:  # radical
            # RO2
            if self.RO2:
                itype += "O2"
            # ROO
            elif "[O+][O-]" in self.smiles:
                itype += "OO"
            # RO
            elif "[O]" in self.smiles:
                itype += "O"

            self.rdc_id = itype
            return

        # Non-radical species
        # Add only functional groups with count > 0 and not basic atoms
        fg_str = ",".join(k for k in self.fgroups if k not in {"C", "H", "N", "O"} and self.fgroups[k] > 0)

        # Optional: SOAP structure inclusion (commented for performance unless needed)
        # if condensable and self.soap_strucs:
        #     s_keys = [k for k, v in self.soap_strucs.items() if v > 0]
        #     if s_keys:
        #         fg_str1 = " " + ",".join(
        #             f"{k}:{self.soap_strucs[k]:.1f}" for k in sorted(s_keys, key=int)
        #         )

        self.rdc_id = itype + fg_str
        return

    def update_with_dict(self, new_dict: Dict[str, Any]) -> None:
        """Check species properties based on input dictionary."""

        if not new_dict:
            return

        # Update from input dict
        for key, val in new_dict.items():
            self.check_and_set_attr(key, val)

        # Update ratios & formula
        self.update_by_functional_group()

        # Update for aerosols
        if self.tg is not None:
            self.update_prop_wgck()

        # Update type info
        self.update_for_reduction()

    def update_condense_status(self, psvoc: float, pnvoc: float) -> None:
        """Update the volatility of a species based on SVOC/NVOC setting.."""

        if self.radical:
            if self.condensable:
                raise ValueError(f"Radical species {self.name} cannot be condensable.")
            return

        if self.psat_atm <= 0.0:
            self.condensable, self.non_volatile = False, False
            return

        icondense, invoc = self.condensable, self.non_volatile

        if psvoc and self.psat_atm > psvoc:  # VOC
            icondense = False

        if pnvoc and self.psat_atm < pnvoc:  # NVOC
            invoc = True

        if self.condensable == icondense and self.non_volatile == invoc:  # No change
            return

        if self.condensable != icondense:
            logger.info("Update condensable for %s: %s -> %s", self.name, self.condensable, icondense)
            logger.info("Psat: %.3e, psvoc: %.3f, pnvoc: %.3f", self.psat_atm, psvoc, pnvoc)
        if self.non_volatile != invoc:
            logger.info("Update non-volatile for %s: %s -> %s", self.name, self.non_volatile, invoc)

        # Update condensable and non-volatile flags
        self.condensable, self.non_volatile = icondense, invoc

    def update_prop_wgck(self) -> None:
        """
        Estimate properties based on the structure using methods rewriten GECKO-A generator subroutines
        wrt_Tg(), diffusion_vol()
        """

        if not self.fgroups:
            raise ValueError(f"Functional groups are not defined for species: {self.name}")
        n_c, n_h, n_o, n_n = self.fgroups["C"], self.fgroups["H"], self.fgroups["O"], self.fgroups["N"]
        nring = self.fgroups.get("Rc", 0) + self.fgroups.get("RN", 0)
        self.tg = compute_tg(n_c, n_h, n_o, n_n)
        self.dvol = compute_vdiffusion(n_c, n_h, n_o, n_n, nring > 0)


def add_basic_to_species_list(species: list, sps_dict: dict, basic_dict: dict) -> None:
    """Add basic species to the species list."""

    for s, imw in basic_dict.items():
        # Add basic species
        if s not in sps_dict:
            logger.info("Add basic species %s to the species list.", s)
            species.append(Species(name=s, mass=imw, status=2))
            sps_dict[s] = species[-1]
        # Update status
        elif sps_dict[s].status != 2:
            logger.info("Update status = 2 for basic species %s.", s)
            sps_dict[s].status = 2


def update_basic_dict(sps_basic: dict, basic_dict: dict) -> int:
    """Update basic species dictionary."""

    n_updated = 0

    for k, v in sps_basic.items():
        if k not in basic_dict:
            logger.info("Basic species %s not found in the basic list. Add it.", k)
            basic_dict[k] = v
            n_updated += 1

    if n_updated:
        logger.info("Updated %d basic species.", n_updated)

    return n_updated


def get_species_category(species: list, wpvocs: Optional[set] = None) -> dict:
    """Get dict of species categories."""
    org, inorg = {}, {}  # {name: index in species list}, species = org + inorg
    # species = rsp + nrsp; org = voc + radical + svoc
    rsp, nrsp, voc, radical, svoc, ro2 = set(), set(), set(), set(), set(), set()

    inorg_wc = {"CO2", "CO", "XCLOST"}
    for s in species:
        if s.status < 1:
            continue
        sname = s.name
        # Reduciable or non-reduciable
        if s.status == 1:
            rsp.add(sname)
        else:
            nrsp.add(sname)
        # Organic or inorganic
        if sname in inorg_wc or s.fgroups.get("C", sname.count("C")) == 0:
            inorg[sname] = s
            continue
        if s.status > 1:
            logger.info("Find organic species %s non-reducible.", sname)
        org[sname] = s
        if s.radical:
            radical.add(sname)
            if s.RO2:
                ro2.add(sname)
        elif s.condensable:
            svoc.add(sname)
        else:
            voc.add(sname)
    sdict = {
        "org": org,
        "inorg": inorg,
        "rsp": rsp,
        "nrsp": nrsp,
        "voc": voc,
        "radical": radical,
        "svoc": svoc,
        "ro2": ro2,
    }

    # Remove pvocs
    if wpvocs:
        pvoc = {s for s in wpvocs if s in rsp}
        if pvoc:
            logger.info("Found # %d pvoc in rsp: %s. Remove them from rsp, voc, svoc, and radical.", len(pvoc), pvoc)
            sdict["pvoc"] = pvoc
            for s in [rsp, voc, svoc, radical]:
                s -= pvoc
        else:
            logger.info("No pvoc found in rsp.")

    logger.info("Catogorized # %d species to org + inorg, or rsp + nrsp. org = voc + radical + svoc", len(species))
    logger.info("==> %s", {k: len(v) for k, v in sdict.items()})

    return sdict


def compute_tg(n_c: float, n_h: float, n_o: float, n_n: float) -> float:
    """
    Estimate glass transition temperature (Tg) in Kelvin based on atom counts.
    This function replicates the logic from the GECKO-A generator subroutine `wrt_Tg`.
    """
    if n_c <= 0 or n_h <= 0:
        return 0.0

    logc = math.log(n_c)
    logh = math.log(n_h)

    if n_o == 0 and n_n == 0:
        tg = (1.96 + logc) * 61.99 + logh * (-113.33) + logc * logh * 28.74
    elif n_n == 0:
        logo = math.log(n_o)
        tg = (12.13 + logc) * 10.95 + logh * (-41.82) + logc * logh * 21.61 + logo * 118.96 + logc * logo * (-24.38)
    else:
        logo = math.log(n_o)
        logn = math.log(n_n)
        tg = (
            551.0146
            + logc * (-61.8525)
            + logh * (-254.1697)
            + logo * 146.0169
            + logn * 136.7473
            + logc * logh * 82.1893
            + logc * logo * (-57.9076)
            + logc * logn * (-43.7932)
        )

    return round(tg, 2)


def compute_vdiffusion(n_c: float, n_h: float, n_o: float, n_n: float, wring: bool) -> float:
    """
    Estimate the Fuller diffusion volume from atom counts and ring flags.
    This function replicates the logic from the GECKO-A generator subroutine `diffusion_vol`.
    """

    v = n_c * 15.9 + n_h * 2.31 + n_o * 6.11 + n_n * 4.54
    if wring:  # Ring correction for aromatic & heterocyclic rings
        v += -18.3
    return round(v, 2)


def update_species_condensibility(species: list, psvoc: float, pnvoc: float) -> None:
    """
    Update the condensibility of species based on saturation vapor pressure.
    This function updates the `condensable` and `non_volatile` flags for each species.
    """
    if not species:
        raise ValueError("Species list is empty. Cannot update condensibility.")
    if psvoc == 0.0 and pnvoc == 0.0:
        return
    if psvoc < 0.0 or pnvoc < 0.0:
        raise ValueError("psvoc and pnvoc must be non-negative values.")

    logger.info("Updating species condensibility with psvoc=%.3f, pnvoc=%.3f ...", psvoc, pnvoc)
    for s in species:
        if s.status > 0:
            s.update_condense_status(psvoc, pnvoc)


def update_species_reduction_info(species: list) -> None:
    """Update rdc_id for species based on their properties."""

    logger.info("Updating species reduction info...")
    for s in species:
        if s.status > 0:
            s.update_for_reduction()


def remove_pvoc_partitioning(species: list, primary_vocs: list) -> None:
    """Remove the gas-particle partitioning of primary vocs."""

    logger.info("Removing gas-particle partitioning for primary VOCs: %s", primary_vocs)
    for s in species:
        if s.name in primary_vocs and s.condensable:
            s.condensable = False  # Remove the condensability
            logger.info("- Removed the condensability of primary VOC: %s", s.name)
