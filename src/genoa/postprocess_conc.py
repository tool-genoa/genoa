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
This moudle contains the main function to run post-processing of the simulation results.
"""

from typing import Optional, Any

import attrs
import numpy as np

from .constants import RCT_DICT, PHASES, NTOP_MAX, EPS
from .gecko_cst import GCODE_DICT
from .record import RecordOption
from .logger import setup_logger
from .mechanism_pack import Mechanism
from .unit_conversion import get_carbon_number, get_unit_conv_factors


# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class PostConc:
    """Time series of concentrations for each species."""

    total: Optional[dict] = None  # Gas + Aerosol
    gas: Optional[dict] = None  # Gas
    aero: Optional[dict] = None  # Aerosol

    steps: Optional[np.ndarray] = None  # Time steps
    mlb: Optional[str] = None  # Label of mechanism

    def __attrs_post_init__(self):
        """Post-initialization"""
        if self.total is None and self.gas and self.aero:
            self.build_total()
        if not self.total:
            raise ValueError("No valid concentration data found")

    def build_total(self) -> None:
        """Build the time series of concentration dicts from gas & aerosol dicts"""
        try:  # Reshape the concentration data
            self.gas, self.aero = reshape_cdata_by_keys([self.gas, self.aero])
        except ValueError as e:
            raise ValueError("Failed to reshape the concentration data") from e
        self.total = {k: v + self.aero[k] for k, v in self.gas.items()}

    def remove_species(self, sps: list) -> dict:
        """Remove species from the dict"""
        if not sps:
            raise ValueError("Empty species list")
        valid_sps = [s for s in sps if s in self.total]
        out = {}
        if not valid_sps:
            logger.warning("No valid species found from the list %s", sps)
            return out
        # Remove species from all phases
        for i in PHASES:
            concs = self.get_phase_concs(i)
            out[i] = {s: concs.pop(s) for s in valid_sps}
        return out

    def get_concs(self, sps: Optional[list] = None, phase: Optional[str] = None) -> dict:
        """Get the time series of concentrations for target species list"""

        if phase:
            concs = {phase: self.get_phase_concs(phase)}
        else:  # All phases
            concs = {s: self.get_phase_concs(s) for s in PHASES}

        if sps:  # Filter species
            concs = {k: {s: v[s] for s in sps if s in v} for k, v in concs.items()}

        return concs

    def get_phase_concs(self, phase: str) -> dict:
        """Get the time series of concentrations for target phase"""
        if phase not in PHASES:
            raise ValueError(f"Phase {phase} not found in the phase list")
        if phase == "gas":
            return self.gas
        if phase == "aero":
            return self.aero
        if phase == "total":
            return self.total
        raise ValueError(f"Phase {phase} not found in the concentration dict")

    def get_sum_concs(self, phase: str, sdict: dict, unit_conv: Optional[dict] = None) -> dict:
        """Return a dict of concentrations by species groups"""
        concs = self.get_phase_concs(phase)
        if unit_conv:
            concs = {k: v * unit_conv[k] for k, v in concs.items() if k in unit_conv}
        return sum_concs_by_groups(concs, sdict)


def reshape_cdata_by_keys(cdata: list, keys: Optional[Any] = None, nt: Optional[int] = None) -> dict:
    """Reshape the concentration data by keys"""
    if not cdata:
        raise ValueError("Empty input concentration data")
    if nt is None:
        nt = len(next(iter(cdata[0].values())))
    if keys is None:
        keys = set().union(*[c.keys() for c in cdata])
    for i, c in enumerate(cdata):
        cdata[i] = {k: c.get(k, np.zeros(nt)) for k in keys}
    return cdata


def sum_concs_by_groups(concs: dict, groups: dict) -> dict:
    """Return a dict of concentrations by species groups"""
    out = {}
    for k, v in groups.items():
        if not v:
            continue
        val_sps = [s for s in v if s in concs]
        if not val_sps:
            continue
        out[k] = np.sum([concs[s] for s in val_sps], axis=0)
    return out


def get_time_series(time_setting: list, boxmodel: str) -> np.ndarray:
    """Get time series for the given time setting and box model"""

    if boxmodel == "SSH":  # In format [t0, t1, dt, hour, year]
        t0, t1, dt, _, _ = time_setting
        a_time = np.arange(t0, t1 + dt, dt)

    elif boxmodel == "GECKO":  # In format [t0, t1, nt, day, month, year]
        t0, t1, nt, _, _, _ = time_setting
        a_time = np.linspace(t0, t1, nt + 1)

    else:
        raise ValueError(f"Error: Box model {boxmodel} not supported.")

    return a_time / 3600  # Convert to hours


def get_sorted_splist(concs: dict, lims: Optional[dict], log: Optional[RecordOption]) -> list:
    """Get the sorted list of species"""
    # Default limits
    if not lims:  # Keys: nlim: int, fac: float, vlim: float, sort: sum/max/min
        lims = {"nlim": NTOP_MAX}

    # Sort species list by sum/ave/max/min
    sort_dict = {"sum": np.sum, "max": np.max, "min": np.min}
    sort_key = lims.get("sort", "sum")
    if sort_key not in sort_dict:
        raise ValueError(f"Sort key {sort_key} not found in the sort dict")
    sort_concs = {s: sort_dict[sort_key](v) for s, v in concs.items()}
    slist = list(concs.keys())
    slist.sort(key=lambda x: sort_concs[x], reverse=True)
    acts = [f"sorted by {sort_key}"]

    # Save sorted species list
    if log:
        log.write(f"\n\n-- Save # {len(slist)} sorted species by {sort_key} and options: {lims}")
        log.list_to_file([f"{i + 1}: {s} {sort_concs[s]}" for i, s in enumerate(slist)])

    # Limit by number threshold
    if lims.get("nlim", -1) > 0:
        vlim = int(lims["nlim"])
        slist = slist[: min(vlim, len(slist))]
        acts.append(f"Top {vlim}")

    # Limit by fraction threshold
    if lims.get("fac", -1) > 0:
        vlim = lims["fac"] * sort_dict[sort_key](sort_concs.values())
        for i, s in enumerate(slist):
            if sort_concs[s] < vlim:
                slist = slist[:i]
                acts.append(f"concs >= {vlim:.2e} ({sort_key} * {lims['fac']})")
                break

    # Limit by value threshold
    if lims.get("vlim", -1) > 0:
        vlim = lims["vlim"]
        for i, s in enumerate(slist):
            if sort_concs[s] < vlim:
                slist = slist[:i]
                acts.append(f"concs >= {vlim:.2e} (lim on {sort_key})")
                break

    logger.info("Sorted # %d species by %s: %s", len(slist), acts, slist)
    return slist


def calc_ekma_ratio(concs: dict, species: list, unit_in: str) -> dict:
    """Return a dict of hc, nox, and hc/nox ratios"""

    # Time steps
    nt = len(next(iter(concs.values())))
    if nt < 1:
        raise ValueError(f"Invalid time steps: {nt}")

    nc_dict = get_carbon_number(species)
    unit_conv = get_unit_conv_factors(species, unit_in=unit_in, unit_out="ppb", nc_dict=nc_dict)

    # Get hydrocarbon mass in ppbC
    hc, nhc = np.zeros(nt), 0
    for s in species:
        if s.status < 1:  # Inactive species
            continue
        sname = s.name
        try:
            ihc = concs[sname] * unit_conv[sname] * nc_dict[sname]
            hc += ihc
            nhc += 1
        except KeyError as e:
            if s.status == 1:
                logger.warning("Species %s not added to HCs mass. Got error: %s", sname, e)
            continue
    logger.info("Read # %d HCs species from # %d species", nhc, len(species))
    # Get NOx mass in ppb
    nox, nnox = np.zeros(nt), 0
    for s in ["NO", "NO2"]:
        try:
            nox += concs[s] * unit_conv[s]
            nnox += 1
        except KeyError as e:
            logger.warning("Species %s not added to NOx mass. Got error: %s", s, e)
    if nnox != 2:
        logger.error("Read # %d != 2 NOx species", nnox)

    # Get HC/NOx ratios
    ratio = np.where(nox > EPS, hc / nox, 0)
    return {"hc": hc, "nox": nox, "hc2nox": ratio}


def calc_reactivity(flag: Optional[str], rct_list: list, tconc: PostConc, mech: Mechanism, unit_in: str) -> dict:
    """Calculate the reactivity in the mechanism in the gas phase"""

    if flag is None:
        flag = "rsp"

    for rct in rct_list:
        if rct not in RCT_DICT:
            raise ValueError(f"Invalid reactivity species: {rct}")

    # Add flag checks
    if flag == "ro2":
        check_sps = {s for s in mech.sinfo["radical"] if mech.sinfo["org"][s].RO2}
    elif flag == "voc":
        check_sps = set(mech.sinfo["org"].keys()) - mech.sinfo["radical"]
    elif flag in mech.sinfo:
        check_sps = mech.sinfo[flag]
    elif flag == "all":  # All species
        check_sps = mech.sinfo["org"].keys() | mech.sinfo["inorg"].keys()
    else:
        raise ValueError(f"Invalid flag {flag} for reactivity calculation")

    logger.info("Calculate reactivity for %s w/ %d species", flag, len(check_sps))

    # Get unit conversion factors
    unit_conv = get_unit_conv_factors(mech.species, unit_in=unit_in, unit_out="molec")

    # Reactivity
    rcts = {s: i for i, s in enumerate(rct_list)}
    cnts = {s: 0 for s in rcts}
    nt = len(tconc.steps)
    concs = tconc.get_phase_concs("gas")
    rct_data = np.zeros((len(rcts), nt))
    for rcn in mech.reactions:
        if rcn.status < 1:
            continue
        r1s = [s for s in rcn.reactants if s in rcts]
        if not r1s:
            continue
        r2s = [s for s in rcn.reactants if s in check_sps]
        if not r2s:
            continue

        # Get the reactivity
        irate = rcn.rate.kvalues[0]
        for s in r1s:
            ircty = np.ones(nt) * irate
            for r in rcn.reactants:
                if s == r:
                    continue
                if r in concs:
                    iconc = concs[r] * unit_conv[r]
                else:
                    logger.warning("Concs not found for species %s for %s in %s", r, rct_list, rcn.to_rcn())
                    iconc = 0
                ircty *= iconc
            rct_data[rcts[s]] += ircty
            cnts[s] += 1  # Counter

    logger.info("Calculated reactivity for %d reactants w/ species counters: %s", len(cnts), cnts)

    return [rct_list, rct_data]


def calc_fgroup(concs: dict, species: list, nt: int, mode: str) -> dict:
    """Calculate the functional group concentrations based on the given mode"""

    if mode == "GECKO":
        return calc_fgroup_gck(concs, species, nt)

    raise ValueError(f"Invalid mode {mode} for functional group calculation")


def calc_fgroup_gck(concs: dict, species: list, nt: int) -> dict:
    """Calculate the functional group concentrations based on GECKO-A model"""
    # Init
    group_dict = {s: np.zeros(nt) for s in GCODE_DICT}  # Mass per functional group
    cmass, nsp = np.zeros(nt), 0  # Total carbon mass, number of species
    nc_dict = get_carbon_number(species)
    unused_note, ngp_sps = {".", "C"}, {}
    for s in species:
        if not (s.status > 0 and s.note and s.name in concs):
            continue
        sname, iconc, n = s.name, concs[s.name], 0
        if nc_dict.get(sname, 0) <= 0:  # No carbon number
            continue
        if s.note.startswith("-"):  # Invalid note
            break
        for g in s.note:
            if g == "_":  # Skip the rest
                break
            if g in unused_note:
                continue
            try:
                group_dict[g] += iconc * nc_dict[sname]
                n += 1
            except KeyError as e:
                raise ValueError(f"Invalid group {g} for species {sname}: {e}") from e
        if n == 0:
            ngp_sps[sname] = s.note
        cmass += iconc * nc_dict[sname]
        nsp += 1

    # Check notes
    if ngp_sps:
        logger.warning("No valid group found for # %d species w/ note: %s", len(ngp_sps), ngp_sps)

    # Normalize the group concentrations
    group_dict = {s: np.where(cmass > EPS, v / cmass, 0) for s, v in group_dict.items()}
    logger.info("Read # %d from # %d species for GECKOA functional groups", nsp, len(species))

    # Update group names and remove empty groups
    group_dict = {GCODE_DICT[k]: v for k, v in group_dict.items() if np.sum(v) > 0}

    return group_dict


def calc_ro2_rct(flag: Optional[str], rct_list: list, tconc: PostConc, mech: Mechanism, unit_in: str) -> dict:
    """Calculate the reactivity in the mechanism in the gas phase"""

    if flag is None:
        flag = "rsp"

    for rct in rct_list:
        if rct not in RCT_DICT:
            raise ValueError(f"Invalid reactivity species: {rct}")

    # Add flag checks
    if flag == "ro2":
        check_sps = {s for s in mech.sinfo["radical"] if mech.sinfo["org"][s].RO2}
    elif flag == "voc":
        check_sps = set(mech.sinfo["org"].keys()) - mech.sinfo["radical"]
    elif flag in mech.sinfo:
        check_sps = mech.sinfo[flag]
    elif flag == "all":  # All species
        check_sps = mech.sinfo["org"].keys() | mech.sinfo["inorg"].keys()
    else:
        raise ValueError(f"Invalid flag {flag} for reactivity calculation")

    logger.info("Calculate reactivity for %s w/ %d species", flag, len(check_sps))

    # Get unit conversion factors
    unit_conv = get_unit_conv_factors(mech.species, unit_in=unit_in, unit_out="molec")

    # Reactivity
    rcts = {s: i for i, s in enumerate(rct_list)}
    cnts = {s: 0 for s in rcts}
    nt = len(tconc.steps)
    concs = tconc.get_phase_concs("gas")
    rct_data = np.zeros((len(rcts), nt))
    for rcn in mech.reactions:
        if rcn.status < 1:
            continue
        r1s = [s for s in rcn.reactants if s in rcts]
        if not r1s:
            continue
        r2s = [s for s in rcn.reactants if s in check_sps]
        if not r2s:
            continue
        # Get the reactivity
        irate = rcn.rate.kvalues[0]
        for s in r1s:
            ircty = np.ones(nt) * irate * concs[s] * unit_conv[s]
            rct_data[rcts[s]] += ircty
            cnts[s] += 1  # Counter

    rct_data[:, 0] = 0.0
    logger.info("Calculated RO2 reactivity for %d reactants w/ species counters: %s", len(cnts), cnts)
    return [rct_list, rct_data]
