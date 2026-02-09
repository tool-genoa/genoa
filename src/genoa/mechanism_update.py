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
This module provides functions to update the mechanism.
"""

from typing import Optional

import attrs

from .constants import SNAME
from .logger import setup_logger
from .mechanism_relation import update_sps_list_with_gen
from .reaction import Reaction, get_reduction_info_from_list, with_norate_reactions
from .setting_global import get_all_settings_map
from .species import add_basic_to_species_list, update_basic_dict
from .utils import calculate_normalized_ratios


# Logger
logger = setup_logger(__name__)


@attrs.define(slots=True)
class MechInfo:
    """Size information of a mechanism"""

    sdict: dict = attrs.field(factory=dict)  # {name: ref} of all species (status >= 0)
    bdict: dict = attrs.field(factory=dict)  # {name: MW} of basic species (status > 1)

    nsps: int = 0  # Number of species (status > 0)
    nrsp: int = 0  # Number of reducible species (status = 1)
    naero: Optional[int] = None  # Number of condensable reducible species
    nrcn: Optional[int] = None  # Number of reactions (status > 0)
    nrrn: Optional[int] = None  # Number of reducible reactions (status = 1)
    rm_rcn: int = 0  # Number of removed reactions
    rm_sps: int = 0  # Number of unsorted species
    pvocs: Optional[list] = None  # List of primary VOCs

    is_valid: bool = True  # Is the mechanism valid?
    is_changed: bool = False  # Has the mechanism been changed?

    def __attrs_post_init__(self):
        """Post initialization."""

        if not (self.sdict and self.bdict):
            raise ValueError("Species dictionary (sdict) and basic species dictionary (bdict) must be provided.")

        # Initialize counts
        self.update_counts_wdict()

    def update_with_used_sps(self, sps_used: set) -> int:
        """Update mechanism information based on used species extracted from reaction list."""

        # Check input used species set
        n_used = len(sps_used)
        if n_used <= len(self.pvocs):
            logger.warning("Only used species: %s containing primary VOCs: %s.", sps_used, self.pvocs)
            self.is_valid = False
            return 0

        # Find non-used reducible species
        unsps = [s for s in self.sdict if s not in sps_used and s not in self.bdict]
        if not unsps:
            return 0  # No unsorted species found

        # Remove unsorted species
        logger.info("Removing # %d unsorted species: %s.", len(unsps), unsps)
        for s in unsps:
            isp = self.sdict.pop(s)
            isp.status = 0  # Set status to 0 (not active)
            if isp.condensable:
                self.naero -= 1

        # Update counts
        self.update_counts_wdict()

        # Check if the number of used species matches the expected number
        if n_used != self.nrsp:
            non_used = [s for s in sps_used if s not in self.sdict]
            logger.warning("Found # %d non-used species: %s.", len(non_used), non_used)
            self.is_valid = False
            raise ValueError(f"{n_used} != {self.nrsp}.")

        return len(unsps)  # Number of removed species

    def update_with_sps_nogen(self, nogens: list) -> None:
        """Update mechanism information based on species not used in generation."""

        if not nogens:
            return

        # Remove species not used in generation
        nogens = [s.name for s in nogens]
        for s in nogens:
            isp = self.sdict.pop(s)
            if isp.condensable:
                self.naero -= 1

        # Update counts
        self.update_counts_wdict()
        logger.info("Removed # %d species not used in generation: %s.", len(nogens), nogens)
        self.rm_sps += len(nogens)

    def update_counts_wdict(self) -> None:
        """Update species counters based on the current species dictionary."""
        self.nsps = len(self.sdict)
        self.nrsp = self.nsps - len(self.bdict)
        if self.nrsp <= 0:
            # logger.warning("Number of reducible species is %s, must be > 0.", self.nrsp)
            self.is_valid = False
        if self.naero is None:  # Count condensable species
            self.naero = sum(1 for s in self.sdict.values() if s.condensable and s.status == 1)
        if self.naero < 0:
            # logger.warning("Number of condensable species is %s, must be >= 0.", self.naero)
            self.is_valid = False

    def update_removed(self, rm_rcn: int = 0, rm_sps: int = 0) -> None:
        """Update the number of removed reactions and species."""
        self.rm_rcn += rm_rcn
        self.rm_sps += rm_sps

    def update_rcn_counts(self, reactions: list) -> None:
        """Update the number of reactions based on the provided reaction list."""

        if not reactions:
            raise ValueError("No reactions provided to update nrcn.")

        nrcn, nrrn = 0, 0  # Reset counts
        for rcn in reactions:
            if rcn.status > 0:
                nrcn += 1
                if rcn.status == 1:
                    nrrn += 1

        if nrrn < 1:
            # logger.warning("Number of reducible reactions is %s, must be > 0.", nrrn)
            self.is_valid = False

        self.nrcn, self.nrrn = nrcn, nrrn

    def read_size(self) -> list:
        """Ouput the size of the mechanism."""

        if self.nrcn is None:
            raise ValueError("Can not read reaction size info.")

        # Output size information: [nrcn, nrsp, naero]
        return [self.nrrn, self.nrsp, self.naero]

    def check_change(self) -> None:
        """Update info if changed"""

        if max(self.rm_rcn, self.rm_sps) > 0:
            self.is_changed = True

    def get_pinfo(self) -> str:
        """Return a string with the mechanism information summary."""

        if self.is_changed:
            pinfo = ["Mechanism after update:"]
        else:
            pinfo = ["Mechanism not changed:"]

        pinfo.append(f"  Species: no. valid / reducible (aerosol): {self.nsps} / {self.nrsp} ({self.naero})")
        pinfo.append(f"  Reactions: no. valid / reducible: {self.nrcn} / {self.nrrn}")

        if self.rm_rcn > 0 or self.rm_sps > 0:
            pinfo.append(f"  Updated w/ removing {self.rm_rcn} reactions and {self.rm_sps} species.")
        if not self.is_valid:
            pinfo.append("!!!Updated mechanism is NOT valid!!!")

        return "\n".join(pinfo)


def update_mechanism(reactions: list, species: list, opts_in: dict) -> list:
    """
    Clean and merge mechanism by removing unused species & reactions and merging reactions.

    Available options:
    - "clean": "wsps", "wpvoc", "wgen"  # Clean modes
    - "check_rate": True  # Check if any reaction has no rate constant
    - "check_kept": set()  # Check if any species in the kept set is missing
    - "merge": True  # Merge reactions
    - "split": True  # Split reactions with multiple products

    Return a list containing:
    - Updated reactions
    - Updated species
    - Mechanism information (MechInfo)
    """

    # Check if empty
    if not (reactions and species):
        raise ValueError("No reactions or species read or no options provided.")

    # Copy options
    opts = opts_in.copy()

    # Clean mechanism
    minfo = clean_mechanism(reactions, species, opts.pop("clean", "wsps"))

    # Check if no rate constants
    if opts.pop("check_rate", False):
        if with_norate_reactions(reactions):
            minfo.is_valid = False
            return [reactions, species, minfo]

    # Check if missing species
    check_kept = opts.pop("check_kept", None)
    if check_kept:
        if _missing_species(minfo.sdict, check_kept):
            minfo.is_valid = False
            return [reactions, species, minfo]

    # Reaction merging
    if opts.pop("merge", False):
        nrcn = len(reactions)
        reactions = reactions_merging(reactions)
        logger.info("Merged reactions, %d -> %d.", nrcn, len(reactions))
        if nrcn != len(reactions):
            minfo.is_changed = True

    # Reaction separation
    if opts.pop("split", False):
        nrcn = len(reactions)
        reactions = reactions_separation(reactions, species)
        logger.info("Separated reactions, %d -> %d.", nrcn, len(reactions))
        if nrcn != len(reactions):
            minfo.is_changed = True

    # Update tracer species if needed
    tracer_to_check = opts.pop("tracer", None)
    if tracer_to_check:
        logger.info("Checking tracer species: %s.", tracer_to_check)
        add_tracer_to_mech(reactions, species, tracer_to_check, minfo)

    # Check if changed at the end after all updates
    minfo.check_change()

    # Update reaction counts
    minfo.update_rcn_counts(reactions)

    # Output mechanism info to check
    logger.info("%s", minfo.get_pinfo())

    # Last check
    if opts:
        logger.warning("Unprocessed mechanism update options: %s.", opts)

    return [reactions, species, minfo]


def clean_mechanism(reactions: list, species: list, clean_mode: str) -> MechInfo:
    """
    Clean mechanism with the given mode.
    Modes:
    - "wsps": Remove inactive species and reactions
    - "wpvoc": Remove reactions not from primary VOCs
    - "wgen": Update generation and remove reactions not from primary VOCs
    """

    check_func = {
        "wsps": valid_scheme,  # Remove inactive species and reactions
        "wpvoc": valid_scheme_from_source,  # 1 + species/reactions only from primary VOCs
        "wgen": valid_scheme_with_new_gen,  # 2 + updated generation
    }
    if clean_mode not in check_func:  # Check mode
        raise ValueError(f"Unknown clean mode {clean_mode}. Use: {check_func.keys()}.")

    return check_func[clean_mode](reactions, species)


def update_species_with_reactions(reactions: list, species: list) -> list:
    """
    Return species dict and basic dict and number of species (all, reducible, condensable)
    NOTE:: basic species dict is updated.
    """

    # Get settings
    gnl = get_all_settings_map()[SNAME.SETUP]()
    basic_dict = gnl.basicspecies_dict

    # Read species info from species list
    sps_dict, new_basic, naero = {}, {}, 0
    for s in species:
        if s.status <= 0:
            continue
        sps_dict[s.name] = s
        if s.status > 1:
            new_basic[s.name] = s.mass
        elif s.condensable:
            naero += 1

    # Check and update basic species dict if needed
    added = update_basic_dict(new_basic, basic_dict)
    ncycle = 0  # Cycle counter for updating basic species

    while added > 0:
        # No.cycle limit
        if ncycle > 5:
            raise ValueError(f"Too many cycles to update basic species, currently {ncycle} and {added} added.")
        ncycle += 1  # Update cycle counter

        # Update species dict
        add_basic_to_species_list(species, sps_dict, basic_dict)

        # Get new & update basic species dict
        new_basic = _get_new_basic_species(reactions, sps_dict, basic_dict)
        if not new_basic:
            break
        added = update_basic_dict(new_basic, basic_dict)

    if ncycle > 0:
        logger.info("Finished update basic species dict with # %s cycle(s).", ncycle)

    # Build species list info for output
    minfo = MechInfo(sdict=sps_dict, bdict=basic_dict, naero=naero, pvocs=gnl.primary_vocs)
    return minfo


def _get_new_basic_species(reactions: list, sps_dict: dict, basic_dict: dict) -> dict:
    """Get new basic species from the reaction list."""
    new_basic = {}  # New basic species
    for rcn in reactions:
        if rcn.status <= 0:
            continue
        tag_basic = True
        # Check reactants
        if rcn.status == 1:
            for s in rcn.reactants:
                if s not in basic_dict:
                    tag_basic = False
                    break

        # Update species and reactions
        if tag_basic:
            rcn.status = 2
            for s in rcn.products:
                if s not in basic_dict:
                    new_basic[s] = sps_dict[s].mass

    if new_basic:
        logger.info("Find # %s new basic species from the species list: %s.", len(new_basic), new_basic)

    return new_basic


def valid_scheme(reactions: list, species: list, minfo: Optional[MechInfo] = None) -> MechInfo:
    """Trim scheme with removing/modifying inactivated species. No check on primary VOCs or generation."""

    # Update species info if not provided
    if minfo is None:
        minfo = update_species_with_reactions(reactions, species)
    sps_dict, sps_basic = minfo.sdict, minfo.bdict

    # Loop over reactions
    nrm_rcn = 0  # Check if mechanism is changed
    for rcn in reactions:

        # Skip non-reducible reactions
        if rcn.status != 1:
            continue

        # Return [0] for inactive reactions or [1, rcts, pdts] for valid
        rinfo = _update_reaction_and_return_info(rcn, sps_dict, sps_basic)

        # Record inactive reactions
        if rinfo[0] == 0:
            nrm_rcn += 1

    # Update removed reactions count
    minfo.update_removed(nrm_rcn)

    logger.info("Updated scheme.")
    return minfo


def valid_scheme_with_new_gen(reactions: list, species: list, minfo: Optional[MechInfo] = None) -> MechInfo:
    """Trim scheme with removing/modifying inactivated species not from primary vocs (check w/ updated generation)."""

    if minfo is None:
        minfo = update_species_with_reactions(reactions, species)

    sps_nogen = update_sps_list_with_gen(reactions, species)
    minfo.update_with_sps_nogen(sps_nogen)

    return valid_scheme(reactions, species, minfo)


def valid_scheme_from_source(reactions: list, species: list, minfo: Optional[MechInfo] = None) -> MechInfo:
    """Trim scheme with removing/modifying inactivated species that can not be derived from primary vocs."""

    # Get species info
    if minfo is None:
        minfo = update_species_with_reactions(reactions, species)
    sps_dict, sps_basic, pvocs = minfo.sdict, minfo.bdict, minfo.pvocs

    if not pvocs:
        logger.warning("No primary VOCs found for valid_scheme_from_source. Use valid_scheme instead.")
        valid_scheme(reactions, species, minfo)

    # Loop over reactions
    nrm_rcn = 0
    species_used, sps_rcn_rel = set(pvocs), {}

    for ircn, rcn in enumerate(reactions):
        if rcn.status != 1:
            continue
        rinfo = _update_reaction_and_return_info(rcn, sps_dict, sps_basic)

        if rinfo[0] == 0:  # Inactive reaction
            nrm_rcn += 1
            continue

        # Record relationship between species and reactions
        rcts, pdts = rinfo[1], rinfo[2]
        if rcts <= species_used:  # All reactants are used species
            check_sps = pdts.copy()  # Record products
            while check_sps:
                s = check_sps.pop()  # Check species s

                if s in sps_rcn_rel:  # Check products from s if recorded
                    to_remove = []
                    for i in sps_rcn_rel[s]:
                        if s == i:  # s is the only reactant
                            check_sps.update(sps_rcn_rel[s][i][0])
                            to_remove.append(i)
                        else:  # need to check all reactants
                            # reactants not in species_used
                            p = [j for j in i.split(" ") if j != s and j not in species_used]
                            if not p:  # All reactants are in species_used
                                check_sps.update(sps_rcn_rel[s][i][0])
                                to_remove.append(i)

                    # Clean sps_rcn_rel: remove processed items
                    for i in to_remove:
                        del sps_rcn_rel[s][i]
                    # Clean sps_rcn_rel: remove empty items
                    if not sps_rcn_rel[s]:
                        del sps_rcn_rel[s]

                # Add s to species_used
                if s not in species_used:
                    species_used.add(s)

        # Record rcts -> pdts relationship in sps_rcn_rel
        else:
            # identifer of this set of reactants
            s = " ".join(sorted(rcts))
            for i in rcts:
                if i not in sps_rcn_rel:
                    sps_rcn_rel[i] = {}
                # shape: {rct: {rcts: [pdts, ircns]}}
                if s in sps_rcn_rel[i]:
                    sps_rcn_rel[i][s][0].update(pdts)  # set of products
                    sps_rcn_rel[i][s][1].append(ircn)  # list of reaction indices
                else:
                    sps_rcn_rel[i][s] = [pdts, [ircn]]  # [set, list]

    # Update species list with used species
    nrm_sps = minfo.update_with_used_sps(species_used)

    # Check if all reactions are derived from primary vOCs
    nrm_unsrc = 0  # No. of reactions not derived from primary VOCs
    for s, vals in sps_rcn_rel.items():
        if not vals:
            continue
        for p, val in vals.items():  # rcts: [pdts, ircns]
            for ircn in val[1]:  # Get reaction index
                # Process reaction
                rcn = reactions[ircn]
                if rcn.status == 1:
                    rcn.status = 0
                    nrm_unsrc += 1

    # Update removed counts
    minfo.update_removed(nrm_rcn + nrm_unsrc, nrm_sps)

    logger.info("Updated scheme from primary VOCs.")
    return minfo


def _missing_species(sps_dict: dict, kept_set: set) -> bool:
    """Check if any species in kept_set is missing from sps_dict."""

    if not sps_dict or not kept_set:
        return False

    unkept = [s for s in kept_set if s not in sps_dict]

    if unkept:
        logger.warning("Found # %d unkept species: %s", len(unkept), unkept)
        return True

    return False


def _update_reaction_and_return_info(rcn: Reaction, sdict: dict, sbasic: dict) -> list:
    """Update a reaction based on the species dictionary and return reactant and product lists if valid."""

    if rcn.status != 1:  # Skip non-reducible reactions
        return [0]

    # Initialize for checking
    rcts, pdts, rm_sps, tag_basic = set(), set(), [], True

    # Check reactants: if not in species list, inactivate reaction
    for s in rcn.reactants:
        if s not in sdict:  # Not found in species list
            rm_sps.append(s)
            break
        # Find reducilbe species
        if s not in sbasic:
            rcts.add(s)
            tag_basic = False
    # Remove reaction as found inactived reactants
    if rm_sps:
        rcn.status = 0
        return [0]

    # Mark reaction as basic if all reactants are basic species
    if tag_basic:
        rcn.status = 2
        rinfo = f"Found basic reaction: {rcn.to_rcn()}"
        for s in rcn.products:  # Check products
            if s not in sbasic:
                rinfo += f"\n  !!! Non-basic product: {s}"
        return [0]

    # Check products: if not in species list, remove the species
    for s in rcn.products:
        if s not in sdict:
            rm_sps.append(s)
        elif s not in sbasic:
            pdts.add(s)
    if rm_sps:  # Remove inactived products
        for s in rm_sps:
            del rcn.products[s]

    rcn.clean_duplicates()  # Clean duplicates in the reaction
    if rcn.status == 0:
        return [0]

    return [1, rcts, pdts]  # Return reactants and products if valid


def balance_carbon_in_reactions(reactions: list, species: list, csname: str = "XCLOST") -> None:
    """Balance carbon loss in the reaction list with the given carbon species."""

    # Get species dict
    sps_dict = {s.name: s for s in species if s.status >= 0}
    if csname not in sps_dict:
        raise ValueError(f"Can not find carbon species: {csname} in the species list.")

    # Get number of carbon atoms in the carbon species
    n = sps_dict[csname].fgroups.get("C", 0)
    if n == 0:
        n = csname.count("C")  # Count C from the name
    if n == 0:
        raise ValueError(f"Can not find carbon atoms in the carbon species: {csname}.")
    logger.info("Balancing carbon loss with species: %s with %s carbon atoms ...", csname, n)

    nrcn = 0
    for rcn in reactions:
        if not rcn.status:
            continue
        loss, gain = 0, 0  # Carbon loss and gain
        # Get carbon loss
        for s in rcn.reactants:
            if s not in sps_dict:
                continue
            c = sps_dict[s].fgroups.get("C", 0)
            if c == 0:
                c = s.count("C")
            loss += c
        if loss == 0:
            continue  # No carbon loss
        # Get carbon gain
        if rcn.rt_fnc is None:
            for s, rt in rcn.products.items():
                if s not in sps_dict:
                    continue
                c = sps_dict[s].fgroups.get("C", 0)
                if c == 0:
                    c = s.count("C")
                gain += c * rt
        else:  # Get carbon gain from the reaction function
            pdts = list(rcn.products.keys())
            rts = rcn.get_ratios(slist=pdts)
            for i, s in enumerate(pdts):
                if s not in sps_dict:
                    continue
                c = sps_dict[s].fgroups.get("C", 0)
                if c == 0:
                    c = s.count("C")
                gain += c * rts[i]

        # Check balance & update product
        loss, gain = round(loss, 5), round(gain, 5)
        if loss == gain:
            continue
        # Add carbon species to reaction
        if loss > gain:
            rcn.products[csname] += (loss - gain) / n
            nrcn += 1  # Update counter
        else:
            logger.warning("Carbon loss %s > gain %s: %s", loss, gain, rcn.to_rcn())

    logger.info("Balanced # %s reactions with species: %s.", nrcn, csname)


def add_tracer_to_mech(reactions: list, species: list, tracer_list: list, minfo: Optional[MechInfo] = None) -> None:
    """Update gain and loss tracers to the reaction list."""

    if not tracer_list or not reactions or not species:
        return  # No tracers or reactions or species to update

    # Update species list
    tracers = add_tracer_to_species(reactions, species, tracer_list, minfo)
    if not tracers:
        logger.info("No tracers to update in reactions, skip.")
        return

    # Update reactions
    add_tracer_to_reactions(reactions, tracer_list[0], tracer_list[1])


def add_tracer_to_species(reactions: list, species: list, tracer_list: list, minfo: Optional[MechInfo]) -> set:
    """Add gain and loss tracers to the species list and return set of all tracers."""

    if minfo is None:
        minfo = update_species_with_reactions(reactions, species)
    sps_dict, sps_basic = minfo.sdict, minfo.bdict

    tracers = set()  # Set of all tracers
    for sdict in tracer_list:
        for k, v in sdict.items():
            if k not in sps_dict:
                continue
            isp = sps_dict[k]  # Target species
            if v not in sps_dict:  # Add new tracer species
                jsp = isp.copy()
                jsp.name, jsp.status = v, 2
                sps_dict[v], sps_basic[v] = jsp, jsp.mass  # Update species dict
                species.append(jsp)  # Add to species list
                logger.info("Added tracer species: %s with status %s.", v, jsp.status)

            elif sps_dict[v].status != 2:  # Update status
                logger.error("Tracer species %s exists with status %s, update to basic!", v, sps_dict[v].status)
                sps_dict[v].status, sps_basic[v] = 2, jsp.mass  # Update status to 2

            tracers.add(v)

    return tracers


def remove_tracer_from_reactions(reactions: list, tracers: Optional[set]) -> None:
    """Remove tracers from reactions."""

    if not tracers or not reactions:
        return  # No tracers or reactions to remove

    # Check and remove tracers from products
    n_removed = 0
    for rcn in reactions:
        for k in list(rcn.products):
            if k in tracers:
                rcn.products.pop(k, None)
                n_removed += 1

    if n_removed > 0:
        logger.info("Removed # %d tracers from # %d reactions.", n_removed, len(reactions))


def add_tracer_to_reactions(reactions: list, gain: dict, loss: dict) -> None:
    """Update reactions with gain and loss tracers."""

    if not (gain or loss):
        return

    # Ensure updating all valid tracers
    sps_to_remove = set(gain.values()).union(loss.values())
    n_removed = 0

    # Update reactions list
    for rcn in reactions:

        if rcn.status != 1:  # Only update reducible reactions
            continue

        # Record tracers
        old = {s for s in rcn.products if s in sps_to_remove}

        # Update loss tracers
        for k in list(rcn.products):
            if k in gain:
                s = gain[k]  # Gain tracer
                rcn.products[s] = rcn.products[k]
                old.discard(s)

        # Update loss tracers
        for k, v in rcn.reactants.items():
            if k in loss:
                s = loss[k]  # Loss tracer
                rcn.products[s] = v
                old.discard(s)

        # Remove old tracers
        if old:
            n_removed += 1
            for s in old:
                rcn.products.pop(s, None)

    if n_removed > 0:
        logger.info("Found and removed invalid tracers from # %d reactions.", n_removed)


def reactions_separation(reactions: list, species: list) -> list:
    """Rewrite reactions with multiple products into elementary-like reactions with single products."""

    new_reactions = []  # For output
    rsps = {s.name for s in species if s.status == 1}  # Get reducible species

    for rcn in reactions:
        if rcn.status != 1 or len(rcn.products) <= 1:
            new_reactions.append(rcn)
            continue

        # Find products need to be separated
        pdts = [s for s in rcn.products if s in rsps]

        # No need or can not separate
        if len(pdts) <= 1:
            new_reactions.append(rcn)
        else:
            new_reactions.extend(split_a_reaction(rcn, pdts))

    return new_reactions


def split_a_reaction(rcn: Reaction, split_products: list) -> list:
    """Split a reaction with multiple products into several reactions with single products."""

    if len(split_products) <= 1:  # No need to split
        return [rcn]

    # Get separation ratios
    rts = rcn.get_ratios(slist=split_products)
    nrts = calculate_normalized_ratios(rts, 5)  # Normalize ratios
    if min(rts) <= 0.0 or min(nrts) <= 0.0:
        logger.warning("No split - Negative ratio %s in reaction: %s", rts, rcn.to_rcn())
        return [rcn]

    # Build new reactions
    new_rcns = []
    # n = len(split_products)

    for s, irt, nirt in zip(split_products, rts, nrts):

        nrcn = rcn.copy()
        # Update products
        for p, prt in rcn.products.items():
            if p == s:
                nrcn.products[p] = irt / nirt
            elif p in split_products:
                nrcn.products.pop(p, None)
            else:
                nrcn.products[p] = prt  # / nirt / n
        # Change kinetic rate with the ratio
        if not nrcn.rate.multiply_rate_by_ratio(nirt):
            logger.warning("Inactive reaction after split: %s, ratio %s.", nrcn.to_rcn(), nirt)
            nrcn.status = 0  # Inactivate reaction if invalid ratio
        new_rcns.append(nrcn)  # Add new reaction

    return new_rcns


def reactions_merging(reactions: list) -> list:
    """Merge reactions with the same reduction information in the reaction list."""

    # Extract reduction info
    rdc_index_info = get_reduction_info_from_list(reactions)

    # Build new reaction list
    new_reactions = []
    # Processed reaction index
    processed = set()

    for i, rcn in enumerate(reactions):

        # Non-reducible
        if rcn.status != 1:
            new_reactions.append(rcn)
            continue

        # Processed reaction
        if i in processed:
            continue

        # Find mergable reactions
        ircns = rdc_index_info[rcn.rcn_id]

        # Only one reaction, no merging
        if len(ircns) <= 1:
            new_reactions.append(rcn)
            continue

        # Merge reactions
        new_rcns = _merge_similar_reactions([reactions[j] for j in ircns])

        # Update new & processed reactions
        new_reactions.extend(new_rcns)
        processed.update(ircns)

    return new_reactions


def _merge_similar_reactions(rcns: list) -> list:
    """Merge reactions based on the first input reaction."""

    if len(rcns) <= 1:  # No merging
        return rcns

    # Get reduction info for the 1st reaction
    rcn_id0 = rcns[0].rcn_id

    # Get products and ratios
    others, pdts, rates = [], {}, []
    for i, rcn in enumerate(rcns):

        # Check type
        if rcn.rcn_id != rcn_id0:
            logger.warning("Reactions with different rinfo: %s != %s, skip merging.", rcn_id0, rcn.rcn_id)
            others.append(rcn)
            continue

        # Get ratio
        rt = rcn.rate.get_ratio()
        rates.append(rt)

        # Record products
        for s in rcn.products:
            if s not in pdts:
                pdts[s] = [[], []]
            pdts[s][0].append(rcn.products[s])  # Product ratio
            pdts[s][1].append(rt)  # Kinetic ratio

    # No merging
    if len(rates) <= 1:
        return rcns

    # Build merged reaction from 1st reaction
    new_rcn = rcns[0].copy()
    sum_rate = sum(rates)
    new_rate = sum_rate / rates[0]

    # Update kinetic rate
    if not new_rcn.rate.multiply_rate_by_ratio(new_rate):
        logger.warning("Inactive reaction after merging: %s, ratio %s.", new_rcn.to_rcn(), new_rate)
        new_rcn.status = 0

    # Update products
    new_rcn.reset_products()
    for s, [prts, rts] in pdts.items():
        # New ratio for product s
        new_rt = sum(prts[i] * v for i, v in enumerate(rts)) / sum_rate
        # Add product
        new_rcn.products[s] += new_rt

    return [new_rcn] + others if others else [new_rcn]
