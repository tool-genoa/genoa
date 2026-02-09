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
This module contains functions to generate reduction information for printing.
"""

import json

from .constants import NOUT_LIM
from .reduction_basic import ReductionListInfo


def get_reduction_pinfos(strategy: str, rdcs: list, reactions: list, species: list, infos: ReductionListInfo) -> list:
    """
    Generate and return a formatted string of reduction information based on the strategy.
    """

    print_dcit = {
        "da": get_da_pinfos,
        "rs": get_rs_pinfos,
        "rm": get_rm_pinfos,
        "lp": get_lp_pinfos,
        "rp": get_rp_pinfos,
        "jp": get_jp_pinfos,
    }

    if strategy in ["rs", "da", "rp", "jp"]:
        return print_dcit[strategy](rdcs)
    if strategy == "rm":
        return print_dcit[strategy](rdcs, reactions)
    if strategy == "lp":
        return print_dcit[strategy](rdcs, species, infos.sps_info)
    if strategy == "none":  # No reduction
        return ["Test mechanism without changes."]

    raise ValueError(f"Unknown reduction strategy: {strategy}")


def get_simple_reduction_pline(strategy: str, rdcs: list) -> str:
    """Generate a simple reduction information string for printing."""

    nrdc = len(rdcs)
    if nrdc > NOUT_LIM:
        return f"# {nrdc}"

    if strategy in {"rs", "rm", "da"}:  # Removal-type strategies
        return f"# {nrdc}: {rdcs}"

    if strategy == "lp":  # Lumping
        sps = (f"{rt} {s}" for rdc in rdcs for rt, s in zip(rdc["ratios"], rdc["species"]))
    elif strategy == "rp":  # Replacement
        sps = (f"{r['species'][0]} -> " + " + ".join(r["species"][1:]) for r in rdcs)
    elif strategy == "jp":  # Jumping
        sps = (f"{r['jumped']} -> " + " + ".join(r["species"]) for r in rdcs)
    else:
        raise ValueError(f"Unknown reduction strategy: {strategy}")

    return f"# {nrdc}: " + ", ".join(sps)


def get_da_pinfos(rdcs: list) -> list:
    """Get print information for removal of gas-particle partitioning species."""

    return [f"Removed # {len(rdcs)} gas-particle partitioning:"] + rdcs


def get_rs_pinfos(rdcs: list) -> list:
    """Get print information for removal of species."""

    return [f"Removed # {len(rdcs)} species:"] + rdcs


def get_rm_pinfos(rdcs: list, reactions: list) -> list:
    """Get print information for removal of reactions."""

    return [f"Removed # {len(rdcs)} reactions:"] + [reactions[rdc].to_rcn() for rdc in rdcs]


def get_lp_pinfos(rdcs: list, species: list, sps_info: dict) -> list:
    """
    Get print information for lumping of species.
    shape of rdc:{
      'species': lumped_species,
      'ratios': ratios,
      'rcns': rcns
    }
    """
    outs = [f"Lumped # {len(rdcs)} species:"]
    for rdc in rdcs:
        header = (
            f"  {rdc['species'][0]} || "
            + " | ".join(f"{rt} {sp}" for rt, sp in zip(rdc["ratios"], rdc["species"]))
            + "\n"
        )

        # Output properties to check
        main_sp = species[sps_info[rdc["species"][0]]]
        keys = ["name", "status", "mass", "fgroups", "precursors"]
        if main_sp.condensable:
            keys.extend(["psat_atm", "dhvap_kj", "soap_strucs"])

        props = []
        for k in keys:
            lines = [f"{k}\t{json.dumps(getattr(main_sp, k))}\t||"]
            for sp in rdc["species"][1:]:
                lines.append(f"\t{json.dumps(getattr(species[sps_info[sp]], k))}")
            props.append("\t".join(lines))

        outs.append(header + "\n".join(props))
    return outs


def get_rp_pinfos(rdcs: list) -> list:
    """Get print information for replacement of species."""
    return [f"Replaced # {len(rdcs)} species:"] + [
        f"  {rdc['species'][0]} -> " + " + ".join(rdc["species"][1:]) for rdc in rdcs
    ]


def get_jp_pinfos(rdcs: list) -> list:
    """
    Get print information for jumping
    shape of rdc_info = {
      "jumped": s,
      "species": children,
      "ratios": children_ratios,
      "crcns": crcns.keys(),
      "prcns": infos.parent_info[s]["ircns"] if s in infos.parent_info else [],
    }
    """
    return [f"Jumped # {len(rdcs)} species:"] + [
        f"{rdc['jumped']} -> "
        + " + ".join(
            f"{rt:6.3E} {sp}" if isinstance(rt, float) else f"({rt}) {sp}"
            for rt, sp in zip(rdc["ratios"], rdc["species"])
        )
        + "\n"
        for rdc in rdcs
    ]
