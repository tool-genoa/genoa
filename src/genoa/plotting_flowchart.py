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
This module contains functions to plot the flowchart of the mechanism.
"""


import os
from typing import Optional

import graphviz
from cairosvg import svg2png
from openbabel import pybel

from .constants import PSAT_BIN_LIMS
from .folder_path import update_output_path
from .logger import setup_logger
from .mechanism_pack import Mechanism
from .reaction import Reaction
from .setting_init import PostProcessOption
from .species import Species


# Logger
logger = setup_logger(__name__)


# Color and shape for flowchart
# More details see: http://www.graphviz.org/doc/info/colors.html
RCOLS = {
    "O3": "blue",
    "OH": "red",
    "HO": "red",
    "NO3": "green",
    "NO2": "yellow3",
    "NO": "orange",
    "RO2": "gray65",
    "HO2": "cyan3",
    "H2O": "orchid3",
    "PHOT": "red4",
    "ISOM": "darkcyan",
    "O2": "deeppink",
    "CO": "black",
    "SO2": "black",
    "CL": "black",
    "VOC": "white",
    "Radical": "white",
    "SVOC": "gray80",
    "LVOC": "gray55",
    "ELVOC": "black",  # gray30
    "NVOC": "black",
}

# Arrowhead
RARRS = {
    "O3": "vee",
    "OH": "odiamond",
    "HO": "odiamond",
    "NO3": "onormal",
    "NO2": "box",
    "NO": "dot",
    "RO2": "diamond",
    "HO2": "odot",
    "H2O": "tee",
    "PHOT": "obox",
    "ISOM": "inv",
    "O2": "invempty",
    "CO": "normal",
    "SO2": "normal",
    "CL": "normal",
}

# Psat_atm upper limits for SVOCs, LVOCs, ELVOCs
RSHPS = {"VOC": "ellipse", "Radical": "plain", "SVOC": "box", "LVOC": "box", "ELVOC": "hexagon", "NVOC": "octagon"}

# Box penwidth
PENS = {"SVOC": "1", "LVOC": "1", "ELVOC": "2", "NVOC": "2"}


def plot_flowchart(mechs: list, opt: PostProcessOption) -> None:
    """Plot the flowchart of the mechanism."""
    if not (opt.plot_flowchart and mechs):  # No flowchart or mechanism
        return

    logger.info("Plotting flowchart for the mechanism w/ input: %s ...", opt.plot_flowchart)
    savpath = update_output_path(os.path.join(opt.savpath, "flowchart"))
    for items in opt.plot_flowchart:
        if len(items) != 3:
            raise ValueError(f"Invalid input for flowchart: {items}. Should be a list of 3 strings.")
        strin, rngs, wfig = items
        strin = strin.strip()
        if strin == "pvocs":
            initsps = opt.plots["pvocs"]
        else:
            initsps = [s.strip() for s in strin.split(",")]
        if not initsps:
            raise ValueError("No initial species for the flowchart.")
        wfdir = os.path.join(opt.savpath, "sfigs") if wfig else None
        dpths = [int(s) for s in rngs.split(",")] if rngs else None
        # Plot for all mechanisms
        for mech in mechs:
            filename = os.path.join(savpath, f"{mech.name}_{strin}")
            logger.info("Plotting flowchart for mech %s w/ settings: %s ...", mech.name, items)
            gens = viz_fig_from_chem(mech, initsps, filename, dpths, wfdir)
            if opt.logt and gens:
                opt.logt.list_to_file(
                    [f"{k} w/ No.{len(v)}: {v}" for k, v in gens.items() if v],
                    f"\n{mech.name}: flowchart w/ settings {items}",
                )


def create_label_with_image(image_path: str, text: str, fontcolor="black", penwidth="0") -> str:
    """Create a label with an image and text in a table format."""
    return f"""<<TABLE BORDER="{penwidth}" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">
    <TR><TD><IMG SRC="{image_path}" SCALE="FALSE"/></TD></TR>
    <TR><TD><FONT COLOR="{fontcolor}">{text}</FONT></TD></TR></TABLE>>"""


def update_radical_smiles(strin: str) -> str:
    """Update the radical SMILES string: C. -> [C], Cn. -> [C]n."""

    if "(.)" in strin:
        strin = strin.replace("(.)", ".")
    if not strin or "." not in strin:  # No radical
        return strin
    if len(strin) < 3 or strin.startswith("."):  # Invalid radical
        raise ValueError(f"Unexpected SMILES string {strin} for radical.")

    # Change the radical SMILES string
    for i, s in enumerate(strin):
        if s != ".":
            continue
        # Check the previous character
        p = strin[i - 1]
        if p.isdigit():
            if i > 1 and strin[i - 2] in ["C", "O"]:
                strin = f"{strin[:i-2]}[{strin[i-2]}]{p}{strin[i+1:]}"
            else:
                raise ValueError(f"Unexpected character {p} before . in {strin}")
        elif p in ["C", "O"]:
            strin = f"{strin[:i-1]}[{p}]{strin[i+1:]}"
        else:
            raise ValueError(f"Unexpected character {p} before . in {strin}")
    return strin


def fetch_chemical_pic(smiles_str: str, savpng: str) -> str:
    """Generate PNG and SVG file for a given SMILES string. Return the PNG file path."""

    if not (smiles_str and savpng and savpng.endswith(".png")):
        logger.warning("Invalid input for chemical picture: %s, %s", smiles_str, savpng)

    if os.path.isfile(savpng):  # Already exist
        return savpng

    # Update radical SMILES
    smiles_str = update_radical_smiles(smiles_str)

    # SMILES -> Pybel molecule
    pybel_mol = pybel.readstring("smi", smiles_str)
    pybel_mol.make2D()

    # Pybel molecule -> SVG
    svg = pybel_mol.write("svg")
    svg_content = svg.split("<!--")[0]  # Removing unnecessary parts
    with open(savpng[:-4] + ".svg", "w", encoding="utf-8") as f:
        f.write(svg_content)

    # SVG -> PNG
    svg2png(bytestring=svg_content, write_to=savpng)
    # Trim PNG
    trim_pic_with_convert(savpng)

    return savpng


def trim_pic_with_convert(image_path: str) -> None:
    """Trim the white space of the image with convert command."""
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"Cannot find the image file {image_path}")
    os.system(f"convert {image_path} -trim {image_path}")


def get_rcn_layout(rcn: Reaction) -> list:
    """Get the arrowhead and color for the input reaction."""

    rcts = [r for r in rcn.reactants if r in RARRS]
    if not rcts:  # Check kinetic rate
        flag = False
        for s in ["RO2", "H2O", "PHOT"]:
            if s in rcn.rate.ssh:
                irct, flag = s, True
                break
        if not flag:  # Check string for ISOM
            for s, v in {"ISOM": "ISOM", "OXYGEN": "O2"}.items():
                if s in rcn.rate.string:
                    irct, flag = v, True
                    break
        if not flag:
            logger.warning("Unknown type of reaction: %s w/ ssh %s", rcn.to_rcn(), rcn.rate.ssh)
            irct = "CO"  # Black
    elif len(rcts) == 1:
        irct = rcts[0]
    else:
        raise ValueError(f"Multiple reactants in reaction {rcn.to_rcn()}")

    # Return arrowhead and color
    return [RARRS[irct], RCOLS[irct]]


def get_sps_layout(sp: Species, sdict: dict) -> list:
    """Get layout for the input species."""

    ipen = "0"  # Default linewidth in string
    if sp.radical:
        return [sdict["Radical"], ipen]
    if not sp.condensable:
        return [sdict["VOC"], ipen]
    for k, v in PSAT_BIN_LIMS.items():
        if sp.psat_atm > v:
            return [sdict[k], PENS.get(k, ipen)]
    raise ValueError(f"Species {sp.name} w/ psat {sp.psat_atm} not found in bin: {PSAT_BIN_LIMS}")


def species_reaction_relations(reactions: list, orgs: dict) -> list:
    """Get the parent and children information for the input reaction list."""
    up_info, down_info = {}, {}  # Species: set(ircns)
    for ircn, rcn in enumerate(reactions):
        if rcn.status != 1:
            continue
        for s in [r for r in rcn.reactants if r in orgs]:
            if s not in down_info:
                down_info[s] = {ircn}
            else:
                down_info[s].add(ircn)
        for s in [p for p in rcn.products if p in orgs]:
            if s not in up_info:
                up_info[s] = {ircn}
            else:
                up_info[s].add(ircn)
    return [up_info, down_info]


def viz_fig_from_chem(mech: Mechanism, initsps: list, svnm: str, dpths: Optional[list], wfdir: Optional[str]) -> dict:
    """Visualizes chemical reactions and species from given data."""

    # Settings for flowchart - adjustable
    nmax_wofig, nmax_wfig = 100, 100  # Maximum number of nodes w/ & w/o figures
    nsep = 0.3  # Node separation: horizontal distance between nodes
    rsep = 0.5  # Rank separation: vertical distance between ranks
    lid = "m"  # notation for lumped species

    if not initsps:
        raise ValueError("No initial species for the flowchart.")

    def _sname(sname: str) -> str:
        """return real name for lumped/non-lumped species"""
        if sname.startswith(lid):  # m for lumped!
            return sname[1:]
        return sname

    def _node_id(sname: str) -> str:
        """return node id"""
        base = _sname(sname)
        if base not in sps_info:
            raise ValueError(f"Not found species {base} in sps_info. Base species may be named after {lid}.")
        isp_base = sps_info.get(base)
        if isp_base.reductions.get("lp"):
            return f"{lid}{base}"
        return base

    # Plot w/ figure saved to wfdir
    if wfdir and isinstance(wfdir, str):
        wfdir = update_output_path(wfdir)
        wfig = True
        logger.info("Flowchart with figures will be saved to %s", wfdir)
    else:
        wfig = False
    nmax = nmax_wfig if wfig else nmax_wofig

    # Get depth in the reaction tree
    if dpths:
        if len(dpths) != 2:
            raise ValueError(f"Invalid depth range: {dpths}. Should be a list of 2 integers.")
        nup, ndown = dpths
        nup = abs(nup)
        logger.info("Get depth range [%d, %d] for initial species: %s", nup, ndown, initsps)
    else:  # Draw all
        nup, ndown = 0, 0

    # Read reaction tree
    reactions, sps_info, nrsps = mech.reactions, mech.sinfo["org"], mech.sinfo["nrsp"]
    up_info, down_info = species_reaction_relations(reactions, sps_info)

    # Generate graph
    dot = graphviz.Digraph(comment="Chemical Reactions.", engine="dot")
    dot.attr(rankdir="TB", nodesep=f"{nsep}", ranksep=f"{rsep}")  # LR #TB

    # Add nodes
    nodes = set(initsps)  # Initial nodes, no lumped species
    nodes_used, paths_used, rcns_used, lumped = set(), set(), set(), set()  # Record used nodes and paths
    ugens, dgens = {s: [s] for s in initsps}, {s: [s] for s in initsps}  # Record paths from/to source species
    sinks = {}  # Record terminal species that sinks {species: set(reaction_type)}

    logger.info("Start to generate flowchart for %d initial species ...", len(initsps))
    while nodes:
        if len(nodes_used) >= nmax:
            logger.info("Reach the maximum number of nodes: %d (adjustable). Stop.", nmax)
            break
        sname = _sname(nodes.pop())  # Real species name
        if sname not in nodes_used:
            isp = sps_info[sname]  # Species obj

            # Check if lumped species: if so, add to lumped set
            nd_name = _node_id(sname)  # Node species name
            if nd_name.startswith(lid):
                lumped.add(sname)

            # Node font color
            if sname in initsps:
                ifcol = "red"
            elif sname in nrsps:
                ifcol = "blue"
            else:
                ifcol = "black"

            if wfig:  # With picture
                # Get figure path w/ name
                sfig = fetch_chemical_pic(isp.smiles, os.path.join(wfdir, f"{sname}.png"))
                iscol, ipen = get_sps_layout(isp, RCOLS)  # Border color
                ilabel = create_label_with_image(sfig, nd_name, ifcol, ipen)  # Figure
                dot.node(nd_name, label=ilabel, color=iscol, shape="plaintext")
            else:
                ishape, ipen = get_sps_layout(isp, RSHPS)  # Node shape
                dot.node(nd_name, shape=ishape, fontcolor=ifcol, penwidth=ipen)

            # Record used species
            nodes_used.add(sname)

        nd_name = _node_id(sname)  # Node id

        # Check relations: down
        if sname in dgens:
            for ircn in down_info.get(sname, []):
                if ircn in rcns_used:
                    continue
                rcn = reactions[ircn]
                iarrow, icol = get_rcn_layout(rcn)
                # Edges
                prods = [s for s in rcn.products if s in sps_info]
                is_draw = False
                for s in prods:
                    if (sname, s, iarrow) in paths_used:
                        continue
                    if ndown and s not in dgens and len(dgens[sname]) >= ndown:  # Reach the depth limit
                        continue
                    if s not in nodes_used:  # Add new node to process
                        nodes.add(s)
                    # Record new path
                    dgens[s] = dgens[sname] + [s]
                    # Add new edge
                    dot.edge(nd_name, _node_id(s), arrowhead=iarrow, color=icol, fontcolor=icol)
                    paths_used.add((sname, s, iarrow))
                    is_draw = True

                if not prods or not is_draw:  # Record sinking
                    sink_id = f"{nd_name}_sink"
                    if sname not in sinks:  # New sink node
                        sinks[sname] = set()
                        # Sink dots shape/color
                        dot.node(sink_id, label="", shape="point", width="0.05", height="0.05", color="gray")
                    if iarrow not in sinks[sname]:  # New sink edge
                        sink_id = f"{nd_name}_sink"
                        dot.edge(nd_name, sink_id, arrowhead=iarrow, color=icol, fontcolor=icol)
                        sinks[sname].add(iarrow)

                rcns_used.add(ircn)  # Record used reaction

        # Check relations: up
        if sname in ugens:
            for ircn in up_info.get(sname, []):
                if ircn in rcns_used:
                    continue
                rcn = reactions[ircn]
                iarrow, icol = get_rcn_layout(rcn)
                # Edges
                for s in [s for s in rcn.reactants if s in sps_info]:
                    if (s, sname, iarrow) in paths_used:
                        continue
                    if nup and s not in ugens and len(ugens[sname]) >= nup:  # Reach the depth limit
                        continue
                    if s not in nodes_used:  # Add new node to process
                        nodes.add(s)
                    # Record new path
                    ugens[s] = ugens[sname] + [s]
                    # Add new edge
                    dot.edge(_node_id(s), nd_name, arrowhead=iarrow, color=icol, fontcolor=icol)
                    paths_used.add((s, sname, iarrow))
                rcns_used.add(ircn)  # Record used reaction

    # Save
    filename = f"{svnm}_wfig" if wfig else svnm
    dot.render(filename, format="pdf")  # Save to DOT file

    if len(nodes_used) > nmax:
        logger.warning("Too many nodes in the flowchart: # %d", len(nodes_used))

    logger.info("Saved flowchart with %d nodes to %s.png", len(nodes_used), filename)
    logger.info("Used # %d up paths and # %d down paths.", len(ugens), len(dgens))
    logger.info("Found # %d lumped species and # %d sinks.", len(lumped), len(sinks))

    # Return generated paths for logging
    gens = {f"Up_{k}": v for k, v in ugens.items() if v}
    gens.update({f"Down_{k}": v for k, v in dgens.items() if v})
    return gens
