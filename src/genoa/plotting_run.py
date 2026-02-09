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
This module provides functions for plotting results.
"""


import os
from typing import Optional

import numpy as np

from .constants import PSAT_BIN_LIMS, PHASES, RATIOS, NTOP_MAX, EPS
from .error_analysis import cmp_error
from .folder_path import update_output_path
from .logger import setup_logger
from .mechanism_pack import Mechanism
from .setting_init import PostProcessOption
from .plotting_basic import FigureOption, plot_dict_to_1fig, plot_ekma, plot_errors  # , plot_gp_distribution
from .plotting_flowchart import fetch_chemical_pic
from .postprocess_conc import PostConc, get_sorted_splist, calc_ekma_ratio, calc_reactivity, calc_ro2_rct, calc_fgroup
from .unit_conversion import get_unit_conv_factors
from .utils import size_difference


# Logger
logger = setup_logger(__name__)


def plot_per_mech(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """Plot for a mechanism and return the plot data"""

    plot_dict = {
        "time_series": plot_time_series,
        "reactivity": plot_reactivity,
        "carbon": plot_carbon,
        "functional_group": plot_functional_group,
        "volatility": plot_psat,
        "ratios": plot_ratios,
        "hc_nox_ratio": plot_ekma_diagram,
        "gp_distribution": None,  # plot_gp_distribution,
    }

    out = {}  # Items for plotting, output data
    for k, v in opt.plots.items():
        if k in ["mech_cmp", "mech_cmp2", "pvocs", "fpre"]:
            continue
        if k not in plot_dict:
            raise ValueError(f"Plot type {k} not found in the plot dictionary")
        if not v:  # Not processing here
            continue
        pdata = plot_dict[k](tconc, mech, opt)
        # Update plotable data
        for k, v in pdata.items():
            if k in out:
                logger.warning("Duplicate key found in the output data: %s, will be overwritten", k)
            out[k] = v
    return out


def plot_among_mechs(data: list, tsteps: np.ndarray, opt: PostProcessOption) -> None:
    """Plot differences among mechanisms"""

    plot_diff_mechs_1d(data, tsteps, opt)  # Time series
    plot_diff_mechs_2d(data, opt)  # 2-mechanism comparison


def plot_among_conditions(conds: list, mechs: list, opt: PostProcessOption) -> None:
    """Plot data among conditions"""

    if not opt.cdata:  # No data to compare among conditions
        return

    # Save path
    save_path = update_output_path(os.path.join(opt.savpath, "errors"))

    # Get size changes
    mech_ids = {m: i for i, m in enumerate(opt.mech_labels)}  # Index
    mref = opt.err_ref_mech
    sizes = {m: [sum(mechs[i].get_size())] for m, i in mech_ids.items()}
    size_diffs = {m: size_difference(sizes[mref], sizes[m], "rel")[0] for m in opt.plots["mech_cmp"] if m != mref}

    # Plot per err type
    logger.info("Plotting error plots among conditions ...")
    for ierr, etype in enumerate(opt.err_types):
        for item in opt.cdata[0]:
            pdata = {s: {} for s in conds}
            for i, s in enumerate(conds):
                for m in size_diffs:
                    pdata[s][m] = opt.cdata[i][item][m][ierr]
            ititle = f"{item} w/ etype {etype} & ref {mref}"
            ipath = os.path.join(save_path, f"{item}_{etype}_{mref}.png".replace("/", "-"))
            plot_errors(size_diffs, pdata, ititle, ipath)

    # Save data
    if opt.tag_csv:

        # Save size data
        csv_file = os.path.join(save_path, "size.csv")
        with open(csv_file, "w", encoding="utf-8") as f:
            f.write("Mech,No.reactions,No.species,No.condensables,Total,Sizediff\n")
            for m, i in mech_ids.items():
                f.write(f"{m},{','.join([str(s) for s in mechs[i].get_size()])},{sizes[m][0]},{size_diffs.get(m, 0)}\n")

        # Save error data
        csv_file = os.path.join(save_path, f"error_{mref}.csv")
        with open(csv_file, "w", encoding="utf-8") as f:
            f.write(f"Item,Mech,Cond,{','.join(opt.err_types)}\n")
            for item in opt.cdata[0]:
                for i, s in enumerate(conds):
                    for m in size_diffs:
                        f.write(f"{item},{m},{s},{','.join([str(e) for e in opt.cdata[i][item][m]])}\n")
        logger.info("Error data among conditions saved to %s.", csv_file)

    logger.info("Plotted %d error plots among conditions.", len(opt.err_types) * len(opt.cdata[0]))


def plot_time_series(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """
    Plot time series of concentrations.
    Example inputs (sgroup_phase/s1,s2,s3+s4_phase): ["org_gas", "pvoc", "NO+NO2_gas"]
    """

    plot_data = {}  # Save later for comparison
    for strin in opt.plots["time_series"]:
        strin = strin.strip()
        if not strin:
            logger.warning("No item found in the input string for plotting time series: %s", strin)
            continue
        item, phase = rsplit_phase(strin, "total")

        # Try stacked plot
        pdata = plot_stack_wrank(item, phase, tconc, mech, opt)
        # Try line plot
        if not pdata:
            pdata = plot_line_series(item, phase, tconc, opt)
        # Save data
        plot_data.update(pdata)

    return plot_data


def plot_reactivity(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """
    Plot reactivity variations.
    Example inputs (s1,s2_mode): ["HO,O3,NO3_voc", "NO,HO2_ro2", "HO,O3_all"]
    """

    plot_data = {}
    for strin in opt.plots["reactivity"]:
        strin = strin.strip()
        if not strin:
            logger.warning("No item found in the input string for plotting reactivity: %s", strin)
            continue
        rct_list, flag = get_splist(strin, sep2=None)
        if strin == "NO,HO2_ro2":
            logger.info("Calculating RO2 rate w/ %s ...", strin)
            pdata = calc_ro2_rct(flag, rct_list, tconc, mech, opt.unit)
        else:
            pdata = calc_reactivity(flag, rct_list, tconc, mech, opt.unit)
        if not pdata:
            logger.warning("No reactivity data found for %s. Skip.", strin)
            continue
        # Plot
        sfx, svnm = f"rct_{strin}", opt.plots["fpre"]
        fopt = FigureOption(sfx=sfx, svnm=svnm, x=tconc.steps, ylists=pdata[1], ykeys=pdata[0], style="line")
        fopt.pltset.update({"title": f"Reactivity of {strin}: {tconc.mlb}", "unit": opt.unit})
        plot_dict_to_1fig(fopt, opt.tag_csv)

        if len(pdata[0]) > 1:  # Plot w/ normalization if multiple species
            fopt.update_to_norm("sum")
            plot_dict_to_1fig(fopt, opt.tag_csv)
        plot_data.update({f"rct_{strin}_{k}": [v] for k, v in zip(pdata[0], pdata[1])})

    return plot_data


def plot_ratios(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """
    Plot ratios of OM/OC, O/C, H/C, N/C.
    Ecample inputs (sgroup_phase): ["org_total", "rsp_total", "rsp", "radical_gas"]
    """
    plot_data = {}
    for strin in opt.plots["ratios"]:
        strin = strin.strip()
        if not strin:
            logger.warning("No item found in the input string for plotting ratios: %s", strin)
            continue
        sps, phase = rsplit_phase(strin, "total")
        if sps not in mech.sinfo:
            raise ValueError(f"Species group not found for ratio plotting: {sps}")
        sps_dict = mech.sinfo["org"]
        if sps != "org":
            sps_dict = {k: sps_dict[k] for k in mech.sinfo[sps] if k in sps_dict}
        # Get ratios
        rts = np.zeros((len(RATIOS) + 1, len(tconc.steps)))
        concs = tconc.get_phase_concs(phase)
        for k, v in sps_dict.items():
            iconcs = concs[v.name]
            for i, s in enumerate(RATIOS):
                rts[i] += iconcs * v.ratios.get(s, 0)
            rts[-1] += iconcs  # Total mass
        rts = {s.replace("/", ":"): np.where(rts[-1] > EPS, rts[i] / rts[-1], 0) for i, s in enumerate(RATIOS)}
        # Plot
        logger.info("Plotting ratios for %s calculated from # %d species ...", sps, len(sps_dict))
        sopt = {"sfx": f"ratios_{strin}", "title": f"Ratios of {sps}: {tconc.mlb}-{phase}"}
        pdata = plot_data_in_one([tconc.steps, rts], sopt, "line", opt)
        plot_data.update({f"rt_{strin}_{k}": [v] for k, v in zip(pdata[0], pdata[1])})

    return plot_data


def plot_carbon(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """Plot carbon distribution. Example inputs (phase): ["gas", "total", "aero"]"""
    # Get species group list
    pvocs, rsp = set(opt.plots["pvocs"]), set(mech.sinfo["rsp"])
    oorgs = set(mech.sinfo["org"]) - pvocs - rsp
    sdict = {
        "Precur.": pvocs,
        "Formed Org.": rsp,
        "Other Org.": oorgs,
        "Lost_C": ["XCLOST"],
        "CO": ["CO"],
        "CO2": ["CO2"],
    }
    # "CO+CO2": ["CO", "CO2"]}

    plot_data = {}  # For output
    unit_conv = get_unit_conv_factors(mech.species, unit_in=opt.unit, unit_out="ppbC")
    for phase in opt.plots["carbon"]:
        phase = phase.strip()
        if not phase or phase not in PHASES:
            logger.warning("Invalid phase for carbon plotting: %s. Skip.", phase)
            continue
        logger.info("Plotting carbon distribution for %s ...", phase)
        concs = tconc.get_sum_concs(phase, sdict, unit_conv)
        sopt = {"sfx": f"carbon_{phase}", "title": f"Carbon Distribution: {tconc.mlb}-{phase}", "unit": "ppbC"}
        plot_data[f"carbon_{phase}"] = plot_data_in_one([tconc.steps, concs], sopt, "stack", opt)
    return plot_data


def plot_psat(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """Plot volatility distribution. Example inputs (phase): ["gas", "total", "aero"]"""

    # Gett svocs
    sdict = {k: set() for k in PSAT_BIN_LIMS}
    svocs = {s: sp.psat_atm for s, sp in mech.sinfo["org"].items() if s in mech.sinfo["svoc"]}
    for s, ipsat in svocs.items():
        if ipsat <= 0:
            raise ValueError(f"Condensable species w/o psat: {s}")
        flag = False
        for k, v in PSAT_BIN_LIMS.items():
            if ipsat > v:
                sdict[k].add(s)  # Add to the bin
                flag = True
                break
        if not flag:
            raise ValueError(f"No bin found for species: {s} w/ psat {ipsat}. Range: {PSAT_BIN_LIMS}")

    sdict.update({"VOC": mech.sinfo["voc"], "Radical": mech.sinfo["radical"]})

    plot_data = {}  # For output
    for phase in opt.plots["volatility"]:
        phase = phase.strip()
        if not phase or phase not in PHASES:
            logger.warning("Invalid phase for psat plotting: %s. Skip.", phase)
            continue
        logger.info("Plotting psat distribution for %s ...", phase)
        concs = tconc.get_sum_concs(phase, sdict)
        sopt = {"sfx": f"psat_{phase}", "title": f"Volatility Distribution: {tconc.mlb}-{phase}"}
        plot_data[f"psat_{phase}"] = plot_data_in_one([tconc.steps, concs], sopt, "stack", opt)

    return plot_data


def plot_functional_group(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """
    Plot functional group distribution.
    Example inputs (phase_mode): ["gas_GECKO", "total_GECKO", "aero_GECKO"]
    """

    plot_data = {}
    for strin in opt.plots["functional_group"]:
        strin = strin.strip()
        if not strin:
            logger.warning("No item found in the input string: %s", strin)
            continue
        phase, mode = rsplit_phase(strin, "GECKO")
        if phase not in PHASES:
            raise ValueError(f"Invalid phase for functional group plotting: {phase}")

        logger.info("Plotting functional group distribution for %s ...", strin)
        # Get data
        concs = calc_fgroup(tconc.get_phase_concs(phase), mech.species, len(tconc.steps), mode)
        sopt = {"sfx": f"fgroup_{strin}", "title": f"Functional Group Distribution: {tconc.mlb}-{phase}"}
        plot_data[f"fgroup_{strin}"] = plot_data_in_one([tconc.steps, concs], sopt, "stack", opt)

    return plot_data


def plot_ekma_diagram(tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """
    Plot an EKMA (Empirical Kineitc Modeling Approach) diagram showing NOx vs. HC concentrations.
    Example inputs (phase_nskip): ["gas_h", total_h", "total_1"]
    """
    plot_data = {}  # Save for output
    for strin in opt.plots["hc_nox_ratio"]:
        strin = strin.strip()
        if not strin:
            logger.warning("No item found in the input string for plotting EKMA: %s", strin)
            continue
        item, nskip = rsplit_phase(strin, "h")
        if item not in ["gas", "total"]:
            raise ValueError(f"Unsupported phase for EKMA plotting: {item}")
        # Get data
        pdata = calc_ekma_ratio(tconc.get_phase_concs(item), mech.species, opt.unit)
        logger.info("Plotting EKMA diagram for %s ...", item)
        # Plot
        plot_ekma([[tconc.mlb, pdata]], tconc.steps, nskip, f"{opt.plots['fpre']}ekma_{item}", tag_csv=opt.tag_csv)
        # Save data
        plot_data.update({f"{k}_{item}": [v] for k, v in pdata.items()})
        plot_data[f"ekma_{item}_{nskip}"] = pdata
    return plot_data


def plot_diff_mechs_1d(plot_data: list, tsteps: np.ndarray, opt: PostProcessOption) -> None:
    """Plot time series among mechanisms. Inputs are a list of mechanism names."""

    if not opt.plots.get("mech_cmp", False):
        return

    mechs, svnm = opt.plots["mech_cmp"], opt.plots["fpre"]
    mech_ids = [opt.mech_labels.index(m) for m in mechs]
    logger.info("Plotting differences among multiple mechanisms: %s ...", mechs)

    ndata = _reorganize_data_1d(plot_data, mech_ids)
    if not ndata:
        raise ValueError(f"No data found for the mechanisms: {mechs}")
    for item, ylist in ndata.items():
        logger.info("Plotting cmp1 for %s ...", item)
        fopt = FigureOption(sfx=item, svnm=svnm, x=tsteps, ylists=ylist, ykeys=mechs, style="line")
        fopt.pltset.update({"title": item, "unit": opt.unit})
        plot_dict_to_1fig(fopt, opt.tag_csv)

    # Error calculation
    compare_errors_among_mechs(ndata, opt)

    # Plot ekma diagram for diff mechs
    if "hc_nox_ratio" in opt.plots:
        items = [s for s in plot_data[mech_ids[0]] if s.startswith("ekma_")]
        if not items:
            logger.info("No EKMA data found for the mechanisms. Skip.")
            return
        logger.info("Plotting EKMA differences for %s ...", items)
        for item in items:
            ndata = [[opt.mech_labels[i], plot_data[i][item]] for i in mech_ids]
            nitem, nskip = rsplit_phase(item, None)
            plot_ekma(ndata, tsteps, nskip, f"{svnm}{nitem}", tag_csv=opt.tag_csv)


def plot_diff_mechs_2d(plot_data: list, opt: PostProcessOption) -> None:
    """Plot 2D differences among mechanisms. Inputs are a list of mechanism names in pairs."""
    if not opt.plots.get("mech_cmp2", None):
        return

    mech_list, svnm = opt.plots["mech_cmp2"], opt.plots["fpre"]
    wdot = [[0, "Start"], [-1, "End"]]  # Add dots
    logger.info("Plotting differences between two mechanisms: %s ...", mech_list)
    for mechs in mech_list:
        if len(mechs) != 2:
            raise ValueError(f"Invalid number of mechanisms for 2D plot: {mechs}")
        ylist_all, items, mlbs = [], None, []
        for m in mechs:
            i = opt.mech_labels.index(m)
            pdata = plot_data[i]
            mlbs.append(opt.mech_labels[i])
            # Get /update items
            if items is None:
                items = [k for k, v in pdata.items() if v and len(v) == 2]
            else:
                items = [k for k in items if k in pdata and pdata[k] and len(pdata[k]) == 2]
            ylist_all.append(pdata)

        # Treat ekma separately
        items = [i for i in items if not i.startswith("ekma_")]

        # Plot per item
        for item in items:
            vals = [ylist_all[0][item], ylist_all[1][item]]
            # Common keys for plotting
            ykeys = [k for k in vals[0][0] if k in vals[1][0]]
            if not ykeys:
                logger.info("No common species found for %s. Skip.", item)
                continue
            # Get data for plotting
            ylists = []
            for i in range(2):
                ylists.append([vals[i][1][vals[i][0].index(k)] for k in ykeys])
            logger.info("Plotting cmp2 differences for %s w/ keys: %s ...", item, ykeys)
            fopt = FigureOption(sfx=item, svnm=svnm, ylists=ylists, ykeys=ykeys, style="2d", mechs=mlbs, wdot=wdot)
            fopt.pltset.update({"xlabel": f"{mlbs[0]}", "ylabel": f"{mlbs[1]}", "title": item})
            plot_dict_to_1fig(fopt, opt.tag_csv)


def plot_stack_wrank(item: str, phase: str, tconc: PostConc, mech: Mechanism, opt: PostProcessOption) -> dict:
    """Plot stacked series for a given item"""

    # Get species group
    if item in mech.sinfo:  # org, inorg, rsp, nrsp, voc, svoc, radical, svoc, ro2
        sgroups = mech.sinfo[item]
    elif item == "pvoc":  # pvoc
        sgroups = opt.plots["pvocs"]
    else:  # No stacked plot
        return {}

    logger.info("Plotting stacked plot for %s w/ # %d species...", item, len(sgroups))
    svnm = opt.plots["fpre"]
    slist_lim = {"sort": "sum", "nlim": NTOP_MAX}

    # Get concs and sorted list
    concs = tconc.get_concs(sgroups, phase)[phase]
    slist = get_sorted_splist(concs, slist_lim, opt.logt)
    ylist = [concs[s] for s in slist]
    ytot = np.sum(list(concs.values()), axis=0)
    strin = f"{item}_{phase}"

    # Plot
    fopt = FigureOption(sfx=strin, svnm=svnm, x=tconc.steps, ylists=ylist, ykeys=slist, style="stack")
    fopt.pltset.update({"title": f"{item} Distribution: {tconc.mlb}-{phase}", "unit": opt.unit})
    if len(slist) < len(concs):
        fopt.wline = [["Total", ytot]]
        fopt.pltset["lgd_title"] = f"Top {len(slist)}"
    if "inorg" not in item:
        plot_species_table(fopt, opt, mech.sinfo["org"])
    plot_dict_to_1fig(fopt, opt.tag_csv)

    if opt.tag_wnorm:
        fopt.update_to_norm("sum")
        plot_dict_to_1fig(fopt, opt.tag_csv)

    return {f"{strin}_tot": [ytot], f"{strin}_top": [slist, ylist]}


def plot_data_in_one(pdata: list, sopt: dict, mode: str, opt: PostProcessOption) -> list:
    """Plot stacked series for a given concs dict"""
    # Get data
    tsteps, concs = pdata
    if not concs:
        logger.warning("No data found for plotting data w/ options: %s", sopt)
        return []
    if mode not in ["stack", "line"]:
        raise ValueError(f"Invalid mode for plotting: {mode}")
    ylists, ykeys = list(concs.values()), list(concs.keys())
    # Plot
    fopt = FigureOption(sfx=sopt["sfx"], svnm=opt.plots["fpre"], x=tsteps, ylists=ylists, ykeys=ykeys, style=mode)
    fopt.pltset.update({"title": sopt["title"], "unit": opt.unit})
    plot_dict_to_1fig(fopt, opt.tag_csv)
    if opt.tag_wnorm:
        fopt.update_to_norm("sum")
        plot_dict_to_1fig(fopt, opt.tag_csv)

    return [ykeys, ylists]


def plot_line_series(item: str, flag: str, tconc: PostConc, opt: PostProcessOption) -> dict:
    """Plot line series for a given item"""

    splist, _ = get_splist(item)
    phase = flag if flag in PHASES else None

    # Get concs
    ylist, ykeys = [], []  # conc list and keys for plotting
    for sps in splist:
        if not sps:
            continue
        if phase is None:  # Read phase from species name
            sdict = get_splist_wphase(sps)
            concs = np.sum([tconc.get_concs(s, p)[p] for p, s in sdict.items()], axis=0)
            phase = ",".join(sdict.keys())
        else:
            concs = tconc.get_sum_concs(phase, {"sps": sps})["sps"]
        ykeys.append("+".join(sps))
        ylist.append(concs)

    if not ylist:
        raise ValueError(f"No species found in the input string: {item}")
    # Plot
    strin = f"{item}_{phase}" if phase else item
    logger.info("Plotting line plot for %s ...", strin)
    fopt = FigureOption(sfx=strin, svnm=opt.plots["fpre"], x=tconc.steps, ylists=ylist, ykeys=ykeys, style="line")
    fopt.pltset.update({"title": f"{item} Distribution: {tconc.mlb}-{phase}", "unit": opt.unit})
    plot_dict_to_1fig(fopt, opt.tag_csv)

    # Outoput
    if len(ykeys) > 1:
        out = {f"{s}_{phase}": [ylist[i]] for i, s in enumerate(ykeys)}
        out.update({f"{strin}_tot": [np.sum(ylist, axis=0)], strin: [ykeys, ylist]})
        return out
    return {strin: ylist}


def plot_species_table(fopt: FigureOption, opt: PostProcessOption, sps_dict: dict) -> None:
    """Plot species figure table"""

    # Species list w/ SMILES
    slist = [(i, s) for i, s in enumerate(fopt.ykeys) if s in sps_dict]
    if not slist:
        logger.warning("No species found for figure table from %s. Skip.", fopt.ykeys)
        return

    # Path for species figures
    spath = update_output_path(os.path.join(opt.savpath, "sfigs"))
    figs = []
    for i, s in slist:
        isp = sps_dict[s]
        if not isp.smiles or isp.name != s:
            logger.warning("No SMILES found for species: %s or name mismatch. Skip.", s)
            continue
        # Generate figure
        sfig = fetch_chemical_pic(isp.smiles, os.path.join(spath, f"{isp.name}.png"))
        iconc = np.mean(fopt.ylists[i])
        sconc = f"{round(iconc, 2)}" if iconc > 0.1 else f"{iconc:.2e}"
        scmd = f'-label "#{i+1} {s} ({sconc})" {sfig}'
        figs.append(scmd)
    if not figs:
        logger.warning("No species figure found for table. Skip.")
        return

    # Save table
    sname = f"{fopt.svnm}{fopt.sfx}_tbl.png"

    # Plot w/ montage
    ncol = 3 if len(figs) > 2 else len(figs)
    cmd_setup = f"-tile {ncol}x -geometry 200x+10+25 -gravity center -pointsize 15"
    cmd = f"montage {' '.join(figs)} {cmd_setup} {sname}"

    # Add title
    title = fopt.pltset.get("lgd_title", "") + " " + fopt.pltset.get("title", "") + f" (mean concs. in {opt.unit})"
    # Saving structure to one table
    cmd += f" && convert {sname} -gravity North -pointsize 20 -annotate +0+0 '{title}' {sname}"
    # logger.info("cmd: %s", cmd)
    os.system(cmd)
    logger.info("Species figure table for %s is saved as %s.", title, sname)


def rsplit_phase(sin: str, default: Optional[str]) -> list:
    """Split species name and phase"""
    if "_" in sin:
        return [s.strip() for s in sin.rsplit("_", 1)]
    return [sin.strip(), default]


def get_splist(spsin: str, sep1: str = ",", sep2: Optional[str] = "+") -> list:
    """Read species list for plotting from input string"""
    splist = []
    spsin, phase = rsplit_phase(spsin, None)
    for sp in spsin.split(sep1):
        sp = sp.strip()
        if not sp:
            continue
        if not sep2:  # No secondary separator
            splist.append(sp)
            continue
        if sep2 in sp:
            splist.append([s.strip() for s in sp.split("+") if s.strip()])
        else:
            splist.append([sp])
    if not splist:
        raise ValueError(f"No species found in the input string: {spsin}")
    return [splist, phase]


def get_splist_wphase(slist: list) -> dict:
    """Read species names and phases from input list."""
    if not slist:
        return {}

    # Separate species w/ phase
    phase_dict = {"G": "gas", "A": "aero", "T": "total"}
    sdict = {s: [] for s in PHASES}
    for s in slist:
        find = False
        for i, p in phase_dict.items():
            if s.startswith(i):
                sdict[p].append(s)
                find = True
                break
        if not find:
            logger.warning("Species w/o phase: %s. Read w/ total. Use G/A/T as prefix to specify phase.", s)
            sdict["total"].append(s)
    return {k: v for k, v in sdict.items() if v}


def _reorganize_data_1d(plot_data: list, mechids: list) -> dict:
    """Reorganize data for plotting among mechanisms"""
    new_pdata = {}
    for mech_id, pdata_in in enumerate(plot_data):
        if mech_id not in mechids:
            continue
        for k, v in pdata_in.items():
            if not v or len(v) != 1:
                continue
            if k not in new_pdata:
                new_pdata[k] = {}
            new_pdata[k][mech_id] = v[0]
    return {k: [v[m] for m in mechids] for k, v in new_pdata.items()}


def compare_errors_among_mechs(data: dict, opt: PostProcessOption) -> None:
    """Compare and record errors among mechanisms"""
    if not (data and opt.err_types and opt.err_ref_mech):
        return

    # Get mechanism indices
    mech_ids = {m: i for i, m in enumerate(opt.mech_labels)}
    mref = opt.err_ref_mech
    mrdcs = [m for m in opt.plots["mech_cmp"] if m != mref and m in mech_ids]
    if not mrdcs:
        logger.warning("No mechanism found for error comparison. Skip.")
        return

    logger.info("Comparing errors among mechanisms %s w/ ref %s. IDs: %s ...", mrdcs, mref, mech_ids)

    # Get all errors {item: {mech: [error1, error2, ...]}}
    err_info, iref = {}, mech_ids[mref]
    for item, ylist in data.items():
        err_info[item] = {}
        for m in mrdcs:
            irdc = mech_ids[m]
            einfo = [cmp_error(ylist[irdc], ylist[iref], s) for s in opt.err_types]
            err_info[item][m] = einfo

    # Save for comparison between conditions
    if opt.cdata is None:
        opt.cdata = []
    opt.cdata.append(err_info)

    # Save to csv
    csv_file = f"{opt.plots['fpre']}error.csv"
    with open(csv_file, "w", encoding="utf-8") as f:
        f.write(f"Item,Mech,Ref,{','.join(opt.err_types)}\n")
        for item, einfos in err_info.items():
            for m, einfo in einfos.items():
                f.write(f"{item},{m},{mref}{','.join([str(e) for e in einfo])}\n")
    logger.info("Errors among mechanisms saved to %s.", csv_file)
