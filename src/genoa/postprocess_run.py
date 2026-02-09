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

import os
import time
from pprint import pformat

import numpy as np
from attrs import asdict

from .constants import SNAME
from .folder_path import update_output_path, get_common_files
from .logger import setup_logger
from .mechanism_pack import load_mechanisms
from .record import build_log
from .setting_init import PostProcessOption
from .setting_global import get_all_settings_map
from .setting_update import setup_4action
from .simulation_refconc import get_refconcs_from_files
from .plotting_flowchart import plot_flowchart
from .plotting_run import plot_per_mech, plot_among_mechs, plot_among_conditions
from .postprocess_conc import PostConc, get_time_series
from .postprocess_rec import plot_rdc_records
from .utils import load_common_conditions


# Logger
logger = setup_logger(__name__)


def run_postprocess() -> None:
    """Post-processing: plot results based on the options"""

    # Get settings
    tbeg_pst = time.perf_counter()
    settings_map = get_all_settings_map()
    gnl = settings_map[SNAME.SETUP]()
    opt = settings_map[SNAME.PST]()

    # Update setttings
    setup_4action(SNAME.PST, opt, gnl)

    if opt.plots.get("no_plot", None):
        logger.info("No plot requested. Stop post-processing.")
        return

    # Options for reading concs for post-processing: {"gas": [], "aero": []}
    copts = {"time": "all", "unit": opt.unit, "winit": True, "use_fake": False, "phase": "both", "total": None}

    # Open log file
    opt.logt = build_log(
        os.path.join(opt.savpath, "postprocess.log"),
        f"Log file for post-processing w/ pid: {os.getpid()}\nOptions:\n{pformat(asdict(opt))}",
    )
    opt.logt.open()

    # Plot reduction records
    plot_rdc_records(opt)

    # Mechanisms
    update_kinetic = bool(opt.plot_reactivity)
    rm_pvocs = gnl.primary_vocs if opt.tag_wopvoc else None
    opt.plots["pvocs"] = gnl.primary_vocs
    mechs_all = load_mechanisms(opt.mech_names, opt.mech_paths, update_kinetic, rm_pvocs)

    # Plot flowchart
    plot_flowchart(mechs_all, opt)

    # Plot simulation results
    if opt.plots.get("no_res_plot", None):
        logger.info("No simulation results plot requested.")
    else:
        # Plot per condition w/ ".T" in the name
        conds = load_common_conditions(opt.mech_respaths)
        for cond in conds:
            # Save path
            savpath = update_output_path(os.path.join(opt.savpath, cond), True)
            logger.info("Now post-processing for condition: %s ==> %s", cond, savpath)

            # Time series
            i = int(cond.rsplit(".T", 1)[1])  # Index of the time setting
            if not gnl.time_settings or i >= len(gnl.time_settings):
                raise ValueError(f"No time settings found for condition: {cond} with index: {i}")
            a_time = get_time_series(gnl.time_settings[i], gnl.boxmodel)

            # Concentrations in desired units
            concs = [get_refconcs_from_files(os.path.join(path, cond), copts) for path in opt.mech_respaths]

            # Process data & plot
            process_and_plot(concs, mechs_all, a_time, savpath, opt)

        # Plot across conditions
        plot_among_conditions(conds, mechs_all, opt)

        # Merge plots for all conditions
        merge_plots(conds, opt)

    logger.info("Post-processing finished. Time used: %.1f s", time.perf_counter() - tbeg_pst)
    opt.logt.close()


def process_and_plot(clist: list, mlist: list, tsteps: np.ndarray, spath: str, opt: PostProcessOption) -> None:
    """Process data and plot results"""

    # Process & plot per mechanism
    pdata = []
    logger.info("Processing and plotting for %d mechanisms w/ options: %s ...", len(mlist), opt.plots)
    for i, (mech, concs_in) in enumerate(zip(mlist, clist)):

        # Save path w/ prefix
        opt.plots["fpre"] = os.path.join(spath, f"{opt.mech_labels[i]}_")

        # Concentrations
        tconcs = PostConc(gas=concs_in["gas"], aero=concs_in["aero"], steps=tsteps, mlb=opt.mech_labels[i])

        # Plot
        pdata.append(plot_per_mech(tconcs, mech, opt))

    # Plot among mechanisms
    opt.plots["fpre"] = os.path.join(spath, "cmp_")
    plot_among_mechs(pdata, tsteps, opt)


def merge_plots(conds: list, opt: PostProcessOption) -> None:
    """Merge plots for all conditions and save to 'merged' folder"""

    if not (opt.tag_merge and conds):  # Skip merging
        return

    fgnames = None  # Common figure names
    for cond in conds:
        fpath = os.path.join(opt.savpath, cond)
        fgnames = get_common_files(fpath, fgnames, ".png")
    if not fgnames:
        logger.info("No figures found for merging.")
        return

    # Save path
    merged_path = update_output_path(os.path.join(opt.savpath, "merged"), True)
    logger.info("Merging # %d figures and saving to %s ...", len(fgnames), merged_path)

    # Add border to figures from certain conditions
    conds_user = ["4remote", "30urban", "60nofix", "61fixhv", "62nohv", "7Wdiurnal"]
    conds_in = [cond for cond in conds if cond.split(".T")[0] in conds_user]
    if conds_in:
        logger.info("Adding border to figures from certain conditions: %s", conds_in)
        for cond in conds_in:
            for fg in fgnames:
                fgpath = os.path.join(opt.savpath, cond, fg)
                if not os.path.isfile(fgpath):
                    continue
                cmd = f"convert {fgpath} -bordercolor red -border 2 {fgpath}"
                os.system(cmd)

    # Merge plots
    cmd_setup = "-tile x2 -geometry +1+1 -gravity North -pointsize 20"
    for fg in fgnames:
        fgwcds = [f'-label "{cond}" {os.path.join(opt.savpath, cond, fg)}' for cond in conds]
        cmd = f"montage {' '.join(fgwcds)} {cmd_setup} {os.path.join(merged_path, fg)}"
        os.system(cmd)
