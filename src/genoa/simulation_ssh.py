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
This module contains functions to run simulations with the SSH-aerosol box model.
"""

import os

from datetime import datetime

from .constants import SNAME
from .ssh_cst import SCOND_ITEMS, SCOND_HOURLY, SMECH_SUFFIX
from .logger import setup_logger
from .setting_global import get_all_settings_map
from .simulation_basic import update_nml_line, get_nline_dict
from .simulation_init import RunSmlSetting
from .utils import add_val_to_str_list, isfloat


# Logger
logger = setup_logger(__name__)


def prepare_ssh_nmls(nml_info: dict, opts: RunSmlSetting) -> list:
    """Prepare SSH-aerosol namelist files and output for writing."""

    # Get settings
    settings_map = get_all_settings_map()
    gnl = settings_map[SNAME.SETUP]()
    ssh_opt = settings_map[SNAME.SSH]()

    # Update default namelist with new settings
    nml_info = update_nml_line(nml_info, _update_ssh_nml_line(opts))

    # Read condition files
    cond_sets = []
    path_cond, conds = opts.conds[0], opts.conds[1:]
    for cond in conds:
        cond_sets.append(_read_from_cond(os.path.join(path_cond, cond)))

    # w/ chem id: updated with id when running simulations
    if opts.wid:
        spath = opts.paths.res
        ipath = opts.paths.nml
        rpath = opts.ref_files_str
        cpath = opts.mech_path
    else:
        cpath = f"{opts.mech_path}/{opts.mech_names[0]}/{opts.mech_names[0]}"
        spath, ipath, rpath = None, None, None

    # If w/ filename
    rfile = "concs.txt" if opts.chemids[0] == "" else ""

    # Get info for writting new namelist
    pool_inputs, items = [], {}  # inputs for pool, items to update
    opts.nmlids, opts.resids, opts.labels = [], [], []  # Reset

    # Check items read from condition files
    cond_items_use = set(SCOND_ITEMS.keys())
    if not ssh_opt.tag_cst_profile:
        cond_items_use -= {"cst_gas", "cst_aero"}

    # Add nout_soa if needed
    if ssh_opt.soa_grps:
        items["nout_soa"] = len(ssh_opt.soa_grps)

    for t, tset in enumerate(gnl.time_settings):

        # Update condition settings
        for cond, cond_set in zip(conds, cond_sets):

            ilabel = f"{cond}.T{t}"  # Label as {conditon}.T{t}
            # Read condition file
            for k, v in cond_set.items():
                if k not in cond_items_use:
                    raise ValueError(f"Invalid info {k}: {v} read from condition file: {cond} in {path_cond}")
                items[SCOND_ITEMS[k]] = v

            # Update hourly initial files
            if ssh_opt.tag_init_hourly:
                for k in SCOND_HOURLY:
                    if k in items:  # Update only if exists - tset[3] is hour
                        items[k] = _add_hour_to_file(items[k], tset[3])

            # Update time settings (if empty: use default namelist settings)
            if tset:
                for k, v in enumerate(["initial_time", "final_time", "delta_t"]):
                    items[v] = tset[k]
                if ssh_opt.tag_time_monthly:  # Add month to time, tset[4] is year
                    k = _get_month_in_seconds(cond, tset[4])
                    for s in ["initial_time", "final_time"]:
                        items[s] = items[s] + k

            # Update id related settings
            if opts.wid:
                inml = f"nml.{len(opts.labels)+1}"  # name of namelist file

                # chem paths w/o suffix
                for k in SMECH_SUFFIX:
                    items[k] = cpath

                # output paths
                items["output_directory"] = spath
                # items["particles_composition_file"] = spath

            else:  # w/o id: real paths
                inml = "nml.ssh"

                # chem file paths w/ suffix
                for k, v in SMECH_SUFFIX.items():
                    items[k] = f"{cpath}{v}"

                # namelist path
                ipath = f"{opts.paths.nml}/{ilabel}"
                if not os.path.exists(ipath):
                    os.makedirs(ipath)  # Create folder

                # result path
                spath = f"{opts.paths.res}/{ilabel}"
                items["output_directory"] = spath
                # items["particles_composition_file"] = f"{spath}/fraction.txt"

                # ref_files_str path
                if opts.ref_files_str:
                    rpath = add_val_to_str_list(opts.ref_files_str, f"/{ilabel}/{rfile}")
                else:
                    rpath = None

            if rpath:  # reference concs.
                items["ref_conc_files_in"] = rpath

            # Write namelist inputa
            snml = os.path.join(ipath, inml)
            new_lines = {nml_info[k]: v for k, v in get_nline_dict(items).items()}
            pool_inputs.append((snml, nml_info["info"], new_lines))

            # Save labels
            opts.labels.append(ilabel)
            opts.nmlids.append(snml)
            opts.resids.append(str(len(opts.labels)) if opts.wid else opts.chemids[0])

    opts.nnml = len(opts.labels)

    return pool_inputs


def get_ssh_sml_inputs(opts: RunSmlSetting) -> list:
    """Get simulation inputs for running SSH-aerosol v2 simulations."""

    # Get inputs for each simulation
    inputs_chem = []
    nerrs = [opts.nerr, opts.nsps] if opts.nerr else None  # Error indicators
    imode = 2 if opts.out_mode >= 2 else 1  # Genoa ssh-aerosol mode: fast 1 or complete 2

    for chemid in opts.chemids:
        inputs = []
        for i, inml in enumerate(opts.labels):
            snml, resid = opts.nmlids[i], opts.resids[i]
            for initid in opts.initids:
                if opts.wid:
                    totof = f"{opts.paths.rec}/toto.{chemid}.{resid}.{initid}"
                    errf = f"{opts.paths.res}/{chemid}/{initid}.{resid}.err"
                elif initid == "":  # no id
                    totof = f"{opts.paths.rec}/{inml}/toto.{inml}"
                    errf = f"{opts.paths.res}/{inml}/errors.txt"
                else:  # - - initid
                    totof = f"{opts.paths.rec}/{inml}/toto.{inml}.{initid}"
                    errf = f"{opts.paths.res}/{inml}/{initid}.err"

                # Command to run
                irun = f"{opts.box_exec} {snml} {imode} {initid} {chemid} {resid}".rstrip()
                inputs.append((irun, totof, errf, nerrs, chemid, opts.err_checks))

        # Add inputs for each chem
        inputs_chem.append(inputs)

    return inputs_chem


def get_ssh_concs_paths(opts: RunSmlSetting) -> dict:
    """Get paths to concentration files for SSH-aerosol."""

    # Get filenames
    fpaths = []
    if opts.wid:
        for chemid in opts.chemids:
            for resid in opts.resids:
                for initid in opts.initids:
                    fpaths.append(f"{opts.soa_path}/{chemid}/{initid}.{resid}.concs")
    else:
        for inml in opts.labels:
            if opts.init_sets:
                for initid in opts.initids:
                    fpaths.append(f"{opts.soa_path}/{inml}/{initid}.concs")
            else:
                fpaths.append(f"{opts.soa_path}/{inml}/concs.txt")

    # Check exist
    for f in fpaths:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Cannot find file: {f}. Check.")

    return fpaths


def _update_ssh_nml_line(opts: RunSmlSetting) -> dict:
    """return a dict containing lines need to be updated in all namelist."""

    # Get settings
    settings_map = get_all_settings_map()
    sml_opt = settings_map[SNAME.SML]()
    env_opt = settings_map[SNAME.ENV]()
    ssh_opt = settings_map[SNAME.SSH]()

    # Initialize dict
    items = {}

    # Initial set
    if sml_opt.init_set_str:
        items["init_species_file"] = sml_opt.init_set_str

    # Error species
    if sml_opt.error_species_str:
        items["err_species_list"] = sml_opt.error_species_str

    # Photolysis file
    if env_opt.photolysis_file:
        items["photolysis_file"] = env_opt.photolysis_file

    # Output options
    items["output_type"] = opts.out_mode  # output mode (3: netcdf, 2: binary 1: text)
    if opts.out_mode == 0:  # With optinoal outputs
        if ssh_opt.out_gas_str:
            items["output_gas_list"] = ssh_opt.out_gas_str
        if ssh_opt.out_aero_str:
            items["output_aero_list"] = ssh_opt.out_aero_str

    # Change formats for new lines
    return get_nline_dict(items)


def _read_from_cond(cond_file: str) -> dict:
    """Return a dict of info read from condition file."""

    # Read content
    with open(cond_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # ... to a dict
    items = {}
    for line in lines:
        line = line.strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        # Check value type
        val = val.strip(" '\"")  # Default: string
        if isfloat(val):  # Number
            val = float(val)
        items[key.strip()] = val

    return items


def _add_hour_to_file(fpath: str, ihour: int) -> str:
    """Add hour to file name."""

    if not fpath:
        return fpath

    # Get file name and extension
    # e,g,, "./init_gas_0h/dat" or "./init_gas.dat" -> "./init_gas_12h.dat"
    fname, fext = os.path.splitext(fpath)
    if fname.endswith("h"):  # Remove hour
        return f"{fname.rsplit('_', 1)[0]}_{ihour}h{fext}"
    return f"{fname}_{ihour}h{fext}"  # Add hour


def _get_month_in_seconds(cond_id: str, year: int) -> int:
    """Extract month info from condition id in the format of m[num]y[num]x[num]"""

    # Settings
    mstr = "m"  # Month indicator

    # Check "m" in the condition id
    i = cond_id.find(mstr) + 1
    if i == 0:
        logger.warning("Invalid condition id: %s - no 'm' found", cond_id)
        return 0

    # Extract month index
    month_str = ""
    while i < len(cond_id) and cond_id[i].isdigit():
        month_str += cond_id[i]
        i += 1
    if not month_str:
        logger.warning("Invalid condition id: %s - no month index found", cond_id)
        return 0
    month = int(month_str)

    # Check month index
    if month == 0:  # January
        return 0
    if month < 0 or month > 11:
        logger.warning("Invalid month index: %s in %s. Should be between 0 and 11", cond_id, month_str)
        return 0
    # Month to seconds
    dt = datetime(year, month + 1, 1) - datetime(year, 1, 1)
    return int(dt.total_seconds())  # Convert to seconds
