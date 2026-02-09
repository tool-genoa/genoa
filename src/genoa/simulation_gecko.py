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
This module contains functions to run simulations using GECKO-A box model.
"""

import os
import subprocess

from .constants import SNAME
from .logger import setup_logger
from .setting_global import get_attrs_from_settings, get_all_settings_map
from .simulation_cmp import set_boxmodel
from .simulation_basic import update_nml_line, get_nline_dict
from .simulation_init import RunSmlSetting
from .utils import add_val_to_str_list


# Logger
logger = setup_logger(__name__)


def run_interp(outputpath: str, cmd_list_add: list, binfile: str) -> bool:
    """Run INTERP30 to convert .mech file to .bin file. Return True if successful."""

    # Get interp execute
    settings = get_attrs_from_settings({SNAME.SETUP: ["path_workplace", "boxmodel"]})
    path_interp = set_boxmodel(settings[0], settings[1], "gecko2box")

    # Run interp
    # ./gecko2box [id_mech] 1
    cmd_run = [path_interp] + cmd_list_add
    with open(os.path.join(outputpath, "toto_interp"), "w", encoding="utf-8") as f:
        subprocess.run(cmd_run, cwd=outputpath, stdout=f, check=True)

    # Check output
    if os.path.exists(os.path.join(outputpath, binfile)):
        return True

    logger.warning("%s not found after running INTERP30.", binfile)
    return False


def prepare_gck_nmls(nml_info: dict, opts: RunSmlSetting) -> list:
    """Prepare namelists and update operation ids (labels, resids, nmlids) in the settings."""

    # Update default namelist with new settings for all simulations
    nml_info = update_nml_line(nml_info, _update_gck_nml_line(opts))

    # w/ chem id: updated with id when running simulations
    if opts.wid:
        spath = opts.paths.res  # result path
        ipath = opts.paths.nml  # namelist path
        rpath = opts.ref_files_str  # ref_files_str path
        cpath = opts.mech_path  # chem path
    else:
        cpath = f"{opts.mech_path}/{opts.mech_names[0]}/{opts.mech_names[0]}"
        spath, ipath, rpath = None, None, None  # Assign later

    if opts.chemids[0] == "":  # If w/ filename
        rfile = "outdat.genoa"
    else:
        rfile = ""

    # Get info for writting new namelist
    pool_inputs, items = [], {}  # inputs for pool, items to update
    opts.nmlids, opts.resids, opts.labels = [], [], []  # Reset

    time_settings = get_attrs_from_settings({SNAME.SETUP: "time_settings"})
    for t, tset in enumerate(time_settings):

        # Update time settings
        if tset:  # if None: use default namelist settings
            for i, s in enumerate(["tstart", "tstop", "ntstep", "iday", "imonth", "iyear"]):
                items[s] = tset[i]

        # Update condition settings
        for cond in opts.conds[1:]:
            # Update condition related settings
            ilabel = f"{cond}.T{t}"  # Label as {conditon}.T{t}
            items["input_cond"] = os.path.join(opts.conds[0], cond)  # Condition

            # Update id related settings
            if opts.wid:
                inml = f"nml.{len(opts.labels)+1}"  # name of namelist file
            else:  # w/o id: real paths
                inml = "nml.gck"

                # namelist path
                ipath = f"{opts.paths.nml}/{ilabel}"
                if not os.path.exists(ipath):
                    os.makedirs(ipath)  # Create folder

                spath = f"{opts.paths.res}/{ilabel}"  # result path
                # ref_files_str path
                if opts.ref_files_str:
                    rpath = add_val_to_str_list(opts.ref_files_str, f"/{ilabel}/{rfile}")
                else:
                    rpath = None
            items["input_dir_chem"] = cpath  # chem path
            items["output_dir"] = spath  # save path
            if rpath:  # reference concs.
                items["ref_conc_list_in"] = rpath

            # Write namelist file
            snml = os.path.join(ipath, inml)
            new_lines = {nml_info[k]: v for k, v in get_nline_dict(items).items()}  # index: new_line
            pool_inputs.append((snml, nml_info["info"], new_lines))

            # Save labels
            opts.labels.append(ilabel)
            opts.nmlids.append(snml)
            opts.resids.append(str(len(opts.labels)) if opts.wid else opts.chemids[0])

    opts.nnml = len(opts.labels)

    return pool_inputs


def get_gck_sml_inputs(opts: RunSmlSetting) -> list:
    """Get simulation inputs for running GECKO-A simulations."""

    # Get inputs for each simulation
    inputs_chem = []
    nerrs = [opts.nerr, opts.nsps] if opts.nerr else None  # Error indicators

    for chemid in opts.chemids:
        inputs = []
        for i, inml in enumerate(opts.labels):
            snml, resid = opts.nmlids[i], opts.resids[i]
            for initid in opts.initids:
                if opts.wid:
                    totof = f"{opts.paths.rec}/toto.{chemid}.{resid}.{initid}"
                    errf = f"{opts.paths.res}/{chemid}/{resid}.{initid}.err"
                elif initid == "":  # no id
                    totof = f"{opts.paths.rec}/{inml}/toto.{inml}"
                    errf = f"{opts.paths.res}/{inml}/outdat.err"
                else:  # - - initid
                    totof = f"{opts.paths.rec}/{inml}/toto.{inml}.{initid}"
                    errf = f"{opts.paths.res}/{inml}/outdat.{initid}.err"

                # Command to run
                irun = f"{opts.box_exec} {snml} {chemid} {resid} {initid}".rstrip()
                inputs.append((irun, totof, errf, nerrs, chemid, opts.err_checks))

        # Add inputs for each chem
        inputs_chem.append(inputs)

    return inputs_chem


def get_gck_concs_paths(opts: RunSmlSetting) -> dict:
    """Get paths to concentration files in the result forlder."""

    # Get file suffix
    if opts.init_sets:
        suffix = [f"{initid}.genoa" for initid in opts.initids]
    else:
        suffix = ["genoa"]

    # Get filenames
    fpaths = []
    if opts.wid:
        for chemid in opts.chemids:
            for resid in opts.resids:
                for s in suffix:
                    fpaths.append(f"{opts.soa_path}/{chemid}/{resid}.{s}")
    else:
        for inml in opts.labels:
            for s in suffix:
                fpaths.append(f"{opts.soa_path}/{inml}/outdat.{s}")
    # Check exist
    for f in fpaths:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Cannot find file: {f}. Check.")

    return fpaths


def _update_gck_nml_line(opts: RunSmlSetting) -> dict:
    """Return a dict containing new lines to be updated in the namelist."""

    # Get settings
    settings_map = get_all_settings_map()
    sml_opt = settings_map[SNAME.SML]()
    env_opt = settings_map[SNAME.ENV]()

    # Initialize dict of {key: new line}
    items = {}

    # Output mode
    items["fg_output"] = opts.out_mode

    # items["fg_interp"] = 0  # product ratio not used for now

    # Initial set
    if sml_opt.init_set_str:
        items["init_conc_list_in"] = sml_opt.init_set_str

    # Error species
    err_sps_str = sml_opt.error_species_str
    if err_sps_str:
        items["err_sps_list_in"] = err_sps_str
        # Get number of error species + 1 for total SOA
        items["fg_out_col"] = _get_num_of_err_sps(err_sps_str) + 1

    # Photolysis file
    if env_opt.photolysis_file:
        items["phot_file"] = env_opt.photolysis_file

    # Change formats for new lines
    return get_nline_dict(items)


def _get_num_of_err_sps(err_sps_str: str) -> int:
    """return the number of error species in the error species string"""

    if ";" in err_sps_str:
        return len(err_sps_str.split(";"))
    return 1
