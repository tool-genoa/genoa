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
This updates reduction parameters and options used for the reduction.
"""

import os
import sys
import configparser
from ast import literal_eval
from typing import Any, Optional
from multiprocessing import cpu_count

from .constants import SNAME, RCT_DICT, UNITS_DICT, TRN_MODE
from .error_analysis import ERR_TYPE_ALL
from .folder_path import update_input_path, update_output_path
from .logger import setup_logger
from .reduction_setting import rdc_checks, log_reduction_settings
from .setting_global import get_all_settings_map
from .setting_init import (
    GlobalSetting,
    PostProcessOption,
    TbrOption,
    TrainingOption,
    TestingOption,
    check_and_set_attr,
)
from .record import setup_loge_4action
from .utils import isfloat, get_condition_list, restart_filename


# Logger
logger = setup_logger(__name__)


def _get_valid_keys(config: configparser.ConfigParser, section_name: str, if_remove: bool) -> set:
    """Return valid keys for processed sections read from config file."""

    # Check if the section exists
    if section_name not in config.sections():
        raise ValueError(f"Cannot find section {section_name} in the config file")

    # Get activated keys
    actived = set()
    for key, val in config[section_name].items():
        if val == "1" or val.lower() == "true":
            actived.add(key)

    # Remove the section if needed
    if if_remove:
        config.remove_section(section_name)

    return actived


def _get_info_to_dict(config: configparser.ConfigParser, section_name: str) -> dict:
    """Output information from a config section as a dictionary."""

    info_dict = {}

    # Get information
    for key, val in config[section_name].items():
        try:
            info_dict[key] = literal_eval(val)
        except (SyntaxError, ValueError):
            info_dict[key] = val

    return info_dict


def _read_cfg_file(cfg_file: str) -> dict:
    """Read information from a configuration (.ini) file"""

    logger.info("Read configuration file: %s ...", cfg_file)

    # Check if the file exists
    if not os.path.exists(cfg_file):
        raise FileNotFoundError(f"File {cfg_file} does not exist!")

    # Read the configuration file
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(cfg_file)

    # Actions
    active_actions = _get_valid_keys(cfg, SNAME.ACN, True)
    active_actions = [s for s in SNAME.ACNS if s in active_actions]
    if not active_actions:
        raise ValueError(f"No activated action in the config file {cfg_file}. Use {SNAME.ACN}")
    logger.info("Find valid actions: %s", active_actions)

    # Read parameters -> {section: {key: value}}
    param_dict = {SNAME.ACN: active_actions}
    for key in cfg.sections():
        if key not in active_actions and key not in SNAME.READ and key in SNAME.ACNS:
            logger.info("Skip section %s due to no activated action ...", key)
            continue
        logger.info("Reading section %s ...", key)
        param_dict[key] = _get_info_to_dict(cfg, key)

    logger.info("Finish reading configuration file. Read sections: %s", param_dict.keys())
    return param_dict


def update_parameters_from_file(cfgfile: str) -> list:
    """Read parameters from the configuration file. Return action dict."""

    # Get user-inputs
    param_dict = _read_cfg_file(cfgfile)

    # Initialize setting map
    settings_map = {k: func() for k, func in get_all_settings_map().items()}

    # Update action set
    settings_map[SNAME.SETUP].actions = param_dict.pop(SNAME.ACN)

    # Update settings with user-inputs
    logger.info("Updating parameters from user inputs ...")
    for dct_id, dct in param_dict.items():

        # Get keys pointed to options
        keys = []

        # Check 1st the specified options
        if dct_id in settings_map:
            keys.append(dct_id)
        # Check 2nd the general settings
        if SNAME.SETUP not in keys:
            keys.append(SNAME.SETUP)

        # Update parameters
        for k, v in dct.items():
            is_updated = False  # Check if updated
            for key in keys:
                is_updated = check_and_set_attr(settings_map[key], k, v)
                if is_updated:
                    logger.info("Updated %s to %s from [%s] for %s.", k, v, dct_id, key)
                    break

            if not is_updated:
                raise ValueError(f"Cannot update {k} to {v} from [{dct_id}] for {keys}")

    logger.info("Finish updating parameters from %s.", cfgfile)


def _get_basicsp_dict(species_files: list, unkept_sps: list) -> dict:
    """Obtain basic species dictionary {name: molar mass} from a set of files."""

    species_dict = RCT_DICT.copy()

    if not species_files:
        logger.warning("No basic species file is provided.")
        # Return default dict
        return species_dict

    # Read species list from files
    for species_file in species_files:
        logger.info("Reading basic species from file %s ...", species_file)
        with open(species_file, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        for line in lines:
            line = line.strip()
            if line == "" or line.startswith("#") or line.startswith("%"):
                continue
            # Get species name and molar mass
            parts = [i for i in line.replace("\t", " ").split(" ") if i != ""]
            if len(parts) != 2:
                raise ValueError(f"Not able to read from line: {line} in {species_file}")
            s, mw = parts[0], float(parts[1])

            # Check
            if mw <= 0.0:
                raise ValueError(f"Species {s} has mass {mw} <= 0.0 in {species_file}")
            if s in unkept_sps:
                continue

            # Load
            if s not in species_dict:
                species_dict[s] = mw
            else:  # already recorded, check value
                if species_dict[s] != mw:
                    raise ValueError(
                        f"Find species {s} with different mass: ",
                        f"{species_dict[s]} and {mw} read from {species_file}",
                    )

    logger.info("Loaded in total # %s basic non-reducible species.", len(species_dict))

    return species_dict


def update_general_settings(gnl: GlobalSetting) -> None:
    """
    Update and check general settings.
    Model-independent and action-independent settings are not updated here.
    """

    logger.info("Update general settings ...")
    # Primary VOCs
    if not gnl.primary_vocs:
        raise ValueError("Check primary_vocs in the config file.")

    # ncpu
    ncpu = cpu_count()
    if gnl.ncpu is None or gnl.ncpu < 1 or gnl.ncpu > ncpu:
        gnl.ncpu = ncpu
    logger.info("Use # %s CPUs for parallel simulation.", gnl.ncpu)

    # Mechanism name
    gnl.mech_name = gnl.prefix + gnl.id_mech
    logger.info("Mechanism name: %s", gnl.mech_name)

    # Reference concentration file if given
    if gnl.refconc_file:
        if not os.path.exists(gnl.refconc_file):
            raise FileNotFoundError(f"Cannot find reference concentration file/folder: {gnl.refconc_file}")
        logger.info("Reference concentration read from: %s & with fake: %s", gnl.refconc_file, gnl.refconc_wfake)

    # Update paths to absolute paths
    gnl.path_init_cond = update_input_path(gnl.path_init_cond)  # Stop if not exist

    # Create new if not exist
    gnl.path_sav_mech = update_output_path(gnl.path_sav_mech)
    gnl.path_sav_res = update_output_path(gnl.path_sav_res)
    gnl.path_workplace = update_output_path(gnl.path_workplace)

    # Setup read mechanism path
    if gnl.path_read_mech:
        gnl.path_read_mech = update_input_path(gnl.path_read_mech)
    else:
        gnl.path_read_mech = gnl.path_sav_mech

    # Simulation date
    for date in gnl.dates:
        year, month, day = date
        if day < 1 or day > 31 or month < 1 or month > 12 or year < 0:
            raise ValueError(f"Find invalid dates in (yy/mm/dd): {year}/{month}/{day}.")

    # Update non-reducible species list
    gnl.basicspecies_dict = _get_basicsp_dict(gnl.basic_species_files, gnl.primary_vocs)

    # Check and build tracer list
    gnl.tracers = _build_tracers(gnl.gain_tracers, gnl.loss_tracers)


def _build_tracers(gain: Optional[dict], loss: Optional[dict]) -> list:
    """Check tracer species and return a list of valid tracers for further use."""

    if not gain and not loss:  # No tracers
        return []

    if not (isinstance(gain, dict) and isinstance(loss, dict)):
        raise ValueError(f"Invalid tracers: {gain}, {loss}. Should be dicts w/ species & tracers as keys & values.")

    all_sps = set(gain.values()) | set(loss.values())
    for s in all_sps:
        if not s or s in gain or s in loss:
            raise ValueError(f"Invalid tracer species: {s}. Should be unique species names.")
    if len(all_sps) != len(gain) + len(loss):
        raise ValueError(f"Duplicate tracer species found in gain: {gain} and loss: {loss}.")

    logger.info("Tracer species from gain: %s and loss: %s.", gain, loss)
    return [gain, loss] if gain or loss else []  # Return empty list if no tracers


def update_gck_settings(settings_map: dict) -> None:
    """Update settings for using GECKO-aerosol model."""

    # Get settings
    gnl = settings_map[SNAME.SETUP]()
    opt = settings_map[SNAME.GCK]()

    # OMP and nsim for GECKO
    if opt.omp_num is None or opt.omp_num < 1 or opt.omp_num > gnl.ncpu:
        opt.omp_num = gnl.ncpu
    opt.nsim = max(gnl.ncpu // opt.omp_num, 1)
    os.environ["OMP_NUM_THREADS"] = str(opt.omp_num)  # Set OMP_NUM_THREADS
    logger.info("Set # %d OMP threads and # %d parallel GECKO-A simulations.", opt.omp_num, opt.nsim)

    # Time_settings
    tsettings = []  # Init
    pinfo = "Simulation time sets in [tstart, tstop, ntstep, iday, imonth, iyear] format for GECKO."
    for tstart in gnl.start_t:  # Starting time in s
        for t in gnl.total_t:
            tstop = tstart + t  # Total time in s
            for dt in gnl.delta_t:
                ntstep = int(t / dt)  # No. time steps
                for idate in gnl.dates:
                    year, month, day = idate
                    tsettings.append([tstart, tstop, ntstep, day, month, year])
    gnl.time_settings = tsettings  # Update
    logger.info("Use # %d %s: %s", len(tsettings), pinfo, tsettings)


def update_ssh_settings(settings_map: dict) -> None:
    """Update settings for using SSH-aerosol model."""

    # Get settings
    gnl = settings_map[SNAME.SETUP]()
    opt = settings_map[SNAME.SSH]()

    # Number of CPUs used for parallel simulation
    opt.nsim = gnl.ncpu
    logger.info("Use # %d CPUs for parallel simulation w/ SSH-aerosol.", opt.nsim)

    # Aerosol species list
    if not (opt.aero_file and os.path.isfile(opt.aero_file)):
        raise FileNotFoundError(f"Cannot find default aerosol species file {opt.aero_file} for SSH-aerosol.")

    # Time_settings
    tsettings = []  # Init
    pinfo = "Simulation time sets in [initial_time, final_time, delta_t, start_hour] format for SSH."
    for dates in gnl.dates:
        _, _, year = dates  # Only year is used
        if not isinstance(year, int) or year < 0:
            raise ValueError(f"Invalid year: {year}. Should be a positive integer.")
        for tstart in gnl.start_t:  # Starting time in s
            # Get hour in [0, 23]
            h = tstart // 3600 % 24
            for t in gnl.total_t:
                tend = tstart + t  # Total time in s
                for dt in gnl.delta_t:
                    tsettings.append([tstart, tend, dt, h, year])
    gnl.time_settings = tsettings  # Update
    logger.info("Use # %d %s: %s", len(tsettings), pinfo, tsettings)

    # Starting hour - used for initial files
    if opt.tag_init_hourly:
        if not tsettings:
            logger.warning("No time settings for hourly initial files. Set tag_init_hourly to 0.")
            opt.tag_init_hourly = 0
    elif gnl.start_t:
        if len(gnl.start_t) > 1:
            logger.warning("Multiple starting hours are given but tag_init_hourly is set to 0.")
            # opt.tag_init_hourly = 1

    # With monthly time settings - used for initial files
    if opt.tag_time_monthly:
        logger.info("Use monthly time settings for generating namelists.")

    # SOA calculation with group id
    if opt.soa_grps:
        if not isinstance(opt.soa_grps, list):
            raise ValueError(f"Invalid SOA groups: {opt.soa_grps}. Should be a list of PVOC groups.")
        soa_grps = []
        for sps in opt.soa_grps:
            # Check species
            for s in sps:
                if s not in gnl.primary_vocs:
                    raise ValueError(f"Species {s} in SOA groups is not a primary VOC. Check: {opt.soa_grps}")
            # Sort and build string
            sps = ",".join(sorted(sps))
            if sps in soa_grps:
                raise ValueError(f"Duplicate SOA group {sps} found. Check: {opt.soa_grps}")
            soa_grps.append(sps)
        if not soa_grps:
            raise ValueError("No valid SOA groups for calculation. Check: {opt.soa_grps}")

        opt.soa_grps = {s: i + 3 for i, s in enumerate(soa_grps)}  # Update as a dict, id default is 2
        logger.info("SOA groups for calculation: %s", opt.soa_grps)


def update_model_related_settings(settings_map: dict) -> None:
    """Update model-related settings."""

    # Check common settings
    box = update_common_model_settings(settings_map)

    # Specific settings for each model
    box_dict = {
        "GECKO": update_gck_settings,
        "SSH": update_ssh_settings,
    }
    if box in box_dict:
        logger.info("Updating specific settings for the %s model ...", box)
        box_dict[box](settings_map)


def update_common_model_settings(settings_map: dict) -> str:
    """Update common settings for all models."""

    # Model type
    box = settings_map[SNAME.SETUP]().boxmodel
    if box not in SNAME.BOX:
        raise ValueError(f"Unsupported model type: {box}. Available models: {SNAME.BOX}")

    # Get model settings
    box_sec = SNAME.BOX[box]
    opt = settings_map[box_sec]()
    logger.info("Update common settings for the %s model ...", box)

    # Path
    if not (opt.path and os.path.isdir(opt.path)):
        raise FileNotFoundError(f"Cannot find path {opt.path} for the {box} model.")
    opt.path = update_input_path(opt.path, True)

    # Namelist
    if not (opt.namelist and os.path.isfile(opt.namelist)):
        raise FileNotFoundError(f"Cannot find namelist file {opt.namelist} for the {box} model.")

    return box


def update_postprocess_settings(opt: PostProcessOption, gnl: GlobalSetting) -> None:
    """Update options for post-processing."""

    # Update plot dictionary
    opt.update_plot_dict()

    # Check if any plot is activated
    if opt.plots.get("no_plot", None):
        return

    # Mechanism names
    if not opt.mech_names:
        opt.mech_names = [gnl.mech_name]
    logger.info("Post-processing mechanisms: %s", opt.mech_names)

    # Mechanism labels
    if not opt.mech_labels:
        opt.mech_labels = opt.mech_names
    elif not isinstance(opt.mech_labels, list):
        raise ValueError(f"Invalid mechanism labels for post-processing. Got {opt.mech_labels}")

    # Mechanism folders
    if not opt.mech_paths or isinstance(opt.mech_paths, str):
        ipath = opt.mech_paths if opt.mech_paths else gnl.path_read_mech
        opt.mech_paths = [ipath for _ in opt.mech_names]
    elif not isinstance(opt.mech_paths, list):
        raise ValueError(f"Invalid mechanism folders for post-processing. Got {opt.mech_paths}")

    # Folder for saving post-process results
    if not opt.savpath:
        raise ValueError("No folder for saving post-process results.")
    if opt.unit:
        opt.savpath += f"_{opt.unit}"
    logger.info("Save post-process results to %s.", opt.savpath)

    if opt.plots.get("no_res_plot", None):  # No need to read results
        logger.info("No need to update settings for post-processing simulation results.")
        return

    # Result folders
    if not opt.mech_respaths or isinstance(opt.mech_respaths, str):
        ipath = opt.mech_respaths if opt.mech_respaths else gnl.path_sav_res
        opt.mech_respaths = [os.path.join(ipath, f"Results_{i}") for i in opt.mech_names]
    elif not isinstance(opt.mech_respaths, list):
        raise ValueError(f"Invalid result folders for post-processing. Got {opt.mech_respaths}")

    # Check
    for i, imech in enumerate(opt.mech_names):
        try:
            if not opt.mech_labels[i]:
                raise ValueError(f"Invalid label for mechanism {imech}")
            opt.mech_paths[i] = update_input_path(opt.mech_paths[i], True)
            opt.mech_respaths[i] = update_input_path(opt.mech_respaths[i], True)
        except ValueError as e:
            raise ValueError(f"Error in mechanism {imech}: {e}") from e

    # Check unit tag
    if not opt.unit:
        unit_dict = {"SSH": "ug", "GECKO": "molec"}
        opt.unit = unit_dict[gnl.boxmodel]
        logger.info("Assign unit tag: %s based on box model type: %s", opt.unit, gnl.boxmodel)
    if opt.unit not in UNITS_DICT:
        raise ValueError(f"Invalid unit tag: {opt.unit}")

    # Reference mechanism for comparison
    if opt.plot_mech_cmp:
        if isinstance(opt.plot_mech_cmp, str) and opt.plot_mech_cmp.lower() == "all":
            opt.plots["mech_cmp"] = opt.mech_labels
        elif not isinstance(opt.plot_mech_cmp, list):
            raise ValueError(f"Invalid mechanism comparison: {opt.plot_mech_cmp}")
        else:
            for s in opt.plot_mech_cmp:
                if s not in opt.mech_labels:
                    raise ValueError(f"Invalid mechanism label: {s} for comparison")
        # Compute errors
        if opt.err_types:
            logger.info("Use error types for mechanism comparison: %s", opt.err_types)
            if isinstance(opt.err_types, str) and opt.err_types.lower() == "all":
                opt.err_types = ERR_TYPE_ALL
            elif isinstance(opt.err_types, list):
                for s in opt.err_types:  # Check error types
                    if "_" in s:
                        s = s.rsplit("_", 1)[0].strip()
                    if s not in ERR_TYPE_ALL:
                        raise ValueError(f"Invalid error type: {s} read from {opt.err_types}")
            else:
                raise ValueError(f"Invalid error types: {opt.err_types}")
            if not opt.err_ref_mech:
                opt.err_ref_mech = gnl.mech_name
            elif opt.err_ref_mech not in opt.mech_labels:
                raise ValueError(f"Invalid reference mechanism: {opt.err_ref_mech} for comparison")
            logger.info("Reference mechanism for comparison: %s", opt.err_ref_mech)

    if opt.plot_mech_cmp2:
        if not isinstance(opt.plot_mech_cmp2, list):
            raise ValueError(f"Invalid mechanism comparison: {opt.plot_mech_cmp2}")
        for s in opt.plot_mech_cmp2:
            if not (isinstance(s, list) and len(s) == 2):
                raise ValueError(f"Invalid mechanism label pair: {s} for comparison")
            if s[0] not in opt.mech_labels or s[1] not in opt.mech_labels:
                raise ValueError(f"Invalid mechanism label pair: {s} for comparison")


def update_tbr_settings(opt: TbrOption, gnl: GlobalSetting) -> None:
    """Update options for threshold-based reduction (TBR)."""

    # Identifier
    if not opt.runid:
        raise ValueError("Threshold-based reduction id is not specified.")

    # Get mechanism name if not given
    if not opt.mech_name:
        opt.mech_name = gnl.prefix + opt.runid
    logger.info("Threshold-reduced mechanism name: %s", opt.mech_name)
    if opt.mech_name == gnl.mech_name:
        raise ValueError("Threshold-reduced mechanism name should be different from the original mechanism.")

    # Get reference mechanism
    if not opt.ref_mech_name:
        opt.ref_mech_name = gnl.mech_name
        opt.ref_mech_path = gnl.path_read_mech
    else:  # Check reference mechanism
        if not os.path.exists(opt.ref_mech_path):
            opt.ref_mech_path = gnl.path_read_mech
        else:
            opt.ref_mech_path = update_input_path(opt.ref_mech_path)
    logger.info("Reference mechanism: %s from %s.", opt.ref_mech_name, opt.ref_mech_path)

    # Get error references
    if opt.tag_sim:
        if not gnl.err_ref_str:  # Get results folder as error references
            gnl.err_ref_str = os.path.join(gnl.path_sav_res, f"Results_{opt.ref_mech_name}")
        logger.info("Reduction error references: %s", gnl.err_ref_str)

        # Check if folders exist and update to absolute path
        abs_dirs = []
        for res_dir in gnl.err_ref_str.split(","):
            if not os.path.isdir(res_dir):
                raise FileNotFoundError(f"Not found {res_dir} folder for threshold-based reduction error references.")
            abs_dirs.append(os.path.abspath(res_dir))
        gnl.err_ref_str = ",".join(abs_dirs)
        logger.info("Absolute path for reduction error references: %s", gnl.err_ref_str)

    # Log file
    setup_loge_4action(SNAME.TBR, gnl, opt)


def update_training_settings(opt: TrainingOption, gnl: GlobalSetting) -> None:
    """Update options for training."""

    # Identifier
    if not opt.runid:
        raise ValueError("Training id is not specified.")
    if not opt.runame:
        opt.runame = f"{gnl.prefix}{opt.runid}"
    logger.info("Training id: %s with mechanism saved name: %s", opt.runid, opt.runame)

    # Parallel or series for writing
    if opt.npara > 1:
        opt.npara = int(min(opt.npara, gnl.ncpu))
        logger.info("Parallelize # %s for training. (CAUTION! Too large may cause Out of Memory!)", opt.npara)
    else:
        logger.info("Use series mode for training.")

    # Paths
    if not opt.path_cond:
        opt.path_cond = gnl.path_init_cond
    opt.path_work = os.path.join(gnl.path_workplace, "training", f"{opt.runid}")
    if not opt.path_mech:
        opt.path_mech = os.path.join(gnl.path_sav_mech, f"{opt.runid}_mechs")
    logger.info("Training work path: %s and mechanism path: %s", opt.path_work, opt.path_mech)

    # Check restart_from file
    if isinstance(opt.restart_from, int):
        if opt.restart_from > 0:  # Set filename if any
            opt.restart_from = restart_filename(opt.path_mech, is_new=False, only_num=False)
    elif not isinstance(opt.restart_from, str):
        raise TypeError(f"Invalid type for restart_from: {opt.restart_from}. Expected str or int.")

    # Check time limiation for training and set restart file
    if opt.tlim:
        if not isfloat(opt.tlim) or opt.tlim <= 0:
            raise ValueError(f"Invalid time limitation for training: {opt.tlim}")
        logger.info("Set time limitation for training: %s seconds.", opt.tlim)
        if opt.restart_to:
            logger.warning("File to save restart training will be reset in training. Got %s", opt.restart_to)

    # Update paths to absolute paths
    opt.path_cond = update_input_path(opt.path_cond, True)
    opt.path_work = update_output_path(opt.path_work, True)
    opt.path_mech = update_output_path(opt.path_mech, not opt.restart_from)

    # Freeze compounds
    if opt.frozen_species:
        logger.info("Read frozen species that do not participate reduction: %s", opt.frozen_species)
    if opt.kept_species:
        logger.info("Read kept species that participate reduction: %s", opt.kept_species)
    opt.kept_all = set(opt.frozen_species + opt.kept_species + gnl.primary_vocs)
    logger.info("Total kept species: %s", opt.kept_all)
    if opt.kept_all:
        opt.trim_opts["check_kept"] = opt.kept_all

    # Starting mechanism
    if opt.mech_name is None:  # Use the default mechanism
        opt.mech_name = gnl.mech_name
    if not opt.mech_path:  # Use the default mechanism path
        opt.mech_path = gnl.path_read_mech
    logger.info("Starting mechanism: %s from %s.", opt.mech_name, opt.mech_path)

    # Reference mechanism
    if not opt.ref_mech_name:  # Use the starting mechanism as reference
        opt.ref_mech_name = opt.mech_name
        opt.ref_mech_path = opt.mech_path
    else:  # Check reference mechanism
        if not os.path.exists(opt.ref_mech_path):
            opt.ref_mech_path = gnl.path_read_mech
        else:
            opt.ref_mech_path = update_input_path(opt.ref_mech_path)
    logger.info("Reference mechanism: %s from %s.", opt.ref_mech_name, opt.ref_mech_path)

    # Check error tolerance for approaching
    if opt.ave_err_approach:  # If used
        if isfloat(opt.ave_err_approach) and 0.0 < opt.ave_err_approach < 1.0:
            logger.info("Error approaching is activated with tolerance: %s", opt.ave_err_approach)
        else:
            raise ValueError(f"Invalid error tolerance for approaching: {opt.ave_err_approach}")

    # Training order
    if opt.group_order_file:  # Read from file
        if not os.path.isfile(opt.group_order_file):
            raise FileNotFoundError(f"Cannot find training order file: {opt.group_order_file}")
        logger.info("Training order read from file: %s", opt.group_order_file)
    if opt.group_order_mode:
        if opt.group_order_mode not in TRN_MODE:
            raise ValueError(f"Invalid training order mode: {opt.group_order_mode}. Should be in {TRN_MODE.keys()}")
    elif not opt.group_order_file:
        raise ValueError("Training order mode is not specified and no training order file is given.")

    # Update mechanism trimming options if tracers are used
    if gnl.tracers:
        opt.trim_opts["tracer"] = gnl.tracers

    # Print reduction settings
    logger.info("%s", log_reduction_settings(rdc_checks))
    logger.info("Mechanism default trimming options: %s", opt.trim_opts)

    # Other parameters will be checked in training_init.py


def update_testing_settings(opt: TestingOption, gnl: GlobalSetting) -> None:
    """Update options for training."""

    # Only used if training is activated
    trn_opt = get_all_settings_map()[SNAME.TRN]() if SNAME.TRN in gnl.actions else None

    # Identifier
    if not opt.runid:
        if not trn_opt:
            raise ValueError("Testing id is not specified and can not be inferred from training.")
        opt.runid = trn_opt.runid
        logger.info("Got Testing id from training: %s", opt.runid)

    # Mechanisms
    if not opt.mech_name:
        opt.mech_name = gnl.mech_name
    if not opt.mech_path:
        opt.mech_path = gnl.path_read_mech
    if not opt.path_work:
        opt.path_work = os.path.join(gnl.path_workplace, "testing", f"{opt.runid}")
    if not opt.ref_mech_name:
        if not trn_opt:
            raise ValueError("No reference mechanism in testing and can not be inferred from training.")
        opt.ref_mech_name = trn_opt.ref_mech_name
        opt.ref_mech_path = trn_opt.ref_mech_path
        logger.info("Got reference mechanism from training: %s from %s", opt.ref_mech_name, opt.ref_mech_path)
    if not opt.ref_mech_path:
        opt.ref_mech_path = gnl.path_read_mech

    # Conditions
    if not opt.path_cond:
        opt.path_cond = gnl.path_init_cond
    if not opt.conds:
        opt.conds = get_condition_list(opt.path_cond)

    # Log file
    setup_loge_4action(SNAME.TST, gnl, opt)


def setup_4action(acn: str, opt: Any, gnl: GlobalSetting) -> None:
    """Update settings for activated actions. Return valid actions as a set."""

    acn_dict = {
        SNAME.PST: update_postprocess_settings,
        SNAME.TBR: update_tbr_settings,
        SNAME.TRN: update_training_settings,
        SNAME.TST: update_testing_settings,
    }

    if acn not in acn_dict:
        raise ValueError(f"Unknown action: {acn} to setup.")

    logger.info("Update settings for %s ...", acn)
    acn_dict[acn](opt, gnl)


def validate_settings() -> None:
    """Update parameters from user-inputs and validate all settings."""

    # Get settings
    settings_map = get_all_settings_map()
    general = settings_map[SNAME.SETUP]()

    # Update general settings
    update_general_settings(general)

    # Update model related settings
    update_model_related_settings(settings_map)

    # Return valid actions
    return general.actions


def set_genoa_parameters():
    """Set parameters from default anf user-defined values."""

    # Use only default values
    if len(sys.argv) <= 1:
        logger.warning("No input config file!")

    # Read configuration and update settings
    else:
        update_parameters_from_file(sys.argv[1])

    # Verify settings
    actions = validate_settings()

    logger.info("Finish setting up parameters.")

    # Return valid action as a set
    return actions
