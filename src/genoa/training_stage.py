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
This module contains functions to read training parameters from a csv file or settings.
"""

import ast
import csv
import os
from typing import Any, List, Optional

from attrs import define, field, asdict

from .constants import RDC_SGY, TRN_MODE
from .logger import setup_logger
from .reduction_setting import rdc_variant_build, group_order_build
from .setting_init import TrainingOption
from .utils import get_condition_list


# Logger
logger = setup_logger(__name__)

# Parameters that need to be read from csv file
CSV_PARAMS = {  # name: type
    "err_max": [[float, list], None],  # Error tolerances, allow list, no default
    "err_ave": [[float, list], None],  # Allow list
    "delta_err": [[float, list], None],  # Allow list
    "strategies": [[str, list], RDC_SGY],  # Reduction strategies, default: RDC_SGY
    "path_cond": [[str], None],  # Path to training conditions, no default
    "ref_cases": [[str], None],  # Reference cases for reduction evaluation
    "nmin_retry": [[int], 0],  # Minimum number of accepted candidates to continue
    "nmax_eval_comb": [[int], rdc_variant_build["ncomb"]],  # Number of full combinations
    "nmax_eval_rdc": [[int], rdc_variant_build["nrdc"]],  # Number of reduction variant
    "nsps_group": [[int], group_order_build["nsmax"]],  # Max number of species in a group
    "flg_split_prod": [[int], 0],  # Elementary-like treatment
    "flg_aerosol": [[int], 0],  # Aerosol-oriented treatment
}

# Valid values for flag
FLG_VAL = {
    "flg_split_prod": [0, 1],
    "flg_aerosol": [0, 1, 2],
}


@define(slots=True)
class StageOption:
    """Stage Parameters that read with CSV_PARAMS"""

    # Error tolerances
    err_max: List[float] = field(factory=list)
    err_ave: List[float] = field(factory=list)
    delta_err: List[float] = field(factory=list)

    # Reduction strategies
    strategies: List[str] = field(factory=list)

    # Training conditions
    path_cond: str = ""
    conds: List[str] = field(factory=list)

    # Reduction evaluation
    ref_cases: List[str] = field(factory=list)

    # Training settings
    flg_split_prod: int = 0  # Elementary-like treatment
    flg_aerosol: int = 0  # Aerosol-oriented treatment
    nmin_retry: int = 0  # Efficient treatment
    nmax_eval_comb: int = 0  # Number of full combinations
    nmax_eval_rdc: int = 0  # Number of reduction variant
    nsps_group: int = 0  # Max number of species in a group

    def init_from_dict(self, istage: int, params_dict: dict) -> None:
        """Initialize the stage parameters from the dictionary."""

        if istage < 0 or istage >= params_dict["nstage"]:
            raise ValueError(f"Invalid stage number: {istage}")
        for k, v in params_dict.items():
            if k == "nstage":
                continue

            # Check if exists
            if hasattr(self, k):
                setattr(self, k, v[istage])
            else:
                raise ValueError(f"Invalid attribute: {k} not found in stage parameters.")

    def get_sginfo(self, istage: int) -> str:
        """Get the stage parameters information."""
        parts = [
            f"Training stage parameters for stage {istage}:",
            f"# {len(self.strategies)} strategies: {self.strategies}",
            f"# {len(self.conds) - 1} conditions from {self.path_cond}",
            f"# {len(self.ref_cases)} reference cases: {self.ref_cases}",
            "\nEvaluation Tolerances\tMaximum\tMean\tIncrease",
        ]

        # Error tolerances
        for merr, aerr, derr in zip(self.err_max, self.err_ave, self.delta_err):
            parts.append(f"{merr}\t{aerr}\t{derr}")

        # Treatment settings
        parts.append("\nTreatment settings:")
        if self.nmax_eval_comb <= 0:
            parts.append("Evaluation w/o reduction variant.")
            self.nmax_eval_rdc = 0  # No reduction variant
        else:
            parts.append(f"Evaluation w/ full reduction variants if candidate size <= {self.nmax_eval_comb}")
            parts.append(f"  Reduction variant size: {self.nmax_eval_rdc} for larger candidate.")
        if self.nsps_group > 0:
            parts.append(f"Max number of {self.nsps_group} species are allowed in a group.")
        if self.flg_aerosol:
            parts.append("Aerosol-oriented treatment")
            if self.flg_aerosol == 2:
                parts.append(" (Rerun cycle w/o aerosol-oriented treatment)")
        if self.flg_split_prod:
            parts.append("Elementary-like treatment")

        parts.append(f"Go to the next settings if No.accepted <= {self.nmin_retry}")

        return "\n".join(parts)


def setup_stage(trn_opt: TrainingOption, iopt: "CycleOption") -> None:
    """Setup training parameters related to stage options."""

    if not iopt.update_stage or iopt.is_stop:  # No update
        return

    # Get stage parameters
    stage_new = get_stage_option(iopt.istage, trn_opt.params_dict)
    if stage_new is None:  # Stop
        iopt.stop("Stop due to stage index exceeds the number of stages.")
        return
    logger.info("Updated stage parameters for stage %d.", iopt.istage)

    # Update stage_new if need
    stage_new.nsps_group = to_update_sgroup(stage_new.nsps_group, trn_opt)

    # Check flags
    iopt.update_ref = to_update_ref(stage_new, iopt.stage)
    iopt.update_kinetic = to_update_kinetic(stage_new, iopt.stage)

    # Update stage after checking flags
    iopt.stage = stage_new
    iopt.update_stage = False  # Reset


def get_stage_option(istage: int, params_dict: dict) -> Optional[StageOption]:
    """Get stage parameters from the dictionary."""
    if istage < 0:
        raise ValueError(f"Invalid stage number: {istage}")
    if istage >= params_dict["nstage"]:
        return None
    # Initialize the stage parameters
    stg = StageOption()
    stg.init_from_dict(istage, params_dict)
    return stg


def to_update_ref(new: StageOption, old: Optional[StageOption] = None) -> int:
    """Check if ref case related parameters are updated compared to the those in the previous stage"""

    if old is None:  # First stage
        return 3

    # Check ref cases
    if old.ref_cases != new.ref_cases:
        return 2
    # Check training conditions
    if os.path.abspath(old.path_cond) != os.path.abspath(new.path_cond):
        return 1

    return 0


def to_update_kinetic(new: StageOption, old: Optional[StageOption] = None) -> int:
    """Check if kinetics related parameters are updated compared to the those in the previous stage"""

    if old is None:  # First stage
        return 3

    # Check elementary-like treatment
    if old.flg_split_prod != new.flg_split_prod:
        return 2
    # Check reduction strategies
    for s in ["jp", "lp", "rp"]:
        if s in new.strategies and s not in old.strategies:
            return 1

    return 0


def to_update_sgroup(nsps_in: int, trn_opt: TrainingOption) -> int:
    """Check and update max number of species for species grouping if needed. Return nsmax for current stage."""

    if nsps_in <= 0 or TRN_MODE.get(trn_opt.group_order_mode, None) != "gprop":  # No limit
        return 0

    trn_opt.nsps_group = nsps_in  # Update
    return nsps_in


def read_training_param_dict(trn_opt: TrainingOption) -> None:
    """Read training parameters per reduction cycle."""

    # Get parameter dictionary
    if trn_opt.training_parameter_table:  # Get parameters from csv file
        param_dict = read_parameters_from_csv(trn_opt.training_parameter_table)
        logger.info("Read training parameters from %s", trn_opt.training_parameter_table)
    else:  # Get parameters from settings
        param_dict = get_training_params_from_dict(asdict(trn_opt))
        logger.info("Read training parameters from settings.")

    # Check if the values are valid
    update_param_dict_values(param_dict)

    trn_opt.params_dict = param_dict


def get_training_params_from_dict(trn_opt_dict: dict) -> dict:
    """Get default reduction stage parameters from settings. Only one reduction stage is considered."""

    # Get parameters from the dictionary
    param_dict = {k: [v] for k, v in trn_opt_dict.items() if k in CSV_PARAMS}

    # Add number of stages
    param_dict["nstage"] = 1

    return param_dict


def read_parameters_from_csv(filename: str, delimiter: str = "&") -> dict:
    """Read parameters for different reduction stages from the input file."""

    # Check if the file exists
    if not (filename and os.path.exists(filename)):
        raise FileNotFoundError(f"Training parameter file not found: {filename}")

    # Read the file
    param_dict, nrow0 = {}, 0
    with open(filename, "r", encoding="utf-8") as csvfile:

        for row in csv.reader(csvfile, delimiter=delimiter):
            row = [r.strip() for r in row]

            if not row or row[0].startswith("#"):
                continue  # Skip empty lines and comments

            # Check length of the row
            nrow = len(row)
            if nrow < 2:
                raise ValueError(f"Row {row} has less than 2 columns")

            if nrow0:
                if nrow != nrow0:
                    raise ValueError(f"Row {row} has inconsistent length with previous rows: {nrow0}")
            else:
                nrow0 = nrow  # Initialize the length of the row

            # Read the parameter name and values
            param_name = row[0]  # First column

            # Check if the parameter name is valid
            if param_name not in CSV_PARAMS:
                raise NameError(f"Invalid parameter name: {param_name}. Defaults: {CSV_PARAMS}")

            # Check for duplicate parameter names
            if param_name in param_dict:
                raise NameError(f"Duplicate parameter found: {param_name}")

            # Assign values to the parameter
            logger.info("Reading parameter %s with values: %s", param_name, row[1:])
            param_dict[param_name] = [ast.literal_eval(r) for r in row[1:]]

    # Add number of stages
    param_dict["nstage"] = nrow0 - 1

    return param_dict


def update_param_dict_values(param_dict: dict) -> None:
    """Check and update if the values in the parameter dictionary are valid."""

    for s, infos in CSV_PARAMS.items():
        ktypes, w_default = infos  # ktypes: expected types, w_default: default value

        if s in param_dict:  # Use read values
            new_values = get_param_values(param_dict[s], ktypes)
            has_none = False
            for k in new_values:
                if k is None or (isinstance(k, list) and None in k):
                    has_none = True
                    break
            if has_none:
                raise ValueError(f"Invalid {s} w/ value {new_values} in {param_dict[s]}.")
            param_dict[s] = new_values

        elif w_default is not None:  # Use default value
            param_dict[s] = [w_default] * param_dict["nstage"]

        else:  # No read or default
            raise ValueError(f"Parameter {s} not found in the dictionary")

    # Check values for each stage
    param_dict["conds"] = []  # Init
    for i in range(param_dict["nstage"]):

        # Check flags
        for k, v in FLG_VAL.items():
            if param_dict[k][i] not in v:
                raise ValueError(f"Invalid value for {k} in stage {i}: {param_dict[k][i]}. Valid values: {v}")

        # Check strategies
        if not isinstance(param_dict["strategies"][i], list):  # Convert to list
            param_dict["strategies"][i] = [param_dict["strategies"][i]]
        if [s for s in param_dict["strategies"][i] if s not in RDC_SGY]:
            raise ValueError(f"Invalid reduction strategies in stage {i}: {param_dict['strategies'][i]}")

        # Check training conditions
        if not os.path.exists(param_dict["path_cond"][i]):
            raise FileNotFoundError(f"Training conditions not found: {param_dict['path_cond'][i]}")
        param_dict["conds"].append(get_condition_list(param_dict["path_cond"][i]))

        # Check reference cases
        if not isinstance(param_dict["ref_cases"][i], list):  # Convert to list
            param_dict["ref_cases"][i] = [param_dict["ref_cases"][i]]
        for s in param_dict["ref_cases"][i]:
            if s.lower() not in ["ref", "pre"]:
                if not os.path.exists(s):
                    raise FileNotFoundError(f"Reference cases not found: {s}")

        # Check error tolerances
        nerr = len(param_dict["ref_cases"][i])  # No.errs should be the same as ref_cases
        for s in ["err_max", "err_ave", "delta_err"]:
            if isinstance(param_dict[s][i], list):
                if len(param_dict[s][i]) != nerr:
                    raise ValueError(f"{s} in stage {i} not has {nerr} of values: {param_dict[s][i]}")
            else:  # Convert to list
                param_dict[s][i] = [param_dict[s][i]] * nerr

            # Check values
            for v in param_dict[s][i]:
                if not 0 < v < 1:
                    raise ValueError(f"Invalid error value: {v} in {s} for stage {i}. Should be in (0, 1)")
        # err_max >= err_ave & delta_err
        for j in range(nerr):
            if param_dict["err_max"][i][j] < param_dict["err_ave"][i][j]:
                raise ValueError(f"err_max < err_ave in stage {i}.\nCheck {param_dict}")


def get_param_values(values: list, ktypes: List[type]) -> list:
    """Check if values are within the expected key type."""

    new_values = []
    for vals in values:
        if isinstance(vals, list):
            new_vals = [check_value_type(v, ktypes) for v in vals]
        else:
            new_vals = check_value_type(vals, ktypes)
        new_values.append(new_vals)

    return new_values


def check_value_type(value: Any, ktypes: List[type]) -> Any:
    """Check if the value is of the expected type."""

    for key_type in ktypes:
        if isinstance(value, key_type):
            return value
        try:
            return key_type(value)
        except (ValueError, TypeError):
            continue

    return None
