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
This module contains the global variables.
"""

from .constants import SNAME
from .setting_init import (
    GlobalSetting,
    GCKOption,
    NewMechOption,
    PostProcessOption,
    SMLOption,
    SSHOption,
    TestingOption,
    TbrOption,
    TrainingOption,
    EnvOption,
)

# Global variables
_GLOBAL_SETTINGS = None
_NEW_MECH_OPTION = None
_SML_OPTION = None
_SSH_OPTION = None
_GCK_OPTION = None
_POST_PROCESS_OPTION = None
_TBR_OPTION = None
_TRAINING_OPTION = None
_TESTING_OPTION = None
_ENV_OPTION = None


def _get_global_settings() -> GlobalSetting:
    """Get the global settings."""
    global _GLOBAL_SETTINGS
    if _GLOBAL_SETTINGS is None:
        _GLOBAL_SETTINGS = GlobalSetting()
    return _GLOBAL_SETTINGS


def _get_new_mech_option() -> NewMechOption:
    """Get the new mechanism options."""
    global _NEW_MECH_OPTION
    if _NEW_MECH_OPTION is None:
        _NEW_MECH_OPTION = NewMechOption()
    return _NEW_MECH_OPTION


def _get_sml_option() -> SMLOption:
    """Get the simulation options."""
    global _SML_OPTION
    if _SML_OPTION is None:
        _SML_OPTION = SMLOption()
    return _SML_OPTION


def _get_ssh_option() -> SSHOption:
    """Get the SSH-aerosol options."""
    global _SSH_OPTION
    if _SSH_OPTION is None:
        _SSH_OPTION = SSHOption()
    return _SSH_OPTION


def _get_gck_option() -> GCKOption:
    """Get the GECKO-A options."""
    global _GCK_OPTION
    if _GCK_OPTION is None:
        _GCK_OPTION = GCKOption()
    return _GCK_OPTION


def _get_post_process_option() -> PostProcessOption:
    """Get the post-processing options."""
    global _POST_PROCESS_OPTION
    if _POST_PROCESS_OPTION is None:
        _POST_PROCESS_OPTION = PostProcessOption()
    return _POST_PROCESS_OPTION


def _get_tbr_option() -> TbrOption:
    """Get threshold-based reduction options."""
    global _TBR_OPTION
    if _TBR_OPTION is None:
        _TBR_OPTION = TbrOption()
    return _TBR_OPTION


def _get_training_option() -> TrainingOption:
    """Get the training options."""
    global _TRAINING_OPTION
    if _TRAINING_OPTION is None:
        _TRAINING_OPTION = TrainingOption()
    return _TRAINING_OPTION


def _get_testing_option() -> TestingOption:
    """Get the testing options."""
    global _TESTING_OPTION
    if _TESTING_OPTION is None:
        _TESTING_OPTION = TestingOption()
    return _TESTING_OPTION


def get_env_option() -> EnvOption:
    """Get the environment options."""
    global _ENV_OPTION
    if _ENV_OPTION is None:
        _ENV_OPTION = EnvOption()
    return _ENV_OPTION


def get_all_settings_map() -> dict:
    """Mapping all settings."""
    return {
        SNAME.SETUP: _get_global_settings,
        SNAME.ENV: get_env_option,
        SNAME.NEW: _get_new_mech_option,
        SNAME.SML: _get_sml_option,
        SNAME.PST: _get_post_process_option,
        SNAME.TBR: _get_tbr_option,
        SNAME.TRN: _get_training_option,
        SNAME.TST: _get_testing_option,
        SNAME.SSH: _get_ssh_option,
        SNAME.GCK: _get_gck_option,
    }


def get_attrs_from_settings(attr_dict: dict):
    """Get the attributes from settings specified in the map."""

    outputs = {}
    for k, v in attr_dict.items():

        # Check setting name
        if k not in get_all_settings_map():
            raise ValueError(f"Key {k} not found in the settings map.")

        # Check input attr names
        if isinstance(v, str):
            attr_names = [v]
        elif isinstance(v, list):
            attr_names = v
        else:
            raise ValueError(f"Input must be a string or a list of strings. Got {v}")

        # Get settings
        settings = get_all_settings_map()[k]()

        # Get attributes
        attr_objs = []
        for s in attr_names:
            if not hasattr(settings, s):
                raise AttributeError(f"Attr. {s} not found in {k} settings.")
            attr_objs.append(getattr(settings, s))

        # Return a single object
        if len(attr_objs) == 1:
            attr_objs = attr_objs[0]

        outputs[k] = attr_objs

    # Return a single object
    if len(outputs) == 1:
        return attr_objs

    return outputs
