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
This module contains the configuration of the logging files.
"""

import os
import sys
import time
import logging
from typing import Optional

from .constants import LOG, NOUT_LIM


# Global variables
_HANDLER = None


def get_handler() -> logging.Handler:
    """Get the handler for the logger."""
    global _HANDLER
    if _HANDLER is None:
        _HANDLER = _create_handlers(LOG.FLAG, LOG.FILE, LOG.LV_CON, LOG.LV_FIL)
    return _HANDLER


def setup_logger(name: str, level: int = LOG.LV_LOG, handler: Optional[logging.Handler] = None) -> logging.Logger:
    """Get logger for different modules."""

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if handler is not None:  # Use input handlers
        if logger.handlers:  # Clear existing handlers
            logger.handlers.clear()
        for h in handler:  # Add new handlers
            logger.addHandler(h)
    elif not logger.handlers:  # Use default handlers
        for h in get_handler():
            logger.addHandler(h)

    return logger


def list_to_log(log: logging.Logger, lst: list, msg: str = "", sep: str = "\n") -> None:
    """Write a list to log file and part to console if needed."""

    if not isout(log, logging.INFO):
        return

    if not lst or len(lst) <= NOUT_LIM:
        log.info("%s%s", msg, sep.join(lst))
    else:
        log.info("%s%s", msg, sep.join(lst[:NOUT_LIM]))
        log.debug("%s", sep.join(lst[NOUT_LIM:]))


def isout(log: logging.Logger, level: int) -> bool:
    """Check if log will be outputted."""
    if not log.handlers:
        return False
    return log.isEnabledFor(level)


def mute_all_loggers() -> None:
    """Set all loggers to error level."""
    for _, log in logging.Logger.manager.loggerDict.items():
        if isinstance(log, logging.Logger):
            log.setLevel(logging.WARNING)


def unmute_all_loggers() -> None:
    """Set all loggers to the default level."""
    for _, log in logging.Logger.manager.loggerDict.items():
        if isinstance(log, logging.Logger):
            log.setLevel(LOG.LV_LOG)


def get_log_suffix() -> str:
    """Get the suffix for the log file."""
    ini_filename = os.path.basename(sys.argv[1]).replace(".ini", "")
    return f".{ini_filename}.{time.strftime('%Y%m%d_%H%M%S')}.{os.getpid()}"


def _create_handlers(is_record: int, log_file_in: str, lv_c: int, lv_f: int) -> list:
    """Create handlers for the logger. Record to console and/or record to file."""

    pinfos = ["\nSetting up logger for running GENOA v3.0 ..."]

    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    handlers = []

    # Check record flag
    if is_record not in [0, 1, 2, 3]:
        raise ValueError(f"Invalid is_record flag. Should in [0, 3]. Got {is_record}.")

    # To console
    if is_record in [0, 1, 3]:
        pinfos.append(f"Logging to console with level {logging.getLevelName(lv_c)}.")
        console_handler = logging.StreamHandler()
        if is_record == 0:  # Print error and above to console regardless of is_record
            lv_c = max(lv_c, logging.ERROR)
        console_handler.setLevel(lv_c)
        console_handler.setFormatter(formatter)
        handlers.append(console_handler)

    # To file
    if is_record in [2, 3]:

        # Set up log file name
        log_file = log_file_in if log_file_in else "log"
        log_file += get_log_suffix()
        try:
            with open(log_file, "w", encoding="utf-8") as f:
                f.close()
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Cannot open log file {log_file}.") from e
        pinfos.append(f"Logging to file {log_file} with level {logging.getLevelName(lv_f)}.")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(lv_f)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    print("\n".join(pinfos), flush=True)
    return handlers
