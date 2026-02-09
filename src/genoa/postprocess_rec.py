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
This module contains functions to post-process the reduction records.
"""

import os
from typing import Any, Optional

import numpy as np
import matplotlib.pyplot as plt

from .constants import LINE_SEP, LINE_SEP2, RDC_SGY
from .folder_path import update_output_path
from .logger import setup_logger
from .record import RecordOption
from .setting_init import PostProcessOption
from .utils import size_difference


# Logger
logger = setup_logger(__name__)


# Plot settings
W_CYCLE = True  # Plot cycle settings
IS_SIZE_LOG = True  # Log scale for size changes


class RECS:
    """Default strings used in the record file."""

    # Cycle information
    TIDR = "Training mechanism:"  # Training mechanism identifier
    CHDR = "Start reduction cycle #"  # Identifier for the cycle header
    MECH = "with mechanism prefix:"  # Mechanism string identifier
    END = "END"  # End cycle string identifier
    PRT = "Pre-testing results:\n"  # Pre-testing string identifier
    STG = "Training stages"  # Training stages string identifier

    # Step information
    SHDR = "Cycle #"  # Identifier for the step header
    ACPT = "Accept:"  # Acceptance string identifier
    BRDC = "with score"  # Brdc string identifier
    RJCT = "No valid reduction found."  # Rejection string identifier
    RDCS = ["Rdc ID", "Stgy", "Size", "Score", "Max Errs", "Ave Errs", "Details"]  # Rdc line items


def plot_rdc_records(opt: PostProcessOption) -> None:
    """Post-process records and plot the results."""

    if not opt.plot_reduction:
        return

    logger.info("Plotting reduction records w/ input: %s ...", opt.plot_reduction)
    for items in opt.plot_reduction:
        if len(items) != 3:
            raise ValueError(f"Invalid record file setting. Got: {items}")
        if not os.path.isfile(items[0]):
            raise FileNotFoundError(f"Record file not found: {items[0]}")
    rfiles, starts, ends = zip(*opt.plot_reduction)
    logger.info("Process # %d record files w/ settings: %s, %s, %s ...", len(rfiles), rfiles, starts, ends)

    train_id = rfiles[0].rsplit("_", 1)[1].split(".")[0]
    savname = os.path.join(opt.savpath, f"rdc_{train_id}_{len(rfiles)}")
    logger.info("Save name: %s", savname)

    # Plot
    plot_record_files(rfiles, starts, ends, savname, opt.logt)


def floats_to_str(nums, is_err=False, sep=" / "):
    """Convert list of floats to 'a/b/c...' format with groups of two di."""
    parts = []
    for i in nums:
        if is_err:  # errors => x100 %
            j = i * 100
            s = str(int(j)) if round(j, 0) == j else f"{j:.1f}"
        else:  # numbers
            s = str(int(i))
            if i >= 1000:
                groups = []
                while s:
                    groups.append(s[-3:])
                    s = s[:-3]
                s = " ".join(reversed(groups))
        parts.append(s)
    return sep.join(parts)


def plot_record_files(record_files: list, starts: list, ends: list, savname: str, log: Optional[RecordOption]) -> None:
    """Plot valid reductions from a list of record files."""

    # Check savname & make output folder
    if not savname:
        raise ValueError("No save name provided.")
    update_output_path(os.path.dirname(os.path.abspath(savname)))

    # w/ keys: sizes, err_tols, steps, pretest
    cdata = read_record_files(record_files, starts, ends)
    if not cdata:
        raise ValueError("No valid cycles found in record files.")

    # e/ keys: steps, sizes, merrs, aerrs, etols, eprts
    pdata = process_record_data(cdata, IS_SIZE_LOG, log)

    # Plot reduction records to one figure
    plot_record_data(pdata, savname, IS_SIZE_LOG, W_CYCLE)


def process_record_data(cycle_data: list, size_log: bool, log: Optional[RecordOption]) -> dict:
    """Process cycle data and return a dictionary for plotting."""

    # Initialize data for plotting & printing
    sizes = [cycle_data[0][1]["sizes"]]  # Initial sizes
    stgys = {i: 0 for i in RDC_SGY}  # Reduction strategies
    etols, eprts = [], [[], []]  # Error tolerances & pre-testing errors (max, ave)
    merrs, aerrs = [], []  # Max and average errors
    nrdcs = 0  # Number of valid reductions
    # Print info w/ header
    pinfos = ["No.cycle\tMech.\tNo.vrdcs\tTrain.Err.Tols\tSize\tTrain.Err.Res\tPre-Test.Res\tStgy\tSize.diff"]

    # Read raw cycle data
    for [imech, cycle] in cycle_data:  # imech, cycle
        n = 0  # Number of valid reductions in current cycle
        stgy = {i: 0 for i in RDC_SGY}
        for step in cycle["steps"]:
            stgy[step[1]] += 1
            sizes.append(step[2])
            merrs.append(np.max(step[4]))
            aerrs.append(np.mean(step[5]))
            n += 1
        if n <= 0:
            continue
        # Record if valid reductions are found
        etols.append([nrdcs, nrdcs + n, [i * 100 for i in cycle["err_tols"]]])
        for i, err in enumerate(cycle["pretest"]):
            eprts[i].append([nrdcs + n, err * 100])
        nrdcs += n
        stgys = {k: stgys[k] + v for k, v in stgy.items()}
        isize = floats_to_str(sizes[-1])
        isize_diff = ",".join([f"{i:.2f}" for i in size_difference(sizes[0], sizes[-1], "rel100")])
        itol = floats_to_str(cycle["err_tols"], True)
        ie_t = floats_to_str([aerrs[-1], merrs[-1]], True, " & ")
        ie_p = floats_to_str(reversed(cycle["pretest"]), True, " & ")
        # Print info
        pinfo = f"{imech}\t{n}\t{itol}\t{isize}\t{ie_t}\t{ie_p}\t{stgy}\t{isize_diff}"
        pinfos.append(f"{len(pinfos)}\t{pinfo}")

    if log:
        log.list_to_file(pinfos, "Reducrtion cycle infos\n")
    else:
        print("\n".join(pinfos))

    # Post-process size info
    sizes = np.array(sizes)
    sizes = sizes / sizes[0]  # Percentage change
    if size_log:
        if np.isin(0, sizes):
            raise ValueError("Log scale not valid for zero values.")
    else:
        sizes *= 100  # Percentage change
    steps = np.arange(len(sizes))

    # Post-process error info
    merrs = np.array(merrs) * 100  # Percentage error
    aerrs = np.array(aerrs) * 100  # Percentage error

    return {"steps": steps, "sizes": sizes, "merrs": merrs, "aerrs": aerrs, "etols": etols, "edots": eprts}


def plot_record_data(plot_data: dict, savname: str, with_log: bool, with_cycle: bool) -> None:
    """Plot the processed record data."""

    # Settings
    ftsize = 12  # Font size

    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    # Plot size changes
    (s1,) = ax1.plot(plot_data["steps"], plot_data["sizes"][:, 0], "b-", label="Rcn.", linewidth=2)
    (s2,) = ax1.plot(plot_data["steps"], plot_data["sizes"][:, 1], "k:", label="Sps.", linewidth=2)
    (s3,) = ax1.plot(plot_data["steps"], plot_data["sizes"][:, 2], "g--", label="SVOCs", linewidth=2)

    # Plot error info
    (e1,) = ax2.plot(plot_data["steps"][1:], plot_data["merrs"], color="m", linestyle="-.", label="Max SBR.Err")
    (e2,) = ax2.plot(plot_data["steps"][1:], plot_data["aerrs"], color="r", linestyle="-", label="Ave SBR.Err")

    lgd = [s1, s2, s3, e1, e2]

    # Plot cycle infos
    if with_cycle:

        # Error tolerances
        for i, (x0, x1, err) in enumerate(plot_data["etols"]):
            e1 = ax2.fill_betweenx([0, err[1]], x0, x1, alpha=0.4, color="r" if i % 2 else "pink", label="Ave SBR.Tol")
            e2 = ax2.fill_betweenx(
                [err[1], err[0]], x0, x1, alpha=0.2, color="gray" if i % 2 else "k", label="Max SBR.Tol"
            )
        lgd.extend([e1, e2])  # Record two legends for all fill plots

        # Pre-testing errors
        lbs = [None, "Ave TBR.Err"]  # ["Max TBR.Err", "Ave TBR.Err"]
        for lb, err, x0, x1 in zip(lbs, plot_data["edots"], ["pink", "orange"], ["^", "o"]):
            if err and lb:
                e1 = ax2.scatter(*zip(*err), color=x0, marker=x1, label=lb, clip_on=False)
                lgd.append(e1)

    # Plot settings
    if with_log:  # Log scale for size changes
        ax1.set_yscale("log")
        ax1.set_ylim(top=1)
        y1label = "Size change in number (factor)"
        ax1.set_ylabel("Size change in number (factor)", color="b", fontsize=ftsize)
    else:
        ax1.set_ylim(0, 100)
        y1label = "Size change in number (%)"

    ax1.set_xlabel("Number of valid reductions (#)", fontsize=ftsize)
    ax1.set_ylabel(y1label, color="b", fontsize=ftsize)
    ax1.set_xlim(plot_data["steps"][0], plot_data["steps"][-1])

    ax2.set_ylabel("Training errors (%)", color="r", fontsize=ftsize)
    ax2.set_ylim(bottom=0)  # Error range

    plt.legend(handles=lgd, scatterpoints=1, framealpha=0.5, frameon=False, loc="center left", fontsize=ftsize)
    plt.tight_layout()

    # Save figure
    savname += "_log" if with_log else ""
    savname += "_wprt" if with_cycle else ""
    fig.savefig(f"{savname}.png", dpi=100)
    logger.info("Record plot is saved as %s.png", savname)


def read_record_files(record_files: list, starts: list, ends: list) -> list:
    """Read valid reductions from a list of record files."""

    if not record_files:
        raise ValueError("No record files provided.")

    if not len(starts) == len(ends) == len(record_files):
        raise ValueError(f"Invalid number. Got: {len(starts)}, {len(ends)} for {len(record_files)} files.")

    cycle_data = []  # Reduction cycle data
    for i, record_file in enumerate(record_files):

        if not os.path.isfile(record_file):
            raise FileNotFoundError(f"Record file not found: {record_file}")

        with open(record_file, "r", encoding="utf-8") as f:
            cycles = f.read().split(LINE_SEP)

        find_start = not bool(starts[i])
        find_end = False

        for cycle in cycles:
            # Get previous mech
            imech = _get_pre_mech(cycle)
            if not imech:  # No previous mech
                continue
            if not find_start and imech == starts[i]:
                find_start = True
            if not find_end and imech == ends[i]:
                find_end = True

            # Check break condition
            if find_end and imech != ends[i]:
                logger.info("Stop reading cycle info as mech: %s > %s ...", imech, ends[i])
                break

            # Read cycle info
            if find_start:
                logger.info("Reading cycle info from mech %s from %s", imech, record_file)
                cinfo = _read_cycle(cycle)
                if cinfo:
                    cycle_data.append([imech, cinfo])

        if not find_start:
            raise ValueError(f"Starting mech {starts[i]} not found in record file: {record_file}")

    logger.info("Found %s valid cycles from # %d files.", len(cycle_data), len(record_files))
    return cycle_data


def _get_pre_mech(intext: str) -> Optional[str]:
    """Read input text block and return the previous mechanism name."""

    # Check if is a cycle block
    if not intext:
        return None

    # Get previous mech name from a cycle block
    if intext.startswith(RECS.CHDR):
        return intext.split("\n", 1)[0].split(RECS.MECH, 1)[1].strip()

    # Get previous mech name from a training block
    if RECS.TIDR in intext:
        return intext.split(RECS.TIDR, 1)[1].split("in", 1)[0].strip()

    logger.warning("No previous mechanism found in input text. Got: %s", intext)
    return None


def _read_cycle(intext: str) -> Optional[dict]:
    """Read input text block and return information about the reduction cycle."""

    # Check if is a cycle block
    if not intext or not intext.startswith(RECS.CHDR):
        if RECS.STG not in intext:
            logger.warning("Input text is not a cycle block. Got: %s", intext)
        return None

    # Split block
    blocks = intext.split(LINE_SEP2)
    if len(blocks) < 3:  # No reduction steps
        logger.warning("No reduction steps found in cycle block. Got: %s", intext)
        return None

    # Read cycle settings from first block
    cycle = _get_cycle_settings(blocks[0])

    # Read pre-testing info
    if not _update_cycle_wpretest(cycle, blocks):
        logger.warning("No pre-testing info found in cycle block.")

    # Read reduction step info
    if not _update_cycle_wrdcs(cycle, blocks):
        logger.warning("No valid reduction steps found in cycle block.")
        return None

    return cycle


def _update_cycle_wrdcs(cycle: dict, blocks: list) -> bool:
    """Update reduction steps info in the cycle settings. Return True if valid steps are found."""

    step_info = []
    for b in blocks[1:-1]:
        info = _get_brdc_info(b)
        if info:
            step_info.append(info)

    if not step_info:  # No valid steps
        return False

    logger.info("Found %s valid reductions in cycle.", len(step_info))
    cycle["steps"] = step_info
    return True


def _update_cycle_wpretest(cycle: dict, blocks: list) -> bool:
    """Update pre-testing info in the cycle settings. Return True if valid pre-testing info is found."""

    for b in [blocks[-2], blocks[-3]]:  # Last two blocks
        info = _get_pretest_info(b)
        if info:
            cycle["pretest"] = info
            logger.info("Found pre-testing info: %s", info)
            return True
    cycle["pretest"] = []  # No pre-testing info
    return False


def _get_pretest_info(intext: str) -> Optional[list]:
    """Read input text block and return information about the pre-testing settings."""
    if not (intext and RECS.PRT in intext):
        return None
    err_strs = ["Max", "Ave"]
    infos = []
    lines = intext.split(RECS.PRT)[1].split("\n")
    for line in lines:
        line = line.strip()
        if not line:
            continue
        for err_str in err_strs:
            if line.startswith(err_str):
                infos.append(float(line.split(":")[1].split("at")[0].strip()))
    if len(infos) != 2:
        raise ValueError(f"Invalid number of error values. Got: {infos} from {lines} in {intext}")
    return infos


def _get_cycle_settings(intext: str) -> dict:
    """Read input text block and return information about the reduction cycle settings."""

    if not intext or not intext.startswith(RECS.CHDR):
        raise ValueError(f"Invalid cycle header block. Got: {intext}")

    infos = {}
    lines = intext.split("\n")
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        # Get size info
        if line.startswith("Mechanism size"):
            sizes = [int(i) for i in line.split(":")[1][:-1].split(",") if i.strip()]
            if len(sizes) != 3:
                raise ValueError(f"Invalid number of sizes. Got: {sizes}")
            infos["sizes"] = sizes
        # Get error tolerance
        elif line.startswith("Evaluation Tolerances"):
            tols = [float(i) for i in lines[i + 1].split("\t") if i.strip()]
            if len(tols) != 3:
                raise ValueError(f"Invalid number of error tolerances. Got: {tols}")
            infos["err_tols"] = tols
    return infos


def _get_brdc_info(intext: str) -> Optional[list]:
    """Read input text block and return information about the accepted reduction brdc."""

    # Check if is a step block
    if not intext or not intext.startswith(RECS.SHDR):
        logger.info("Input text is not a step block. Got: %s", intext)
        return None

    # Check if is rejected
    if RECS.RJCT in intext:
        return None

    # Split block
    blocks = intext.split(RECS.ACPT)
    if len(blocks) != 2:
        raise ValueError(f"Invalid number of blocks for accepted record. Got: {intext}")

    # Get brdc index
    brdc = int(blocks[1].split(RECS.BRDC)[0].strip())

    # Find brdc info
    bitems = None
    for line in blocks[0].split("\n"):
        if not (line and "\t" in line):
            continue
        items = line.split("\t")
        if items[0].strip().isdigit() and int(items[0]) == brdc:
            bitems = items
            break
    if not bitems:
        raise ValueError(f"Accepted brdc not found. Got: {intext}")

    # Check if is a valid rdc line
    if len(items) != len(RECS.RDCS):
        raise ValueError(f"Invalid number of items in rdc line. Got: {items}")
    return [_parse_rdc_str(i) for i in items]


def _parse_rdc_str(item: str) -> Any:
    """Parse input item to the correct type."""
    item = item.strip()
    if not item:  # Empty item
        logger.error("Empty item found in rdc line.")
        return None
    if item.isdigit():  # Integer
        return int(item)
    if item.startswith("["):  # List of floats
        if "," in item:  # Comma separated
            return [float(i) for i in item[1:-1].split(",")]
        return [float(i) for i in item[1:-1].split(" ") if i.strip()]  # Space separated
    return item
