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
This module contains functions for plotting figures.
"""

import os
from typing import Optional

import attrs
import numpy as np
import matplotlib.pyplot as plt

from .constants import UNITS_DICT, EPS
from .logger import setup_logger


# Logger
logger = setup_logger(__name__)

# Plot settings
STYLES = ["-", "--", "-.", ":"]
COLORS = ["b", "g", "r", "c", "m", "y", "k"]
MARKERS = ["o", "s", "D", "v", "^", "<", ">", "p", "P", "*", "X", "d", "H", "h"]
COLOR_MAP = [
    "Blues_r",
    "Oranges_r",
    "Greens_r",
    "Purples_r",
    "YlOrBr_r",
    "hot_r",
    "Reds_r",
    "spring_r",
    "summer_r",
    "autumn_r",
    "winter_r",
    "cool_r",
]
NCMP = len(COLOR_MAP)  # Number of color maps
FGSIZE_L = (8, 6)  # Large figure size
FGSIZE_S = (5, 4)  # Small figure size
DPI = 100  # Resolution
SDOT = 50  # Dot size

# Set font size
plt.rcParams.update({"font.size": 12})


@attrs.define(slots=True)
class FigureOption:
    """Options for figure plotting."""

    # Plot data
    x: Optional[np.ndarray] = None  # Time steps
    ylists: Optional[list] = None
    ykeys: Optional[list] = None

    # Save option
    sfx: str = "test"  # Suffix for the figure name
    svnm: str = "./figure_"
    mechs: Optional[list] = None  # Mechanism names

    # Figure settings
    style: str = "line"  # line/stack/2d
    is_norm: bool = False  # Normalize ylists
    pltset: dict = attrs.field(factory=dict)  # Plot settings

    # For stacked plot
    wline: Optional[list] = None  # Lines for stacked plot
    wdot: Optional[list] = None  # Dots for 2D plot

    def update_to_norm(self, mode: str) -> None:
        """Update ylists to normalized values."""
        if self.is_norm:
            raise ValueError("ylists are already normalized")

        if self.style == "2d":
            self.ylists = [normalize_list(s, mode) for s in self.ylists]
        else:
            self.ylists = normalize_list(self.ylists, mode)
        if self.wline:  # Reset wline
            self.wline = None
        self.is_norm = True
        self.sfx += "_norm"
        self.pltset["ylabel"] = "Fraction"
        self.pltset["title"] = f"Normalized {self.pltset['title']}"
        self.pltset.pop("unit", None)


def normalize_list(alist: list, mode: str) -> list:
    """Normalize a list of arrays."""
    if mode == "sum":
        val = np.sum(alist, axis=0)
    elif mode == "max":
        val = np.max(alist, axis=1)
    elif mode == "init":
        val = np.array([a[0] for a in alist])
    else:
        raise ValueError(f"Invalid mode for normalization: {mode}")
    return np.where(val > EPS, alist / val, 0)


def plot_dict_to_1fig(fopt: FigureOption, tag_csv: bool = False) -> None:
    """Plot one stacked/line figure from a dictionary of data."""

    plot_dict = {
        "stack": plot_stack_fig,
        "line": plot_line_fig,
        "2d": plt_2d_fig,
    }

    if fopt.style not in plot_dict:
        raise ValueError(f"Invalid style for plotting: {fopt.style}. Available: {plot_dict.keys()}")
    plt.figure(figsize=FGSIZE_L if fopt.style == "stack" else FGSIZE_S)
    plot_dict[fopt.style](fopt)
    _setup_and_save_figure(fopt)
    if tag_csv:
        _save_csv(fopt)


def plot_stack_fig(fopt: FigureOption) -> None:
    """Plot stacked figure."""

    plt.stackplot(fopt.x, *fopt.ylists, labels=fopt.ykeys, alpha=0.5)
    # Add lines
    if fopt.wline:
        for i, (k, v) in enumerate(fopt.wline):
            plt.plot(fopt.x, v, label=k, color=COLORS[i % len(COLORS)], linestyle=STYLES[i % len(STYLES)])


def plot_line_fig(fopt: FigureOption) -> None:
    """Plot line figure."""
    for i, (k, v) in enumerate(zip(fopt.ykeys, fopt.ylists)):
        plt.plot(fopt.x, v, label=k, color=COLORS[i % len(COLORS)], linestyle=STYLES[i % len(STYLES)])


def plt_2d_fig(fopt: FigureOption) -> None:
    """Plot 2D figure."""

    if not len(fopt.ylists) == len(fopt.mechs) == 2:
        raise ValueError("2D plot requires 2 lists of data and 2 mechanisms")
    if not len(fopt.ykeys) == len(fopt.ylists[0]) == len(fopt.ylists[1]):
        logger.error("ykeys: %s, lens: %s", fopt.ykeys, [len(i) for i in fopt.ylists])
        raise ValueError("2D plot requires the same number of keys and data points")
    for s, x, y in zip(fopt.ykeys, fopt.ylists[0], fopt.ylists[1]):
        (l,) = plt.plot(x, y, label=s, alpha=0.7)
        # Add dots
        if fopt.wdot:
            c = l.get_color()
            for i, (j, _) in enumerate(fopt.wdot):
                plt.scatter(x[j], y[j], c=c, marker=MARKERS[i % len(MARKERS)], s=SDOT, alpha=0.7)

    # Add legend for dots
    if fopt.wdot:
        for i, (_, s) in enumerate(fopt.wdot):
            plt.scatter([], [], c="grey", marker=MARKERS[i % len(MARKERS)], s=SDOT, label=s)

    # Plot 1:1 line
    min_xy, max_xy = np.min(fopt.ylists), np.max(fopt.ylists)
    plt.plot([min_xy, max_xy], [min_xy, max_xy], "k--", label="1:1", alpha=0.5)


def _setup_and_save_figure(fopt: FigureOption) -> None:
    """Set up figure settings."""

    # Based on style of plot
    if fopt.style == "2d":
        # Set figure to square
        plt.gca().set_aspect("equal", adjustable="box")
    else:
        plt.xlim(fopt.x[0], fopt.x[-1])
        if "xlabel" not in fopt.pltset:
            fopt.pltset["xlabel"] = "Time (hour)"
        if "unit" in fopt.pltset:
            unit = fopt.pltset["unit"]
            fopt.pltset["ylabel"] = "Mixing ratio (ppb)" if unit == "ppb" else f"Concentration ({UNITS_DICT[unit]})"

    # Based on input dict: fopt.pltset
    # Labels
    plt.xlabel(fopt.pltset.get("xlabel", None))
    plt.ylabel(fopt.pltset.get("ylabel", None))
    # Title
    if "title" in fopt.pltset:
        plt.title(fopt.pltset["title"].replace("_", " "), pad=10)
    # Scale
    if fopt.pltset.get("log", False):
        plt.yscale("log")
        plt.xscale("log")
    # Grid
    if fopt.pltset.get("grid", False):
        plt.grid(True, linestyle="--", alpha=0.5)
    # Legend
    ltitle = fopt.pltset.get("lgd_title", None)
    if len(fopt.ykeys) <= 3:
        plt.legend(loc="best", framealpha=0.5, title=ltitle)
    else:
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), framealpha=0.5, title=ltitle)

    plt.tight_layout()
    # plt.show()

    # Save figure
    savname = f"{fopt.svnm}{fopt.sfx}_{fopt.style}.png"
    plt.savefig(savname, dpi=DPI)
    logger.info("Figure saved as %s", savname)
    plt.close()


def _save_csv(fopt: FigureOption) -> None:
    """Save data to csv file."""
    if fopt.style == "2d":
        keys = [f"{i}_{j}" for i in fopt.mechs for j in fopt.ykeys]
        vals = [j for i in fopt.ylists for j in i]
    else:
        keys, vals = fopt.ykeys, fopt.ylists
    hdr = f"{fopt.sfx}: Time(s)," + ", ".join(keys)
    savname = f"{fopt.svnm}{fopt.sfx}_{fopt.style}.csv"
    np.savetxt(savname, np.vstack((fopt.x, *vals)).T, delimiter=",", header=hdr)


def plot_ekma(ratios: list, steps: list, nskip_in: str, savname: str, tag_csv: bool = False) -> None:
    """Plot EKMA diagram for NOx vs. HC ratios."""

    # Settings
    clim = 0.8  # Color upper limit

    # Time points
    if nskip_in.isdigit():
        nskip = max(1, int(nskip_in))
        clb = f"Time Steps with {nskip} skip" if nskip > 1 else "Time Step"
    elif nskip_in == "h":
        nskip = int(1.0 / (steps[1] - steps[0]))
        clb = "Time (h)"
        logger.info("EKMA plotting with hourly time steps: %d", nskip)
    else:
        raise ValueError(f"Invalid time step skip: {nskip_in} for EKMA plotting")
    nsteps = range(len(steps))[::nskip]
    nt = len(nsteps)
    if nt < 2:
        raise ValueError(f"Invalid time step skip: {nskip} for EKMA plotting. nsteps: {len(steps)}")
    nmech = len(ratios)

    # Colors of time variation
    if nmech > NCMP:
        raise ValueError(f"Number of keys {nmech} is larger than the number of color maps {NCMP}")
    colors = {k: plt.get_cmap(COLOR_MAP[i])(np.linspace(0, clim, nt)) for i, (k, _) in enumerate(ratios)}

    # Plot
    plt.figure(figsize=FGSIZE_S if nmech == 1 else FGSIZE_L)
    for i, (k, pdata) in enumerate(ratios):
        x, y = pdata["hc"][nsteps], pdata["nox"][nsteps]
        plt.scatter(x, y, c=colors[k], marker=MARKERS[i % len(MARKERS)], s=SDOT, label=k, alpha=0.8)

    # Plot the limit line: [HCs] / [NOx] = 8
    range_x, range_y = plt.xlim(), plt.ylim()
    min_xy, max_xy = min(range_x[0] / 8, range_y[0]), max(range_x[1], range_y[1] * 8)
    plt.plot([min_xy, max_xy], [min_xy, max_xy / 8], "b--", label="HCs/NOx = 8")

    # Colorbar legend for time
    sm = plt.cm.ScalarMappable(cmap="Greys_r" if nmech > 1 else COLOR_MAP[0], norm=plt.Normalize(vmin=0, vmax=nt))
    sm.set_array([])
    sm.set_clim(0, nt * clim)
    plt.gca().figure.colorbar(sm, ax=plt.gca(), label=clb, orientation="vertical", aspect=40, pad=0.02)

    # Layout
    plt.xlabel("HCs (ppbC)")
    plt.ylabel("NOx (ppb)")
    title = os.path.basename(savname).replace("_", " ")
    plt.title(title if title else "EKMA Diagram", pad=10)

    # Legend
    # plt.legend(loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=nmech + 1)
    plt.legend(loc="best", framealpha=0.5)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"{savname}_{nskip_in}.png")

    # Save zoom-in figure
    plt.xlim(range_x)
    plt.ylim(range_y)
    plt.savefig(f"{savname}_{nskip_in}_zoom.png")
    plt.close()

    # Save data
    if tag_csv:
        csv_data, hdr = [steps], "Time(s),"  # Data and header
        for k, pdata in ratios:
            hdr += f"{k}_NOx(ppb),{k}_HCs(ppbC),"
            csv_data.extend([pdata["nox"], pdata["hc"]])
        np.savetxt(f"{savname}_{nskip_in}.csv", np.asarray(csv_data).T, delimiter=",", header=hdr)


def plot_errors(sdiff_dict: dict, err_dict: dict, title: str, savname: str) -> None:
    """Plot errors and size changes for different conditions."""

    plt.figure(figsize=FGSIZE_S if len(sdiff_dict) <= 3 else FGSIZE_L)
    colors = {m: COLORS[i % len(COLORS)] for i, m in enumerate(list(sdiff_dict.keys()))}
    shapes = {s: MARKERS[i % len(MARKERS)] for i, s in enumerate(list(err_dict.keys()))}

    # Plot data
    mechs = list(sdiff_dict.keys())
    for cond, einfo in err_dict.items():
        for i, m in enumerate(mechs):
            plt.plot(i + 1, einfo[m], color=colors[m], marker=shapes[cond], alpha=0.7)

    # Set legend
    lbs = []
    for s in err_dict:
        lbs.append(plt.Line2D([], [], color="grey", marker=shapes[s], linestyle="None", label=s.split(".T")[0]))

    if len(lbs) <= 3:
        plt.legend(handles=lbs, loc="best", framealpha=0.5)
    else:
        plt.legend(handles=lbs, loc="center left", bbox_to_anchor=(1, 0.5), framealpha=0.5)

    # Add x-axis mechanism names
    ax = plt.gca()
    mech_locs = np.arange(1, len(mechs) + 1)
    sizes = [(1 - sdiff_dict[m]) * 100.0 for m in mechs]
    ax.set_xticks(mech_locs)
    ax.set_xticklabels([f"{s:.1f}" if s > 0.1 else f"{s:.3f}" for s in sizes])
    m_ax = ax.secondary_xaxis("top")
    m_ax.set_xticks(mech_locs)
    m_ax.set_xticklabels(mechs)
    m_ax.set_xlabel("Mechanisms")

    # Layout
    plt.xlabel("Relative Size (%)")
    plt.ylabel("Error")
    plt.title(title, pad=10)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()

    # Save figure
    plt.savefig(savname)
    plt.close()
