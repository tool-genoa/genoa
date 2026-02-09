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
GENOA v3: Package entry point.
"""

from .genoa import run_genoa


def main():
    """Entry point for GENOA v3."""

    run_genoa()


if __name__ == "__main__":
    main()
