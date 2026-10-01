"""Thin line graphs of the obs and every model for one site and time window.

Run from the repo root:
    python scripts/05_line_plots.py AU-Cum 2013-12-11 2013-12-21
    python scripts/05_line_plots.py AU-Cum 2013-01-01 2013-12-31 --freq 1D
    python scripts/05_line_plots.py AU-How 2012-06-01 2012-10-31 --freq 1D --vars Qle GPP

Arguments:
    site            e.g. AU-Cum (reads data/harmonised/<site>.nc)
    start end       window, e.g. 2013-12-11 2013-12-21 (use "all all" for the whole record)
    --freq F        average to this step first: 1D (daily), 7D (weekly), 1MS (monthly).
                    Default is half-hourly, which only reads well for a few days to
                    a few weeks.
    --vars V ...    only these fluxes (Qle Qh Rnet Qg NEE GPP). Default: all six.

Writes outputs/lines/<site>_<start>to<end>.png and prints the path.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # find ccrc/ from scripts/
from ccrc.lineplots import line_plot


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("site")
    p.add_argument("start")
    p.add_argument("end")
    p.add_argument("--freq", default=None)
    p.add_argument("--vars", nargs="+", default=None)
    a = p.parse_args()

    start, end = (None, None) if a.start == "all" else (a.start, a.end)
    out = line_plot(Path("data/harmonised") / f"{a.site}.nc", start, end,
                    variables=a.vars, freq=a.freq)
    print(out)


if __name__ == "__main__":
    main()
