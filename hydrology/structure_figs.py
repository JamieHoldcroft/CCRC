# must run from root

import os
from pathlib import Path
import shutil

here = os.getcwd()

sites = sorted(Path(f"{here}/data/harmonised/").glob("*.nc"))
ws_timeseries = sorted(Path("hydrology/figures/ws_timeseries/").glob("*.png"))

destination = Path("hydrology/sites/")         # every figure is saved here
destination.mkdir(exist_ok=True)

for source in sites:
    site_folder = destination / source
    site_folder.mkdir(exist_ok=True)

    
