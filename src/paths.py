"""Project paths for use in other modules"""

from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SRC_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
FIGURE_DIR = PROJECT_ROOT / "figures"

CMIP6_ARCHIVE_DIR = DATA_DIR / "cmip6"
CMIP6_CATALOG_FILE = DATA_DIR / "cmip6_catalog.parquet"
CMIP6_HEADERS_FILE = DATA_DIR / "cmip6_headers.parquet"

# Configure this archive separately on each machine with data/noaa_olr.
NOAA_OLR_FILE = DATA_DIR / "noaa_olr" / "OLR-Daily_v02r00_s19790101_e20251231.nc"
NOAA_OLR_JOB_DIR = DATA_DIR / "noaa_olr_event_jobs"

# The regridded historical rlut archive is configured with a local symlink.
CMIP6_RLUT_DIR = DATA_DIR / "cmip6_rlut"
CMIP6_EVENT_DIR = DATA_DIR / "cmip6_event_catalogs"
CMIP6_EVENT_JOB_DIR = DATA_DIR / "cmip6_event_jobs"

# New NOAA-grid products are kept separate from the older IMERG-grid files.
# Configure the output symlink locally before processing climate fields.
CMIP6_NOAA_GRID_DIR = DATA_DIR / "cmip6_rlut_noaa1deg"
CMIP6_REGRID_JOB_DIR = DATA_DIR / "cmip6_regrid_jobs"

# data/archived points to the project's existing archived/ folder on Levante.
ARCHIVED_EVENT_DIR = DATA_DIR / "archived" / "event_catalogs"
