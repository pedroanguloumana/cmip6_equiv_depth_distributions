"""Project paths for use in other modules"""

from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SRC_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
FIGURE_DIR = PROJECT_ROOT / "figures"

CMIP6_ARCHIVE_DIR = DATA_DIR / "cmip6"
CMIP6_CATALOG_FILE = DATA_DIR / "cmip6_catalog.parquet"
CMIP6_HEADERS_FILE = DATA_DIR / "cmip6_headers.parquet"
