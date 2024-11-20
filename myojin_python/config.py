from pathlib import Path


# Paths
PROJ_ROOT = Path(__file__).resolve().parents[1]
# logger.info(f"PROJ_ROOT path is: {PROJ_ROOT}")

# PATH TO MAIN CONDUIT SWEEP OUTPUT FILES
DATA_DIR = Path("/Users/crrowell/Kahuna/data/myojin/mainSweep2")

SECONDARY_DATA_DIR = PROJ_ROOT / "data"

OUTCOME_CSV = PROJ_ROOT / "references/OutComeCodeTable.csv"
SIMPLIFIED_OUTCOME_CSV = PROJ_ROOT / "references/simplifiedOutComeCodeTable.csv"


# RAW_DATA_DIR = SECONDARY_DATA_DIR / "raw"
# INTERIM_DATA_DIR = SECONDARY_DATA_DIR / "interim"
# PROCESSED_DATA_DIR = SECONDARY_DATA_DIR / "processed"
# EXTERNAL_DATA_DIR = SECONDARY_DATA_DIR / "external"


REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = PROJ_ROOT / "figures"

# # If tqdm is installed, configure loguru with tqdm.write
# # https://github.com/Delgan/loguru/issues/135
# try:
#     from tqdm import tqdm

#     logger.remove(0)
#     logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True)
# except ModuleNotFoundError:
#     pass
