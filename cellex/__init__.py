name = "CELLEX"

### GLOBAL CONSTANTS
# These must be defined first, as they are required in sub-packages.
ES_METRICS = ["det", "ep", "ges", "nsi"]


__author__ = ", ".join([
    "Tobias O. Stannius",
    ])

# Load CELLEX modules
from . import metrics
from . import plotting
from . import preprocessing
from . import utils
from .esobject import ESObject
from .summarydata import SummaryData
