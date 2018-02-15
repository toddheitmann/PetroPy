from .datasets import log_data
from .download import ul_lands_download, kgs_download, create_log_inventory_table
from .electrofacies import electrofacies
from .graphs import LogViewer
from .log import Log

__version__ = '0.1.6'

def version():
    print(__version__)
