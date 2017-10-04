from .datasets import log_data
from .download import ul_lands_download, kgs_download
from .electrofacies import electrofacies
from .graphs import LogViewer
from .log import Log

__version__ = '0.1.1'

def version():
    print(__version__)
