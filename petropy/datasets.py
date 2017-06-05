# -*- coding: utf-8 -*-
"""
Datasets is a way to retrieve included logs with petropy

"""

import os
from las import Las

def log_data(source):
    """
    retrieves log data for a formation

    Parameters
    ----------
    source : str {'WFMP'}
        source location for log data

    Returns
    -------
    log : Log
        Log object of data source

    Raises
    ------
    ValueError
        If source is not in dictionary key

    """

    file_dir = os.path.dirname(__file__)

    paths = {
        'WFMP': os.path.join(file_dir, '..', 'data', 'sample', '42303347740000.las')
    }

    if source in paths:
        las_path = paths[source]
    else:
        raise ValueError('%s is not valid source' % source)

    log = Las(las_path)

    return log
