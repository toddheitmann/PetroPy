# -*- coding: utf-8 -*-
"""
Datasets is a way to retrieve included logs with petropy. It currently
supports reading a sample log from the Permain Basin in Reagan County.

"""

import os
from .log import Log

def log_data(source):
    """
    retrieves log data for a formation

    Parameters
    ----------
    source : str {'WFMP'}
        source location for log data

    Returns
    -------
    :class:`petropy.Log`
        Log object of data source

    Raises
    ------
    ValueError
        If source is not in dictionary key

    Example
    -------
    >>> import petropy as ptr
    # reads sample Wolfcamp Log from las file
    >>> log = ptr.log_data('WFMP')

    """

    file_dir = os.path.dirname(__file__)

    paths = {
        'WFMP': os.path.join(file_dir, 'data', '42303347740000.las')
    }

    p = os.path.join(file_dir, 'data', 'tops.csv')

    if source in paths:
        las_path = paths[source]
    else:
        raise ValueError('%s is not valid source' % source)

    log = Log(las_path)
    log.tops_from_csv()

    return log
