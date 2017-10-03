# -*- coding: utf-8 -*-
"""
Download

This module downloads files from different public datasets. Each
function downloads the specific dataset to parse and unzip.

"""


import os
import sys
from ftplib import FTP
from zipfile import ZipFile
from io import BytesIO

if sys.version_info[0] < 3:
    from urllib2 import urlopen
else:
    from urllib.request import urlopen


def ul_lands_download(save_dir = None):
    """
    Downloads las files from University Lands Texas

    This function downloads files from the university lands ftp
    website located at publiftp.utlands.utsystem.edu.

    The bulk of this script is provided courtesy of Jon Reynolds,
    with Glacier Geosciences at:

    http://www.glaciergeosciences.com/


    Parameters
    ----------
    save_dir : str (default None)
        path to directory to save data. defaults to data folder
        within petropy

    Returns
    -------
    bool
        True if successful

    Examples
    --------
    >>> import petropy as ptr
    >>> ptr.ul_lands_download()

    >>> import petropy as ptr
    >>> p = 'path/to/my/folder/'
    >>> ptr.ul_lands_download(p)

    """

    ftp_site = 'publicftp.utlands.utsystem.edu'

    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(__file__), '..',
                                'data', 'ul')

    ftp = FTP(ftp_site)
    ftp.login()
    root_dir = 'ScannedLogs'

    for direct in ftp.nlst(root_dir):
        for subdir in ftp.nlst(os.path.join(root_dir, direct)):
            files = ftp.nlst(os.path.join(root_dir, direct, subdir))
            for f in files:
                if '.las' in f or '.LAS' in f:
                    src = os.path.join(root_dir, direct, subdir, f)
                    dest = os.path.join(save_dir, direct, subdir, f)
                    if os.path.exists(os.path.join(save_dir,
                                      direct, subdir)) == False:
                        os.makedirs(os.path.join(save_dir, direct,
                                                 subdir))
                    if os.path.exists(dest) == False:
                        with open(dest, 'wb') as las_file:
                            ftp.retrbinary('RETR ' + src,
                                           las_file.write)
    ftp.quit()

    return True

def kgs_download(save_dir = None):
    """
    Downloads las files from Kansas Geologic Society

    This function downloads files from the Kansas Geologic Society.
    These are zip files inside zip files, so the function parses out
    all las files and saves them in the folder data/kgs.


    Parameters
    ----------
    save_dir : str (default None)
        path to directory to save data. defaults to data folder within
        petropy

    Returns
    -------
    bool
        True if successful

    Examples
    --------
    >>> import petropy as ptr
    >>> ptr.kgs_download()

    >>> import petropy as ptr
    >>> p = 'path/to/my/folder/'
    >>> ptr.kgs_download(p)

    """

    urls = [
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2016.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2015.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2014.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2013.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2012.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2006_2011.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2001_2005.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/1999.zip'
    ]

    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(__file__), '..',
                                'data', 'kgs')

    for url in urls:

        year_dir = os.path.join(save_dir,
                                url.split('/')[-1].split('.')[0])
        if not os.path.isdir(year_dir):
            os.makedirs(year_dir)

        url = urlopen(url)
        zip_file = ZipFile(BytesIO(url.read()))

        zip_file.extractall(save_dir)
        for zip_name in zip_file.namelist():
            zip_path = os.path.join(save_dir, zip_name)
            las_zip = ZipFile(zip_path)
            las_zip.extractall(year_dir)
            las_zip.close()
            os.remove(zip_path)

    return True
