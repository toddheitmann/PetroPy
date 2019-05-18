# -*- coding: utf-8 -*-
"""
Download

This module downloads files from different public datasets. Each
function downloads the specific dataset to parse and unzip.

"""


import os
import sys
import time
import fnmatch
from ftplib import FTP
from zipfile import ZipFile
from io import BytesIO
import pandas as pd

if sys.version_info[0] < 3:
    from urllib2 import urlopen
else:
    from urllib.request import urlopen

from .log import Log


def ul_lands_download(save_dir = None):
    """
    Downloads las files from University Lands Texas

    This function downloads files from the university lands ftp
    website located at publiftp.utlands.utsystem.edu.  It inventories
    readable logs into a csv file containing header data in the
    save_dir. This inventory data is also returned as a DataFrame.

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
    df : :class:`pandas.DataFrame`
        DataFrame of header data for all logs downloaded and read.

    Examples
    --------
    >>> import petropy as ptr
    >>> ptr.ul_lands_download()

    >>> import petropy as ptr
    >>> p = r'path/to/my/folder/'
    >>> ptr.ul_lands_download(p)

    Note
    ----
    Function takes approximately twelve hours to scan ftp site,
    download, and inventory 30 GB of log data. YMMV depending on
    internet and processor speed.

    """
    start_time = time.time()
    ftp_site = 'publicftp.utlands.utsystem.edu'

    if save_dir is None:
        save_dir = os.path.join(os.path.dirname(__file__), 'data',
                                'ul')

    ftp = FTP(ftp_site)
    ftp.login()
    root_dir = 'ScannedLogs'

    ftp_source_files = []
    for direct in ftp.nlst(root_dir):
        print(direct)
        for subdir in ftp.nlst(os.path.join(root_dir, direct)):
            files = ftp.nlst(os.path.join(root_dir, direct, subdir))
            ftp_source_files.append((direct, subdir, files))
    ftp.quit()

    time.sleep(2)

    total_time = (time.time() - start_time) / 60.0
    print('Find Files Time: %.2f minutes' % total_time)

    process_time = time.time()

    paths = []
    for direct, subdir, files in ftp_source_files:
        las_files = [f for f in files if '.las' in f.lower()]
        for f in las_files:
            src = os.path.join(root_dir, direct, subdir, f)
            dest = os.path.join(save_dir, direct, subdir, f)
            paths.append((src, dest))

    chunk_time = time.time()
    # chunck paths into groups of 1000 las files
    print('Number of las files: %i' % len(paths))
    total_time = (time.time() - process_time) / 60.0
    print('Prcoess Filenames Time: %.2f minutes' % total_time)

    n = 400
    for i in range(0, len(paths), n):
        ftp = FTP(ftp_site)
        ftp.login()
        files = paths[i:i + n]
        for src, dest in files:
            direct = os.path.dirname(dest)
            if os.path.exists(direct) == False:
                os.makedirs(direct)
            if os.path.exists(dest) == False:
                with open(dest, 'wb') as las_file:
                    ftp.retrbinary('RETR ' + src, las_file.write)
        ftp.quit()
        time.sleep(5)

    total_time = (time.time() - chunk_time) / 60.0
    print('Chunk Files Time: %.2f minutes' % total_time)

    inventory_time = time.time()
    df = create_log_inventory_table(save_dir)
    total_time = (time.time() - chunk_time) / 60.0
    print('Inventory Time: %.2f minutes' % total_time)

    return df


def kgs_download(save_dir = None):
    """
    Downloads las files from Kansas Geologic Society

    This function downloads files from the Kansas Geologic Society.
    These are zip files inside zip files, so the function parses out
    all las files and saves them in the folder input save_dir or with
    package data in the folder data/kgs. It inventories readable logs
    into a csv file containing header data in the save_dir. This
    inventory data is also returned as a DataFrame.

    Parameters
    ----------
    save_dir : str (default None)
        path to directory to save data. defaults to data folder within
        petropy

    Returns
    -------
    df : :class:`pandas.DataFrame`
        DataFrame of header data for all logs downloaded and read.

    Examples
    --------
    >>> import petropy as ptr
    >>> ptr.kgs_download()

    >>> import petropy as ptr
    >>> p = r'path/to/my/folder/'
    >>> ptr.kgs_download(p)

    Note
    ----
    Function takes approximately one hour to download, unzip, and
    inventory 20GB of log data. YMMV depending on internet and
    processor speed.

    """

    urls = [
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2018.zip',
        'http://www.kgs.ku.edu/PRS/Scans/Log_Summary/2017.zip',
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
        save_dir = os.path.join(os.path.dirname(__file__), 'data',
                                'kgs')

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

    df = create_log_inventory_table(save_dir)

    return df

def create_log_inventory_table(save_dir, folder_copy = None):
    """
    Scans all folders and subfolders (recursive scan) for las files,
    and opens them as a :class:`petropy.Log` object. Extracts header
    data and curve names. Returns DataFrame of data after saving to a
    csv file in the save_dir folder.

    Parameters
    ----------
    save_dir : str
        path to folder for recusive scan

    Returns
    -------
    df : :class:`pandas.DataFrame`
        DataFrame of header data for all logs downloaded and read.

    Example
    -------
    >>> import petropy as ptr
    >>> p = r'path/to/folder/'
    >>> df = ptr.create_log_inventory_table(p)
    >>> # filter logs with triple-combo for processing
    >>> tc_df = df[df.GR_N == 'Y' & df.RESDEEP_N == 'Y' &
    ...            df.NPHI_N == 'Y' & df.RHOB_N == 'Y']
    >>> # print count of useable logs
    >>> print(len(tc_df))


    """

    if folder_copy is not None:
        if not os.path.isdir(folder_copy):
            os.makedirs(folder_copy)

    error_log = os.path.join(save_dir, 'error_log.txt')
    with open(error_log, 'w') as f:
        f.write('LAS INVENTORY ERROR LOG\n')
        f.write('-----------------------\n')

    log_data = []
    for root, dirnames, filenames in os.walk(save_dir):
        for filename in fnmatch.filter(filenames, '*.las'):
            try:
                las_path = os.path.join(root, filename)
                log = Log(las_path)

                if folder_copy is not None:

                    keys = ['WELL', 'UWI', 'API']
                    new_name = ''
                    for k in keys:
                        if k in log.well:
                            if len(log.well[k].value) > 0:
                                new_name = log.well[k].value + '.las'
                                break

                    if new_name == '':
                        new_name = os.path.basename(path)

                    new_name = os.path.join(folder_copy, new_name)
                    log.write(new_name)

            except:
                with open(error_log, 'a') as f:
                    f.write(las_path)
                    f.write('\n')
                continue

            data = {'PATH': las_path}

            for w in log.well:
                data[w.mnemonic] = w.value

            for c in log.curves:
                data[c.mnemonic] = 'Y'

            log_data.append(data)

    df = pd.DataFrame(log_data)
    log_inventory_path = os.path.join(save_dir, 'log_inventory.csv')
    df.to_csv(log_inventory_path, index = False)

    return df
