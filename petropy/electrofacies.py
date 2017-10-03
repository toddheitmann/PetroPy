# -*- coding: utf-8 -*-
"""
Electrofacies is a model to calculate numerical facies from log data.
It uses sckit-learn for standardization and clustering.

"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

def electrofacies(logs, formations, curves, n_clusters, log_scale = [],
                  n_components = 0.85, curve_name = 'FACIES'):
    """
    Electrofacies function to group intervals by rock type. Also
    referred to as heterogenous rock analysis.

    Parameters
    ----------
    logs : list of :class:`ptr.Log` objects
        List of Log objects
    formations: list of formation names
        List of strings containg formation names which should be
        previously loaded into Log objects
    curves : list of curve names
        List of strings containing curve names as inputs in the
        electrofacies calculations
    n_clusters : int
        Number of clusters to group intervals. Number of electrofacies.
    log_scale : list of curve names
        List of string containing curve names which are preprocessed
        on a log scale. For example, deep resistivity separates better
        on a log scale, and is graph logarithmically when viewing data
        in a log viewer.
    n_components : int, float, None or string (default 0.85)
        Number of principal components to keep. If value is less than
        one, the number of principal components be the number required
        to exceed the explained variance.
    curve_name : str (default 'FACIES')
        Name of the new electrofacies curve.

    Examples
    --------
    >>> # loads sample Wolfcamp calculates electrofacies for that well
    >>> import petropy as ptr
    # reads sample Wolfcamp Log from las file
    >>> log = ptr.log_data('WFMP')
    >>> logs = [log]
    >>> f = ['WFMPA', 'WFMPB', 'WFMPC']
    >>> c = ['GR_N', 'RESDEEP_N', 'RHOB_N', 'NPHI_N', 'PE_N']
    >>> scale = ['RESDEEP_N']
    >>> logs = electrofacies(logs, f, c, 8, log_scale = scale)

    >>> import petropy as ptr
    # loads logs from a list of paths and
    # calculates electrofacies across the wells
    #
    # defin file_paths for las files to analyze
    >>> file_paths = ['path/to/log1.las', 'path/to/log2.las',
    ... 'path/to/log3.las', 'path/to/log4.las']
    # create list of Log objects
    >>> logs = [ptr.Log(x) for x in file_paths]
    # define csv with tops for all wells
    >>> tops_csv = 'path/to/tops.csv'
    # add formation tops to wells
    >>> for log in logs:
    ...     log.tops_from_csv(tops_csv)
    # define list of formation tops. If single formation, f = ['FORM']
    >>> f = ['FORM1', 'FORM2']
    # list of curves to use for classification
    >>> c = ['GR_N', 'RESDEEP_N', 'RHOB_N', 'NPHI_N', 'PE_N']
    >>> scale = ['RESDEEP_N']
    # run electrofacies across logs in list
    >>> logs = electrofacies(logs, f, c, 8, log_scale = scale)
    # save las in renamed file
    >>> for i, log in enumerate(logs):
    ...     new_file_name = file_paths[i].split('.')[0]+'_with_HRA.las'
    ...     log.write(new_file_name)

    """

    df = pd.DataFrame()

    for log in logs:

        if log.well['UWI'] is None:
            raise ValueError('UWI required for log identification.')

        log_df = log.df()
        log_df['UWI'] = log.well['UWI'].value
        log_df['DEPTH_INDEX'] = np.arange(0, len(log[0]))

        for formation in formations:
            top = log.tops[formation]
            bottom = log.next_formation_depth(formation)
            depth_index = np.intersect1d(np.where(log[0] >= top)[0],
                                         np.where(log[0] < bottom)[0])
            df = df.append(log_df.iloc[depth_index])

    for s in log_scale:
        df[s] = np.log(df[s])

    not_null_rows = pd.notnull(df[curves]).any(axis = 1)

    X = StandardScaler().fit_transform(df.loc[not_null_rows, curves])

    pc = PCA(n_components = n_components).fit(X)

    components = pd.DataFrame(data = pc.transform(X),
                              index = df[not_null_rows].index)

    components.columns = \
                   ['PC%i' % x for x in range(1, pc.n_components_ + 1)]

    components['UWI'] = df.loc[not_null_rows, 'UWI']
    components['DEPTH_INDEX'] = df.loc[not_null_rows, 'DEPTH_INDEX']

    size = len(components) // 20
    if size > 10000:
        size = 10000
    elif size < 100:
        size = 100

    df.loc[not_null_rows, curve_name] = \
                            MiniBatchKMeans(n_clusters = n_clusters,
                            batch_size = size).fit_predict(components)

    for log in logs:

        uwi = log.well['UWI'].value

        for v, vector in enumerate(pc.components_):
            v += 1
            pc_curve = 'PC%i' % v

            ### add eigenvector data to header ###

            if pc_curve in log.keys():
                data = log[pc_curve]
                depth_index = components.loc[components.UWI == uwi,
                                             'DEPTH_INDEX']
                data[depth_index] = \
                          np.copy(components.loc[components.UWI == uwi,
                                                pc_curve])
            else:
                data = np.empty(len(log[0]))
                data[:] = np.nan
                depth_index = components.loc[components.UWI == uwi,
                                             'DEPTH_INDEX']
                data[depth_index] = \
                          np.copy(components.loc[components.UWI == uwi,
                          pc_curve])

                log.add_curve(pc_curve, np.copy(data),
                                       descr = 'Pricipal Component %i \
                                       from electrofacies' % v)

        if curve_name in log.keys():
            data = log[curve_name]
            depth_index = df.loc[df.UWI == uwi, 'DEPTH_INDEX']
            data[depth_index] = df.loc[df.UWI == uwi, curve_name]
        else:
            data = np.empty(len(log[0]))
            data[:] = np.nan
            depth_index = components.loc[components.UWI == uwi,
                                         'DEPTH_INDEX']

            data[depth_index] = \
               np.copy(components.loc[components.UWI == uwi, pc_curve])

            log.add_curve(curve_name, np.copy(data),
                          descr = 'Electrofacies')

    return logs
