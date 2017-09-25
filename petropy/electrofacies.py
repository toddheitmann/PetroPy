# -*- coding: utf-8 -*-
"""
Electrofacies is a model to calculate numerical facies from log data. It uses sckit-learn for
standardization and clustering.

"""
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans

from log import Parameter

def electrofacies(logs, curves, formations, n_clusters, log_scale = [], n_components = 0.85):

    df = pd.DataFrame()

    for log in logs:

        if log.uwi == '':
            raise ValueError('UWI required for log identification.')

        log.curve_df['UWI'] = log.uwi

        for formation in formations:
            top = log.tops[formation]
            bottom = log.next_formation_depth(formation)
            df = df.append(log.curve_df[(log.curve_df.DEPTH >= top) & (log.curve_df.DEPTH < bottom)])

        log.curve_df.drop(['UWI'], 1, inplace = True)

    for s in log_scale:
        df[s] = np.log(df[s])

    not_null_rows = pd.notnull(df[curves]).any(axis = 1)

    X = StandardScaler().fit_transform(df.loc[not_null_rows, curves])

    pc = PCA(n_components = n_components).fit(X)

    components = pd.DataFrame(data = pc.transform(X), index = df[not_null_rows].index)
    components.columns = ['PC%i' % x for x in range(1, pc.n_components_ + 1)]
    components['UWI'] = df.loc[not_null_rows, 'UWI']

    size = len(components) // 20
    if size > 10000:
        size = 10000
    elif size < 100:
        size = 100

    facies = MiniBatchKMeans(n_clusters = n_clusters, batch_size = size).fit_predict(components)

    df.loc[not_null_rows, 'FACIES'] = facies

    for log in logs:

        uwi = log.uwi

        for v, vector in enumerate(pc.components_):
            v += 1
            column = 'PC%i' % v

            ### add eigenvector data to header ###

            pc_parameter = Parameter(name = column, unit = '', des = 'Principal Component from electrofacies')
            log.curve_parameters.append(pc_parameter)
            log.curve_values.append(column)
            log.curve_df = log.curve_df.join(components[components.UWI == uwi][column])

        facies_parameter = Parameter(name = 'FACIES', unit = '', des = 'Electrofacies')
        log.curve_parameters.append(facies_parameter)
        log.curve_values.append('FACIES')
        log.curve_df = log.curve_df.join(df[df.UWI == uwi]['FACIES'])

    return logs
