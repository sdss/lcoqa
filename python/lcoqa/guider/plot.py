#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Mar 3, 2018
# @Filename: plot.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import itertools

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from scipy import stats


def add_series(df, header):

    dateobs = pd.to_datetime(header.reindex(df.index)['date_obs'])
    timediff = dateobs.diff().dt.seconds

    chunks = np.split(np.arange(len(timediff)), np.where(timediff > 100)[0])
    series_idx = [[ii + 1] * len(qq) for ii, qq in enumerate(chunks)]
    series_idx = list(itertools.chain.from_iterable(series_idx))

    df['serie'] = series_idx
    df['dateobs'] = dateobs

    return df


def plot_residuals(h5_file, table, compare_table='guider'):
    """Plots the residuals between two sets of fit data.

    Parameters:
        h5_file (str):
            An HDF file containing the data to plot.
        table (str):
            The name of the table containing the first set of data.
        compare_table(str):
            The name of the table containing the second set of data.

    """

    if not table.endswith('_fit'):
        table += '_fit'

    if not compare_table.endswith('_fit'):
        compare_table += '_fit'

    data = pd.read_hdf(h5_file, table)
    header = pd.read_hdf(h5_file, 'header')

    compare_data = pd.read_hdf(h5_file, compare_table).loc[data.index]

    diff = data - compare_data
    if 'mScale' in data and 'mScale' in compare_data:
        data['mScale'] = data.mScale / compare_data.mScale

    diff = diff.dropna(how='all', axis=1)
    diff = diff.dropna().drop(columns=['n_points'])

    diff = diff[(np.abs(stats.zscore(diff)) < 10).all(axis=1)]

    diff.loc[:, ['mRA', 'mDec', 'mRot', 'pos_error']] *= 3600.

    diff = diff[['mRA', 'mDec', 'mRot'] + (['mScale'] if 'mScale' in diff else [])]
    diff = add_series(diff, header)

    melted = diff.melt(['dateobs', 'serie'], var_name='axis', value_name='value')

    fg = sns.FacetGrid(melted, col='axis', sharey=False, col_wrap=2, hue='serie')
    fg.map(plt.plot, 'dateobs', 'value')

    return diff


def plot_pos_error(data, x, mode='kde', **kwargs):

    assert mode in ['kde', 'scatter']

    if mode == 'kde':
        sns.set_style('white')
        ax = sns.jointplot(x=x, y='pos_error', data=data, kind='kde', space=0, color='g',
                           stat_func=None, n_levels=50, **kwargs)
    elif mode == 'scatter':
        ax = data.plot.scatter(x, 'pos_error', s=0.5, alpha=0.4)

    return ax
