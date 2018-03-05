#!/usr/bin/env python
# encoding: utf-8
#
# @Author: José Sánchez-Gallego
# @Date: Feb 22, 2018
# @Filename: guider.py
# @License: BSD 3-Clause
# @Copyright: José Sánchez-Gallego


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import functools

import numpy as np
import pandas as pd
import sqlalchemy

from tqdm import tqdm

from ..utils.fit_data import ModelFit, TransRotModel, TransRotScaleModel
from ..utils.umeyama import umeyama


__ALL__ = ['db_to_hdf', 'get_position_error', 'umeyama_fit',
           'guider_fit', 'transrotscale_fit', 'transrot_fit',
           'fit_data']


# GProbe bits
GOOD = 0x00        # A happy, working, in-use probe has all bits set to 0.
BROKEN = 0x01      # A known broken probe, labeled as such in plPlugMap
NOSTAR = 0x02      # A probe with no star in plPlugMap (e.g. tritium)
DISABLED = 0x04    # A probe that has been disabled by the observers
ABOVEFOCUS = 0x08  # A fiber that is above the focal plane
BELOWFOCUS = 0x10  # A fiber that is below the focal plane
TOOFAINT = 0x20    # The observed star is too faint to be reliably used for guiding
UNKNOWN = 0xff     # This shouldn't ever happen

PLUGPLATESCALE = {'APO': 217.7358, 'LCO': 330.275}


def db_to_hdf(db_file, min_mjd=None):
    """Converts a guider database to a Pandas table."""

    engine = sqlalchemy.create_engine(f'sqlite:///{db_file}')

    header = pd.read_sql_table('header', engine)
    frame = pd.read_sql_table('frame', engine)

    bintable = pd.read_sql_table('bintable', engine)
    bintable = bintable.drop(columns=['pk']).set_index(['frame_pk', 'fiberid'])

    frame_header = pd.merge(frame, header, left_on=['pk'], right_on=['frame_pk'])
    frame_header = frame_header.drop(columns=['pk_x', 'pk_y', 'header_blob']).set_index('frame_pk')

    if min_mjd:
        frame_header = frame_header.loc[frame_header.mjd >= min_mjd]

    bintable = bintable.reindex(frame_header.index, level=0)

    dateobs = pd.to_datetime(frame_header.date_obs)
    frame_header['dateobs'] = dateobs

    basename = db_file.split('.')[0]
    hdf_name = basename + '.h5'

    frame_header.to_hdf(hdf_name, 'header', complib='zlib', mode='w')
    bintable.to_hdf(hdf_name, 'bintable', complib='zlib')

    return


def get_position_error(X, Y, trans, rot, scale):
    """Position error between two sets of points given a transformation."""

    theta = np.deg2rad(rot)
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta), np.cos(theta)]])

    pos_error = Y - (((scale * rot_matrix) @ X.T).T + trans)
    n_points = X.shape[0]

    return np.sqrt((pos_error ** 2).sum(axis=1).sum()) / n_points


def guider_fit(fibres, observatory='APO'):
    """Fits a set of fibres using the guider algorithm.

    Parameters:
        fibres (`~pandas.DataFrame`):
            A `~pandas.DataFrame` object with the fibres to fit. Must contain
            columns for ``dRA, dDec, xFocal, yFocal, dx, dy``.

    Returns:
        series (`~pandas.Series`):
            A `~pandas.Series` with the best fit for the measured errors,
            as well as the position error. All errors are in units of arcsec,
            except `mScale` which is the relative scale change.

    """

    dRA = fibres.dRA
    dDec = fibres.dDec
    raCenter = fibres.xFocal
    decCenter = fibres.yFocal

    arcsec_mm = 3600 / PLUGPLATESCALE[observatory]
    guideRMS = np.sqrt((fibres.dx**2 + fibres.dy**2).sum() / len(fibres)) * arcsec_mm

    bb = np.zeros(3, dtype=np.float64)
    AA = np.zeros((3, 3), dtype=np.float64)

    bb[0] = dRA.sum()
    bb[1] = dDec.sum()
    bb[2] = (raCenter * dDec - decCenter * dRA).sum()

    AA[0, 0] = len(fibres)
    AA[0, 1] = 0
    AA[0, 2] = -decCenter.sum()

    AA[1, 0] = 0
    AA[1, 1] = len(fibres)
    AA[1, 2] = raCenter.sum()

    AA[2, 2] = (raCenter * raCenter + decCenter * decCenter).sum()

    AA[2, 0] = AA[0, 2]
    AA[2, 1] = AA[1, 2]

    b3 = (raCenter * dRA + decCenter * dDec).sum()

    xx = np.linalg.solve(AA, bb)

    mRA = xx[0] / PLUGPLATESCALE[observatory]
    mDec = xx[1] / PLUGPLATESCALE[observatory]
    mRot = -np.rad2deg(xx[2])

    mScale = b3 / AA[2, 2] + 1.

    p0 = fibres[['xCenter', 'yCenter']].as_matrix()
    p1 = p0 - fibres[['dRA', 'dDec']].as_matrix()

    pos_error = get_position_error(p0, p1, [xx[0], xx[1]], -mRot, mScale)
    pos_error /= PLUGPLATESCALE[observatory]

    n_points = len(fibres)

    return pd.Series({'mRA': mRA, 'mDec': mDec, 'mRot': mRot,
                      'mScale': mScale, 'guideRMS': guideRMS,
                      'pos_error': pos_error, 'n_points': n_points})


def umeyama_fit(fibres, observatory='APO'):
    """Fits a set of fibres using the Umeyama (1991) algorithm.

    Parameters:
        fibres (`~pandas.DataFrame`):
            A `~pandas.DataFrame` object with the fibres to fit. Must contain
            columns for ``dRA, dDec, xFocal, yFocal``.

    Returns:
        series (`~pandas.Series`):
            A `~pandas.Series` with the best fit for the measured errors,
            as well as the position error. All errors are in units of arcsec,
            except `mScale` which is the relative scale change.

    """

    # xFocal/yFocal have the same direction as RA/Dec at APO and oposite at LCO.
    focal_to_sky = 1 if observatory == 'APO' else -1

    p0 = fibres[['xCenter', 'yCenter']].as_matrix().T
    p1 = p0 + focal_to_sky * fibres[['dRA', 'dDec']].as_matrix().T

    try:
        cc, rot, tt = umeyama(p0, p1)
    except ValueError:
        return pd.Series({'mRA': np.nan, 'mDec': np.nan, 'mRot': np.nan,
                          'mScale': np.nan, 'pos_error': np.nan})

    mRA = tt[0] / PLUGPLATESCALE[observatory]
    mDec = tt[1] / PLUGPLATESCALE[observatory]
    mRot = -np.rad2deg(np.arctan2(rot[1, 0], rot[0, 0]))
    mScale = cc

    pos_error = get_position_error(p0.T, p1.T, [tt[0], tt[1]], -mRot, mScale)
    pos_error /= PLUGPLATESCALE[observatory]

    n_points = len(fibres)

    return pd.Series({'mRA': mRA, 'mDec': mDec, 'mRot': mRot,
                      'mScale': mScale, 'pos_error': pos_error,
                      'n_points': n_points})


def transrotscale_fit(fibres, observatory='APO'):
    """Fits a set of fibres using `.TransRotScaleModel`.

    Parameters:
        fibres (`~pandas.DataFrame`):
            A `~pandas.DataFrame` object with the fibres to fit. Must contain
            columns for ``dRA, dDec, xFocal, yFocal``.

    Returns:
        series (`~pandas.Series`):
            A `~pandas.Series` with the best fit for the measured errors,
            as well as the position error. All errors are in units of arcsec,
            except `mScale` which is the relative scale change.

    """

    # xFocal/yFocal have the same direction as RA/Dec at APO and oposite at LCO.
    focal_to_sky = 1 if observatory == 'APO' else -1

    p0 = fibres[['xCenter', 'yCenter']].as_matrix()
    p1 = p0 + focal_to_sky * fibres[['dRA', 'dDec']].as_matrix()

    mfit = ModelFit(TransRotScaleModel(), p1, p0)

    [mRA, mDec], mRot, mScale = mfit.model.getTransRotScale()

    pos_error = get_position_error(p0, p1, [mRA, mDec], mRot, mScale)
    pos_error /= PLUGPLATESCALE[observatory]

    mRA /= PLUGPLATESCALE[observatory]
    mDec /= PLUGPLATESCALE[observatory]
    mRot = -mRot

    n_points = len(fibres)

    return pd.Series({'mRA': mRA, 'mDec': mDec, 'mRot': mRot,
                      'mScale': mScale, 'pos_error': pos_error,
                      'n_points': n_points})


def transrot_fit(fibres, observatory='APO'):
    """Fits a set of fibres using `.TransRotModel`.

    This fitter does not try to fit scale iteratively. Instead, the scale
    is calculated in the same way as in `.guider_fit`.

    Parameters:
        fibres (`~pandas.DataFrame`):
            A `~pandas.DataFrame` object with the fibres to fit. Must contain
            columns for ``dRA, dDec, xFocal, yFocal``.

    Returns:
        series (`~pandas.Series`):
            A `~pandas.Series` with the best fit for the measured errors,
            as well as the position error. All errors are in units of arcsec,
            except `mScale` which is the relative scale change.

    """

    dRA = fibres.dRA
    dDec = fibres.dDec
    raCenter = fibres.xFocal
    decCenter = fibres.yFocal

    # xFocal/yFocal have the same direction as RA/Dec at APO and oposite at LCO.
    focal_to_sky = 1 if observatory == 'APO' else -1

    p0 = fibres[['xCenter', 'yCenter']].as_matrix()
    p1 = p0 + focal_to_sky * fibres[['dRA', 'dDec']].as_matrix()

    mfit = ModelFit(TransRotModel(), p1, p0)

    [mRA, mDec], mRot = mfit.model.getTransRot()

    AA22 = (raCenter * raCenter + decCenter * decCenter).sum()
    b3 = (raCenter * dRA + decCenter * dDec).sum()
    mScale = b3 / AA22 + 1.

    pos_error = get_position_error(p0, p1, [mRA, mDec], mRot, mScale)
    pos_error /= PLUGPLATESCALE[observatory]

    mRA /= PLUGPLATESCALE[observatory]
    mDec /= PLUGPLATESCALE[observatory]
    mRot = -mRot

    n_points = len(fibres)

    return pd.Series({'mRA': mRA, 'mDec': mDec, 'mRot': mRot,
                      'mScale': mScale, 'pos_error': pos_error,
                      'n_points': n_points})


def fit_data(h5_file, min_fibres=5, min_mjd=None, observatory='LCO'):
    """Runs all the fitters in sequence and saves the data to a file.

    The parameter ``min_fibres`` defines how many good fibres (i.e., in which
    the position of the centre of the fibre and the measured star are not NaN)
    are needed for a frame to be included.

    """

    bintable = pd.read_hdf(h5_file, 'bintable')

    if min_mjd is not None:
        frames = pd.read_hdf(h5_file, 'header')
        bintable = bintable.loc[frames.mjd > min_mjd]

    invalid_fibre = BROKEN | NOSTAR | DISABLED | TOOFAINT

    filt = bintable.gprobe_exists == 1
    filt &= bintable.gprobebits & invalid_fibre == 0
    filt &= ~bintable['dRA'].isna() & ~bintable['dDec'].isna()

    bintable_filtered = bintable[filt]
    bintable_min_size = bintable_filtered.groupby(level=0).filter(lambda gg: len(gg) > min_fibres)
    grouped_min_size = bintable_min_size.groupby(level=0)

    tqdm.pandas()

    print(f'Observatory: {observatory}')
    print('# groups: {}'.format(len(grouped_min_size)))

    print('Running guider algorithm ...')

    guider_fit_obs = functools.partial(guider_fit, observatory=observatory)
    guider_df = grouped_min_size.progress_apply(guider_fit_obs)
    guider_df.to_hdf(h5_file, 'guider_fit')

    print('Running Umeyama algorithm ...')

    umeyama_fit_obs = functools.partial(umeyama_fit, observatory=observatory)
    umeyama_df = grouped_min_size.progress_apply(umeyama_fit_obs)
    umeyama_df.to_hdf(h5_file, 'umeyama_fit')

    print('Running TransRotScaleModel algorithm ...')

    transrotscale_fit_obs = functools.partial(transrotscale_fit, observatory=observatory)
    transrotscale_df = grouped_min_size.progress_apply(transrotscale_fit_obs)
    transrotscale_df.to_hdf(h5_file, 'transrotscale_fit')

    print('Running TransRotModel algorithm ...')

    transrot_fit_obs = functools.partial(transrot_fit, observatory=observatory)
    transrot_df = grouped_min_size.progress_apply(transrot_fit_obs)
    transrot_df.to_hdf(h5_file, 'transrot_fit')

    return
