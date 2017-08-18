from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Process guide camera images into a single roll-up file with statistics.

import numpy as np
import os
import fitsio
import sdss_access
import pydl.pydlutils.yanny as yanny
import pydl.pydlutils.spheregroup as spheregroup

spath = sdss_access.path.Path()


def signal_dtype():
    return([('expno', np.dtype(np.int32), 1),
            ('mediansky', np.dtype(np.float32), 1),
            ('sn2ref', np.dtype(np.float32), 1),
            ('signalref', np.dtype(np.float32), 1)])


def signal_fibers_dtype():
    return([('fiberid', np.dtype(np.int32), 1),
            ('tmass_h', np.dtype(np.float32), 1),
            ('xfocal', np.dtype(np.float32), 1),
            ('yfocal', np.dtype(np.float32), 1),
            ('sn2', np.dtype(np.float32), 1),
            ('sn2_model', np.dtype(np.float32), 1),
            ('signal', np.dtype(np.float32), 1),
            ('signal_model', np.dtype(np.float32), 1)])


def mapfile(name=None, mjd=None):
    mdir = os.path.join(os.getenv('MAPPER_DATA_2S'), str(mjd))
    basename = 'plPlugMapM-{name}.par'.format(name=name)
    return(os.path.join(mdir, basename))


def fit_sn2(tmass_h=None, sn2=None):
    pivot = 12.1
    slope = - 0.5
    offset = np.log10(400.)
    x_fit = tmass_h - pivot
    y_fit = np.log10(sn2) - offset - x_fit * slope
    coeff = np.polyfit(x_fit, y_fit, 0)
    sn2ref = 10.**(offset + coeff)
    model = 10.**(coeff + offset + x_fit * slope)
    return (sn2ref, model)


def fit_signal(tmass_h=None, signal=None):
    pivot = 12.1
    slope = - 0.4
    offset = np.log10(400.)
    x_fit = tmass_h - pivot
    y_fit = np.log10(signal) - offset - x_fit * slope
    coeff = np.polyfit(x_fit, y_fit, 0)
    signalref = 10.**(offset + coeff)
    model = 10.**(coeff + offset + x_fit * slope)
    return (signalref, model)


def signal(exposure=None):
    """Evaluate signal for an exposure

    Parameters:
    -----------

    exposure : ndarray
      information from exposures-[mjd].fits file for a signal exposure
    telescope : str
      which telescope ('lco' or 'apo')

    Comments:
    --------

    """

    # Read in flux and errors
    flux = fitsio.read(exposure['expfile'], ext=1)
    err = fitsio.read(exposure['expfile'], ext=2)

    # Replace bad error values
    ibad = np.nonzero(err <= 0)
    err[ibad] = 1

    plugfile = mapfile(name=exposure['name'].strip(), mjd=exposure['mjd'])
    if(plugfile == ''):
        print(exposure['name'])
        print(exposure['mjd'])
        return None, None
    pmap = yanny.yanny(plugfile)

    # Not obvious what this fix does?
    try:
        ibad = np.nonzero(pmap['PLUGMAPOBJ']['fiberId'] < 1)[0]
        pmap['PLUGMAPOBJ']['fiberId'][ibad] = 1
    except KeyError:
        print("Bad plug file: {plugfile}".format(plugfile=plugfile))
        return (None, None)

    # Look up science and sky holes
    iobj = np.nonzero((pmap['PLUGMAPOBJ']['objType'] == b'STAR_BHB') |
                      (pmap['PLUGMAPOBJ']['objType'] == b'SPECTROPHOTO_STD'))[0]
    fiobj = 299 - (pmap['PLUGMAPOBJ']['fiberId'][iobj] - 1)
    isky = np.nonzero(pmap['PLUGMAPOBJ']['objType'] == b'SKY')[0]
    fisky = 299 - (pmap['PLUGMAPOBJ']['fiberId'][isky] - 1)

    # Calculate sky
    fsky = np.median(flux[fisky, :], 0)

    # Calculate median S/N and median signal
    sn = flux * 0.
    signal = flux * 0.
    for i in np.arange(flux.shape[0]):
        signal[i, :] = (flux[i, :] - fsky)
        sn[i, :] = (signal[i, :] / err[i, :])
    mediansn = np.median(sn, 1)
    mediansky = np.median(fsky)
    mediansignal = np.median(signal, 1)

    # Get more information from plateHolesSorted
    holefile = spath.full('plateHolesSorted', plateid=exposure['plateid'])
    holes = yanny.yanny(holefile)
    (m1, m2, d12) = spheregroup.spherematch(pmap['PLUGMAPOBJ']['ra'][iobj],
                                            pmap['PLUGMAPOBJ']['dec'][iobj],
                                            holes['STRUCT1']['target_ra'],
                                            holes['STRUCT1']['target_dec'],
                                            1. / 3600.)
    fiberid = pmap['PLUGMAPOBJ']['fiberId'][iobj[m1]]
    xf = holes['STRUCT1']['xfocal'][m2]
    yf = holes['STRUCT1']['yfocal'][m2]
    tmass_h = holes['STRUCT1']['tmass_h'][m2]
    sn2 = (mediansn[fiobj[m1]])**2
    signal = mediansignal[fiobj[m1]]

    ifit = np.where((tmass_h > 11.0) & (tmass_h < 12.5) &
                    (sn2 > 30.) & (sn2 < 3000.) &
                    (signal > 0.) & (signal < 10000.))[0]
    if(len(ifit) == 0):
        return(None, None)

    (sn2ref, sn2_model) = fit_sn2(tmass_h[ifit], sn2[ifit])
    # sn2ref_ql = fit_sn2_ql(tmass_h, sn2)
    (signalref, signal_model) = fit_signal(tmass_h[ifit], signal[ifit])

    signal_summary = np.zeros(1, dtype=signal_dtype())
    signal_summary['expno'] = exposure['expno']
    signal_summary['mediansky'] = mediansky
    signal_summary['sn2ref'] = sn2ref
    signal_summary['signalref'] = signalref

    signal_fibers = np.zeros(len(ifit), dtype=signal_fibers_dtype())
    signal_fibers['fiberid'] = fiberid[ifit]
    signal_fibers['tmass_h'] = tmass_h[ifit]
    signal_fibers['xfocal'] = xf[ifit]
    signal_fibers['yfocal'] = yf[ifit]
    signal_fibers['sn2'] = sn2[ifit]
    signal_fibers['sn2_model'] = sn2_model
    signal_fibers['signal'] = signal[ifit]
    signal_fibers['signal_model'] = signal_model

    return(signal_summary, signal_fibers)
