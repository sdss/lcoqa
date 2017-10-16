from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

sn2range = 0.25
signalrange = 0.125


class Expected(object):
    """Expected signal given seeing

    Attributes
    ----------
    signal : interpolate.interpolate.interp1d instance
      function returning expected signal given FWHM seeing, arcsec)

    """
    def __init__(self):
        coeffs = np.array([0.16691948, -0.12886663, -0.15449739, 0.30526956, -0.41347824, 0.41610774])
        deg = 5
        seeings = 0.5 + np.arange(20) * 0.1
        smodel = np.zeros(len(seeings))
        for indx in np.arange(deg + 1):
            smodel = smodel + coeffs[indx] * (seeings - 1.5)**(deg - indx)
        seeings = np.append([0.], seeings)
        smodel = np.append([1.], smodel)
        seeings = np.append(seeings, [5.])
        smodel = np.append(smodel, [0.])
        oseeings = (16000. / 7600.)**0.2 * seeings
        scale = 470.
        self.signal = interpolate.interp1d(oseeings, smodel * scale,
                                           fill_value='extrapolate')

expected = Expected()


def _title(exposure=None):
    title_template = "mjd={mjd}; plate={plate}; exp={expno}"
    title = title_template.format(mjd=exposure['mjd'],
                                  plate=exposure['plateid'],
                                  expno=exposure['expno'])
    return(title)


def seeing(outfile=None, gcam=None, exposure=None):
    """Plot seeing during an exposure

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    gcam : ndarray
      structured array of guider camera summary for MJD
    """
    indx = np.where((gcam['indx'] >= exposure['gstart']) &
                    (gcam['indx'] <= exposure['gend']))[0]
    plt.plot(gcam['indx'][indx], gcam['seeing'][indx], label='seeing')
    plt.plot(gcam['indx'][indx], gcam['fwhm_median'][indx],
             label='FWHM in focus')
    plt.xlim((gcam['indx'][indx].min() - 3, gcam['indx'][indx].max() + 3))
    plt.ylim((0.5, 2.2))
    plt.xlabel('Guide camera file number')
    plt.ylabel('Seeing (arcsec)')
    plt.title(_title(exposure=exposure))
    plt.legend()
    plt.savefig(outfile)
    plt.clf()


def guider_rms(outfile=None, gcam=None, exposure=None):
    """Plot guider RMS during an exposure

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    gcam : ndarray
      structured array of guider camera summary for MJD
    """
    indx = np.where((gcam['indx'] >= exposure['gstart']) &
                    (gcam['indx'] <= exposure['gend']))[0]
    plt.plot(gcam['indx'][indx], gcam['gdrms'][indx])
    plt.xlim((gcam['indx'][indx].min() - 3, gcam['indx'][indx].max() + 3))
    plt.ylim((0., 1.5))
    plt.xlabel('Guide camera file number')
    plt.ylabel('Guider RMS (arcsec)')
    plt.title(_title(exposure=exposure))
    plt.legend()
    plt.savefig(outfile)
    plt.clf()


def signal(outfile=None, exposure=None, signal=None,
           signal_fibers=None):
    """Plot signal vs. H magnitude

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    signal : ndarray
      structured array of signal information for this exposure
    signal_fibers : ndarray
      structured array of signal information for fibers in this exposure
    """
    cm = plt.cm.get_cmap('RdYlBu')
    isort = np.argsort(signal_fibers['tmass_h'])
    rf = np.sqrt(signal_fibers['xfocal']**2 + signal_fibers['yfocal']**2)
    sc = plt.scatter(signal_fibers['tmass_h'][isort],
                     np.log10(signal_fibers['signal'][isort]),
                     c=rf[isort], s=50, alpha=0.8, vmin=0.,
                     vmax=330., cmap=cm)
    plt.xlabel('H magnitude')
    plt.ylabel('log$_{10}$ signal (median green chip)')
    plt.colorbar(sc, label='Radius (mm)')
    plt.plot(signal_fibers['tmass_h'][isort],
             np.log10(signal_fibers['signal_model'][isort]),
             color='black', linewidth=2, label='Model')
    scale = (expected.signal(exposure['fwhm_median_median']) /
             signal['signalref'])
    plt.plot(signal_fibers['tmass_h'][isort],
             np.log10(scale * signal_fibers['signal_model'][isort]),
             color='blue', linewidth=2, label='Expected given seeing')
    plt.xlim((10.5, 13.5))
    plt.ylim((2., np.log10(3000.)))
    plt.legend(loc=1)
    plt.savefig(outfile)
    plt.title(_title(exposure=exposure))
    plt.clf()


def sn2(outfile=None, exposure=None, signal=None, signal_fibers=None):
    """Plot signal-to-noise ratio vs. H magnitude

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    signal : ndarray
      structured array of signal information for this exposure
    signal_fibers : ndarray
      structured array of signal information for fibers in this exposure
    """
    cm = plt.cm.get_cmap('RdYlBu')
    isort = np.argsort(signal_fibers['tmass_h'])
    rf = np.sqrt(signal_fibers['xfocal']**2 + signal_fibers['yfocal']**2)
    sc = plt.scatter(signal_fibers['tmass_h'][isort],
                     np.log10(signal_fibers['sn2'][isort]),
                     c=rf[isort], s=50, alpha=0.8, vmin=0.,
                     vmax=330., cmap=cm)
    plt.xlabel('H magnitude')
    plt.ylabel('log$_{10}$ (S/N)$^2$ (median green chip)')
    plt.colorbar(sc, label='Focal radius (mm)')
    plt.plot(signal_fibers['tmass_h'][isort],
             np.log10(signal_fibers['sn2_model'][isort]),
             color='black', linewidth=2, label='Model')
    plt.plot(np.array([10.5, 13.5]),
             np.log10(np.array([signal['sn2ref'], signal['sn2ref']])),
             color='black', linestyle='dotted', linewidth=2,
             label='Reference (S/N)')
    plt.xlim((10.5, 13.5))
    plt.ylim((2., np.log10(3000.)))
    plt.legend(loc=1)
    plt.title(_title(exposure=exposure))
    plt.savefig(outfile)
    plt.clf()


def signal_focal(outfile=None, exposure=None, signal=None,
                 signal_fibers=None):
    """Plot signal vs. focal plane location

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    signal : ndarray
      structured array of signal information for this exposure
    signal_fibers : ndarray
      structured array of signal information for fibers in this exposure
    """
    cm = plt.cm.get_cmap('RdYlBu')
    diff = (np.log10(signal_fibers['signal']) -
            np.log10(signal_fibers['signal_model']))
    sc = plt.scatter(signal_fibers['xfocal'],
                     signal_fibers['yfocal'],
                     c=diff, s=50, alpha=0.8,
                     vmin=- signalrange, vmax=signalrange, cmap=cm)
    plt.ylim((-330., 330))
    plt.xlim((-330., 330))
    plt.xlabel("X (West)")
    plt.ylabel("Y (South)")
    plt.colorbar(sc, label='log$_{10}$ Signal / Model')
    xpa = np.array([- 1., 1.]) * 300. * np.sin(exposure['pa'] * np.pi / 180.)
    ypa = np.array([- 1., 1.]) * 300. * np.cos(exposure['pa'] * np.pi / 180.)
    plt.plot(xpa, ypa, color='blue', label='Parallactic direction')
    plt.legend(loc=1)
    plt.title(_title(exposure))
    plt.savefig(outfile)
    plt.clf()


def sn2_focal(outfile=None, exposure=None, signal=None,
              signal_fibers=None):
    """Plot signal-to-noise ratio vs. focal plane location

    Parameters
    ----------
    outfile : string
      file name to save figure to
    exposure : ndarray
      structured array of for this science exposure
    signal : ndarray
      structured array of signal information for this exposure
    signal_fibers : ndarray
      structured array of signal information for fibers in this exposure
    """
    cm = plt.cm.get_cmap('RdYlBu')
    diff = (np.log10(signal_fibers['sn2']) -
            np.log10(signal_fibers['sn2_model']))
    sc = plt.scatter(signal_fibers['xfocal'],
                     signal_fibers['yfocal'],
                     c=diff, s=50, alpha=0.8,
                     vmin=- sn2range, vmax=sn2range, cmap=cm)
    plt.ylim((-330., 330))
    plt.xlim((-330., 330))
    plt.xlabel("X (West)")
    plt.ylabel("Y (South)")
    plt.colorbar(sc, label='log$_{10}$ (S/N)$^2$ / Model')
    xpa = np.array([- 1., 1.]) * 300. * np.sin(exposure['pa'] * np.pi / 180.)
    ypa = np.array([- 1., 1.]) * 300. * np.cos(exposure['pa'] * np.pi / 180.)
    plt.plot(xpa, ypa, color='blue', label='Parallactic direction')
    plt.legend(loc=1)
    plt.title(_title(exposure))
    plt.savefig(outfile)
    plt.clf()


def signal_summary(outfile=None, summary=None, title=None):
    """Plot signal summary vs seeing

    Parameters
    ----------
    outfile : string
      file name to save figure to
    summary : ndarray
      structured array with summary information for exposures
    title : string
      title for plot
    """
    nfwhm = 200
    fwhm = 0.1 + 3. * (np.arange(nfwhm) / np.float32(nfwhm - 1))
    signal = expected.signal(fwhm)
    cm = plt.cm.get_cmap('RdYlBu')
    sc = plt.scatter(summary['fwhm_median_median'],
                     summary['signalref'],
                     c=summary['gdrms_median'], s=30, alpha=0.6,
                     vmin=0.2, vmax=0.7, cmap=cm)
    plt.plot(fwhm, signal, linewidth=2)
    plt.xlim((0.2, 2.6))
    plt.ylim((-10., 550.))
    plt.xlabel("FWHM in focus (arcsec)")
    plt.ylabel("Signal (counts) at H=12.1")
    plt.colorbar(sc, label='Guider RMS (arcsec)')
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()


def signal_summary_gdrms(outfile=None, summary=None, title=None):
    """Plot signal summary vs guide RMS

    Parameters
    ----------
    outfile : string
      file name to save figure to
    summary : ndarray
      structured array with summary information for exposures
    title : string
      title for plot
    """
    signal = expected.signal(summary['fwhm_median_median'])
    cm = plt.cm.get_cmap('RdYlBu')
    diff = summary['signalref'] / signal
    sc = plt.scatter(summary['gdrms_median'], diff,
                     c=summary['fwhm_median_median'], s=30, alpha=0.6,
                     vmin=0.7, vmax=2.4, cmap=cm)
    plt.xlim((0.15, 1.2))
    plt.ylim((-0.05, 1.4))
    plt.xlabel("Guider RMS (arcsec)")
    plt.ylabel("Reference signal / Expected signal given seeing")
    plt.colorbar(sc, label='FWHM in focus (arcsec)')
    plt.title(title)
    plt.savefig(outfile)
    plt.clf()
