from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import os
import re
import glob
import fitsio
import astropy.time as time
import sdss_access
import pydl.pydlutils.yanny as yanny

spath = sdss_access.path.Path()


def parangle(ha=None, dec=None, lat=-29.0):
    """Return parallactic angle of observation

    Parameters
    ----------
    ha : np.float32
      hour angle of observation (deg)
    dec : np.float32
      declination of observation (deg)
    lat : np.float32
      latitude of observatory (deg)

    Returns:
    -------
    pa : np.float32
      parallactic angle (deg East of North)
    """
    r2d = 180. / np.pi
    d2r = np.pi / 180.
    pa = - r2d * np.arctan2(- np.sin(d2r * ha),
                            np.cos(d2r * dec) * np.tan(d2r * lat) -
                            np.sin(d2r * dec) * np.cos(d2r * ha))
    return(pa)


def exposures(mjd=None, telescope=None, gcam=None):
    """Gather exposure information from ap1D or as1D files

    Parameters
    ----------
    mjd : np.int32, or int
      MJD of observation
    telescope : string
      telescope of observations ('apo' or 'lco')
    gcam : ndarray
      structured array of guider camera summary for MJD

    Returns:
    -------
    exposures : ndarray
      structured array of science exposures for MJD

    Comments:
    --------

    Calculates some basic information from the processed guide images;
      guider_rms_median : median guideRMS value (calculated)
      gdrms_median  : median guideRMS value (from proc-gimg headers)
      gstart  : first guider image
      gend  : last guider image
      seeing_median  : median seeing value
      fwhm_median_median : median in-focus guide FWHM, median over exposure

    """

    # Set location based on telescope
    if(telescope == 'apo'):
        prefix = 'ap'
        longitude = np.float64(254.1803)
        latitude = np.float64(32.7797556)
        instrument = 'apogee-n'
    else:
        prefix = 'as'
        longitude = np.float64(289.3075)
        latitude = np.float64(-29.01472)
        instrument = 'apogee-s'
    location = (str(longitude) + 'd', str(latitude) + 'd')

    # Read in full list of plates
    plist = os.path.join(os.getenv('PLATELIST_DIR'),
                         'platePlans.par')
    plans = yanny.yanny(plist)['PLATEPLANS']

    # Find list of exposures to analyze
    ap1ddir = os.path.dirname(spath.full('ap1D', apred='apogee2',
                                         prefix=prefix,
                                         instrument=instrument,
                                         chip='a', mjd=mjd, num=0))
    efiles = glob.glob(os.path.join(ap1ddir, prefix + '1D-b-*.fits'))

    if(len(efiles) == 0):
        print("No exposures on MJD {mjd}".format(mjd=mjd))
        return None

    # Initialize array of exposures
    e_dtype = [('raCen', np.dtype(np.float64), 1),
               ('decCen', np.dtype(np.float64), 1),
               ('expno', np.dtype(np.int32), 1),
               ('expfile', np.dtype(np.str_), 200),
               ('name', np.dtype(np.str_), 200),
               ('imagetyp', np.dtype(np.str_), 8),
               ('date-obs', np.dtype(np.str_), 21),
               ('exptime', np.dtype(np.float32), 1),
               ('lst', np.dtype(np.float32), 1),
               ('ha', np.dtype(np.float32), 1),
               ('programname', np.dtype(np.str_), 200),
               ('ha_design', np.dtype(np.float32), 1),
               ('ha_observable_min', np.dtype(np.float32), 1),
               ('ha_observable_max', np.dtype(np.float32), 1),
               ('pa', np.dtype(np.float32), 1),
               ('plateid', np.dtype(np.int32), 1),
               ('dithpix', np.dtype(np.int32), 1),
               ('seeing_median', np.dtype(np.float32), 1),
               ('fwhm_median_median', np.dtype(np.float32), 1),
               ('guider_rms_median', np.dtype(np.float32), 1),
               ('gdrms_median', np.dtype(np.float32), 1),
               ('gstart', np.dtype(np.int32), 1),
               ('gend', np.dtype(np.int32), 1),
               ('mjd_start', np.dtype(np.float64), 1),
               ('mjd_end', np.dtype(np.float64), 1),
               ('mjd_mid', np.dtype(np.float64), 1),
               ('mjd', np.dtype(np.int32), 1)]
    exps = np.zeros(len(efiles), dtype=e_dtype)
    exps['expfile'] = efiles
    exps['mjd'] = mjd

    # Loop through each exposure
    for exp in exps:
        # Basic info
        header = fitsio.read_header(exp['expfile'])
        expno = re.search(".*" + prefix + "1D-b-([0-9]{8})\.fits$",
                          exp['expfile'])
        exp['expno'] = expno.group(1)
        exp['imagetyp'] = header['IMAGETYP']
        exp['name'] = header['NAME']
        exp['date-obs'] = header['DATE-OBS']
        exp['exptime'] = header['EXPTIME']
        exp['dithpix'] = header['DITHPIX']
        try:
            exp['plateid'] = header['PLATEID']
        except ValueError:
            exp['plateid'] = 0

        # What is the time of the exposure?
        dt = str(header['DATE-OBS'].decode())
        tt = time.Time(dt, location=location)
        tt.format = 'mjd'
        exp['mjd_start'] = tt.value
        exp['mjd_end'] = exp['mjd_start'] + exp['exptime'] / 3600. / 24.
        exp['mjd_mid'] = 0.5 * (exp['mjd_start'] + exp['mjd_end'])
        tt_mid = time.Time(exp['mjd_mid'], format='mjd',
                           location=location)
        exp['lst'] = tt_mid.sidereal_time('apparent').value * 15.

        # Summarize guider information for exposure
        igcam = np.where((gcam['mjd'] > exp['mjd_start']) &
                         (gcam['mjd'] < exp['mjd_end']))[0]
        if(len(igcam) > 0):
            exp['seeing_median'] = np.nanmedian(gcam['seeing'][igcam])
            exp['fwhm_median_median'] = np.nanmedian(gcam['fwhm_median'][igcam])
            exp['guider_rms_median'] = np.nanmedian(gcam['guider_rms'][igcam])
            exp['gdrms_median'] = np.nanmedian(gcam['gdrms'][igcam])
            exp['gstart'] = gcam['indx'][igcam].min()
            exp['gend'] = gcam['indx'][igcam].max()
        else:
            exp['seeing_median'] = -1.
            exp['fwhm_median_median'] = -1.
            exp['guider_rms_median'] = -1.
            exp['gstart'] = -1
            exp['gend'] = -1

        # Now store some information about the plate
        if(exp['plateid'] > 0):
            iplans = np.where(plans['plateid'] == exp['plateid'])[0]
            plugfile = spath.full('plPlugMapP', plateid=exp['plateid'])
            plug = yanny.yanny(plugfile)
            exp['raCen'] = plans['raCen'][iplans]
            exp['decCen'] = plans['decCen'][iplans]
            exp['ha'] = (exp['lst'] - exp['raCen'] + 360.) % 360.
            exp['programname'] = plans['programname'][iplans[0]]
            exp['ha_design'] = plans['ha'][iplans, 0]
            exp['ha_observable_min'] = np.float32(plug['ha_observable_min'].split()[0])
            exp['ha_observable_max'] = np.float32(plug['ha_observable_max'].split()[0])
            exp['pa'] = parangle(ha=exp['ha'], dec=exp['decCen'],
                                 lat=latitude)

    # Keep only object exposures
    ikeep = np.where(exps['imagetyp'] == 'Object  ')[0]
    exps = exps[ikeep]

    return(exps)
