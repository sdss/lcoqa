from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Process guide camera images into a single roll-up file with statistics.

import numpy as np
import os
import fitsio
import astropy.time as time
import re


def gcam(mjd=None, telescope=None):
    """Process an MJD of processed guide camera files into a single file

    Parameters:
    -----------

    mjd : integer
      MJD (days) of observation
    telescope : str
      which telescope ('lco' or 'apo')

    Comments:
    --------

    Return array containing a row for each guider image.

    Note that 'date-obs' from each guider image header is translated into
    a floating point 'mjd'.

    This copies all of the information from extension 6 of each gcam file.
    It also calculates:

      gdrms - Guider rms offset of in-focus fibers (arcsec)
      seeing - Seeing reported by guider (arcsec)
      fwhm_median - Median of FWHM of in-focus fibers (arcsec)
      fwhm_mean - Mean of FWHM of in-focus fibers (arcsec)
    """

    if(telescope == 'lco'):
        gdir = os.path.join(os.getenv('GCAM_DATA_2S'), str(mjd))
    else:
        gdir = os.path.join(os.getenv('GCAM_DATA'), str(mjd))

    if(os.path.isdir(gdir) is False):
        print("No MJD: {mjd}".format(mjd=mjd))
        return

    files = os.listdir(gdir)
    gcam0_dtype = [('guider_rms', np.dtype(np.float32), 1),
                   ('gdrms', np.dtype(np.float32), 1),
                   ('seeing', np.dtype(np.float32), 1),
                   ('fwhm_median', np.dtype(np.float32), 1),
                   ('fwhm_mean', np.dtype(np.float32), 1),
                   ('indx', np.dtype(np.int32), 1),
                   ('date-obs', np.dtype(np.str_), 21),
                   ('mjd', np.dtype(np.float64), 1),
                   ('enabled', np.dtype(np.int32), 17),
                   ('ra', np.dtype(np.float64), 17),
                   ('dec', np.dtype(np.float64), 17),
                   ('xFocal', np.dtype(np.float32), 17),
                   ('yFocal', np.dtype(np.float32), 17),
                   ('focusOffset', np.dtype(np.float32), 17),
                   ('xstar', np.dtype(np.float32), 17),
                   ('ystar', np.dtype(np.float32), 17),
                   ('xCenter', np.dtype(np.float32), 17),
                   ('yCenter', np.dtype(np.float32), 17),
                   ('dx', np.dtype(np.float32), 17),
                   ('dy', np.dtype(np.float32), 17),
                   ('dRA', np.dtype(np.float32), 17),
                   ('dDec', np.dtype(np.float32), 17),
                   ('fwhm', np.dtype(np.float32), 17),
                   ('flux', np.dtype(np.float32), 17),
                   ('mag', np.dtype(np.float32), 17),
                   ('sky', np.dtype(np.float32), 17),
                   ('skymag', np.dtype(np.float32), 17)]

    count = 0
    for file in sorted(files):
        if(file.startswith('proc-')):
            count = count + 1
    gcam = np.zeros(count, dtype=gcam0_dtype)

    count = 0
    for file in sorted(files):
        if(file.startswith('proc-')):
            indx_search = re.search("proc-gimg-([0-9]{4}).fits", file)
            indx = np.int32(indx_search.group(1))
            header = fitsio.read_header(os.path.join(gdir, file), ext=0)
            if('object' in header['IMAGETYP']):
                rescale = (3600. /
                           np.float32(header['PLATSCAL']))
                data = fitsio.read(os.path.join(gdir, file), ext=6)
                ii = np.nonzero((data['enabled']) &
                                (data['gprobebits'] == 0) &
                                (data['dx'] == data['dx']))[0]
                if(len(ii) > 0):
                    soff = (data['dx'][ii]**2 + data['dy'][ii]**2)
                    rms = np.sqrt(soff.mean())
                    if(rms == rms):
                        gcam['guider_rms'][count] = rms * rescale
                ii = np.nonzero((data['enabled']) &
                                (data['gprobebits'] == 0) &
                                (np.abs(data['focusOffset']) < 1.) &
                                (data['fwhm'] == data['fwhm']))[0]
                if(len(ii) > 0):
                    fwhm_median = np.median(data['fwhm'][ii])
                    if(fwhm_median == fwhm_median):
                        gcam['fwhm_median'][count] = fwhm_median
                    fwhm_mean = np.mean(data['fwhm'][ii])
                    if(fwhm_mean == fwhm_mean):
                        gcam['fwhm_mean'][count] = fwhm_mean
                for (name, pos) in zip(gcam.dtype.names,
                                       range(len(gcam.dtype.names))):
                    if(pos >= 8):
                        try:
                            gcam[name][count] = data[name]
                        except ValueError:
                            pass
                gcam['gdrms'][count] = np.float32(header['GDRMS'])
                gcam['seeing'][count] = np.float32(header['SEEING'])
                dt = str(header['DATE-OBS'].decode())
                gcam['date-obs'][count] = dt
                gcam['indx'][count] = indx
                tt = time.Time(dt)
                tt.format = 'mjd'
                gcam['mjd'][count] = tt.value
                count = count + 1

    gcam = gcam[0:count - 1]
    return(gcam)
