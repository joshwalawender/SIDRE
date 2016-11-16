import os
import sys
import numpy as np

import ccdproc as ccd
import sep

import astropy.units as u

from astropy import log
log.setLevel('WARNING')
log.disable_warnings_logging()

import matplotlib.pyplot as plt
from matplotlib import rcParams


def solve_astrometry(im):
    pass



if __name__ == '__main__':
    filepath = '/Users/vysosuser/astrometry_test'
    filename = 'V5_Taurus-1-PSr3-20161023at095814.fts'
    file = os.path.join(filepath, filename)
    im = ccd.fits_ccddata_reader(file, unit='adu')
    print(im.wcs.to_header())
    

    m, s = np.mean(im.data), np.std(im.data)
    plt.imshow(im.data, interpolation='nearest', cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
    plt.colorbar()

    bkg = sep.Background(im.data, bw=64, bh=64, fw=3, fh=3)
    bkg_rms = bkg.rms()
    data_sub = im.data - bkg
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
#     for o in objects:
#         print(o['x'], o['y'])
    print(len(objects))   
#     flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
#                                      3.0, err=bkg.globalrms, gain=1.6)

