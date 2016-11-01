import os
import sys
import re
import numpy as np
from scipy.ndimage.filters import median_filter

import ccdproc as ccd

import astropy.units as u
from astropy import table

from astropy import log
log.setLevel('WARNING')
log.disable_warnings_logging()




if __name__ == '__main__':
    import sort
    filepath = '/Users/vysosuser/ShutterMap/V5/20160723UT/'
    image_table = sort.sort(filepath)
    bycat = image_table.group_by('CATEGORY')

    bias_table = bycat.groups[bycat.groups.keys['CATEGORY'] == b'Bias']
    dark_table = bycat.groups[bycat.groups.keys['CATEGORY'] == b'Dark']
    flat_table = bycat.groups[bycat.groups.keys['CATEGORY'] == b'Flat']
    flat_table = flat_table.group_by('EXPTIME')

    bias_files = [os.path.join(x['DIR'].decode('utf8'), x['file'])
                  for x in bias_table]
    bias_images = [ccd.fits_ccddata_reader(f, unit='adu') for f in bias_files]
    master_bias = ccd.combine(bias_images, combine='median')

    flats = []
    for exptime in flat_table.groups.keys['EXPTIME']:
        thisexptime = flat_table.groups[flat_table.groups.keys['EXPTIME'] == exptime]
        print('Found {} images with exposure time {:.1f}'.format(
              len(thisexptime), float(exptime) ) )
        ## Combine each group of flats with same exposure times
        files = [os.path.join(x['DIR'].decode('utf8'), x['file'])
                 for x in thisexptime]
        images = [ccd.fits_ccddata_reader(f, unit='adu') for f in files]
        images = [ccd.subtract_bias(im, master_bias) for im in images]
        combined_flat = ccd.combine(images, method='median', scale=lambda x: 1./np.median(x))
        flats.append(combined_flat)

    ratios = []
    for i,flat in enumerate(flats):
        if i > 0:
            for j in range(i):
                ratio_ij = flats[i].divide(flats[j])
                ratios.append(ratio_ij)

    masks = np.array([ccd.ccdmask(ratio, ncsig=25, nlsig=25, lsigma=7, hsigma=7, findbadcolumns=True)
                      for ratio in ratios])
    from astropy.io import fits
    for i,m in enumerate(masks):
        maskint = np.array(m, dtype=np.int)
        hdu = fits.PrimaryHDU(maskint)
        hdu.writeto('mask{}.fits'.format(i))
    mask = np.any(masks, axis=0)
    n_bad_pix = np.sum(np.array(mask, dtype=np.int))
    pct_bad_pix = n_bad_pix/len(mask.ravel())*100.
    print('Found {:d} ({:.1f}%) bad pixels.'.format(n_bad_pix, pct_bad_pix))

    if os.path.exists('mask.fits'): os.remove('mask.fits')
    from astropy.io import fits
    maskint = np.array(mask, dtype=np.int)
    hdu = fits.PrimaryHDU(maskint)
    hdu.writeto('mask.fits')
