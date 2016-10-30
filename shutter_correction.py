import os
import re
import numpy as np
from scipy.ndimage.filters import median_filter

import ccdproc as ccd

import astropy.units as u
from astropy import table
from astropy import log
log.setLevel('WARNING')
log.disable_warnings_logging()

import sort

def determine_shutter_correction_map(filepath):
    image_table = sort.sort(filepath)
    bycat = image_table.group_by('CATEGORY')

    bias_table = bycat.groups[bycat.groups.keys['CATEGORY'] == b'Bias']
    flat_table = bycat.groups[bycat.groups.keys['CATEGORY'] == b'Flat']

    ## Bias Subtract all flats
    bias_files = [os.path.join(x['DIR'].decode('utf8'), x['file'])
                  for x in bias_table]
    bias_images = [ccd.fits_ccddata_reader(f, unit='adu') for f in bias_files]
    master_bias = ccd.combine(bias_images, combine='median')

    flat_table = flat_table.group_by('EXPTIME')
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
        combined_flat = ccd.combine(images, method='average', sigma_clip=True,
                                    sigma_clip_low_thresh=5,
                                    sigma_clip_high_thresh=5)
        flats.append(combined_flat)
        assert combined_flat.shape == flats[0].shape

    def nflat(flat, delta=25):
        nl, nc = flat.data.shape
        normalization_factor = np.median(flat.data[nl-delta:nl+delta,nc-delta:nc+delta])
        nflat = flat.data / normalization_factor
        return nflat

    gain = 1.6
    a = gain*np.sum([nflat(f) for f in flats])
    b = gain*np.sum([nflat(f)/float(f.header['EXPTIME']) for f in flats])
    c = gain*np.sum([nflat(f)/float(f.header['EXPTIME'])**2 for f in flats])

    Dijcomb = ccd.combine(flats, method='average')
    Dij = gain*Dijcomb.data*len(flats)

    Eijm = [ccd.CCDData(data=f.data/float(f.header['EXPTIME']), unit='adu')
            for f in flats]
    Eijcomb = ccd.combine(Eijm, method='average')
    Eij = gain*Eijcomb.data*len(flats)

    alpha_ij = 1.0 / (a*c - b**2) * (c*Dij - b*Eij)
    delta_alpha = np.sqrt( c / (a*c - b**2) )
    beta_ij = 1.0 / (a*c - b**2) * (a*Eij - b*Dij)
    delta_beta = np.sqrt( a / (a*c - b**2) )

    SF = beta_ij / alpha_ij  # SF = -t_SH * (1-SH_ij) in the paper terminology
    delta_SF = np.sqrt( (alpha_ij**-1 * delta_beta)**2 + 
                        (beta_ij * alpha_ij**-2 * delta_alpha)**2 )
    SFF = median_filter(SF, size=(25, 25))

    from astropy.io import fits
    fits_file = 'SFF.fits'
    if os.path.exists(fits_file):
        os.remove(fits_file)
    hdup = fits.PrimaryHDU(SF)
    hdus = fits.ImageHDU(delta_SF)
    hdul = fits.HDUList([hdup, hdus])
    hdul.writeto(fits_file)


if __name__ == '__main__':
    filepath = '/Users/vysosuser/ShutterMap/V5'
    determine_shutter_correction_map(filepath)
