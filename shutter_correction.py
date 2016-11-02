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

def determine_shutter_correction_map(filepath, output='ShutterMap.fits',
              gain=1.6):
    image_table = sort.sort(filepath)
    bycat = image_table.group_by('CATEGORY')

    assert 'Bias' in set(image_table['CATEGORY'])
    assert 'Flat' in set(image_table['CATEGORY'])
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
        if len(files) > 1:
            images = [ccd.fits_ccddata_reader(f, unit='adu') for f in files]
            images = [ccd.subtract_bias(im, master_bias) for im in images]
            sigma_clip = len(files) > 5
            combined_flat = ccd.combine(images, method='average',
                                        sigma_clip=sigma_clip,
                                        sigma_clip_low_thresh=5,
                                        sigma_clip_high_thresh=5)
            flats.append(combined_flat)
            assert combined_flat.shape == flats[0].shape
        else:
            flats.append(ccd.fits_ccddata_reader(files[0], unit='adu'))

    def nflat(flat, delta=25):
        nl, nc = flat.data.shape
        normalization_factor = np.median(flat.data[nl-delta:nl+delta,nc-delta:nc+delta])
        nflat = flat.data / normalization_factor
        return nflat

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
    SNR = np.mean(abs(SF)/delta_SF)
    print('Typical pixel SNR of correction = {:.1f}'.format(SNR))
    ShutterMap = ccd.CCDData(data=SF, uncertainty=delta_SF, unit=u.second,
                     meta={'SNR': (SNR, 'Mean signal to noise per pixel')})

    if output:
        if os.path.exists(output):
            os.remove(output)
        ccd.fits_ccddata_writer(ShutterMap, output)

    return ShutterMap


if __name__ == '__main__':
#     filepath = '/Users/vysosuser/ShutterMap/V5/20161027UT'
#     determine_shutter_correction_map(filepath)

    filepath = ['/Volumes/Drobo/V5/Images/20161025UT/AutoFlat',
                '/Volumes/Drobo/V5/Images/20161025UT/Calibration']
    determine_shutter_correction_map(filepath, output='ShutterMap_twi.fits')
