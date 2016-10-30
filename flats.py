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


def ccdmask(ratio, ncmed=7, nlmed=7, ncsig=15, nlsig=15, lsigma=6, hsigma=6,
            ngood=5, linterp=2, cinterp=3, eqinterp=2):
    '''
    Uses method from the IRAF ccdmask task to generate a mask based on the given
    input.

    The Following documentation is copied directly from:
    http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?ccdmask

    The input image is first subtracted by a moving box median. The median is
    unaffected by bad pixels provided the median size is larger that twice the
    size of a bad region. Thus, if 3 pixel wide bad columns are present then the
    column median box size should be at least 7 pixels. The median box can be a
    single pixel wide along one dimension if needed. This may be appropriate for
    spectroscopic long slit data.

    The median subtracted image is then divided into blocks of size nclsig by
    nlsig. In each block the pixel values are sorted and the pixels nearest the
    30.9 and 69.1 percentile points are found; this would be the one sigma
    points in a Gaussian noise distribution. The difference between the two
    count levels divided by two is then the local sigma estimate. This algorithm
    is used to avoid contamination by the bad pixel values. The block size must
    be at least 10 pixels in each dimension to provide sufficient pixels for a
    good estimate of the percentile sigma. The sigma uncertainty estimate of
    each pixel in the image is then the sigma from the nearest block.

    The deviant pixels are found by comparing the median subtracted residual to
    a specified sigma threshold factor times the local sigma above and below
    zero (the lsigma and hsigma parameters). This is done for individual pixels
    and then for column sums of pixels (excluding previously flagged bad pixels)
    from two to the number of lines in the image. The sigma of the sums is
    scaled by the square root of the number of pixels summed so that
    statistically low or high column regions may be detected even though
    individual pixels may not be statistically deviant. For the purpose of this
    task one would normally select large sigma threshold factors such as six or
    greater to detect only true bad pixels and not the extremes of the noise
    distribution.

    As a final step each column is examined to see if there are small segments
    of unflagged pixels between bad pixels. If the length of a segment is less
    than that given by the ngood parameter all the pixels in the segment are
    also marked as bad.
    
    ncmed = 7, nlmed = 7
    The column and line size of a moving median rectangle used to estimate the
    uncontaminated local signal. The column median size should be at least 3
    pixels to span single bad columns.

    ncsig = 15, nlsig = 15
    The column and line size of regions used to estimate the uncontaminated
    local sigma using a percentile. The size of the box should contain of order
    100 pixels or more.

    lsigma = 6, hsigma = 6
    Positive sigma factors to use for selecting pixels below and above the
    median level based on the local percentile sigma.

    ngood = 5
    Gaps of undetected pixels along the column direction of length less than
    this amount are also flagged as bad pixels.

    linterp = 2
    Mask code for pixels having a bounding good pixel separation which is
    smaller along lines; i.e. to use line interpolation along the narrower
    dimension.

    cinterp = 3
    Mask code for pixels having a bounding good pixel separation which is
    smaller along columns; i.e. to use columns interpolation along the narrower
    dimension.

    eqinterp = 2
    Mask code for pixels having a bounding good pixel separation which is equal
    along lines and columns.
    '''
    medsub = ratio.data - median_filter(ratio.data, size=(nlsig, ncsig))

    nbl = int(np.floor(ratio.data.shape[0] / nlsig))
    nbc = int(np.floor(ratio.data.shape[1] / ncsig))
    mask = np.zeros(ratio.data.shape, dtype=np.bool)
    for i in range(nbl):
        for j in range(nbc):
            block = medsub[i*nlsig:(i+1)*nlsig,j*ncsig:(j+1)*ncsig]
            high = np.percentile(block.ravel(), 69.1)
            low = np.percentile(block.ravel(), 30.9)
            block_sigma = (high-low)/2.0
            block_mask = ( (block > hsigma*block_sigma) |
                           (block < -lsigma*block_sigma) )
            mblock = np.ma.MaskedArray(block, mask=block_mask)
            csum = np.ma.sum(mblock, axis=0)
            csum_sigma = np.array([np.sqrt(ncsig-x)*block_sigma
                                   for x in np.ma.sum(mblock.mask, axis=0)])
            colmask = ( (csum > hsigma*csum_sigma) |
                        (csum < -lsigma*csum_sigma) )
            colmask = colmask.filled(True)
            for c,masked in enumerate(colmask):
                block_mask[:,c] = block_mask[:,c] | np.array([masked]*nlsig)
            mask[i*nlsig:(i+1)*nlsig,j*ncsig:(j+1)*ncsig] = block_mask
    return mask


if __name__ == '__main__':
    import sort
    filepath = '/Users/vysosuser/ShutterMap/V5'
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

    masks = np.array([ccdmask(ratio, ncsig=50, nlsig=50, lsigma=7, hsigma=7)
                      for ratio in ratios])
#     from astropy.io import fits
#     for i,m in enumerate(masks):
#         maskint = np.array(m, dtype=np.int)
#         hdu = fits.PrimaryHDU(maskint)
#         hdu.writeto('mask{}.fits'.format(i))
    mask = np.any(masks, axis=0)
    n_bad_pix = np.sum(np.array(mask, dtype=np.int))
    pct_bad_pix = n_bad_pix/len(mask.ravel())*100.
    print('Found {:d} ({:.1f}%) bad pixels.'.format(n_bad_pix, pct_bad_pix))

    if os.path.exists('mask.fits'): os.remove('mask.fits')
    from astropy.io import fits
    maskint = np.array(mask, dtype=np.int)
    hdu = fits.PrimaryHDU(maskint)
    hdu.writeto('mask.fits')
