import os

from datetime import datetime as dt
from datetime import timedelta as tdelta

import numpy as np
from astropy.io import fits
import astropy.units as u
import ccdproc
from ccdproc import fits_ccddata_reader as ccdread

from .config import get_config
from .sort import *

def make_master_bias(date):
    config = get_config()
    date_dto = dt.strptime(date, '%Y%m%dUT')
    one_day = tdelta(1)

    bias_root = config.get('BiasPath', os.path.abspath('.'))
    bias_root = os.path.expanduser(bias_root)
    bias_root = os.path.abspath(bias_root)

    bias_files = []
    for i in np.arange(config['CalsWindow']):
        this_day = dt.strftime((date_dto - i*one_day), '%Y%m%dUT')
        bias_path = bias_root.replace('[YYYYMMDDUT]', this_day)
        if os.path.exists(bias_path):
            image_table = get_image_table(bias_path, 'Bias')
            new_files = [os.path.join(bias_path, fn) for fn in image_table['file']]
            bias_files.extend(new_files)

    print('Combining {:d} files in to master bias'.format(len(bias_files)))
    bias_images = []
    for bias_file in bias_files:
        try:
            bias_image = ccdread(bias_file)
        except ValueError:
            bias_image = ccdread(bias_file, unit='adu')
        # Estimate uncertainty as read noise
        readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
        try:
            readnoise.value_from(bias_image.header)
        except KeyError:
            readnoise.value = config.get('RN', 10.0)
        bias_image.uncertainty = np.ones(bias_image.data.shape) * readnoise.value

    master_bias = ccdproc.combine(bias_images, combine='median')

    # Update header
    master_bias.header.add_history(
           '{:d} files combined to make master bias'.format(len(bias_files)))
    for f in bias_files:
        master_bias.header.add_history(' {}'.format(f))

    # Write master bias to file
    mbfn = '{}_{}.fits'.format(config.get('MasterBiasRootName', 'MasterBias'),
                               date)
    print('Writing {}'.format(mbfn))
    ccdproc.fits_ccddata_writer(master_bias, mbfn)

    return master_bias


def make_master_dark(date, master_bias=None):
    config = get_config()
    date_dto = dt.strptime(date, '%Y%m%dUT')
    one_day = tdelta(1)

    dark_root = config.get('DarkPath', os.path.abspath('.'))
    dark_root = os.path.expanduser(dark_root)
    dark_root = os.path.abspath(dark_root)

    dark_files = []
    for i in np.arange(config['CalsWindow']):
        this_day = dt.strftime((date_dto - i*one_day), '%Y%m%dUT')
        dark_path = dark_root.replace('[YYYYMMDDUT]', this_day)
        if os.path.exists(dark_path):
            image_table = get_image_table(dark_path, 'Dark')
            new_files = [os.path.join(dark_path, fn) for fn in image_table['file']]
            dark_files.extend(new_files)

    print('Combining {:d} files in to master dark'.format(len(dark_files)))
    if master_bias:
        try:
            dark_images = [ccdread(f).subtract_bias(master_bias)
                           for f in dark_files]
        except ValueError:
            dark_images = [ccdread(f, unit='adu').subtract_bias(master_bias)
                           for f in dark_files]
    else:
        try:
            dark_images = [ccdread(f) for f in dark_files]
        except ValueError:
            dark_images = [ccdread(f, unit='adu') for f in dark_files]
    
    master_dark = ccdproc.combine(dark_images, combine='median')

    # Update header
    master_dark.header.add_history(
           '{:d} files combined to make master dark'.format(len(dark_files)))
    for f in dark_files:
        master_dark.header.add_history(' {}'.format(f))

    # Write master dark to file
    mbfn = '{}_{}.fits'.format(config.get('MasterDarkRootName', 'MasterDark'),
                               date)
    print('Writing {}'.format(mbfn))
    ccdproc.fits.ccddata_writer(master_dark, mbfn)
