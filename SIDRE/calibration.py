import os

from datetime import datetime as dt
from datetime import timedelta as tdelta

import numpy as np
from astropy.io import fits
import astropy.units as u
import ccdproc

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
    try:
        bias_images = [ccdproc.fits_ccddata_reader(bf) for bf in bias_files]
    except ValueError:
        bias_images = [ccdproc.fits_ccddata_reader(bf, unit='adu')
                       for bf in bias_files]
    master_bias = ccdproc.combine(bias_images, combine='median')

    # Update header
    master_bias.header.add_history(
           '{:d} files combined to make master bias'.format(len(bias_files)))
    for f in bias_files:
        master_bias.header.add_history(' {}'.format(f))

    # Write master bias to file
    mbfn = '{}_{}.fits'.format(config.get('MasterBiasRootName', 'MasterBias_'),
                               date)
    print('Writing {}'.format(mbfn))
    ccdproc.fits.ccddata_writer(master_bias, mbfn)
