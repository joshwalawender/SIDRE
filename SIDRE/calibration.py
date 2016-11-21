import os

from datetime import datetime as dt
from datetime import timedelta as tdelta

import numpy as np
from astropy.io import fits
import astropy.units as u
import ccdproc
from ccdproc import fits_ccddata_reader as ccdread
from astropy.nddata import StdDevUncertainty

from .config import get_config
from .sort import *

def make_master_bias(date, clobber=True):
    '''
    Given a date string (in YYYYMMDDUT format), collect files which should be
    used to make a master bias file and then combine them to make the master
    bias.
    '''
    config = get_config()
    date_dto = dt.strptime(date, '%Y%m%dUT')
    one_day = tdelta(1)

    bias_root = config.get('BiasPath', os.path.abspath('.'))
    bias_files = []
    # Collect filenames of bias files within time window
    for i in np.arange(config['CalsWindow']):
        this_day = dt.strftime((date_dto - i*one_day), '%Y%m%dUT')
        bias_path = bias_root.replace('[YYYYMMDDUT]', this_day)
        if os.path.exists(bias_path):
            image_table = get_image_table(bias_path, 'Bias')
            new_files = [os.path.join(bias_path, fn) for fn in image_table['file']]
            bias_files.extend(new_files)

    nbiases = len(bias_files)
    print('Combining {:d} files in to master bias'.format(nbiases))
    bias_images = []
    for i,bias_file in enumerate(bias_files):
        print('  Reading {:d}/{:d}: {}'.format(i+1, nbiases, bias_file))
        try:
            bias_image = ccdread(bias_file, verify=True)
        except ValueError:
            bias_image = ccdread(bias_file, unit='adu', verify=True)
        # Estimate uncertainty as read noise
        readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
        try:
            readnoise.value_from(bias_image.header)
        except KeyError:
            readnoise.value = config.get('RN')
        bias_image.uncertainty = StdDevUncertainty(np.ones(bias_image.data.shape) * readnoise.value)
        bias_images.append(bias_image)

    master_bias = ccdproc.combine(bias_images, combine='median')

    # Update header
    master_bias.header.add_history(
           '{:d} files combined to make master bias'.format(nbiases))
    for i,f in enumerate(bias_files):
        master_bias.header.add_history('{:d}/{:d}: {}'.format(i, nbiases, os.path.basename(f)))

    # Write master bias to file
    mbfn = '{}_{}.fits'.format(config.get('MasterBiasRootName', 'MasterBias'),
                               date)
    mbf = os.path.join(config.get('MasterPath', '/'), mbfn)
    print('Writing {}'.format(mbf))
    if clobber and os.path.exists(mbf):
        os.remove(mbf)
    ccdproc.fits_ccddata_writer(master_bias, mbf, checksum=True)

    return master_bias


def make_master_dark(date, master_bias=None, clobber=True):
    '''
    Given a date string (in YYYYMMDDUT format), collect files which should be
    used to make a master dark file and then combine them to make the master
    dark.
    '''
    config = get_config()
    date_dto = dt.strptime(date, '%Y%m%dUT')
    one_day = tdelta(1)

    dark_root = config.get('DarkPath', os.path.abspath('.'))

    dark_files = []
    for i in np.arange(config['CalsWindow']):
        this_day = dt.strftime((date_dto - i*one_day), '%Y%m%dUT')
        dark_path = dark_root.replace('[YYYYMMDDUT]', this_day)
        if os.path.exists(dark_path):
            image_table = get_image_table(dark_path, 'Dark')
            new_files = [os.path.join(dark_path, fn) for fn in image_table['file']]
            dark_files.extend(new_files)

    ndarks = len(dark_files)
    print('Combining {:d} files in to master dark'.format(ndarks))
    dark_images = []
    for i,dark_file in enumerate(dark_files):
        print('  Reading {:d}/{:d}: {}'.format(i+1, ndarks, dark_file))
        try:
            dark_image = ccdread(dark_file, verify=True)
        except ValueError:
            dark_image = ccdread(dark_file, unit='adu', verify=True)
        if master_bias:
            dark_image = ccdproc.subtract_bias(dark_image, master_bias)

        # Estimate uncertainty
        readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
        try:
            readnoise.value_from(dark_image.header)
        except KeyError:
            readnoise.value = config.get('RN')
        gain = ccdproc.Keyword('GAIN', unit=u.electron/u.adu)
        try:
            gain.value_from(dark_image.header)
        except KeyError:
            gain.value = config.get('Gain')

        dark_image = ccdproc.gain_correct(dark_image, gain.value)
        dark_image = ccdproc.create_deviation(dark_image, readnoise=readnoise.value)
        dark_images.append(dark_image)

    master_dark = ccdproc.combine(dark_images, combine='median')

    # Update header
    master_dark.header.add_history(
           '{:d} files combined to make master dark'.format(ndarks))
    for i,f in enumerate(dark_files):
        master_dark.header.add_history('{:d}/{:d}: {}'.format(i, ndarks, os.path.basename(f)))

    # Write master dark to file
    mdfn = '{}_{}.fits'.format(config.get('MasterDarkRootName', 'MasterDark'),
                               date)
    mdf = os.path.join(config.get('MasterPath', '/'), mdfn)
    print('Writing {}'.format(mdf))
    if clobber and os.path.exists(mdf):
        os.remove(mdf)
    ccdproc.fits_ccddata_writer(master_dark, mdf, checksum=True)
    
    return master_dark