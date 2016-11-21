import os

from datetime import datetime as dt

from astropy.io import fits
import astropy.units as u
import ccdproc

from .config import get_config
from .calibration import *


def get_image_date(im, datekw='DATE-OBS', datefmt='%Y-%m-%dT%H:%M:%S'):
    '''
    Return string with the UT date of the image in YYYYMMDDUT format
    '''
    dto = dt.strptime(im.header[datekw], datefmt)
    return dto.strftime('%Y%m%dUT')


def get_master_file(date, type='Bias'):
    '''
    
    '''
    assert type in ['Bias', 'Dark', 'Flat']
    config = get_config()
    mp = config.get('MasterPath', os.path.abspath('.'))
    assert os.path.exists(mp)
    mfileroot = config.get('Master{}RootName'.format(type),
                           'Master{}_'.format(type))
    mfile = '{}{}.fits'.format(mfileroot, date)
    if os.path.exists(os.path.join(mp, mfile)):
        master = ccdproc.fits_ccddata_reader(os.path.join(mp, mfile), verify=True)
        master.header.set('FILENAME', value=mfile, comment='File name')
    else:
        master = None
    return master


def analyze_image(file,
                  datekw='DATE-OBS', datefmt='%Y-%m-%dT%H:%M:%S'):
    '''
    Run a single file through the processing and analysis steps.
    '''
    start_DR = dt.utcnow()
    metadata = fits.Header()
    metadata.set('DRSTART', value=start_DR.iso(),
                 comment='UT time of start of analysis')
    
    config = get_config()
    try:
        im = ccdproc.fits_ccddata_reader(file, ext=config.get('RawFileExt', 0))
    except ValueError:
        im = ccdproc.fits_ccddata_reader(file, ext=config.get('RawFileExt', 0),
                     unit='adu')

    date = get_image_date(im, datekw=datekw, datefmt=datefmt)

    # Bias correct the image
    master_bias = get_master_file(date, type='Bias')
    if master_bias:
        ccdproc.subtract_bias(im, master_bias, add_keyword=None)
        metadata.set('BIASFILE', value=master_bias.header['FILENAME'],
                 comment='Filename of master bias file')
        metadata.set('BIASCSUM', value=master_bias.header['CHECKSUM'],
                 comment='CHECKSUM of master bias file')
        metadata.set('BIASDSUM', value=master_bias.header['DATASUM'],
                 comment='DATASUM of master bias file')

    # Dark correct the image
    master_dark = get_master_file(date, type='Dark')
    if master_dark:
        ccdproc.subtract_dark(im, master_dark, add_keyword=None)
        metadata.set('DARKFILE', value=master_dark.header['FILENAME'],
                 comment='Filename of master dark file')
        metadata.set('DARKCSUM', value=master_dark.header['CHECKSUM'],
                 comment='CHECKSUM of master dark file')
        metadata.set('DARKDSUM', value=master_dark.header['DATASUM'],
                 comment='DATASUM of master dark file')

    # Gain correct the image
    gain = ccdproc.Keyword('GAIN', unit=u.electron / u.adu)
    try:
        gain.value_from(im.header)
    except KeyError:
        gain.value = config.get('Gain')
    readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
    try:
        readnoise.value_from(im.header)
    except KeyError:
        readnoise.value = config.get('RN')
    im = ccdproc.gain_correct(im, gain)
    im.uncertainty = ccdproc.create_deviation(im, gain=gain, readnoise=readnoise)

    # Shutter correct the image
    shutter_map = get_master_file(date, type='ShutterMap')
    if shutter_map:
        ccdproc.apply_shutter_map(image, shutter_map)

    # Flat correct the image
    master_flat = get_master_file(date, type='Flat')
    if master_flat:
        ccdproc.flat_correct(im, master_flat, add_keyword=None)
        metadata.set('FLATFILE', value=master_flat.header['FILENAME'],
                 comment='Filename of master flat file')
        metadata.set('FLATCSUM', value=master_flat.header['CHECKSUM'],
                 comment='CHECKSUM of master flat file')
        metadata.set('FLATDSUM', value=master_flat.header['DATASUM'],
                 comment='DATASUM of master flat file')



    end_DR = dt.utcnow()
    metadata.set('DREND', value=end_DR.iso(),
                 comment='UT time of start of analysis')

    for hdu in photometry:
        if 'CHECKSUM' not in hdu.header.keys():
            hdu.add_checksum()
    return photometry
