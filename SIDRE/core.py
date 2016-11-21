import os
import re
from glob import glob
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


def get_master(date, type='Bias'):
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


def get_master_shutter_map(date):
    '''
    '''
    date_dto = dt.strptime(date, '%Y%m%dUT')
    config = get_config()
    mp = config.get('MasterPath', os.path.abspath('.'))
    assert os.path.exists(mp)
    mfileroot = config.get('MasterShutterMapRootName', 'ShutterMap_')
    # Look for files with date name nearest to date being analyzed
    shutter_map_files = glob(os.path.join(mp, '{}*.fits'.format(mfileroot)))
    if len(shutter_map_files) < 1:
        mfile = None
    elif len(shutter_map_files) == 1:
        mfile = shutter_map_files[0]
    else:
        dates = []
        for file in shutter_map_files:
            match = re.match('{}(\d{8}UT)\.fits', file)
            if match:
                filedate = dt.strptime(match.group(1), '%Y%m%dUT')
                timediff = abs((date_dto - filedate).total_seconds())
                dates.append((timediff, file))
        dates.sort()
        mfile = dates[0][1]
    if mfile:
        shutter_map = ccdproc.fits_ccddata_reader(os.path.join(mp, mfile), verify=True)
        shutter_map.header.set('FILENAME', value=mfile, comment='File name')
    return shutter_map


def analyze_image(file,
                  datekw='DATE-OBS', datefmt='%Y-%m-%dT%H:%M:%S'):
    '''
    Run a single file through the processing and analysis steps.
    '''
    start_DR = dt.utcnow()
    metadata = fits.Header()
    metadata.set('DRSTART', value=start_DR.isoformat(),
                 comment='UT time of start of analysis')
    
    config = get_config()
    try:
        im = ccdproc.fits_ccddata_reader(file, ext=config.get('RawFileExt', 0))
    except ValueError:
        im = ccdproc.fits_ccddata_reader(file, ext=config.get('RawFileExt', 0),
                     unit='adu')

    date = get_image_date(im, datekw=datekw, datefmt=datefmt)

    # Bias correct the image
    master_bias = get_master(date, type='Bias')
    if master_bias:
        ccdproc.subtract_bias(im, master_bias, add_keyword=None)
        metadata.set('BIASFILE', value=master_bias.header['FILENAME'],
                 comment='Filename of master bias file')
        metadata.set('BIASCSUM', value=master_bias.header['CHECKSUM'],
                 comment='CHECKSUM of master bias file')
        metadata.set('BIASDSUM', value=master_bias.header['DATASUM'],
                 comment='DATASUM of master bias file')

    # Dark correct the image
    master_dark = get_master(date, type='Dark')
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
    im = ccdproc.gain_correct(im, gain.value)
    im = ccdproc.create_deviation(im, readnoise=readnoise.value)

    # Shutter correct the image
    shutter_map = get_master_shutter_map(date)
    if shutter_map:
        ccdproc.apply_shutter_map(im, shutter_map)

    # Flat correct the image
    master_flat = get_master(date, type='Flat')
    if master_flat:
        ccdproc.flat_correct(im, master_flat, add_keyword=None)
        metadata.set('FLATFILE', value=master_flat.header['FILENAME'],
                 comment='Filename of master flat file')
        metadata.set('FLATCSUM', value=master_flat.header['CHECKSUM'],
                 comment='CHECKSUM of master flat file')
        metadata.set('FLATDSUM', value=master_flat.header['DATASUM'],
                 comment='DATASUM of master flat file')





    end_DR = dt.utcnow()
    metadata.set('DREND', value=end_DR.isoformat(),
                 comment='UT time of start of analysis')

    print(metadata)