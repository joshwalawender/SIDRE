import os
import re
from datetime import datetime as dt
from glob import glob

import astropy.units as u
from astropy import stats
import ccdproc

from .config import get_config


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


def mode(data):
    '''
    Return mode of image.  Assumes int values (ADU), so uses binsize of one.
    '''
    bmin = np.floor(min(data.ravel())) - 1./2.
    bmax = np.ceil(max(data.ravel())) + 1./2.
    bins = np.arange(bmin,bmax,1)
    hist, bins = np.histogram(data.ravel(), bins=bins)
    centers = (bins[:-1] + bins[1:]) / 2
    w = np.argmax(hist)
    mode = int(centers[w])
    return mode


def imstat(file, ext=0, sigma=5, iters=1):
    '''
    '''
    assert os.path.exists(file)
    with fits.open(file, 'readonly') as hdul:
        mean, median, stddev = stats.sigma_clipped_stats(hdul[ext].data,
                                     sigma=sigma,
                                     iters=iters)
        mode = get_mode(hdul[ext].data)
        print(mean, median, mode, stddev)

