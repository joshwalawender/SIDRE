import os
import sys
import re
import tempfile
import logging
from glob import glob
from datetime import datetime as dt
if sys.version_info.major == 2:
    import subprocess32 as subprocess
elif sys.version_info.major == 3:
    import subprocess

from astropy.io import fits
from astropy import wcs
import astropy.units as u
import astropy.coordinates as c
import ccdproc

from .config import get_config
from .calibration import *
from .utils import *


def preprocess_header(file):
    with fits.open(file, 'update', verify=True) as hdul:
        updated = False
        if not 'RADESYSa' in hdul[0].header.keys()\
           and 'RADECSYS' in hdul[0].header.keys():
            hdul[0].header.set('RADESYSa', hdul[0].header['RADECSYS'])
            updated = True
        if not 'BUNIT' in hdul[0].header.keys():
            hdul[0].header.set('BUNIT', 'adu')
            updated = True
        else:
            if hdul[0].header['BUNIT'] == 'ADU':
                hdul[0].header.set('BUNIT', 'adu')
                updated = True
        if updated:
            hdul.flush()


class ScienceImage(object):
    '''
    '''
    def __init__(self, file, preprocess_header=preprocess_header,
                 logfile=None, verbose=False, unit='adu'):
        assert os.path.exists(file)
        self.file = file
        self.filename = os.path.basename(file)
        self.fileroot = os.path.splitext(self.filename)[0]
        if preprocess_header:
            preprocess_header
        try:
            self.ccd = ccdproc.fits_ccddata_reader(file, verify=True)
        except ValueError:
            self.ccd = ccdproc.fits_ccddata_reader(file, verify=True, unit=unit)

        self.config = get_config()
        self.datekw = self.config.get('DATE-OBS', 'DATE-OBS')
        self.datefmt = self.config.get('DATEFMT', '%Y-%m-%dT%H:%M:%S')
        self.image_date = None
        self.date = self.get_image_date()
        self.add_logger()
        self.log.info('Processing File: {}'.format(self.filename))

        self.header_pointing = None


    def add_logger(self, logfile=None, verbose=False):
        '''Create a logger object to use with this image.  The logger object
        will be available as self.log
        
        Parameters
        ----------
        logfile : file to write log to
        
        verbose : Defaults to False.  If verbose is true, it sets the logging
            level to DEBUG (otherwise level is INFO).
        '''
        self.log = logging.getLogger(self.filename.replace('.', '_'))
        if len(self.log.handlers) == 0:
            self.log.setLevel(logging.DEBUG)
            LogFormat = logging.Formatter('%(asctime)23s %(levelname)8s: %(message)s')
            ## Log to a file
            if logfile:
                self.logfile = logfile
                self.logfilename = os.path.split(self.logfile)[1]
                if os.path.exists(logfile): os.remove(logfile)
                LogFileHandler = logging.FileHandler(logfile)
                LogFileHandler.setLevel(logging.DEBUG)
                LogFileHandler.setFormatter(LogFormat)
                self.log.addHandler(LogFileHandler)
            ## Log to console
            LogConsoleHandler = logging.StreamHandler(stream=sys.stdout)
            if verbose:
                LogConsoleHandler.setLevel(logging.DEBUG)
            else:
                LogConsoleHandler.setLevel(logging.INFO)
            LogConsoleHandler.setFormatter(LogFormat)
            self.log.addHandler(LogConsoleHandler)


    def get_image_date(self):
        '''
        Return string with the UT date of the image in YYYYMMDDUT format
        '''
        self.image_date = dt.strptime(self.ccd.header[self.datekw], self.datefmt)
        self.date = self.image_date.strftime('%Y%m%dUT')
        return self.date



    def analyze_image(self):
        '''
        Run the image through the processing and analysis steps.
        '''
        start_DR = dt.utcnow()
        metadata = fits.Header()
        metadata.set('DRSTART', value=start_DR.isoformat(),
                     comment='UT time of start of analysis')
    
    
        # Bias correct the image
        master_bias = get_master(self.date, type='Bias')
        if master_bias:
            self.log.info('Subtracting master bias')
            self.ccd = ccdproc.subtract_bias(self.ccd, master_bias, add_keyword=None)
            metadata.set('BIASFILE', value=master_bias.header['FILENAME'],
                     comment='Filename of master bias file')
            metadata.set('BIASCSUM', value=master_bias.header['CHECKSUM'],
                     comment='CHECKSUM of master bias file')
            metadata.set('BIASDSUM', value=master_bias.header['DATASUM'],
                     comment='DATASUM of master bias file')

        # Gain correct the image
        gain = ccdproc.Keyword('GAIN', unit=u.electron / u.adu)
        try:
            gain.value_from(self.ccd.header)
        except KeyError:
            gain.value = self.config.get('Gain')

        self.log.info('Gain correcting image')
        self.ccd = ccdproc.gain_correct(self.ccd, gain.value)

        # Dark correct the image
        exptime = ccdproc.Keyword(self.config.get('EXPTIME', 'EXPTIME'),
                                  unit=u.second)
        exptime.value_from(self.ccd.header)
        master_dark = get_master(self.date, type='Dark')
        if master_dark:
            self.log.info('Dark correcting image')
            self.ccd = ccdproc.subtract_dark(self.ccd, master_dark,
                           data_exposure=exptime.value,
                           dark_exposure=1.0*u.second,
                           scale=True,
                           add_keyword=None)
            metadata.set('DARKFILE', value=master_dark.header['FILENAME'],
                     comment='Filename of master dark file')
            metadata.set('DARKCSUM', value=master_dark.header['CHECKSUM'],
                     comment='CHECKSUM of master dark file')
            metadata.set('DARKDSUM', value=master_dark.header['DATASUM'],
                     comment='DATASUM of master dark file')

#         readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
#         try:
#             readnoise.value_from(self.ccd.header)
#         except KeyError:
#             readnoise.value = self.config.get('RN')
#         self.log.info('Estimating uncertainty')
#         self.ccd = ccdproc.create_deviation(self.ccd, readnoise=readnoise.value)

        # Shutter correct the image
        shutter_map = get_master_shutter_map(self.date)
        if shutter_map:
            self.log.info('Applying shutter map')
            self.ccd = ccdproc.apply_shutter_map(self.ccd, shutter_map)

        # Flat correct the image
        master_flat = get_master(self.date, type='Flat')
        if master_flat:
            self.log.info('Flat fielding image')
            self.ccd = ccdproc.flat_correct(im, master_flat, add_keyword=None)
            metadata.set('FLATFILE', value=master_flat.header['FILENAME'],
                     comment='Filename of master flat file')
            metadata.set('FLATCSUM', value=master_flat.header['CHECKSUM'],
                     comment='CHECKSUM of master flat file')
            metadata.set('FLATDSUM', value=master_flat.header['DATASUM'],
                     comment='DATASUM of master flat file')

        end_DR = dt.utcnow()
        metadata.set('DREND', value=end_DR.isoformat(),
                     comment='UT time of start of analysis')

        for key in metadata.keys():
            print(key, metadata[key])

    def get_header_pointing(self):
        '''
        Read the pointing coordinate from the header RA and DEC keywords.
        '''
        RAkwd = self.config.get('RA', 'RA')
        DECkwd = self.config.get('DEC', 'DEC')
        coord_frame = self.config.get('CoordinateFrame', 'Fk5')
        coord_format = self.config.get('CoordinateFormat', 'HMSDMS')
        RA = self.ccd.header[RAkwd]
        DEC = self.ccd.header[DECkwd]

        if coord_format == 'HMSDMS':
            self.header_pointing = c.SkyCoord(RA, DEC, frame=coord_frame,
                                     unit=(u.hourangle, u.deg))
        elif coord_format == 'decimal degrees':
            self.header_pointing = c.SkyCoord(float(RA), float(DEC),
                                     frame=coord_frame,
                                     unit=(u.deg, u.deg))
        elif coord_format == 'decimal hours degrees':
            self.header_pointing = c.SkyCoord(float(RA), float(DEC),
                                     frame=coord_frame,
                                     unit=(u.hourangle, u.deg))
        else:
            self.header_pointing = None
            self.log.warning('CoordinateFormat not understood.')


    def render_jpeg(self, jpegfilename=None):
        '''
        '''
        if not jpegfilename:
            jpegfilename = '{}.jpg'.format(os.path.basename(self.filename))
        
        
    def solve_astrometry(self, downsample=4, SIPorder=2):
        '''
        
        '''
        solvefield_args = self.config.get('SolveFieldArgs', [])
        solvefield_args.extend(['-z', '{:d}'.format(downsample)])
        solvefield_args.extend(['-t', '{:d}'.format(SIPorder)])
        
        # Create temporary directory
        tdir = tempfile.mkdtemp()
        tfile = os.path.join(tdir, self.filename)
        ccdproc.fits_ccddata_writer(self.ccd, tfile)
        
        # run astrometry.net on the temporary fits file
        cmd = ['solve-field', '-p']
        cmd.extend(solvefield_args)
        if not self.header_pointing:
            self.get_header_pointing()
        if self.header_pointing:
            cmd.extend(['-3', '{:.4f}'.format(self.header_pointing.ra.degree)])
            cmd.extend(['-4', '{:.4f}'.format(self.header_pointing.dec.degree)])
        cmd.append(tfile)
        fail = subprocess.call(cmd)
        solved_file = os.path.join(tdir, '{}.solved'.format(self.fileroot))
        wcs_file = os.path.join(tdir, '{}.wcs'.format(self.fileroot))
        if not fail and os.path.exists(solved_file) and os.path.exists(wcs_file):
            new_wcs = wcs.WCS(wcs_file)
            if new_wcs.is_celestial:
                self.ccd.wcs = new_wcs
                self.ccd.header.add_history('New WCS solved by astrometry.net')
                new_header = new_wcs.to_header(relax=True)
                for key in new_header.keys():
                    self.ccd.header.set(key, new_header[key], new_header.comments[key])
        # Cleanup temporary directory
        tfiles = glob(os.path.join(tdir, '*'))
        for tf in tfiles:
            os.remove(tf)
        os.rmdir(tdir)
        
        
