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

import astropy.units as u
import astropy.coordinates as c
from astropy.time import Time
from astropy.io import fits
from astropy import wcs
import ccdproc
import sep

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

from .config import get_config
from .calibration import *
from .utils import *


def preprocess_header(file):
    print('Fixing header')
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
        # Do splitext twice to hanfle .fits.fz extensions
        self.fileroot = os.path.splitext(os.path.splitext(self.filename)[0])[0]

        # Initialize metadata for analysis
        start_DR = dt.utcnow()
        self.metadata = fits.Header()
        self.metadata.set('DRSTART', value=start_DR.isoformat(),
                          comment='UT time of start of analysis')

#         if preprocess_header:
#             preprocess_header(self.file)
        try:
            self.ccd = ccdproc.fits_ccddata_reader(file, verify=True)
        except ValueError:
            self.ccd = ccdproc.fits_ccddata_reader(file, verify=True, unit=unit)

        self.config = get_config()
        self.datekw = self.config.get('DATE-OBS', 'DATE-OBS')
        self.datefmt = self.config.get('DATEFMT', '%Y-%m-%dT%H:%M:%S')
        self.obstime = None
        self.date = self.get_date()
        self.add_logger(verbose=verbose)
        self.log.info('Opened File: {}'.format(self.filename))

        self.header_pointing = None
        self.wcs_pointing = None
        self.altaz = None
        try:
            lat = c.Latitude(self.config.get('Latitude') * u.degree)
            lon = c.Longitude(self.config.get('Longitude') * u.degree)
            height = self.config.get('Elevation') * u.meter
            self.loc = c.EarthLocation(lon, lat, height)
            self.altazframe = c.AltAz(location=self.loc, obstime=self.obstime,
                                temperature=self.config.get('Temperature')*u.Celsius,
                                pressure=self.config.get('Pressure')/1000.*u.bar)
#             self.moon = c.get_moon(self.obstime, location=self.loc)
        except:
            self.loc = None
            self.altazframe = None
            self.moon = None
        self.back = None
    
        self.gain = ccdproc.Keyword('GAIN', unit=u.electron / u.adu)
        try:
            self.gain.value_from(self.ccd.header)
        except KeyError:
            self.gain.value = self.config.get('Gain')
    
        self.exptime = ccdproc.Keyword(self.config.get('EXPTIME', 'EXPTIME'),
                                       unit=u.second)
        self.exptime.value_from(self.ccd.header)
    
        self.readnoise = ccdproc.Keyword('RDNOISE', unit=u.electron)
        try:
            self.readnoise.value_from(self.ccd.header)
        except KeyError:
            self.readnoise.value = self.config.get('RN')

    
    
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
            LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                                          datefmt='%Y%m%d %H:%M:%S')
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
    
    
    def get_date(self):
        '''
        Return string with the UT date of the image in YYYYMMDDUT format.
        Also populates the `obstime` property with an `astropy.time.Time`
        instance.
        '''
        dto = dt.strptime(self.ccd.header[self.datekw], self.datefmt)
        self.obstime = Time(dto)
        self.date = dto.strftime('%Y%m%dUT')
        return self.date
    
    
    def get_header_pointing(self):
        '''
        Read the pointing coordinate from the header RA and DEC keywords
        and populate the `header_pointing` property with an 
        `astropy.coordinates.SkyCoord` instance.
        '''
        self.log.info('Reading pointing from header')
        RAkwd = self.config.get('RA', 'RA')
        DECkwd = self.config.get('DEC', 'DEC')
        equinox = self.config.get('Equinox', 2000.0)
        if type(equinox) is str:
            equinox = self.ccd.header[equinox]
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
                                     equinox=equinox,
                                     obstime=self.obstime,
                                     unit=(u.deg, u.deg))
        elif coord_format == 'decimal hours degrees':
            self.header_pointing = c.SkyCoord(float(RA), float(DEC),
                                     frame=coord_frame,
                                     equinox=equinox,
                                     obstime=self.obstime,
                                     unit=(u.hourangle, u.deg))
        else:
            self.header_pointing = None
            self.log.warning('CoordinateFormat not understood.')

        if self.loc and self.header_pointing:
            self.header_altaz = self.header_pointing.transform_to(self.altazframe)

        if self.header_pointing:
            self.log.info('Header pointing: {}'.format(self.header_pointing.to_string('hmsdms')))
        else:
            self.log.warning('Failed to parse pointing from header')

        return self.header_pointing
    
    
    def solve_astrometry(self, downsample=4, SIPorder=3):
        '''
        Use a local install of astrometry.net to solve the image
        WCS and populate the `ccd.wcs` and `ccd.header` with the
        updated WCS info.
        '''
        self.log.info('Calculating astrometric solution')
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
        self.log.info('  Calling: {}'.format(' '.join(cmd)))
        cmd.append(tfile)
        with open(os.path.join(tdir, 'output.txt'), 'w') as output:
            fail = subprocess.call(cmd, stdout=output, stderr=output)
        with open(os.path.join(tdir, 'output.txt'), 'r') as output:
            lines = output.readlines()
        for line in lines:
            self.log.debug('  {}'.format(line.strip('\n')))
        rootname = os.path.splitext(self.filename)[0]
        solved_file = os.path.join(tdir, '{}.solved'.format(rootname))
        wcs_file = os.path.join(tdir, '{}.wcs'.format(rootname))
        new_wcs = None
        if not bool(fail):
            new_wcs = wcs.WCS(wcs_file)
            if new_wcs.is_celestial:
                self.ccd.wcs = new_wcs
                self.ccd.header.add_history('New WCS solved by astrometry.net')
                new_header = new_wcs.to_header(relax=True)
                for key in new_header.keys():
                    self.ccd.header.set(key, new_header[key], new_header.comments[key])
                self.get_wcs_pointing()
                self.log.info('  Pointing center: {} ({})'.format(
                              self.wcs_pointing.to_string('hmsdms', precision=1),
                              self.wcs_pointing.to_string('decimal', precision=4) ) )
            else:
                new_wcs = None
        # Cleanup temporary directory
        tfiles = glob(os.path.join(tdir, '*'))
        for tf in tfiles:
            os.remove(tf)
        os.rmdir(tdir)
        return new_wcs
    
    
    def get_wcs_pointing(self):
        '''
        Populate the `wcs_pointing` property with and 
        `astropy.coordinates.SkyCoord` instance based on the information
        in the `ccd.wcs` property.
        '''
        equinox = self.config.get('Equinox', 2000.0)
        if equinox in self.ccd.header.keys():
            equinox = self.ccd.header[equinox]
        if not self.ccd.wcs:
            self.wcs_pointing = None
        else:
            coord_frame = self.config.get('CoordinateFrame', 'Fk5')
            nx, ny = self.ccd.data.shape
            r, d = self.ccd.wcs.all_pix2world([nx/2.], [ny/2.], 1)
            self.wcs_pointing = c.SkyCoord(r[0], d[0], frame=coord_frame,
                                           unit=(u.deg, u.deg),
                                           equinox='J2000',
                                           obstime=self.obstime)
        if self.loc and self.wcs_pointing:
            self.wcs_altaz = self.wcs_pointing.transform_to(self.altazframe)
        return self.wcs_pointing
    
    
    def calculate_pointing_error(self):
        '''
        Assming that the `header_pointing` extracted from the RA and DEC
        FITS header keywords represents the intended pointing of the
        image and that the `wcs_pointing` property represents the actual
        pointing, calculate the pointing error magnitude.
        '''
        if not self.header_pointing:
            self.get_header_pointing()
        if not self.wcs_pointing:
            self.get_wcs_pointing()
        if self.header_pointing and self.wcs_pointing:
            sep = self.header_pointing.separation(self.wcs_pointing)
        self.log.info('Pointing error = {:.1f}'.format(sep.to(u.arcmin)))
        return sep
    
    
    def get_UCAC4(self):
        '''
        '''
        self.log.info('Get UCAC4 stars in field of view')
        if not self.header_pointing:
            self.get_header_pointing()
        if not self.wcs_pointing:
            self.get_wcs_pointing()
        if self.wcs_pointing:
            pointing = self.wcs_pointing
        elif im.header_pointing:
            pointing = self.wcs_pointing
        else:
            return None
        from astroquery.vizier import Vizier
        from astropy.coordinates import Angle

        v = Vizier(columns=['_RAJ2000', '_DEJ2000','imag'],
                   column_filters={"imag":">0"})
        v.ROW_LIMIT = 1e5

        fp = self.ccd.wcs.calc_footprint(axes=self.ccd.data.shape)
        dra = fp[:,0].max() - fp[:,0].min()
        ddec = fp[:,1].max() - fp[:,1].min()
        radius = np.sqrt( (dra*np.cos(fp[:,1].mean()*np.pi/180.))**2 + ddec**2 )/2.
        catalog = v.query_region(self.header_pointing, catalog='I/322A',
                                 radius=Angle(radius, "deg") )[0]

        if self.ccd.wcs:
            x, y = self.ccd.wcs.all_world2pix(catalog['_RAJ2000'], catalog['_DEJ2000'], 1)
            w = (x > 0) & (x < 4096) & (y > 0) & (y < 4096)
            catalog = catalog[w]
        self.log.info('  Found {:d} catalog stars'.format(len(catalog)))
        return catalog
    
    
    def subtract_background(self):
        '''
        '''
        self.log.info('Subtracting Background')
        sbc = self.config.get('Background', {})
        bkg = sep.Background(self.ccd.data, mask=self.ccd.mask,
                             bw=sbc.get('bw', 64),
                             bh=sbc.get('bh', 64),
                             fw=sbc.get('fw', 3),
                             fh=sbc.get('fh', 3),
                             fthresh=sbc.get('fthresh', 0))
        self.back = bkg
        self.ccd.data -= bkg.back()
        self.ccd.header.add_history('Background subtracted using SEP')
        self.ccd.header.add_history('  bw={:d}, bh={:d}, fw={:d}, fh={:d}'.format(
                                    sbc.get('bw', 64), sbc.get('bh', 64),
                                    sbc.get('fw', 3), sbc.get('fh', 3) ))
        self.ccd.header.add_history('  fthresh={:f}'.format(sbc.get('fthresh', 0)))
    
    
    def extract(self):
        '''
        '''
        self.log.info('Extracting sources')
        extract_config = self.config.get('Extract', {})
        thresh = extract_config.get('thresh', 5)
        minarea = extract_config.get('minarea', 5)
        objects = sep.extract(self.ccd.data, err=self.ccd.uncertainty, mask=self.ccd.mask,
                              thresh=thresh, minarea=minarea)
        self.log.info('  Found {:d} sources'.format(len(objects)))
        return objects
    
    
    def render_jpeg(self, jpegfilename=None, binning=1,
                    overplot_UCAC4=False):
        '''
        '''
        self.log.info('Render JPEG of image to {}'.format(jpegfilename))
        if not jpegfilename:
            jpegfilename = '{}.jpg'.format(self.fileroot)
        import matplotlib as mpl
        from matplotlib import pyplot as plt
        vmin = np.percentile(self.ccd.data, 0.5)
        vmax = np.percentile(self.ccd.data, 99.0)
        dpi=72
        nx, ny = self.ccd.data.shape
        sx = nx/dpi/binning
        sy = ny/dpi/binning
        fig = plt.figure(figsize=(sx, sy), dpi=dpi)
        ax = fig.gca()
        plt.imshow(self.ccd.data, cmap='gray', vmin=vmin, vmax=vmax)
        plt.xticks([])
        plt.yticks([])

        if overplot_UCAC4:
            self.log.info('  Overlaying UCAC4 catalog stars')
            UCAC4 = self.get_UCAC4()
            x, y = self.ccd.wcs.all_world2pix(UCAC4['_RAJ2000'], UCAC4['_DEJ2000'], 1)
            for xy in zip(x, y):
                c = plt.Circle(xy, radius=5, edgecolor='r', facecolor='none')
                ax.add_artist(c)

        plt.savefig(jpegfilename, dpi=dpi)
    
    
    def bias_correct(self):
        '''
        '''
        master_bias = get_master(self.date, type='Bias')
        if master_bias:
            self.log.info('Subtracting master bias')
            self.ccd = ccdproc.subtract_bias(self.ccd, master_bias,
                               add_keyword=None)
            self.metadata.set('BIASFILE', value=master_bias.header['FILENAME'],
                              comment='Filename of master bias file')
            self.metadata.set('BIASCSUM', value=master_bias.header['CHECKSUM'],
                              comment='CHECKSUM of master bias file')
            self.metadata.set('BIASDSUM', value=master_bias.header['DATASUM'],
                              comment='DATASUM of master bias file')
    
    
    def gain_correct(self):
        '''
        '''
        self.log.info('Gain correcting image')
        self.ccd = ccdproc.gain_correct(self.ccd, self.gain.value)
    
    
    def dark_correct(self):
        '''
        '''
        master_dark = get_master(self.date, type='Dark')
        if master_dark:
            self.log.info('Dark correcting image')
            self.ccd = ccdproc.subtract_dark(self.ccd, master_dark,
                           data_exposure=self.exptime.value,
                           dark_exposure=1.0*u.second,
                           scale=True,
                           add_keyword=None)
            self.metadata.set('DARKFILE', value=master_dark.header['FILENAME'],
                              comment='Filename of master dark file')
            self.metadata.set('DARKCSUM', value=master_dark.header['CHECKSUM'],
                              comment='CHECKSUM of master dark file')
            self.metadata.set('DARKDSUM', value=master_dark.header['DATASUM'],
                              comment='DATASUM of master dark file')
    
    
    def shutter_correct(self):
        '''
        '''
        shutter_map = get_master_shutter_map(self.date)
        if shutter_map:
            self.log.info('Applying shutter map')
            self.ccd = ccdproc.apply_shutter_map(self.ccd, shutter_map)
            self.metadata.set('SHUTFILE', value=shutter_map.header['FILENAME'],
                          comment='Filename of master shutter correction file')
            self.metadata.set('SHUTCSUM', value=shutter_map.header['CHECKSUM'],
                          comment='CHECKSUM of master shutter correction file')
            self.metadata.set('SHUTDSUM', value=shutter_map.header['DATASUM'],
                          comment='DATASUM of master shutter correction file')
    
    
    def flat_correct(self):
        '''
        '''
        master_flat = get_master(self.date, type='Flat')
        if master_flat:
            self.log.info('Flat fielding image')
            self.ccd = ccdproc.flat_correct(im, master_flat, add_keyword=None)
            self.metadata.set('FLATFILE', value=master_flat.header['FILENAME'],
                              comment='Filename of master flat file')
            self.metadata.set('FLATCSUM', value=master_flat.header['CHECKSUM'],
                              comment='CHECKSUM of master flat file')
            self.metadata.set('FLATDSUM', value=master_flat.header['DATASUM'],
                              comment='DATASUM of master flat file')


