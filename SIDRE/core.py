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

import matplotlib as mpl
from matplotlib import pyplot as plt

import astropy.units as u
import astropy.coordinates as c
from astropy.time import Time
from astropy.io import fits
from astropy import wcs
from astropy.table import Table, Column
from astropy import stats
import ccdproc
import sep
import photutils as phot

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
        # Convert RADECSYS to RADESYSa
        if not 'RADESYSa' in hdul[0].header.keys()\
           and 'RADECSYS' in hdul[0].header.keys():
            hdul[0].header.set('RADESYSa', hdul[0].header['RADECSYS'])
            updated = True
        # Add BUNIT
        if not 'BUNIT' in hdul[0].header.keys():
            hdul[0].header.set('BUNIT', 'adu')
            updated = True
        # If BUNIT is ADU convert to adu
        elif hdul[0].header['BUNIT'] == 'ADU':
            hdul[0].header.set('BUNIT', 'adu')
            updated = True

        if updated:
            hdul.flush()


class ScienceImage(object):
    '''
    '''
    def __init__(self, file, preprocess_header=preprocess_header,
                 logfile=None, verbose=False, unit='adu'):
        self.log = None
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

        # Preprocess header if needed
        if preprocess_header:
            preprocess_header(self.file)

        # Read FITS file
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
        self.source_mask = None
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
        self.UCAC4 = None
        self.extracted = None
        self.assoc = None
    
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


    def __del__(self):
        if self.log:
            self.log.info('Done.')


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

    #--------------------------------------------------------------------------
    # Image Reduction Methods
    #--------------------------------------------------------------------------
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
        Wrapper around `ccdproc.gain_correct` function.
        '''
        self.log.info('Gain correcting image')
        self.ccd = ccdproc.gain_correct(self.ccd, self.gain.value)


    def create_deviation(self):
        '''
        Wrapper around `ccdproc.create_deviation` function.
        '''
        self.log.info('Estimating Uncertainty in Image')
        if self.ccd.unit == u.electron:
            self.ccd = ccdproc.create_deviation(self.ccd, #gain=self.gain.value,
                               readnoise=self.readnoise.value)
        elif self.ccd.unit == u.adu:
            self.ccd = ccdproc.create_deviation(self.ccd, gain=self.gain.value,
                               readnoise=self.readnoise.value)
        else:
            self.log.error('Coud not estimate uncertainty')


    def dark_correct(self):
        '''
        Wrapper around `ccdproc.subtract_dark` function.
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
        Wrapper around `ccdproc.apply_shutter_map` function.
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
        Wrapper around `ccdproc.flat_correct` function.
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

    #--------------------------------------------------------------------------
    # WCS Methods
    #--------------------------------------------------------------------------
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
            self.log.warning('  CoordinateFormat not understood.')

        if self.loc and self.header_pointing:
            self.header_altaz = self.header_pointing.transform_to(self.altazframe)

        if self.header_pointing:
            self.log.info('  Header pointing: {}'.format(
                          self.header_pointing.to_string('hmsdms')))
        else:
            self.log.warning('  Failed to parse pointing from header')

        return self.header_pointing


    def solve_astrometry(self, downsample=2, SIPorder=4):
        '''
        Use a local install of astrometry.net to solve the image
        WCS and populate the `ccd.wcs` and `ccd.header` with the
        updated WCS info.
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
        self.log.info('Calculating astrometric solution')
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
                    self.ccd.header.set(key, new_header[key],
                                        new_header.comments[key])
                self.get_wcs_pointing()
                self.log.info('  Pointing center: {} ({})'.format(
                              self.wcs_pointing.to_string('hmsdms', precision=1),
                              self.wcs_pointing.to_string('decimal', precision=4)))
            else:
                new_wcs = None
        else:
            self.log.warning('solve-field failed')
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


    def test_astrometry(self, nradii=10, min_best_match=0.10, min_regional_match=0.50):
        '''
        '''
        assert self.assoc
        pass_astrometry = True
        ny, nx = self.ccd.shape
        radius = np.sqrt((nx/2.)**2 + (ny/2.)**2)
        radii = np.linspace(0, radius, nradii+1)
        frac = np.zeros(nradii)
        for i in range(nradii):
            n_assoc = len(self.assoc[(self.assoc['r'] > radii[i])\
                                   & (self.assoc['r'] <= radii[i+1])])
            n_detected = len(self.extracted[(self.extracted['r'] > radii[i])\
                                          & (self.extracted['r'] <= radii[i+1])])
            frac[i] = float(n_assoc)/float(n_detected)
        
        AQM = frac/max(frac)
        self.log.info('Astrometery verificaton metrics:')
        best_match = max(frac)
        self.log.info('  Best match = {:.1f} %'.format(best_match*100.))
        if best_match < min_best_match:
            pass_astrometry = False
            self.log.warning('Astrometry test failed.')
            self.log.warning('Best region absolute success rate < {:.1f} %'.format(min_best_match*100.))
        
        for i,val in enumerate(AQM):
            self.log.info(' AQM ({:d}) = {:.2f}'.format(i, val))
        if min(AQM) < min_regional_match:
            pass_astrometry = False
            self.log.warning('Astrometry test failed.')
            self.log.warning('Worst region relative success rate < {:.1f} %'.format(min_regional_match*100.))

        return pass_astrometry


    #--------------------------------------------------------------------------
    # Source Extractor Methods
    #--------------------------------------------------------------------------
    def sx_subtract_background(self):
        '''
        Use `sep.Background` to generate a background model and then
        subtract it from the image.
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
        self.ccd.header.add_history(' bw={:d},bh={:d},fw={:d},fh={:d}'.format(
                                    sbc.get('bw', 64), sbc.get('bh', 64),
                                    sbc.get('fw', 3), sbc.get('fh', 3) ))
        self.ccd.header.add_history(' fthresh={:f}'.format(sbc.get('fthresh',0)))


    def extract(self):
        '''
        Wrapper around the `sep.extract` function.
        '''
        self.log.info('Extracting sources')
        extract_config = self.config.get('Extract', {})
        thresh = extract_config.get('thresh', 5)
        minarea = extract_config.get('minarea', 5)
        if self.ccd.uncertainty:
            objects = sep.extract(self.ccd.data, err=self.ccd.uncertainty.array,
                                  mask=self.ccd.mask,
                                  thresh=float(thresh), minarea=minarea)
        else:
            objects = sep.extract(self.ccd.data,
                                  mask=self.ccd.mask,
                                  thresh=float(thresh), minarea=minarea)

        self.log.info('  Found {:d} sources'.format(len(objects)))
        self.extracted = Table(objects)

        next = len(self.extracted)
        self.extracted.add_column(Column(data=[None]*next, name='RA', dtype=np.float))
        self.extracted.add_column(Column(data=[None]*next, name='Dec', dtype=np.float))
        self.extracted.add_column(Column(data=[None]*next, name='mag', dtype=np.float))

        ny, nx = self.ccd.shape
        r = np.sqrt((self.extracted['x']-nx/2.)**2 + (self.extracted['y']-ny/2.)**2)
        self.extracted.add_column(Column(data=r.data, name='r', dtype=np.float))

        return self.extracted


    def associate(self, input, magkey='rmag'):
        '''
        Associate entries from the input stellar catalog (e.g. UCAC4) with
        the results from the `extract` method using a simple nearest neighbor
        algorithm.
        '''
        if '_RAJ2000' in input.keys():
            rakey = '_RAJ2000'
        elif 'RA' in input.keys():
            rakey = 'RA'
        else:
            rakey = None

        if '_DEJ2000' in input.keys():
            deckey = '_DEJ2000'
        elif 'DEC' in input.keys():
            deckey = 'DEC'
        else:
            deckey = None
    
        if not rakey or not deckey:
            self.log.warning('Could not parse input catalog table')
            return None
    
        assocconf = self.config.get('Assoc')
        assoc_r = assocconf.get('radius', 3) # pixels
        
        ny, nx = self.ccd.shape
        for i,cstar in enumerate(input):
            x, y = self.ccd.wcs.all_world2pix(cstar[rakey], cstar[deckey], 1)
            dist = np.sqrt( (self.extracted['x']-x)**2 + (self.extracted['y']-y)**2 )
            ncandidates = len(dist[dist < assoc_r])
            if ncandidates > 0:
                id = dist.argmin()
                self.extracted[id]['RA'] = cstar[rakey]
                self.extracted[id]['Dec'] = cstar[deckey]
                self.extracted[id]['mag'] = cstar[magkey]
        self.assoc = self.extracted[~np.isnan(self.extracted['RA'])]
        self.log.info('  Associated {:d} extracted sources with catalog'.format(
                      len(self.assoc)))
        return input


    #--------------------------------------------------------------------------
    # Photometry Methods
    #--------------------------------------------------------------------------
    def make_source_mask(self, snr=5, npixels=5, **kwargs):
        self.log.info('Masking sources in image')
        self.source_mask = phot.make_source_mask(self.ccd.data, snr, npixels, **kwargs)
        return self.source_mask


    def subtract_background(self, box_size=128):
        '''
        '''
        self.log.info('Subtracting Background')
        if self.source_mask is not None and self.ccd.mask is not None:
            mask = (self.source_mask | self.ccd.mask)
        elif self.source_mask is not None:
            mask = self.source_mask
        elif self.ccd.mask is not None:
            mask = self.ccd.mask
        else:
            mask = None
        bkg = phot.Background2D(self.ccd.data, box_size=box_size, mask=mask,
                                sigma_clip=phot.SigmaClip())
        self.back = bkg
        self.ccd.data -= bkg.background
        self.ccd.header.add_history('Background subtracted using photutils')


    def find_stars(self):
        '''
        '''
        self.log.info('Finding stars in image')
        mean, median, std = stats.sigma_clipped_stats(self.ccd.data, sigma=3.0, iters=5)
        star_finder = phot.DAOStarFinder(fwhm=2.0, threshold=5.*std)
        assert self.back
        sources = star_finder(self.ccd.data)
        return sources


    def measure_stars(self, x, y, radius=5):
        '''
        Iterate through a set of catalog stars and perform aperture photometry
        on each star.
        '''
        assert self.back
        self.log.info('Performing circular aperture photometry (r={:.1f} pix) on {:d} stars'.format(
                      radius, len(x)))
        positions = zip(x, y)
        capertures = phot.CircularAperture(positions, r=radius)
        phot_table = phot.aperture_photometry(self.ccd.data, capertures,
                          error=self.ccd.uncertainty.array)
        return phot_table





    def calculate_zero_point(self, plot=None):
        '''
        Estimate the photometric zero point of the image using the associated
        catalog.  Find the mean difference between instrumental magnitude and
        catalog magnitude.
        '''
        if not self.assoc:
            return None
        instmag = -2.512*np.log10(self.assoc['flux'])
        diffs = instmag - self.assoc['mag']
        zp, zp_med, zp_std = stats.sigma_clipped_stats(diffs, sigma=5.0, iters=3)
        print(zp, zp_med, zp_std, zp_std/len(diffs)**0.5)

        ## Make plots
        if plot:
            fig = plt.figure(figsize=(10,10), dpi=72)
            plt.plot(self.assoc['mag'], diffs, 'bo')
            plt.xlabel('Catalog Magnitude')
            plt.ylabel('Zero Point')
            plt.ylim([-31, -22])
            plt.xlim([7, 16])
            
            plt.savefig(plot, dpi=72)




    #--------------------------------------------------------------------------
    # Other Methods
    #--------------------------------------------------------------------------
    def get_UCAC4(self):
        '''
        Use `astroquery` to get the UCAC4 catalog for this image from the
        Vizier service.
        '''
        self.log.info('Get UCAC4 stars in field of view')
        if not self.header_pointing:
            self.get_header_pointing()
        if not self.wcs_pointing:
            self.get_wcs_pointing()
        if self.wcs_pointing:
            pointing = self.wcs_pointing
        elif self.header_pointing:
            pointing = self.wcs_pointing
        else:
            return None
        from astroquery.vizier import Vizier
        from astropy.coordinates import Angle

        v = Vizier(columns=['_RAJ2000', '_DEJ2000', 'rmag'],
                   column_filters={"rmag":">0"})
        v.ROW_LIMIT = 1e5

        fp = self.ccd.wcs.calc_footprint(axes=self.ccd.data.shape)
        dra = fp[:,0].max() - fp[:,0].min()
        ddec = fp[:,1].max() - fp[:,1].min()
        radius = np.sqrt((dra*np.cos(fp[:,1].mean()*np.pi/180.))**2 + ddec**2)/2.
        catalog = v.query_region(self.header_pointing, catalog='I/322A',
                                 radius=Angle(radius, "deg") )[0]

        if self.ccd.wcs:
            x, y = self.ccd.wcs.all_world2pix(catalog['_RAJ2000'],
                                              catalog['_DEJ2000'], 1)
            w = (x > 0) & (x < 4096) & (y > 0) & (y < 4096)
            catalog = catalog[w]
        self.log.info('  Found {:d} catalog stars'.format(len(catalog)))
        catalog.add_column(Column(catalog['_RAJ2000'], name='RA'))
        catalog.add_column(Column(catalog['_DEJ2000'], name='DEC'))
        self.UCAC4 = catalog
        return catalog


    def render_jpeg(self, jpegfilename=None, binning=1, radius=6,
                    overplot_UCAC4=False, overplot_extracted=False,
                    overplot_assoc=False, overplot_pointing=False,
                    overplot_catalog=None):
        '''
        Render a jpeg of the image with optional overlays.
        '''
        self.log.info('Render JPEG of image to {}'.format(jpegfilename))
        if not jpegfilename:
            jpegfilename = '{}.jpg'.format(self.fileroot)
        vmin = np.percentile(self.ccd.data, 0.5)
        vmax = np.percentile(self.ccd.data, 99.0)
        dpi=72
        nx, ny = self.ccd.data.shape
        sx = nx/dpi/binning
        sy = ny/dpi/binning
        fig = plt.figure(figsize=(sx, sy), dpi=dpi)
        ax = fig.gca()
        mdata = np.ma.MaskedArray(self.ccd.data, mask=self.ccd.mask)
        palette = plt.cm.gray
        palette.set_bad('r', 1.0)
        plt.imshow(mdata, cmap=palette, vmin=vmin, vmax=vmax)
        plt.xticks([])
        plt.yticks([])

        if overplot_catalog:
            self.log.info('  Overlaying catalog stars')
            x, y = self.ccd.wcs.all_world2pix(overplot_catalog['RA'],
                                              overplot_catalog['DEC'], 1)
            for xy in zip(x, y):
                c = plt.Circle(xy, radius=radius, edgecolor='b', facecolor='none')
                ax.add_artist(c)

        if overplot_UCAC4:
            self.log.info('  Overlaying UCAC4 catalog stars')
            if not self.UCAC4:
                self.get_UCAC4()
            x, y = self.ccd.wcs.all_world2pix(self.UCAC4['_RAJ2000'],
                                              self.UCAC4['_DEJ2000'], 1)
            for xy in zip(x, y):
                c = plt.Circle(xy, radius=radius, edgecolor='b', facecolor='none')
                ax.add_artist(c)

        if overplot_extracted:
            self.log.info('  Overlaying extracted stars')
            x, y = self.ccd.wcs.all_world2pix(self.extracted['RA'],
                                              self.extracted['Dec'], 1)
            for xy in zip(x, y):
                c = plt.Circle(xy, radius=radius, edgecolor='r', facecolor='none')
                ax.add_artist(c)

        if overplot_assoc:
            self.log.info('  Overlaying associated stars')
            x, y = self.ccd.wcs.all_world2pix(self.assoc['RA'],
                                              self.assoc['Dec'], 1)
            for xy in zip(x, y):
                c = plt.Circle(xy, radius=radius, edgecolor='g', facecolor='none')
                ax.add_artist(c)

        if overplot_pointing:
            if not self.ccd.wcs.is_celestial:
                return None
            if not self.header_pointing:
                self.get_header_pointing()
            if self.header_pointing:
                x, y = self.ccd.wcs.all_world2pix(self.header_pointing.ra.degree,
                                                  self.header_pointing.dec.degree, 1)
                plt.plot([0,nx], [ny/2,ny/2], 'y-', alpha=0.7)
                plt.plot([nx/2, nx/2], [0,ny], 'y-', alpha=0.7)
                # Draw crosshair on target
                ms = radius*6
                c = plt.Circle((x, y), radius=ms, edgecolor='g', alpha=0.7,
                               facecolor='none')
                ax.add_artist(c)
                plt.plot([x, x], [y+0.6*ms, y+1.4*ms], 'g', alpha=0.7)
                plt.plot([x, x], [y-0.6*ms, y-1.4*ms], 'g', alpha=0.7)
                plt.plot([x-0.6*ms, x-1.4*ms], [y, y], 'g', alpha=0.7)
                plt.plot([x+0.6*ms, x+1.4*ms], [y, y], 'g', alpha=0.7)
        plt.savefig(jpegfilename, dpi=dpi)
    
    
