#!/usr/env/python

## Import General Tools
from pathlib import Path
import argparse
import logging


##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description='''
''')
## add flags
p.add_argument("-v", "--verbose", dest="verbose",
    default=False, action="store_true",
    help="Be verbose! (default = False)")
## add options
p.add_argument("--input", dest="input", type=str,
    help="The input.")
## add arguments
p.add_argument('argument', type=int,
               help="A single argument")
p.add_argument('allothers', nargs='*',
               help="All other arguments")
args = p.parse_args()


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('MyLogger')
log.setLevel(logging.DEBUG)
## Set up console output
LogConsoleHandler = logging.StreamHandler()
LogConsoleHandler.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
LogConsoleHandler.setFormatter(LogFormat)
log.addHandler(LogConsoleHandler)
## Set up file output
# LogFileName = None
# LogFileHandler = logging.FileHandler(LogFileName)
# LogFileHandler.setLevel(logging.DEBUG)
# LogFileHandler.setFormatter(LogFormat)
# log.addHandler(LogFileHandler)


##-------------------------------------------------------------------------
## Exceptions
##-------------------------------------------------------------------------
class MEFDataError(Exception):
    """Base class for exceptions in this module."""
    pass


class IncompatiblePixelData(MEFDataError):
    """Raise when trying to operate on multiple MEFData
    objects which have incompatible pixeldata.
    """
    def __init__(self, message):
        super().__init__(f"MEFData objects have incompatible pixeldata. {message}")


class IncorrectNumberOfExtensions(MEFDataError):
    """Raise when verify method fails for a specific instrument.
    """
    def __init__(self, datatype, expected, kd):
        msg = f"Incorrect number of {datatype} entries.  Expected {expected} for {type(kd)}"
        print(msg)
        super().__init__(msg)


##-------------------------------------------------------------------------
## MEFData Classes
##-------------------------------------------------------------------------
class MEFData(object):
    """Our data model.
    
    Attributes:
    pixeldata -- a list of CCDData objects containing pixel values.
    tabledata -- a list of astropy.table.Table objects
    headers -- a list of astropy.io.fits.Header objects
    """
    def __init__(self, *args, **kwargs):
        self.pixeldata = []
        self.tabledata = []
        self.headers = []

    def verify(self):
        """Method to check the data against expectations. For the 
        MEFData class this simply passes and does nothing, but
        subclasses for specific instruments can populate this
        with appropriate tests.
        """
        pass

    def add(self, kd2):
        """Method to add another MEFData object to this one and return
        the result.  This uses the CCDData object's add method and
        simply loops over all elements of the pixeldata list.
        """
        if len(self.pixeldata) != len(kd2.pixeldata):
            raise IncompatiblePixelData
        for i,pd in enumerate(self.pixeldata):
            self.pixeldata[i] = pd.add(kd2.pixeldata[i])

    def subtract(self, kd2):
        """Method to subtract another MEFData object to this one
        and return the result.  This uses the CCDData object's
        subtract method and simply loops over all elements of
        the pixeldata list.
        """
        if len(self.pixeldata) != len(kd2.pixeldata):
            raise IncompatiblePixelData
        for i,pd in enumerate(self.pixeldata):
            self.pixeldata[i] = pd.subtract(kd2.pixeldata[i])

    def multiply(self, kd2):
        """Method to multiply another MEFData object by this one
        and return the result.  This uses the CCDData object's
        multiply method and simply loops over all elements of
        the pixeldata list.
        """
        if len(self.pixeldata) != len(kd2.pixeldata):
            raise IncompatiblePixelData
        for i,pd in enumerate(self.pixeldata):
            self.pixeldata[i] = pd.multiply(kd2.pixeldata[i])

    def get(self, kw):
        """Method to loop over all headers and get the specified keyword value.
        Returns the first result it finds and doe not check for duplicate
        instances of the keyword in subsequent headers.
        """
        for hdr in self.headers:
            val = hdr.get(kw, None)
            if val is not None:
                return val


##-------------------------------------------------------------------------
## Get HDU Type
##-------------------------------------------------------------------------
def get_hdu_type(hdu):
    """Function to examine a FITS HDU object and determine its type.  Returns
    one of the following strings:
    
    'header' -- This is a PrimaryHDU or ImageHDU with no pixel data.
    'pixeldata' -- This is a PrimaryHDU or ImageHDU containing pixel data.
    'uncertainty' -- This is a pixeldata HDU which is associated with the
                     uncertainty data written by either CCDData or MEFData.
    'mask' -- This is a pixeldata HDU which is associated with the mask
              data written by either CCDData or MEFData.
    'tabledata' -- This is a TableHDU type HDU.
    """
    if type(hdu) in [fits.PrimaryHDU, fits.ImageHDU] and hdu.data is None:
        # This is a header only HDU
        return 'header'
    elif type(hdu) in [fits.PrimaryHDU, fits.ImageHDU] and hdu.data is not None:
        # This is a pixel data HDU
        extname = hdu.header.get('EXTNAME', '').strip()
        if extname == 'MASK':
            # This is a mask HDU
            return 'mask'
        elif extname == 'UNCERT':
            # This is an uncertainty HDU
            return 'uncertainty'
        else:
            # This must be pixel data
            return 'pixeldata'
    elif type(hdu) == fits.TableHDU:
            # This is table data
            return 'tabledata'


##-------------------------------------------------------------------------
## MEFData Reader
##-------------------------------------------------------------------------
def fits_MEFdata_reader(file, defaultunit='adu', datatype=MEFData):
    """A reader for MEFData objects.
    
    Currently this is a separate function, but should probably be
    registered as a reader similar to fits_ccddata_reader.
    
    Arguments:
    file -- The filename (or pathlib.Path) of the FITS file to open.
    Keyword arguments:
    defaultunit -- If the BUNIT keyword is unable to be located or
                   parsed, the reader will assume this unit.  Defaults
                   to "adu".
    datatype -- The output datatype.  Defaults to MEFData, but could
                be a subclass such as MOSFIREData.  The main effect of
                this is that it runs the appropriate verify method on
                the data.
    """
    try:
        hdul = fits.open(file, 'readonly')
    except FileNotFoundError as e:
        print(e.msg)
        raise e
    except OSError as e:
        print(e.msg)
        raise e
    # Loop though HDUs and read them in as pixel data or table data
    md = datatype()
    while len(hdul) > 0:
        print('Extracting HDU')
        hdu = hdul.pop(0)
        md.headers.append(hdu.header)
        hdu_type = get_hdu_type(hdu)
        print(f'  Got HDU type = {hdu_type}')
        if hdu_type == 'header':
            pass
        elif hdu_type == 'tabledata':
            md.tabledata.append(Table(hdu.data))
        elif hdu_type == 'pixeldata':
            # Check the next HDU
            mask = None
            uncertainty = None
            if len(hdul) > 0:
                next_type = get_hdu_type(hdul[0])
                if next_type == 'mask':
                    mask = hdul[0].data
                elif next_type == 'uncertainty':
                    uncertainty = hdul[0].data
            if len(hdul) > 1:
                next_type2 = get_hdu_type(hdul[1])
                if next_type2 == 'mask':
                    mask = hdul[1].data
                elif next_type2 == 'uncertainty':
                    uncertainty = hdul[1].data               
            # Sanitize "ADU per coadd" BUNIT value
            if hdu.header.get('BUNIT') == "ADU per coadd":
                hdu.header.set('BUNIT', 'adu')
            # Populate the CCDData object
            c = CCDData(hdu.data, mask=mask, uncertainty=uncertainty,
                        meta=hdu.header,
                        unit=hdu.header.get('BUNIT', defaultunit),
                       )
            md.pixeldata.append(c)
    print(f'Read in {len(md.headers)} headers, '
          f'{len(md.pixeldata)} sets of pixel data, '
          f'and {len(md.tabledata)} tables')
    md.verify()
    return md
