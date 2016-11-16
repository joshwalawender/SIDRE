import os
import re

import ccdproc as ccd

import astropy.units as u
from astropy import table

from SIDRE.config import get_config

def get_ImageFileCollection(filepath):
    '''
    Given a directory path with FITS files in it, use the header keywords (hard
    coded in this function) to categorize each file as one of:
    
    Science:       A science exposure
    Bias:          A bias frame
    Dark:          A dark frame
    Flat:          A flat field frame (twilight or dome)
    Rejected:      A file that has been rejection for any reason.
    Uncategorized: A file which was not categorized as one of the above.

    A column called "CATEGORY" is added to the `ImageFileCollection.summary`
    table and populated with a string of the above category.
    
    This method can be replaced to customize the code to any particular header
    or metadata convention.
    '''
    assert os.path.exists(os.path.abspath(filepath))
    temperature_deadband = get_config().get('TemperatureDeadband', 1.0)

    keywords = ['EXPTIME', 'SET-TEMP', 'CCD-TEMP', 'XBINNING', 'YBINNING',
                'IMAGETYP', 'OBJECT', 'DATE-OBS']
    ifc = ccd.ImageFileCollection(filepath, keywords=keywords)
    ifc.summary.add_column(table.Column(data=['']*len(ifc.summary),
                name='CATEGORY', dtype='a12'))
    for i,entry in enumerate(ifc.summary):
        tempdiff = float(entry['SET-TEMP']) - float(entry['CCD-TEMP'])
        if abs(tempdiff) > temperature_deadband:
            ifc.summary[i]['CATEGORY'] = b'Rejected'
        elif re.search('Light Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
            ifc.summary[i]['CATEGORY'] = b'Science'
        elif re.search('Bias Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
            ifc.summary[i]['CATEGORY'] = b'Bias'
        elif re.search('Dark Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
            ifc.summary[i]['CATEGORY'] = b'Dark'
        elif re.search('Flat', entry['IMAGETYP'], flags=re.IGNORECASE):
            ifc.summary[i]['CATEGORY'] = b'Flat'
        else:
            ifc.summary[i]['CATEGORY'] = b'Uncategorized'

    return ifc




if __name__ == '__main__':
    filepath = '/Volumes/Drobo/V5/Images/20161023UT'
    files = get_ImageFileCollection(filepath)
    print(files)
