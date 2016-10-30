import os
import re

import ccdproc as ccd

import astropy.units as u
from astropy import table

def sort(filepath):
    '''Given a list of FITS files
    '''
    temp_threshold = 1.0

    if type(filepath) is str:
        filepaths = [filepath]
    elif type(filepath) is list:
        filepaths = filepath

    keywords = ['EXPTIME', 'SET-TEMP', 'CCD-TEMP', 'XBINNING', 'YBINNING', 
                'IMAGETYP', 'OBJECT', 'DATE-OBS']
    image_table = None
    for fp in filepaths:
        if os.path.exists(os.path.abspath(fp)):
            ifc = ccd.ImageFileCollection(fp, keywords=keywords)
            ifc.summary.add_column(table.Column(data=['']*len(ifc.summary),
                        name='CATEGORY', dtype='a12'))
            ifc.summary.add_column(table.Column(data=['']*len(ifc.summary),
                        name='DIR', dtype='a80'))
            print('Examining {} files in {}'.format(len(ifc.summary), fp))
            for i,entry in enumerate(ifc.summary):
                ifc.summary[i]['DIR'] = ifc.location
                tempdiff = float(entry['SET-TEMP']) - float(entry['CCD-TEMP'])

                if abs(tempdiff) > temp_threshold:
                    ifc.summary[i]['CATEGORY'] = b'Rejected'
                elif re.search('Light Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
                    ifc.summary[i]['CATEGORY'] = b'Object'
                elif re.search('Bias Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
                    ifc.summary[i]['CATEGORY'] = b'Bias'
                elif re.search('Dark Frame', entry['IMAGETYP'], flags=re.IGNORECASE):
                    ifc.summary[i]['CATEGORY'] = b'Dark'
                elif re.search('Flat', entry['IMAGETYP'], flags=re.IGNORECASE):
                    ifc.summary[i]['CATEGORY'] = b'Flat'
                else:
                    ifc.summary[i]['CATEGORY'] = b'Uncategorized'
            if image_table is None:
                image_table = ifc.summary
            else:
                image_table = table.vstack([image_table, ifc.summary])

    bycat = image_table.group_by('CATEGORY')
    for cat in bycat.groups.keys['CATEGORY']:
        thiscat = bycat.groups[bycat.groups.keys['CATEGORY'] == cat]
        print('Found {} {} files.'.format(len(thiscat), cat.decode('utf8').upper()))

    return image_table


if __name__ == '__main__':
#     filepath = '/Volumes/Drobo/V5/Images/20161023UT'
    filepath = ['/Volumes/Drobo/V5/Images/20161023UT',
                '/Volumes/Drobo/V5/Images/20161023UT/Calibration',
                '/Volumes/Drobo/V5/Images/20161023UT/AutoFlat']
    files = sort(filepath)
    print(files)
