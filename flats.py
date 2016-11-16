import os
import sys
import re
import numpy as np
from scipy.ndimage.filters import median_filter

import ccdproc as ccd

import astropy.units as u
from astropy import table

from astropy import log
log.setLevel('WARNING')
log.disable_warnings_logging()

from SIDRE.config import get_config
import SIDRE




if __name__ == '__main__':
    filepath = '/Users/vysosuser/ShutterMap/V5/20160723UT/'
    files = SIDRE.sort.get_ImageFileCollection(filepath)
