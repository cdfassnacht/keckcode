#from esi import esi_pipeline
from esi import bgsub

dir = '../Raw/'
out = 'calib'

import glob
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

# Try to just select the science files
# NOTE: the line below that looks at the first letter in the TARGNAME
#  and only processes objects where the name begins with E or B is
#  specific to this data set.  All of the science targets in this 2013_05
#  observing run were either EELs or had B1950-coordinate names.
#
#files = glob.glob('%s/e*_00[345]?.fits'%dir)
files = glob.glob('%s/e*.fits'%dir)
files.sort()
for f in files:
    print f
    h = pyfits.open(f)[0].header
    if h['TARGNAME'][0:2] not in ['HS']:
        continue
    #if h['TARGNAME'][0] not in ['1']:
    #    continue
    #if h['TARGNAME'][0] not in ['E', 'A', 'J', '1', 'F']:
    #    continue
    print h['TARGNAME']
    o = '%s_%s'%(h['TARGNAME'],f.split('_')[1].split('.')[0])
    bgsub.bgsub(dir,f.split('/')[-1],o,out)


