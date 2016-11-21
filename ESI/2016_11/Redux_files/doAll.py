#from esi import esi_pipeline
from esi import bgsub

rawdir  = '../../Raw/'
out     = 'calib'
obsdate = '161121'

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
files = glob.glob('%s/e%s_005?.fits'%(rawdir,obsdate))
#files = glob.glob('%s/e%s_0061.fits'%(rawdir,obsdate))
#files = glob.glob('%s/e%s*.fits'%(rawdir,obsdate))
files.sort()
for f in files:
    print f
    h = pyfits.open(f)[0].header
    #if h['TARGNAME'][0:5] not in ['Feige']:
    #    print '  Skipping target %s' % h['targname']
    #    continue
    #if h['TARGNAME'][0:5] not in ['S2252','S2356','D0053']:
    #    print '  Skipping target %s' % h['targname']
    #    continue
    #if h['TARGNAME'][0] not in ['1']:
    #    continue
    if h['TARGNAME'][0].upper() not in ['A', 'J', 'S', 'D']:
        continue
    print h['TARGNAME']
    o = '%s_%s'%(h['TARGNAME'],f.split('_')[1].split('.')[0])
    bgsub.bgsub(rawdir,f.split('/')[-1],o,out)


