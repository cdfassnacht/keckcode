"""
Code to extract some key info from the zresults*fits file that gets produced
after running zspec on calibrated DEIMOS data
"""

import sys
from astropy.io import fits as pf
from astropy.table import Table

maskname = sys.argv[1]

hdu = pf.open(maskname)
tdat = hdu[1].data
nobj = len(tdat)
for i in range(nobj):
    i = tdat[i]
    print '%-15s %-3s %-7s %1d %7.4f %g %s' % \
        (i['objname'],i['slitname'],i['maskname'],i['zquality'],i['z'], \
        i['z_err'],i['comment'])
    
