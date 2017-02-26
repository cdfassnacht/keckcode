"""
A set of functions related to preparing files related to using the DSIMULATOR
package for preparing DEIMOS slitmasks
"""

import astropy
if astropy.__version__[:3] == '0.3':
   from astropy.coordinates import ICRS as SkyCoord
else:
   from astropy.coordinates import SkyCoord
from astropy import units as u

def write_dsim(outcat, outradec, outid, pri, theta, selband, magname, outfile,
               verbose=False):
    """
    Writes the selected objects in a catalog out in the proper format to
    be used as inputs to DSIMULATOR
    """

    ngal = len(outradec)
    f = open(outfile,'w')
    tmpra = outradec.ra.to_string(unit=u.hourangle,decimal=False,sep=':',
                                  precision=3,pad=True)
    tmpdec = outradec.dec.to_string(decimal=False,sep=':',precision=3,
                                    alwayssign=True,pad=True)
    for i in range(ngal):
        g = outcat[i]
        tmpstr = '%-16s %s %s 2000.0 %5.2f %s %d 1 0 %.1f' % \
            (outid[i],tmpra[i],tmpdec[i],g[magname],selband,pri[i],theta[i])
        f.write('%s\n' % tmpstr)
        if verbose:
            print tmpstr

    f.close()
