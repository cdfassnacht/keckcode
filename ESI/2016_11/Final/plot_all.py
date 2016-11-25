"""
Does a summary plot of all the spectra
"""

from astropy.io import ascii
from matplotlib import pyplot as plt
import spec_simple as ss

plotinfo = ascii.read('plotinfo.txt')
for i in plotinfo:
    infile = '%s_spec_%s.fits' % (i['obj'],i['aper'])
    tmpspec = ss.Spec1d(infile,informat='fits',logwav=True)
    if i['smo']>1:
        tmpspec.smooth_boxcar(i['smo'])
        smo = True
    else:
        tmpspec.plot()
        smo = False
    plt.xlim(4180.,10260.)
    plt.ylim(-0.03,i['ymax'])
    tmpspec.mark_speclines(i['type'],i['z'],usesmooth=smo)
    plt.show()
