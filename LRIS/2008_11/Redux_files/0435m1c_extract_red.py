"""
A series of commands to extract 1D spectra from the 2D spectra in the
0435m1c red-side slitmask.
"""

import spec_simple as ss
import numpy as np

""" Start by reading in the slit boundary definitions """
slitnum,x1,x2,y1,y2,y1b,y2b = np.loadtxt('0435m1c_slits_red.txt',unpack=True)

""" Loop through the slits """
for i in range(slitnum.size):
    print ''
    print '====================================================================='
    print ''
    print 'Extracting slit %.0f' % (slitnum[i])
    """ 
    Start by reading in the coadded background-subtracted slit.  The data
    in this slit will be used to define the trace in the individual exposures
    """
    bgslit = ss.Spec2d('0435m1c_red_bgsub.fits',xtrim=[x1[i],x2[i]],
                       ytrim=[y1b[i],y2b[i]])
    bgslit.find_and_trace(doplot=False)
    aper = [-1.5*bgslit.sig0, 1.5*bgslit.sig0]

    """ Now loop through the individual exposures and extract the spectra """
    for j in range(5):
        inname  = '0435m1c_red_straight_%d.fits' % (j+1)
        outname = '0435m1c_slit%02d_%d.spec' % (slitnum[i],(j+1))
        slitj = ss.Spec2d(inname,xtrim=[x1[i],x2[i]],ytrim=[y1[i],y2[i]])
        slitj.aper = aper
        slitj.find_and_trace()
        slitj.mu = bgslit.mu + slitj.mu0 - bgslit.mu0
        slitj.extract_spectrum(outfile=outname,doplot=False)
        del slitj
    del bgslit


