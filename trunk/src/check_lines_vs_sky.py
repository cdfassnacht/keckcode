"""
Given the redshift of an object, plots the locations of typical emission
lines with respect to spectra of the NIR sky emission lines and transmission.
This code can, thus, be used to select targets for NIR spectroscopy.
"""

import os
try:
    from astropy.io import fits as pf
except:
    try:
        import pyfits as pf
    except:
        print ''
        print 'ERROR. Could not load either astropy.io.fits or pyfits'
        print 'Please make sure one of these is installed and in your'
        print '  PYTHONPATH'
        print ''
        exit()
import numpy as np
import matplotlib.pyplot as plt
import spec_simple as ss

""" Check the command-line syntax """

""" Get the object redshift from the command line """

w1 = None
w2 = None

""" Set up a fake wavelength vector, to run from 9000 Ang to 2 microns """
w = np.arange(9100.,25000.)
if w1 is not None:
    xmin = w1
else:
    xmin = w.min() - 100.
if w2 is not None:
    xmax = w2
else:
    xmax = w.max() + 100.

""" Get rid of the space between the subplots"""
plt.subplots_adjust(hspace=0.001)

""" Plot the atmospheric transmission spectrum """
ax1 = plt.subplot(211)
ss.plot_atm_trans(w,scale=1.)
plt.xlim(xmin,xmax)
plt.ylim(-0.15,1.1)

""" Plot the night-sky emission lines """
ax2 = plt.subplot(212,sharex=ax1)
skymod = ss.make_sky_model(w)
skyspec = ss.Spec1d(wav=w,flux=skymod)
skyspec.plot(title=None)
plt.xlim(xmin,xmax)
dy = 0.05 * skymod.max()
plt.ylim(-dy,(skymod.max()+dy))

""" Adjust the plot """
plt.show()
