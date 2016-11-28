"""
Code to take over after the calibration, etc., that is done by the code
in the esi directory.
"""

import numpy as np
from astropy.io import fits as pf
import spec_simple as ss
from matplotlib import pyplot as plt

class Esi2d(ss.Spec2d):
    """
    A class for ESI 2D spectra, for which Matt's / Lindsay's calibration
    code produces a multi-extension fits files containing 10 2D spectra,
    one for each of the 10 orders produced by the spectrograph
    """

    def __init__(self, infile):
        """

        Create an instance of this class by loading the data from the input
        file into 10 Spec2d instances.

        """

        """ Open the multiextension fits file that contains the 10 orders """
        self.infile = infile
        self.hdu = pf.open(infile)
        print ''
        print self.infile

        """ Load each order into its own Spec2d container """
        self.order = []
        print ''
        print 'Order  Shape    Dispaxis'
        print '----- --------- --------'
        for i in range(10):
            tmpspec = ss.Spec2d(None,hdulist=self.hdu,hext=i+1,verbose=False,
                                logwav=True,fixnans=False)
            print ' %2d   %dx%d     %s' % \
                ((i+1),tmpspec.data.shape[1],tmpspec.data.shape[0],
                 tmpspec.dispaxis)
            self.order.append(tmpspec)

        """ Each order for ESI has a different plate scale in arcsec/pix """
        self.arcsecperpix = np.array([
                0.120, # order 1
                0.127, # order 2
                0.134, # order 3
                0.137, # order 4
                0.144, # order 5
                0.149, # order 6
                0.153, # order 7
                0.158, # order 8
                0.163, # order 9
                0.168  # order 10
                ])

        """ 
        Set the range of valid pixels for each order, expressed as 
        blue (start of good pixels) to red (end of good pixels) 
        """
        self.blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
        self.red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]

    #-----------------------------------------------------------------------

    def plot_profiles(self):
        """

        Plots, in one figure, the spatial profiles for all the 10 orders

        """

        for i in range(10):
            plt.subplot(2,5,(i+1))
            self.order[i].spatial_profile()

    #-----------------------------------------------------------------------

    def extract_all(self):
        """
        Goes through each of the 10 orders on the ESI spectrograph and
        does the two-step extraction procedure from the Spec2d class in
        spec_simple.py: find_and_trace and extract_spectrum
        """
