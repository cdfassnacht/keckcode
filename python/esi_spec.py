"""
Code to take over after the calibration, etc., that is done by the code
in the esi directory.
"""

import spec_simple as ss

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
        
        self.infile = infile
        self.order = []
        for i in range(10):
            if i==0:
                verbose = True
            else:
                verbose = False
            tmpspec = ss.Spec2d(infile,hext=i+1,verbose=verbose)
            self.order.append(tmpspec)

