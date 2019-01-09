from astropy.io import fits as pf
from specim import specfuncs as ss

"""
============================== Esi1d class ==============================
"""


class Esi1d(ss.Spec1d):
    """
    A class for ESI 1D spectra, which have been extracted by the Esi2d
    methods, but have not yet been combined into one final output spectrum.
    Therefore, there are 10 extracted 1d spectra, one for each order.
    These 10 extracted spectra will be stored in an array of Spec1d instances.

    The main purpose of this class is to combine the 10 orders into one
    output spectrum.  This functionality is split out from the Spec2d
    class because in some cases, e.g., co-adding several spectra of the same
    object, it may be easier to deal with the orders separately rather than
    after they have been combined into one Spec1d spectrum (however, this
    assertione may not be correc)
    """

    def __init__(self, infile):
        """

        Create an instance of this class by loading the data from the input
        file into 10 Spec1d instances.

        """

        """ Open the multiextension fits file that contains the 10 orders """
        self.infile = infile
        self.hdu = pf.open(infile)
        print('')
        print('Science file:  %s' % self.infile)

