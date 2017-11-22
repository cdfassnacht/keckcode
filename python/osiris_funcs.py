"""

osiris_funcs.py

"""

import numpy as np
import imfuncs as imf

class osCube(imf.Image):
    """
    A class used to visualize and analyze OSIRIS data
    """

    def __init__(self, infile, verbose=True):
        """
        Loads in an OSIRIS data cube that was processed through the
        standard data reduction pipeline and, possibly, has had some additional
        processing done.
        """

        """ Read the data into an Image structure """
        imf.Image.__init__(self, infile, verbose=verbose)

    # -----------------------------------------------------------------------

    def set_imslice(self, imslice=0, hext=0):
        """
        Sets the 2-dimension slice to use for the display functions.

        Inputs:
          imslice - which image slice to use.
        """

        """ Get the number of dimensions in the input image """
        hdr = self.hdu[hext].header
        if 'NAXIS' in hdr.keys():
            ndim = hdr['naxis']
        else:
            raise KeyError

        """
        The OSIRIS data-reduction pipeline produces a 3-dimensional data cube.
        Check this
        """
        if ndim != 3:
            print('')
            print('ERROR: Expected a 3-dimensional data cube but found '
                  '%d dimensions' % ndim)
            print('')
            raise ValueError

        """
        Select the image slice to use.  Right now this assumes the standard
        axis order produced by the pipeline, namely that CTYPE1 is WAVE,
        CTYPE2 is RA, and CTYPE3 is DEC.  This means that the data array
        structure is hard-wired to this assumption
        """

        nwave = hdr['naxis1']
        if imslice >= nwave:
            print('')
            print('ERROR: Requested an image slice outside the available '
                  'range')
            print('Maximum available slice value is %d' % (nwave-1))
            print('')
            raise IndexError

        self.data = self.hdu[hext].data[:, :, imslice]


