"""

esiset.py - Code to perform actions on multiple 2d ESI data sets

"""

from os import path
from matplotlib import pyplot as plt
from .esi2d import Esi2d


class EsiSet(list):
    """

    The EsiSet class is effectively just a list of Esi2d instances that
    includes operations that act on the whole list (e.g., coadds)

    """

    # ------------------------------------------------------------------------

    def __init__(self, inlist, usevar=True, indir='.', prefix=None,
                 suffix='bgsub'):

        """
        Reads the input data and possibly its associated variance.
        The input files can be designated in one of two ways:

          1. A list of filenames - used if the prefix parameter is None
             e.g., ['J0745+271_0035_bgsub.fits', 'J0745+271_0036_bgsub.fits']

          2. A list of frame numbers - used if the prefix parameter is set.
             e.g., [35, 36]
             In this case, the input files would be created from the
             frame numbers, the prefix, and the suffix, as
              [prefix]_[frame]_[suffix].fits
             NOTE: the code assumes that the frame number in the filename
              is a 4-digit number with leading zeros.  Thus, having an
              inlist of [35, 36], a prefix of 'J0745', and a suffix of
              'bgsub' would produce filenames of:
                'J0745_0035_bgsub.fits', etc.
        """

        """ Read in the data """
        for i in range(len(inlist)):
            if prefix is not None:
                filename = '%s_%04d_%s.fits' % (prefix, inlist[i], suffix)
            else:
                filename = inlist[i]
            if isinstance(indir, list):
                datadir = indir[i]
            else:
                datadir = indir
            indat = path.join(datadir, filename)

            if usevar:
                varname = indat.replace(suffix, 'var')
            else:
                varname = None
            d = Esi2d(indat, varfile=varname)

            self.append(d)

    # ------------------------------------------------------------------------

    def extract(self, doplot=True, **kwargs):
        """

        Loops through each 2d spectrum in the list and extracts a 1d spectrum
         from each of the 10 orders in the 2d spectrum

        """

        for i in self:

            if doplot:
                plt.figure()

            """ Extract the 1d spectrum """
            print('')
            print(i.infile)
            i.extract_all(**kwargs)

            if doplot:
                plt.show()
