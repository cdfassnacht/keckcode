"""

esiset.py - Code to perform actions on multiple 2d ESI data sets

"""

from os import path
from matplotlib import pyplot as plt
from specim.specfuncs.specset1d import SpecSet1d
from specim.specfuncs.ech1dset import Ech1dSet
from .esi2d import Esi2d


class EsiSet(list):
    """

    The EsiSet class is effectively just a list of Esi2d instances that
    includes operations that act on the whole list (e.g., coadds)

    """

    # ------------------------------------------------------------------------

    def __init__(self, inlist, usevar=True, indir='.', prefix=None,
                 suffix='bgsub', verbose=True):

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
            if datadir == '.':
                indat = filename
            else:
                indat = path.join(datadir, filename)

            if usevar:
                varname = indat.replace(suffix, 'var')
            else:
                varname = None

            if verbose:
                print('')
                print(filename)
            d = Esi2d(indat, varfile=varname, verbose=verbose)

            self.append(d)

    # ------------------------------------------------------------------------

    def extract(self, doplot=True, verbose=True, debug=False, **kwargs):
        """

        Loops through each 2d spectrum in the list and, from each of the
         2d spectra extracts a the individual spectra from each order
         and returns them as an Esi1d object
        The final output from this method is thus a list of Esi1d objects

        """

        """ Set up the container for the extracted spectra """
        extract_list = []

        """
        Loop through the input 2d spectra and extract all of the spectral
        orders from each one with the extract_all method
        """
        for i in self:

            """ Extract the 1d spectra """
            if verbose:
                print('')
                print(i.infile)
            tmpspec = i.extract_all(plot_profiles=doplot,
                                    plot_extracted=doplot, **kwargs)
            extract_list.append(tmpspec)
            del(tmpspec)

            if doplot:
                plt.show()

        """
        Convert the list of extracted spectra to an Ech1dSet object
        and return
        """
        if debug:
            print(type(extract_list[0]))
            print(type(extract_list[0][0]))
        return Ech1dSet(extract_list)

    # ------------------------------------------------------------------------

    def coadd1d(self, doplot=True, outfile=None):
        """

        Takes a set of Esi1d objects

        """

        """ Loop through the orders """
        coaddlist = []
        for i in range(10):

            """ Create a list of Spec1d objects """
            speclist = []
            for espec in self:
                if espec[i].spec1d is not None:
                    speclist.append(espec[i].spec1d)
            if len(speclist) == 0:
                print('')
                print('ERROR: Called coadd1d but inputs do not have '
                      'extracted spectra yet')
                print('')
                raise ValueError

            """ Coadd the spectra in the list """
            specall = SpecSet1d(spec1dlist=speclist)
            coaddlist.append(specall.coadd(doplot=doplot))
            plt.show()
            del(speclist)

        return coaddlist
