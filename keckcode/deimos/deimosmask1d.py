"""

deimosmask1d  - Code to define a DeimosMask1d class, which is used to
                visualize and process the spec1d*fits files that come
                out of the PypeIt data reduction pipeline for DEIMOS data

"""

from astropy.io import fits as pf
from specim.specfuncs import spec1d


class DeimosMask1d(list):
    """

    Class to visualize and perhaps process the 1d spectra that are contained
    in the spec1d*fits that come out of the PypeIt data reduction pipeline
    for DEIMOS data

    """

    def __init__(self, infile, verbose=True):
        """

        Instantiate a DeimosMask1d object

        Inputs:
         infile   - a spec1d*fits file that has been produced by the PypeIt
                    data reduction pipeline

        """

        """
        Load the data from the input file and get information from the
        primary header
        """
        hdu = pf.open(infile)
        hdr0 = hdu[0].header
        self.nspec = hdr0['nspec']

        """ Load the spectra into a list of Spec1d objects """
        colnames = ['opt_wave', 'opt_counts', 'opt_counts_ivar',
                    'opt_counts_sky']
        if verbose:
            print('')
            print('Loading %d spectra from %s' % (self.nspec, infile))
        for i in range(1, self.nspec+1):
            spec = spec1d.Spec1d(hdu[i], colnames=colnames, verbose=False)
            mask = spec['var'] > 0.
            spec['var'][mask] = 1. / spec['var'][mask]
            spec['var'][~mask] = 25. * spec['var'].max()
            self.append(spec)

    # -----------------------------------------------------------------------

    def plot(self, index, **kwargs):
        """

        Plots the spectrum designated by index, where index runs from
        0 to nspec-1

        """

        self[index].plot(**kwargs)

    # -----------------------------------------------------------------------

    def smooth(self, index, filtwidth, **kwargs):
        """

        Plots a smoothed version of the spectrum designated by index,
         where index runs from 0 to nspec-1

        Inputs:
         index     - index of spectrum to plot (between 0 and nspec-1)
         filtwidth - width of kernel to be used for smoothing

        """

        self[index].smooth(filtwidth, **kwargs)
