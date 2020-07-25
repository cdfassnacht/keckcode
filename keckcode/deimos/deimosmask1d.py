"""

deimosmask1d  - Code to define a DeimosMask1d class, which is used to
                visualize and process the spec1d*fits files that come
                out of the PypeIt data reduction pipeline for DEIMOS data

"""

import numpy as np
from collections import OrderedDict
from astropy.io import fits as pf
from specim.specfuncs import spec1d

import sys
pyversion = sys.version_info.major

# ---------------------------------------------------------------------------


class DeimosMask1d(OrderedDict):
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

        """ Set up the empty container by calling the superclass """
        if pyversion == 2:
            super(DeimosMask1d, self).__init__()
        else:
            super().__init__()

        """
        Load the data from the input file and get information from the
        primary header
        """
        hdu = pf.open(infile)
        hdr0 = hdu[0].header
        self.nspec = hdr0['nspec']
        self.slitid = np.zeros(self.nspec)

        """ Load the spectra into a list of Spec1d objects """
        colnames = ['opt_wave', 'opt_counts', 'opt_counts_ivar',
                    'opt_counts_sky']
        if verbose:
            print('Read %d spectra from %s' % (self.nspec, infile))
        for i in range(1, self.nspec+1):
            spec = spec1d.Spec1d(hdu[i], colnames=colnames, verbose=False)
            mask = spec['var'] > 0.
            spec['var'][mask] = 1. / spec['var'][mask]
            spec['var'][~mask] = 25. * spec['var'].max()
            # self.append(spec)
            slitid = hdu[i].header['slitid']
            self[slitid] = spec
            self.slitid[(i-1)] = slitid

    # -----------------------------------------------------------------------

    def __str__(self):
        """

        Gives a cleaner output when doing a print statement for the object

        """

        return('DeimosMask1d with %d spectra' % self.nspec)
    
    # -----------------------------------------------------------------------

    def __repr__(self):
        """

        Gives a cleaner representation of the object

        """

        return('DeimosMask1d with %d spectra' % self.nspec)
    
    # -----------------------------------------------------------------------

    def plot(self, slitid, **kwargs):
        """

        Plots the spectrum for the slit designated by the slitid paramters

        Inputs:
           slitid    - slit ID (an integer value)
        """
        return self[slitid].plot(**kwargs)

    # -----------------------------------------------------------------------

    def smooth(self, slitid, filtwidth, **kwargs):
        """

        Plots a smoothed version of the spectrum designated by the
         slitid parameter

        Inputs:
         slitid    - slit ID (an integer value)
         filtwidth - width of kernel to be used for smoothing

        """

        self[index].smooth(filtwidth, **kwargs)

    # -----------------------------------------------------------------------

    def mark_lines(self, lines, z, index, **kwargs):
        """
        Draws a plot with spectral line marks on a given spectra.

        Inputs:
          index - index of spectrum to plot (between 0 and nspec-1)
          z     - redshift estimate
          lines - dictionary. Key should be the line name, value should be
                  boolean

        """
        if lines['em'] and lines['strongem']:
            lines.pop('strongem')

        for k,v in lines.items():
            if v:
                self[0].mark_lines(k, z, **kwargs)
