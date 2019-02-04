"""

nsxspec.py 

Primarily a definition of the NsxSpec class, which is written to handle
files produced by Tom Barlow's nsx (NIRES Spectral Extraction) code.

"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from specim import specfuncs as ss
from specim.specfuncs import echelle1d

# ===========================================================================


class NsxSpec(echelle1d.Ech1d):
    """
    A class used to visualize and analyze NIRES data that have been processed
    by Tom Barlow's nsx code
    """

    def __init__(self, root, frame, frame2=None, hasspec=True):
        """
        Loads, at minimum, the spatial profiles associated with either
        a single frame (if frame2 is None) or the result of a AB or BA
        subtraction
        """

        """ Set up some standard information """
        self.hasspec = False
        dtype = [('order', int), ('pixmin', int), ('pixmax', int)]
        oinfo = np.array([
                (3, 13, 2044),
                (4, 13, 2044),
                (5, 13, 2044),
                (6, 13, 2044),
                (7, 13, 1023),
                ], dtype=dtype)
        self.ordinfo = Table(oinfo)

        """ Set the input file name """
        if frame2 is not None:
            self.inroot = '%s_%04d-%04d' % (root, frame, frame2)
        else:
            self.inroot = '%s_%04d' % (root, frame)

        """ Load the total profile """
        intot = '%s-pro.tbl' % self.inroot
        self.totprof = ascii.read(intot)

        """
        If nsx was run in single-spectrum mode, load the other profiles
        """
        if frame2 is None:
            self.prof = []
            for info in self.ordinfo:
                pname = '%s-pro%d.tbl' % (self.inroot, info['order'])
                self.prof.append(ascii.read(pname))

        """
        Load the spectra unless the user has requested only the profiles
        """
        if hasspec:
            for info in self.ordinfo:
                sname = '%s-sp%d.tbl' % (self.inroot, info['order'])
                ospec = ss.Spec1d(sname, informat='nsx')

                """ Trim to the good region """
                w = ospec['wav'][info['pixmin']:info['pixmax']]
                f = ospec['flux'][info['pixmin']:info['pixmax']]
                v = ospec['var'][info['pixmin']:info['pixmax']]

                spec = ss.Spec1d(wav=w, flux=f, var=v)
                self.append(spec)

            self.hasspec = True
            # self.read_spec()

    # -----------------------------------------------------------------------

    def plot_profnsx(self, mode='median', order=None):
        """
        Plots either the total profile (default) or the profile for an
        individual order [NOT YET IMPLEMENTED]
        """

        """ Select the profile to be plotted """
        if order is not None:
            print('Order selection is not yet implemented')
        else:
            profile = self.totprof

        """ Plot the profile """
        if mode == 'average' or mode == 'mean':
            y = profile['average']
        else:
            y = profile['median']
        plt.plot(profile['arcsec'], y)

    # -----------------------------------------------------------------------

    def _make_resp345(self, smo=51, outfile=None, doplot=False):
        """
        Makes the response function for the three reddest orders, i.e.,
          Order 3 - K
          Order 4 - H
          Order 5 - J
        This method takes advantage of the fact that, when observing a
         standard star, the extracted spectra in these three orders
         look fairly similar, except for regions of strong atmospheric
         absorption (especially in order 5) and intrinsic spectral lines
         in the standard star (primarily hydrogen absorption).

        """

        """
        First do a rough correction for the atmospheric absorption.
        For now, this uses the model NIR atmospheric model in the specim code

        Once the correction has been done, normalize the corrected spectra
        by dividing by their maximum value
        """

        nspec = []
        for i in range(3):
            self[i].atm_corr()
            flux = self[i].atmcorr['flux']
            ns = flux / flux.max()
            nspec.append(ns)

        """
        Set up the unsmoothed spectrum to use as the basis for the
        response curve
        """
        resp = nspec[0]
        resp[:760] = nspec[2][:760]
        resp[915:999] = nspec[2][915:999]
        rspec = ss.Spec1d(wav=np.arange(resp.size), flux=resp)

        """ Smooth the curve and save it """
        rspec.smooth(smo, doplot=False)

        """ Plot if requested """
        if doplot:
            for spec in nspec:
                plt.plot(spec)
            plt.plot(resp, 'r')
            rspec.plot(mode='smooth', color='k')

        """ Save the smoothed response file if requested """
        if outfile is not None:
            rspec.smospec.save(outfile)

        """ Return the smoothed output file """
        return rspec.smospec
