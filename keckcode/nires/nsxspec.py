"""

nsxspec.py 

Primarily a definition of the NsxSpec class, which is written to handle
files produced by Tom Barlow's nsx (NIRES Spectral Extraction) code.

"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from specim import specfuncs as ss
from specim.specfuncs import echelle1d

# ===========================================================================


class NsxSpec(echelle1d.Ech1d):
    """
    A class used to visualize and analyze NIRES data that have been processed
    by Tom Barlow's nsx code
    """

    def __init__(self, root, frame, frame2=None, profileonly=False):
        """
        Loads, at minimum, the spatial profiles associated with either
        a single frame (if frame2 is None) or the result of a AB or BA
        subtraction
        """

        """ Set up some standard information """
        self.hasspec = False
        self.orders = [3, 4, 5, 6, 7]
        if frame2 is not None:
            self.inroot = '%s_%04d-%04d' % (root, frame, frame2)
        else:
            self.inroot = '%s_%04d' % (root, frame)

        """ Load the profiles """
        intot = '%s-pro.tbl' % self.inroot
        self.totprof = ascii.read(intot)

        """
        Load the spectra unless the user has requested only the profiles
        """
        if not profileonly:
            self.read_spec()

    # -----------------------------------------------------------------------

    def read_spec(self):
        """
        Reads in the extracted spectra for each order that have been produced
        by running nsx, and convert them into Spec1d format
        """

        """ Read in the spectra and convert them """
        for order in self.orders:
            infile = '%s-sp%d.tbl' % (self.inroot, order)
            nsxspec = ascii.read(infile)
            wav = nsxspec['angstrom']
            flux = nsxspec['object']
            var = nsxspec['error']**2

            """
            Often the starting and ending pixels have no signal, so
            mask them out
            """
            mask = var > 0.
            wav = wav[mask]
            flux = flux[mask]
            var = var[mask]

            """ Convert to a Spec1d instance and save """
            spec = ss.Spec1d(wav=wav, flux=flux, var=var)
            self.append(spec)

        self.hasspec = True

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

    def plot_spec(self, smo=None, z=None, title=None, **kwargs):
        """
        Plots all of the extracted spectra in one plot
        """

        """ Load the spectra if they have not already been loaded """
        if self.hasspec is False:
            self.read_spec()

        """ Set starting values """
        i = 0
        wmin = 30000.
        wmax = 0.

        """ Plot the spectra """
        for order in self.orders:
            """ Set the spectrum """
            spec = self[i]

            """ Set plot labeling """
            if order == 7:
                if title is None:
                    title = 'Extracted Spectrum'
                showz = True
            else:
                title = False
                showz = False

            """ Plot the spectrum """
            if smo is not None:
                spec.smooth(smo)
            else:
                spec.plot(**kwargs)

            """ Mark spectral lines if desired """
            if z is not None:
                if smo is not None:
                    spec.mark_lines('strongem', z, showz=showz, usesmooth=True)
                else:
                    spec.mark_lines('strongem', z, showz=showz)

            """ Adjust the plot limits """
            if spec['wav'].min() < wmin:
                wmin = spec['wav'].min()
            if spec['wav'].max() > wmax:
                wmax = spec['wav'].max()
            i += 1

        """ Set final plot limits """
        plt.xlim(wmin, wmax)
        
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
