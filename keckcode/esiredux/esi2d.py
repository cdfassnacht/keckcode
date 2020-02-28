"""
Code to take over after the calibration, etc., that is done by the code
in the esi directory.  There are two new classes, plus some code that deals
with combining multiple input files.  The two classes are:

* Esi2d - used to find and trace spectra, and to convert the 2d spectra into
          an output 1d spectrum in the form of an Esi1d instance.
          This class is effectively an array of Spec2d instances, plus
          some additional code that handles the complexity of a 10-member
          structure.
* Esi1d - Similar to Esi2d, but in this case effectively an array of
          Spec1d instances, plus code to deal with the 10-order structure
          of ESI data.  There is also code to combine these ten orders into
          a single Spec1d output.
          UNDER CONSTRUCTION
"""

import numpy as np
from math import log10
from scipy import ndimage, interpolate
from matplotlib import pyplot as plt

from astropy.io import fits as pf
from astropy.table import Table

import special_functions as sf

from cdfutils import datafuncs as df
from specim import specfuncs as ss
from specim.specfuncs import echelle2d as ech2d

"""
============================== Esi2d class ==============================
"""


class Esi2d(ech2d.Ech2d):
    """
    A class for ESI 2D spectra, for which Matt's / Lindsay's calibration
    code produces a multi-extension fits files containing 10 2D spectra,
    one for each of the 10 orders produced by the spectrograph.
    Essentially this class is just a list of Spec2d objects
    """

    def __init__(self, infile, varfile=None):
        """

        Create an instance of this class by loading the data from the input
        file into 10 Spec2d instances.

        """

        """
        Start by setting up some default values and put it all into a
         recarray structure for ease in extracting values.
         pixscale     - each order for ESI has a different plate scale in
                        arcsec/pix
         blue and red - define the range of good pixels that the Oldham
                        extraction uses to define its aperture
        """
        dtype = [('order', int), ('name', 'S8'), ('pixscale', float),
                 ('pixmin', int), ('pixmax', int)]
        orderinfo = np.array([
                (1, 'Order_01', 0.120, 1500, 3000),
                (2, 'Order_02', 0.127, 1400, 3400),
                (3, 'Order_03', 0.134, 1300, 3700),
                (4, 'Order_04', 0.137, 1200,   -1),
                (5, 'Order_05', 0.144, 1100,   -1),
                (6, 'Order_06', 0.149,  900,   -1),
                (7, 'Order_07', 0.153,  600,   -1),
                (8, 'Order_08', 0.158,  200,   -1),
                (9, 'Order_09', 0.163,    0,   -1),
                (10, 'Order_10', 0.168,   0,   -1),
                ], dtype=dtype)
        self.ordinfo = Table(orderinfo)

        """ Read in the spectra by a call to the superclass """
        super(Esi2d, self).__init__(infile, varspec=varfile, logwav=True,
                                    ordinfo=self.ordinfo, fixnans=False)

        """ Set the default gridding for multi-order plots """
        self.plotgrid = (2, 5)

    # ------------------------------------------------------------------------

    def get_ap_oldham(self, slit, apcent, nsig, ordinfo, doplot=True):
        """
        Defines a uniform aperture as in the example ESI extraction scripts
        from Lindsay
        """

        B = ordinfo['pixmin']
        R = ordinfo['pixmax']
        xproj = np.median(slit[:, B:R], 1)
        m, s = df.sigclip(xproj)

        smooth = ndimage.gaussian_filter(xproj, 1)
        if ordinfo['name'] == 'Order 3':
            smooth = ndimage.gaussian_filter(xproj[:-30], 1)
        x = np.arange(xproj.size) * 1.

        """
        The four parameters immediately below are the initial guesses
        bkgd, amplitude, mean location, and sigma for a Gaussian fit
        """
        fit = np.array([0., smooth.max(), smooth.argmax(), 1.])
        fit = sf.ngaussfit(xproj, fit)[0]

        cent = fit[2] + apcent / ordinfo['pixscale']
        apymax = 0.1 * xproj.max()
        ap = np.where(abs(x-cent) < nsig / ordinfo['pixscale'], 1., 0.)
        # slit.apmin = cent - nsig
        # slit.apmax = cent + nsig

        if doplot:
            plt.subplot(2, 5, (ordinfo['order']))
            plt.plot(x, apymax*ap)  # Scale the aperture to easily see it
            plt.plot(x, xproj)
            plt.ylim(-apymax, 1.1*xproj.max())
            plt.axvline(cent, color='k', ls='dotted')

        ap = ap.repeat(slit.shape[1]).reshape(slit.shape)
        return ap, fit

    # --------------------------------------------------------------------

    def _extract_oldham(self, spec2d, ordinfo, apcent, nsig, normap=False):
        """
        Implements Lindsay Oldham's (or perhaps Matt Auger's) method
         for extracting the spectrum from an individual spectral order
         on ESI.

        Inputs:
          spec2d  - the requested spectral order
        """

        """
        Make temporary copies of the science and variance spectra, and
        get information about their properties
        """
        spec = spec2d
        slit = spec.data.copy()
        vslit = spec.vardata.copy()
        vslit[vslit <= 0.] = 1e9
        vslit[np.isnan(vslit)] = 1e9
        vslit[np.isnan(slit)] = 1e9
        h = spec.header
        x = np.arange(slit.shape[1])*1.
        w = 10**(h['CRVAL1'] + x * h['CD1_1'])

        """ Make the apertures """
        ap, fit = self.get_ap_oldham(slit, apcent, nsig, ordinfo)
        cent = fit[2] + apcent / ordinfo['pixscale']
        print(cent)
        spec2d.apmin = cent - nsig
        spec2d.apmax = cent + nsig

        """
        Set up to do the extraction, including normalizing the aperture
        if requested
        """
        ap[vslit >= 1e8] = 0.
        if normap:
            ap = ap / ap.sum(0)
        ap[np.isnan(ap)] = 0.
        slit[np.isnan(slit)] = 0.

        """
        Extract the spectrum
        If the aperture is requested to be normalized (normap is True) then
         first normalize ap (above) and then use the weighting (which I don't
         understand) that Lindsay used in her extraction code.
        Otherwise, use the weighting (which makes more sense to me) that
         was used in the old makeOrderCorr.py files.
        """
        if normap:
            spec = (slit*ap**2).sum(0)
            vspec = (vslit*ap**4).sum(0)
        else:
            spec = (slit*ap).sum(0)
            vspec = (vslit*ap**2).sum(0)

        """
        Normalize the spectrum and its associated variance
        Need to do the variance first, or else you would incorrectly
        be using the normalized spectrum to normalize the variance
        """
        vspec /= np.median(spec)**2
        spec /= np.median(spec)

        """ Save the output in a Spec1d container """
        # newspec = ss.Spec1d(wav=w,flux=spec,var=vspec)
        spec2d.spec1d = ss.Spec1d(wav=w, flux=spec, var=vspec)

    # --------------------------------------------------------------------

    def _extract_cdf(self, spec, info, muorder=-1, sigorder=-1,
                     apmin=-1., apmax=1., weight='gauss', normalize=False,
                     plot_traces=False):
        """

        Extracts the spectrum from an individual order using the
         Spec2d class functionality from the specim.specfuncs module

        For now assume that the data reduction has properly rectified
        the 2D spectra, so that it is not necessary to trace a varying
        position and width as the position changes
        NOTE: This is almost certainly not true, but keep for now
        """
        
        """ Clear the figure container if we're going to plot the trace  """
        if plot_traces:
            plt.figure(info['order'] + 2)
            plt.clf()
            
        """ Set the aperture boundaries, IN ARCSECONDS """
        spec.apmin = apmin / info['pixscale']
        spec.apmax = apmax / info['pixscale']

        """ Set the range of valid pixels for fitting the trace """
        B = info['pixmin']
        R = info['pixmax']

        """ Trace and then extract the spectrum using the Spec2d methods  """
        print('')
        print('==================================================')
        print('%s' % info['name'])
        print('')
        spec.find_and_trace(doplot=plot_traces, muorder=muorder,
                            sigorder=sigorder, fitrange=[B, R],
                            verbose=False)
        spec.extract(extrange=[B, R], weight=weight, doplot=False,
                     verbose=False)

        """
        Also normalize the flux and variance of the extracted
        spectrum in the same way that the Oldham extraction does
        """
        if normalize:
            medflux = np.median(spec.spec1d['flux'])
            spec.spec1d['flux'] /= medflux
            spec.spec1d['var'] /= medflux**2

   # --------------------------------------------------------------------

    def plot_extracted(self, method='1x1', xmin=3840., xmax=10910.,
                       ymin=None, ymax=None, color='b'):
        """

        Plots the 10 extracted 1d spectra in one plot.  There are two
        ways to do this, which are set by the method parameter:
          method='oneplot' - (default) Have one plot window showing the
                              full wavelength range and the overlapping orders
          method='tenplot' - The plot has 10 separate windows, one for each
                              order

        """
        if method == 'tenplot':
            print('NOTE: tenplot is not yet implemented')
            return
        else:
            for i in self:
                i.spec1d.plot(color=color)
            plt.xlim(xmin, xmax)
            if ymin is not None and ymax is not None:
                plt.ylim(ymin, ymax)

    # --------------------------------------------------------------------

    def extract_all(self, method='oldham', apcent=0., nsig=1.0,
                    apmin=-1., apmax=1., muorder=-1, sigorder=-1,
                    normap=False, weight='gauss', plot_profiles=True,
                    plot_traces=False, plot_extracted=True,
                    xmin=3840., xmax=10910., ymin=-0.2, ymax=5.,
                    apnum=None, showfit=False):
        """
        Goes through each of the 10 orders on the ESI spectrograph and
        extracts the spectrum via one of two procedures:
        1. Using the Spec2d class methods from spec_simple.py (method='cdf')
           This approach does the two-step extraction procedure using the
           relevant methods in the Spec2d class in spec_simple.py, namely
           (1) find_and_trace and (2) extract_spectrum.
          xxxxx
        2. Using the code from Lindsay Oldham (or perhaps Matt Auger)
           implemented above (method='oldham')
           In this approach, the two step procedure is (1) get_ap_oldham, and
           (2) _extract_oldham
          xxxxx
        """

        """
        Set up for plotting profiles, which needs to be done here since the
        'oldham' and 'cdf' techniques do the plotting at different points
        """
        if plot_profiles:
            plt.figure()

        """
        Extract the spectra
        """
        # for i in range(10):
        for spec, info in zip(self, self.ordinfo):
            if method == 'cdf':
                self._extract_cdf(spec, info, plot_traces=plot_traces,
                                  muorder=muorder, sigorder=sigorder,
                                  apmin=apmin, apmax=apmax, weight=weight)
            elif method == 'oldham':
                self._extract_oldham(spec, info, apcent, nsig, normap=normap)

        """ Plot all the spatial profiles, along with the profile fits """
        if plot_profiles and method == 'cdf':
            self.plot_profiles(showfit=showfit)

        """
        Plot the extracted spectra
        """
        if plot_extracted:
            plt.figure()
            self.plot_extracted()

    # --------------------------------------------------------------------

    def respcorr_oldham(self, respfile):
        """
        Does the response correction using the technique that Lindsay Oldham /
        Matt Auger have used.
        """
        corr = np.load(respfile)
        right = None
        rb = None
        rr = None

        """ Set up the wavelength vector for the corrected 1-d spectrum """
        scale = 1.7e-5
        w0 = log10(self[0].spec1d['wav'][0])
        w1 = log10(self[9].spec1d['wav'][-1])
        outwav = np.arange(w0, w1, scale)
        outflux = np.ones((outwav.size, 10)) * np.nan
        outvar = outflux.copy()

        for i in range(1, 10):
            """
            Create the response correction for this order and then correct
            the flux and the variance
            """
            w = self[i].spec1d['wav']
            w0, w1, mod = corr[i+1]
            mod = sf.genfunc(w, 0., mod)
            spec = self[i].spec1d['flux'] / mod
            var = self[i].spec1d['var'] / mod**2

            """ Mask out NaNs """
            mask = np.isnan(spec)
            spec[mask] = 0.
            var[mask] = 1e9

            """
            Only correct for the ranges over which the correction is valid
            """
            mask = (w > w0) & (w < w1)
            w = w[mask]
            spec = spec[mask]
            var = var[mask]
            if right is not None:
                left = np.median(spec[(w > rb) & (w < rr)])
                spec *= right / left
                var *= (right / left)**2
            try:
                """
                Find the overlap region
                 - blue end (rb) is the start of the next order
                 - red end (rr) is the end of this current spectrum
                """
                rb = self[i+1].spec1d['wav'][0]
                rr = w[-1]
                right = np.median(spec[(w > rb) & (w < rr)])
            except:
                pass

            lw = np.log10(w)
            c = (outwav >= lw[0]) & (outwav <= lw[-1])
            mod = interpolate.splrep(lw, spec, k=1)
            outflux[c, i-1] = interpolate.splev(outwav[c], mod)
            mod = interpolate.splrep(lw, var, k=1)
            outvar[c, i-1] = interpolate.splev(outwav[c], mod)
            # self[i].spec1d['flux'] = spec
            # self[i].spec1d['var'] = var

        outvar[outvar == 0.] = 1.e9
        flux = np.nansum(outflux/outvar, 1) / np.nansum(1./outvar, 1)
        var = np.nansum(1./outvar, 1)**-1
        self.spec1d = ss.Spec1d(wav=outwav, flux=flux, var=var, logwav=True)

    # --------------------------------------------------------------------

    def stitch(self, respfile, method='oldham', respfunc='divide'):
        """

        This method does the two final steps needed to create a 1d spectrum:

         1. Corrects each order for the spectral response function.  If the
            observations were done in a specific manner, then this would be
            flux calibration, but most observations of the standard star
            are taken with a fairly narrow slit and, therefore, true flux
            calibration cannot be done.  Instead, the response correction will
            produce a spectrum for which the relative flux densities are
            correct even though the absolute values are not.

         2. Stitches the 10 spectral orders together into a single spectrum,
            doing the necessary small corrections in the overlap regions to
            make the final spectrum relatively smooth.

        Inputs:
           respfile - File containing the response correction
        Optional inputs:
           method   - Technique for doing the response correction.  Only two
                      possibilities: 'oldham' or 'cdf
                      ('cdf' NOT functional yet)
           respfunc - How to apply the response correction.  Two possibilities:
                      1. 'divide':   divide orders by the response fn (default)
                      2. 'multiply': multiply orders by the response fn (NOT
                          yet implemented)
        """


