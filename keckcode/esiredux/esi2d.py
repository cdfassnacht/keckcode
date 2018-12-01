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
from scipy import ndimage,interpolate
from matplotlib import pyplot as plt
from astropy.io import fits as pf
import special_functions as sf
from ccdredux import sigma_clip
from specim import specfuncs as ss

"""
============================== Esi2d class ==============================
"""

class Esi2d(ss.Spec2d):
    """
    A class for ESI 2D spectra, for which Matt's / Lindsay's calibration
    code produces a multi-extension fits files containing 10 2D spectra,
    one for each of the 10 orders produced by the spectrograph
    """

    def __init__(self, infile, varfile=None, startord=None):
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
        dtype = [('name','S8'), ('pixscale',float), ('blue',int), ('red',int)]
        self.orderinfo = np.array([
                ('order 1', 0.120, 1500, 3000),
                ('order 2', 0.127, 1400, 3400),
                ('order 3', 0.134, 1300, 3700),
                ('order 4', 0.137, 1200,   -1),
                ('order 5', 0.144, 1100,   -1),
                ('order 6', 0.149,  900,   -1),
                ('order 7', 0.153,  600,   -1),
                ('order 8', 0.158,  200,   -1),
                ('order 9', 0.163,    0,   -1),
                ('order 10',0.168,    0,   -1),
                ], dtype=dtype)
        self.orderinfo = self.orderinfo.view(np.recarray)

        """ Open the multiextension fits file that contains the 10 orders """
        self.infile = infile
        self.hdu = pf.open(infile)
        print ''
        print 'Science file:  %s' % self.infile

        """ If there is an external variance file, open that """
        self.extvar = None
        if varfile is not None:
            self.extvar = pf.open(varfile)
            print 'Variance file: %s' % varfile

        """ Load each order into its own Spec2d container """
        self.order = []
        if startord is not None:
            ostart = startord
        else:
            ostart = 0
        print ''
        print 'Order  Shape    Dispaxis'
        print '----- --------- --------'
        for i in range(ostart, 10):
            hext = i + 1 - ostart
            if self.extvar is not None:
                tmpspec = ss.Spec2d(self.hdu, hext=hext,
                                    extvar=self.extvar,verbose=False,
                                    logwav=True, fixnans=False)
            else:
                tmpspec = ss.Spec2d(self.hdu, hext=hext, verbose=False,
                                    logwav=True, fixnans=False)
            #if self.extvar is not None:
            #    tmpspec = ss.Spec2d(None,hdulist=self.hdu,hext=i+1,
            #                        extvar=self.extvar,verbose=False,
            #                        logwav=True,fixnans=False)
            #else:
            #    tmpspec = ss.Spec2d(None,hdulist=self.hdu,hext=i+1,verbose=False,
            #                        logwav=True,fixnans=False)
            print ' %2d   %dx%d     %s' % \
                ((i+1),tmpspec.data.shape[1],tmpspec.data.shape[0],
                 tmpspec.dispaxis)
            self.order.append(tmpspec)

    #-----------------------------------------------------------------------

    def plot_profiles(self, bgsub=True, showfit=False, fitrange=None):
        """

        Plots, in one figure, the spatial profiles for all the 10 orders

        Optional inputs
          bgsub   - Set to true (the default) if the data have had the sky
                     subtracted already at this point.  If this is the case,
                     then the "sky" level in each profile should be close to
                     zero and, therefore, the subplots can be displayed in a
                     way that does not require y-axis labels for each profile
          showfit - Show the fit to the profiles? Default=False
        """

        if bgsub:
            normspec = True
            plt.subplots_adjust(wspace=0.001)
        else:
            normspec = False
        plt.subplots_adjust(hspace=0.001)

        """ Set the common widths of the plots """
        #xmax = 0
        #for i in range(1,11):
        #    if self.hdu[i].data.shape[0] > xmax:
        #        xmax = self.hdu[i].data.shape[0]
            

        """ Set up the figure and the full-sized frame for the final labels """
        fig = plt.gcf()
        ax = fig.add_subplot(111)
        plt.setp(ax.get_xticklabels(),visible=False)
        plt.setp(ax.get_yticklabels(),visible=False)

        for i in range(10):
            axi = fig.add_subplot(2,5,(i+1))
            B = self.orderinfo.blue[i]
            R = self.orderinfo.red[i]
            if showfit:
                self.order[i].spatial_profile(normalize=normspec,title=None,
                                              fit=self.order[i].p0,
                                              pixrange=[B,R])
            else:
                self.order[i].spatial_profile(normalize=normspec,title=None,
                                              pixrange=[B,R])
            plt.xlim(-1,120)
            if normspec:
                plt.ylim(-0.1,1.1)
            if i==0 or i==5:
                pass
            else:
                plt.setp(axi.get_yticklabels(),visible=False)
                axi.set_ylabel('',visible=False)
            if i<5:
                plt.setp(axi.get_xticklabels(),visible=False)
            axi.set_xlabel('',visible=False)
            axi.annotate('%d'%(i+1),(10.,0.95))

        ax.set_title('Spatial Profiles')
        ax.set_xlabel('Spatial Direction')
        ax.set_ylabel('Relative Flux')
        ax.xaxis.set_label_coords(0.5, -0.05)
        ax.yaxis.set_label_coords(-0.03,0.5)

    #---------------------------------------------------------------------------

    def get_ap_oldham(self, slit, B, R, apcent, apnum, wid, order, 
                      doplot=True):
        """
        Defines a uniform aperture as in the example ESI extraction scripts
        from Lindsay
        """

        xproj = np.median(slit[:,B:R],1) 
        m,s = sigma_clip(xproj)
           
        smooth = ndimage.gaussian_filter(xproj,1)
        if order == 2: # NB: This was 3 in the original file, see note below
            smooth = ndimage.gaussian_filter(xproj[:-30],1)
        x = np.arange(xproj.size)*1. 

        """ 
        The four parameters immediately below are the initial guesses
        bkgd, amplitude, mean location, and sigma for a Gaussian fit
        """
        fit = np.array([0.,smooth.max(),smooth.argmax(),1.])
        fit = sf.ngaussfit(xproj,fit)[0] 

        cent = fit[2] + apcent[apnum]/self.orderinfo.pixscale[order]
        print cent
        apmax = 0.1 * xproj.max()
        ap = np.where(abs(x-cent)<wid/self.orderinfo.pixscale[order],1.,0.)

        if doplot & order<60:
            plt.subplot(2,5,(order+1))
            plt.plot(x,apmax*ap) # Scale the aperture to easily see it
            plt.plot(x,xproj)
            plt.ylim(-apmax,1.1*xproj.max())
            plt.axvline(cent, color='k', ls='dotted')
    
        ap = ap.repeat(slit.shape[1]).reshape(slit.shape)
        return ap,fit

    #-----------------------------------------------------------------------

    def extract_oldham(self, order, apcent, apnum, wid):
        """
        Implements Lindsay Oldham's (or perhaps Matt Auger's) method
         for extracting the spectrum from an individual spectral order
         on ESI.

        Inputs:
          order  - the requested spectral order
        """

        """ 
        Make temporary copies of the science and variance spectra, and
        get information about their properties
        """
        spec  = self.order[order]
        slit  = spec.data.copy()
        vslit = spec.vardata.copy()
        vslit[vslit<=0.] = 1e9
        vslit[np.isnan(vslit)] = 1e9
        vslit[np.isnan(slit)] = 1e9
        h = spec.hdr
        x = np.arange(slit.shape[1])*1.
        w = 10**(h['CRVAL1']+x*h['CD1_1']) 

        """ Make the apertures """
        B = self.orderinfo.blue[order]
        R = self.orderinfo.red[order]
        ap, fit = self.get_ap_oldham(slit, B, R, apcent, apnum, wid, order)

        """ Set up to do the extraction, including by normalizing ap """
        ap[vslit>=1e8] = 0.
        ap = ap/ap.sum(0)
        ap[np.isnan(ap)] = 0.
        slit[np.isnan(slit)] = 0.

        """ Extract the spectrum (I do not understand the weighting here) """
        spec = (slit*ap**2).sum(0) 
        vspec = (vslit*ap**4).sum(0)

        """ 
        Normalize the spectrum and its associated variance 
        Need to do the variance first, or else you would incorrectly
        be using the normalized spectrum to normalize the variance
        """
        vspec /= np.median(spec)**2
        spec /= np.median(spec)

        """ Save the output in a Spec1d container """
        #newspec = ss.Spec1d(wav=w,flux=spec,var=vspec)
        self.order[order].spec1d = ss.Spec1d(wav=w,flux=spec,var=vspec)

    #-----------------------------------------------------------------------

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
        w0 = log10(self.order[0].spec1d['wav'][0])
        w1 = log10(self.order[9].spec1d['wav'][-1])
        outwav  = np.arange(w0,w1,scale)
        outflux = np.ones((outwav.size,10))*np.nan
        outvar  = outflux.copy()

        for i in range(1,10):
            """ 
            Create the response correction for this order and then correct
            the flux and the variance
            """
            w    = self.order[i].spec1d['wav']
            w0,w1,mod = corr[i+1]
            mod = sf.genfunc(w,0.,mod)
            spec = self.order[i].spec1d['flux'] / mod
            var  = self.order[i].spec1d['var'] / mod**2 

            """ Mask out NaNs """
            mask = np.isnan(spec)
            spec[mask] = 0.
            var[mask] = 1e9

            """ 
            Only correct for the ranges over which the correction is valid 
            """
            mask = (w>w0)&(w<w1)
            w    = w[mask]
            spec = spec[mask]
            var  = var[mask]
            if right is not None:
                left = np.median(spec[(w>rb)&(w<rr)])
                spec *= right/left
                var *= (right/left)**2
            try:
                """ 
                Find the overlap region
                 - blue end (rb) is the start of the next order
                 - red end (rr) is the end of this current spectrum
                """
                rb = self.order[i+1].spec1d['wav'][0]
                rr = w[-1]
                right = np.median(spec[(w>rb)&(w<rr)]) 
            except:
                pass
            
            lw = np.log10(w)
            c = (outwav>=lw[0])&(outwav<=lw[-1])
            mod = interpolate.splrep(lw,spec,k=1)
            outflux[c,i-1] = interpolate.splev(outwav[c],mod)
            mod = interpolate.splrep(lw,var,k=1)
            outvar[c,i-1] = interpolate.splev(outwav[c],mod)
            #self.order[i].spec1d['flux'] = spec
            #self.order[i].spec1d['var'] = var

        outvar[outvar==0.] = 1.e9
        flux = np.nansum(outflux/outvar,1)/np.nansum(1./outvar,1)
        var  = np.nansum(1./outvar,1)**-1
        self.spec1d = ss.Spec1d(wav=outwav,flux=flux,var=var,logwav=True)
            

    #-----------------------------------------------------------------------

    def plot_spec1d(self, xmin=3840., xmax=10910., ymin=-0.2, ymax=5.,
                    color='b'):
        """
        Plots the 10 extracted 1d spectra on a single plot
        """
        for i in range(10):
            self.order[i].spec1d.plot(color=color)
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)

    #-----------------------------------------------------------------------

    def extract_all(self, method='oldham', apnum=0, apcent=[0.,], wid=1.0, 
                    apmin=-1., apmax=1.,
                    plot_profiles=True, plot_traces=False, plot_extracted=True, 
                    xmin=3840., xmax=10910., ymin=-0.2, ymax=5.,
                    startord=None):
        """
        Goes through each of the 10 orders on the ESI spectrograph and
        extracts the spectrum via one of two procedures:
        1. Using the Spec2d class methods from spec_simple.py (method='cdf')
          This approach does the two-step extraction procedure using the
          relevant methods in the Spec2d class in spec_simple.py, namely
          (1) find_and_trace and (2) extract_spectrum.
          xxxxx
        2. Using the code from Lindsay Oldham (or perhaps Matt Auger) implemented
          above (method='oldham')
          In this approach, the two step procedure is (1) get_ap_oldham, and
          (2) extract_oldham
          xxxxx
        """

        """ 
        Set up for plotting profiles, which needs to be done here since the
        'oldham' and 'cdf' techniques do the plotting at different points
        """
        if plot_profiles:
            plt.figure(1)
            plt.clf()

        """ 
        Extract the spectra 
        """
        if startord is not None:
            ostart = startord
        else:
            ostart = 0
        for i in range(ostart, 10):
            if method == 'cdf':
                if plot_traces:
                    plt.figure(i+3)
                    plt.clf()
                """ 
                For now assume that the data reduction has properly rectified
                the 2D spectra, so that it is not necessary to trace a varying
                position and width as the position changes
                """
                muorder = -1
                sigorder = -1
                print ''
                print '=================================================='
                print 'Order: %d' % (i+1)
                print ''
                self.order[i].apmin = apmin / self.orderinfo.pixscale[i]
                self.order[i].apmax = apmax / self.orderinfo.pixscale[i]
                B = self.orderinfo.blue[i]
                R = self.orderinfo.red[i]
                self.order[i].find_and_trace(doplot=plot_traces,muorder=muorder,
                                             sigorder=sigorder,fitrange=[B,R],
                                             verbose=False)
                self.order[i].extract()
                print self.order[i].p0

                """
                Also normalize the flux and variance of the extracted
                spectrum in the same way that the Oldham extraction does
                """
                medflux = np.median(self.order[i].spec1d['flux'])
                self.order[i].spec1d['flux'] /= medflux
                self.order[i].spec1d['var'] /= medflux**2


            elif method == 'oldham':
                self.extract_oldham(i,apcent,apnum,wid)

        """ Plot all the spatial profiles, along with the profile fits """
        if plot_profiles and method == 'cdf':
            self.plot_profiles(showfit=True)

        """
        Plot the extracted spectra
        """
        if plot_extracted:
            plt.figure(2)
            plt.clf()
            self.plot_spec1d()

    #-----------------------------------------------------------------------    
        
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
                      possibilities: 'oldham' or 'cdf' ('cdf' NOT functional yet)
           respfunc - How to apply the response correction.  Two possibilities:
                      1. 'divide':   divide orders by the response fn (default)
                      2. 'multiply': multiply orders by the response fn (NOT
                          yet implemented)
        """

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
        print ''
        print 'Science file:  %s' % self.infile

        # """ If there is an external variance file, open that """
        # self.extvar = None
        # if varfile is not None:
        #     self.extvar = pf.open(varfile)
        #     print 'Variance file: %s' % varfile
        # 
        # """ Load each order into its own Spec2d container """
        # self.order = []
        # print ''
        # print 'Order  Shape    Dispaxis'
        # print '----- --------- --------'
        # for i in range(10):
        #     if self.extvar is not None:
        #         tmpspec = ss.Spec2d(None,hdulist=self.hdu,hext=i+1,
        #                             extvar=self.extvar,verbose=False,
        #                             logwav=True,fixnans=False)
        #     else:
        #         tmpspec = ss.Spec2d(None,hdulist=self.hdu,hext=i+1,verbose=False,
        #                             logwav=True,fixnans=False)
        #     print ' %2d   %dx%d     %s' % \
        #         ((i+1),tmpspec.data.shape[1],tmpspec.data.shape[0],
        #          tmpspec.dispaxis)
        #     self.order.append(tmpspec)

