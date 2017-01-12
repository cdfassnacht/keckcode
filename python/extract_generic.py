"""
This version of the extraction code has started with Lindsay's example
but has been modified to work with the ESI code in esi_spec.py
(in the KeckCDF github repo).
"""

try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import numpy as np,pylab as plt
from scipy import ndimage,interpolate
import special_functions as sf
from spectra import spectools as st
import indexTricks as iT
import spec_simple as ss
import esi_spec as esi
import sys


apsize = []

def extract(pref, name, frames, apnum, apcent, aplab, stdOrderCorr,
            indir='.', wid=1., wht=False, method='oldham', plot_extracted=False,
            apmin=-4., apmax=4.):
    ''' 
    frames = input frame numbers - give a list
    wht    = True gives a Gaussian aperture 
    wid    = how many sigmas wide your aperture is 
    '''

    """ Create list to hold the Spec2d instance for each frame """
    speclist = []

    for numIndx in range(len(frames)):
        num = frames[numIndx]
        print pref,num
        if type(indir)==list:
            idir = indir[numIndx]
        else:
            idir = indir
        specname = '%s/%s_%04d_bgsub.fits'%(idir,pref,num)
        varname  = specname.replace('bgsub','var')
        d = esi.Esi2d(specname,varfile=varname)

        plt.figure()

        """
        The esi_spec code has now been re-written to loop over the orders
        to do the extraction
        """
        d.extract_all(method,apnum,apcent,wid,apmin=apmin,apmax=apmax,
                      plot_extracted=plot_extracted)
        plt.show()

        """
        This scales variable does not seem to get used anywhere below,
        but keep it for now for legacy purposes
        """
        scales = []
        for order in range(10):
            h = d.order[order].hdr
            scales.append(h['CD1_1'])

        """ Add to the lists of Spec2d containers """
        speclist.append(d)

    """ Coadd the spectra """
    print 'Finished the loop'
    coadd(speclist,stdOrderCorr,name,aplab,apnum,pref)

#---------------------------------------------------------------------------

""" 
-----------------------------------------------------------------------
Start of coadd
-----------------------------------------------------------------------
"""

def coadd(speclist, stdOrderCorr, name, aplab, apnum, pref):

    """ Transfer the information into the expected structures """
    ospex = {} # spectrum
    ovars = {} # variance
    owave = {} # wavelength (one for each order of the echelle)
    print ''
    for order in range(1,11):
        ospex[order] = []
        ovars[order] = []
    for i in range(len(speclist)):
        for j in range(1,11):
            ospex[j].append(speclist[i].order[j-1].spec1d.flux)
            ovars[j].append(speclist[i].order[j-1].spec1d.var)
            if i==0:
                owave[j] = speclist[i].order[j-1].spec1d.wav

    """
    Now we have a spectrum for each order of the echelle, covering different 
    wavelength ranges...
    """
    scale = 1.7e-5 # of wavelengths. NB. cd1_1 = 1.65e-5, which is about 0.06 arcseconds/pixel
    w0 = np.log10(owave[1][0])
    w1 = np.log10(owave[10][-1]) # total wavelength coverage
    outwave = np.arange(w0,w1,scale)
    outspec = np.zeros((outwave.size,10))*np.nan
    outvar = outspec.copy()

    corr = np.load(stdOrderCorr)
    right = None
    rb = None
    rr = None

    """
    Sum the different exposures, one order at a time.
    """
    for order in range(2,11):
        w = owave[order]
        s = w*0.
        v = w*0.
        """ This is an inverse-variance weighted sum"""
        for i in range(len(speclist)):
            tmp = ndimage.median_filter(ospex[order][i],7)
            s += tmp/ovars[order][i]
            v += 1./ovars[order][i]

        os = s/v
        ov = 1./v
        r = ov.max()*100

        """ 
        This second pass rejects pixels in individual images that are outliers
         compared to the smoothed summed image.  
        If an individual pixel differs by more than 5 sigma from the smoothed
         summed image then it is rejected.
        """
        for j in range(1):
            os = ndimage.median_filter(os,5)
            s = np.empty((os.size,len(speclist)))
            v = s.copy()
            spec = w*0.
            var = w*0.
            for i in range(len(speclist)):
                s[:,i] = ospex[order][i]
                v[:,i] = ovars[order][i]

            S2N = (os-s.T).T/(v.T+ov).T**0.5 

            c = abs(S2N)<5. 

            s[~c] = np.nan
            v[~c] = np.nan
            os = np.nansum(s/v,1)/np.nansum(1./v,1)
            ov = np.nansum(1./v,1)**-1
            
        """
        ------------------------------------------------------------
        Do response correction
        ------------------------------------------------------------
        """
        spec = os
        var = ov

        w0,w1,mod = corr[order]
        mod = sf.genfunc(w,0.,mod)
        spec /= mod
        var /= mod**2 

        c = np.isnan(spec)
        spec[c] = 0.
        var[c] = 1e9

        c = (w>w0)&(w<w1)

        w = w[c]
        spec = spec[c]
        var = var[c]
        if right is not None:
            left = np.median(spec[(w>rb)&(w<rr)])
            spec *= right/left
            var *= (right/left)**2
        try:
            rb = owave[order+1][0] # blue end is start of next order
            rr = w[-1] # red end is end of this spectrum
            right = np.median(spec[(w>rb)&(w<rr)]) 
        except:
            pass

        lw = np.log10(w)
        c = (outwave>=lw[0])&(outwave<=lw[-1])
        mod = interpolate.splrep(lw,spec,k=1)
        outspec[c,order-1] = interpolate.splev(outwave[c],mod)
        mod = interpolate.splrep(lw,var,k=1)
        outvar[c,order-1] = interpolate.splev(outwave[c],mod)
        
    spec = np.nansum(outspec/outvar,1)/np.nansum(1./outvar,1)

    var = np.nansum(1./outvar,1)**-1
    ow,s,v = outwave,spec,var
    plt.figure()
    tmpspec = ss.Spec1d(wav=10.**ow,flux=s,var=v)
    tmpspec.plot()
    plt.xlim(4150.,10300.)
    plt.ylim(-0.05,0.6)
    plt.xlabel('Observed wavelength')
    plt.suptitle(name)
    #outplt = '%s_%s.png' % (name,aplab[apnum])
    #plt.savefig(outplt)
    plt.show()
    outname = '%s_spec_%s.fits' % (name,aplab[apnum])
    hdu  = pf.HDUList()
    phdu = pf.PrimaryHDU()
    hdr = phdu.header
    hdr['object'] = pref
    outwv   = pf.ImageHDU(ow,name='wavelength')
    outflux = pf.ImageHDU(s,name='flux')
    outvar  = pf.ImageHDU(v,name='variance')
    hdu.append(phdu)
    hdu.append(outwv)
    hdu.append(outflux)
    hdu.append(outvar)
    hdu.writeto(outname,clobber=True)
    
    #return outwave,spec,var
