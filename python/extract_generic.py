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


""" 
Set the range of valid pixels for each order.
 blue - sets the starting pixel to use for the extraction
 red  - sets the ending pixel to use for the extraction
"""
blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]
arcsecperpix = [0.120,0.127,0.134,0.137,0.144,0.149,0.153,0.158,0.163,0.168]
apsize = []

def clip(arr,nsig=3.):
    a = arr.flatten()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s

#---------------------------------------------------------------------------

def get_ap(slit, B, R, apcent, apnum, wid, order):
    xproj = np.median(slit[:,B:R],1) 
    m,s = clip(xproj)
           
    smooth = ndimage.gaussian_filter(xproj,1)
    if order == 2: # NB: This was 3 in the original file, see note below
        smooth = ndimage.gaussian_filter(xproj[:-30],1)
    x = np.arange(xproj.size)*1. 
    """ 
    The four parameters immediately below are 
    bkgd, amplitude, mean location, and sigma for a Gaussian fit
    """
    fit = np.array([0.,smooth.max(),smooth.argmax(),1.])
    fit = sf.ngaussfit(xproj,fit)[0] 

    cent = fit[2] + apcent[apnum]/arcsecperpix[order-1]
    print cent
    apmax = 0.1 * xproj.max()
    ap = np.where(abs(x-cent)<wid/arcsecperpix[order-1],1.,0.)

    if order<60.:
        plt.subplot(2,5,order)
        plt.plot(x,apmax*ap) # Scale the aperture to easily see it
        plt.plot(x,xproj)
        plt.ylim(-apmax,1.1*xproj.max())
    
    ap = ap.repeat(slit.shape[1]).reshape(slit.shape)
    return ap,fit

#---------------------------------------------------------------------------

def extract(pref, name, frames, apnum, apcent, aplab, stdOrderCorr,
            indir='.', wid=1., wht=False):
    ''' 
    frames = input frame numbers - give a list
    wht    = True gives a Gaussian aperture 
    wid    = how many sigmas wide your aperture is 
    '''

    """ Create lists to hold the HDULists produced for each frame """
    slist = []
    vlist = []
    wlist = []

    for numIndx in range(len(frames)):
        num = frames[numIndx]
        print pref,num
        if type(indir)==list:
            idir = indir[numIndx]
        else:
            idir = indir
        specname = '%s/%s_%04d_bgsub.fits'%(idir,pref,num)
        varname  = specname.replace('bgsub','var')
        d = esi.Esi2d(specname)
        v = esi.Esi2d(varname)

        """ 
        Set up HDUList containers for the extracted spectra (one for each order)
        and the associated variance and wavelength files
        """
        sphdu = pf.PrimaryHDU()
        vphdu = pf.PrimaryHDU()
        wphdu = pf.PrimaryHDU()
        shdu = pf.HDUList([sphdu])
        vhdu = pf.HDUList([vphdu])
        whdu = pf.HDUList([wphdu])

        scales = []
        plt.figure()
        #plt.subplot(111)
        #plt.title('%s Frame %d' % (pref,num))
        """
        In the original code, the spectral order in the following loop ran
        from 1 to 10, since the spectral orders in the HDUList were in 
        HDU1 to HDU10, and the PHDU was in HDU0.  However, when reading the
        data in into a Esi2d structure, the spectral orders will run from
        0 to 9.
        """
        for order in range(10):
            B = blue[order]
            R = red[order]
            slit = d.order[order].data.copy()
            vslit = v.order[order].data.copy()
            vslit[vslit<=0.] = 1e9
            vslit[np.isnan(vslit)] = 1e9
            vslit[np.isnan(slit)] = 1e9
            h = d.order[order].hdr
            x = np.arange(slit.shape[1])*1.
            w = 10**(h['CRVAL1']+x*h['CD1_1']) 

            """ Make the apertures """
            ap,fit = get_ap(slit,B,R,apcent,apnum,wid,order)
            #apsize.append(2.355*fit[3]*arcsecperpix[order])

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

            #ospex[order].append(spec)
            #ovars[order].append(vspec)
            #owave[order] = w
            shdu.append(pf.ImageHDU(spec,name='order%02d'%num))
            vhdu.append(pf.ImageHDU(vspec,name='order%02d'%num))
            if numIndx == 0:
                whdu.append(pf.ImageHDU(w,name='order%02d'%num))
            scales.append(h['CD1_1'])
            
        plt.show()

        """ Add to the lists of HDULists """
        slist.append(shdu)
        vlist.append(vhdu)
        wlist.append(whdu)

    """ Coadd the spectra """
    coadd(slist,vlist,wlist,stdOrderCorr,name)

#---------------------------------------------------------------------------

""" 
-----------------------------------------------------------------------
Start of coadd
-----------------------------------------------------------------------
"""

def coadd(speclist, varlist, wlist, stdOrderCorr, name):

    """ Transfer the information into the expected structures """
    ospex = {} # spectrum
    ovars = {} # variance
    owave = [] # wavelength (one for each order of the echelle)
    print ''
    for order in range(10):
        ospex[order] = []
        ovars[order] = []
    for i in range(len(speclist)):
        for j in range(10):
            ospex[j].append((speclist[i])[j+1].data)
            ovars[j].append((varlist[i])[j+1].data)
            if i==0:
                owave.append((wlist[i])[j+1].data)

    return ospex,ovars,owave

    """
    Now we have a spectrum for each order of the echelle, covering different 
    wavelength ranges...
    """
    scale = 1.7e-5 # of wavelengths. NB. cd1_1 = 1.65e-5, which is about 0.06 arcseconds/pixel
    w0 = np.log10(owave[0][0])
    w1 = np.log10(owave[9][-1]) # total wavelength coverage
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
    for order in range(1,10):
        w = owave[order]
        s = w*0.
        v = w*0.
        for i in range(len(speclist)):
            tmp = ndimage.median_filter(ospex[order][i],7)
            s += tmp/ovars[order][i]
            v += 1./ovars[order][i]

        os = s/v
        ov = 1./v
        r = ov.max()*100

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
    
    return outwave,spec,var
