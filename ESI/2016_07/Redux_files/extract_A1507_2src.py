try:
    import pyfits as py
except:
    from astropy.io import fits as py
import numpy as np,pylab as pl
from scipy import ndimage,interpolate
import special_functions as sf
from spectra import spectools as st
import indexTricks as iT

""" 
Set the range of valid pixels for each order.
 blue - sets the starting pixel to use for the extraction
 red  - sets the ending pixel to use for the extraction
"""
blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
red  = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]
arcsecperpix = [0.120,0.127,0.134,0.137,0.144,0.149,0.153,0.158,0.163,0.168]
apsize = []
stdOrderCorr = 'orderCorr_HZ44.dat'

def clip(arr,nsig=3.):
    a = arr.flatten()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s

def extract(pref,nums,apnum,wid=1.,wht=False):
    ''' 
    nums = input frame numbers - give a list
    wht  = True gives a Gaussian aperture 
    wid  = how many sigmas wide your aperture is 
    '''
    ospex = {} # spectrum
    ovars = {} # variance
    owave = {} # wavelength (one for each order of the echelle)
    for order in range(1,11):
        ospex[order] = []
        ovars[order] = []

    for numIndx in range(len(nums)):
        num = nums[numIndx]
        print pref,num
        d = py.open('%s_%04d_bgsub.fits'%(pref,num))
        v = py.open('%s_%04d_var.fits'%(pref,num))

        scales = []
        pl.figure()
        #pl.subplot(111)
        #pl.title('%s Frame %d' % (pref,num))
        for order in range(1,11):
            B = blue[order-1]
            R = red[order-1]
            slit = d[order].data.copy()
            vslit = v[order].data.copy()
            vslit[vslit<=0.] = 1e9
            vslit[np.isnan(vslit)] = 1e9
            vslit[np.isnan(slit)] = 1e9
            h = d[order].header
            x = np.arange(slit.shape[1])*1.
            w = 10**(h['CRVAL1']+x*h['CD1_1']) 

            slice = np.median(slit[:,B:R],1) 
            m,s = clip(slice)
           
            smooth = ndimage.gaussian_filter(slice,1)
            if order == 3.:
                smooth = ndimage.gaussian_filter(slice[:-30],1)
            x = np.arange(slice.size)*1. 
            fit = np.array([0.,smooth.max(),smooth.argmax(),1.])
            fit = sf.ngaussfit(slice,fit)[0] 

            cent = fit[2] + apcent[apnum]/arcsecperpix[order-1]
            print cent
            ap = np.where(abs(x-cent)<wid/arcsecperpix[order-1],1.,0.)

            if order<60.:
                pl.subplot(2,5,order)
                pl.plot(x,ap)
                pl.plot(x,slice)
            ap = ap.repeat(slit.shape[1]).reshape(slit.shape) 
                

            ap[vslit>=1e8] = 0.
            ap = ap/ap.sum(0)
            ap[np.isnan(ap)] = 0.
            slit[np.isnan(slit)] = 0.
            apsize.append(2.355*fit[3]*arcsecperpix[order-1])
            spec = (slit*ap**2).sum(0) 
            vspec = (vslit*ap**4).sum(0)
            vspec /= np.median(spec)**2
            spec /= np.median(spec)

            ospex[order].append(spec)
            ovars[order].append(vspec)
            owave[order] = w
            scales.append(h['CD1_1'])

        pl.show()

    # ok. So now we have a spectrum for each order of the echelle, covering different wavelength ranges...
    scale = 1.7e-5 # of wavelengths. NB. cd1_1 = 165e-5, which is about 0.06 arcseconds/pixel
    w0 = np.log10(owave[1][0])
    w1 = np.log10(owave[10][-1]) # total wavelength coverage
    outwave = np.arange(w0,w1,scale)
    outspec = np.zeros((outwave.size,10))*np.nan
    outvar = outspec.copy()

    corr = np.load(stdOrderCorr)
    right = None
    rb = None
    rr = None
    for order in range(2,11): # purpose is to sum exposures -- here, for each order on the echelle
        w = owave[order]
        s = w*0.
        v = w*0.
        for i in range(len(nums)):
            tmp = ndimage.median_filter(ospex[order][i],7)
            s += tmp/ovars[order][i]
            v += 1./ovars[order][i]

        os = s/v
        ov = 1./v
        r = ov.max()*100

        for j in range(1):
            os = ndimage.median_filter(os,5)
            s = np.empty((os.size,len(nums)))
            v = s.copy()
            spec = w*0.
            var = w*0.
            for i in range(len(nums)):
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
    pl.figure()
    pl.subplot(211)
    pl.plot(10**ow,s,'k')
    pl.xlim(4000.,9000.)
    pl.ylim([-0.5,1])
    pl.subplot(212)
    pl.plot(10**ow,v,'r')
    pl.xlim(4000.,9000.)
    pl.ylim([0,0.05])
    pl.xlabel('Observed wavelength')
    pl.suptitle(name)
    outplt = '%s_%s.png' % (name,aplab[apnum])
    pl.savefig(outplt)
    pl.show()
    outname = '%s_spec_ap_%s.fits' % (name,aplab[apnum])
    hdu  = py.HDUList()
    phdu = py.PrimaryHDU()
    hdr = phdu.header
    hdr['object'] = pref
    outwv   = py.ImageHDU(ow,name='wavelength')
    outflux = py.ImageHDU(s,name='flux')
    outvar  = py.ImageHDU(v,name='variance')
    hdu.append(phdu)
    hdu.append(outwv)
    hdu.append(outflux)
    hdu.append(outvar)
    hdu.writeto(outname,clobber=True)
    
    return outwave,spec,var

fullname = 'A1507-1442'
name='A1507'
inframes = [44,45,46]
apcent = [0.,-4.25]
aplab = ['A','B']

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,inframes,jj,wid=1.0)
