import pyfits,numpy,pylab
from scipy import ndimage,interpolate
import special_functions as sf


blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]

def clip(arr,nsig=3.):
    a = arr.flatten()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s

def extract(pref,nums,wid=1.,loc=None,wht=False):
    ospex = {}
    ovars = {}
    owave = {}
    for order in range(1,11):
        ospex[order] = []
        ovars[order] = []

    for numIndx in range(len(nums)):
        num = nums[numIndx]
        d = pyfits.open('%s_%04d_bgsub.fits'%(pref,num))
        v = pyfits.open('%s_%04d_var.fits'%(pref,num))

        scales = []
        for order in range(1,11):
            B = blue[order-1]
            R = red[order-1]
            slit = d[order].data.copy()
            vslit = v[order].data.copy()
            vslit[vslit<=0.] = 1e9
            vslit[numpy.isnan(vslit)] = 1e9
            vslit[numpy.isnan(slit)] = 1e9
            h = d[order].header
            x = numpy.arange(slit.shape[1])*1.
            w = 10**(h['CRVAL1']+x*h['CD1_1'])

            slice = numpy.median(slit[:,B:R],1)
            m,s = clip(slice)

            smooth = ndimage.gaussian_filter(slice,1)
            x = numpy.arange(slice.size)*1.
            if loc is not None:
                pos = loc[numIndx]+9.47933-2.1065*(order-1)
                fit = numpy.array([0.,smooth[pos],pos,1.])
            else:
                fit = numpy.array([0.,smooth.max(),smooth.argmax(),1.])
            fit = sf.ngaussfit(slice,fit)[0]
            
            if wht:
                fit[0] = 0.
                fit[3] *= wid 
                ap = sf.ngauss(x,fit)
            else:
                ap = numpy.where(abs(x-fit[2])<wid*2.355*fit[3],1.,0.)
            ap = ap.repeat(slit.shape[1]).reshape(slit.shape)
            print fit[2]

            ap[vslit>=1e8] = 0.
            ap = ap/ap.sum(0)


#            spec = (slit.T*ap).sum(1)
#            vspec = (vslit.T*ap**2).sum(1)
            spec = (slit*ap**2).sum(0)
            vspec = (vslit*ap**4).sum(0)
            vspec /= numpy.median(spec)**2


            if order>=25:
                pylab.imshow(vslit[:,2270:2420],origin='bottom')
                pylab.colorbar()
                pylab.figure()
                pylab.imshow((vslit*ap)[:,2270:2420],origin='bottom')
                pylab.colorbar()
                pylab.figure()
                pylab.plot(vspec)
                pylab.figure()
                pylab.plot(spec)
                pylab.show()


            spec /= numpy.median(spec)
            ospex[order].append(spec)
            ovars[order].append(vspec)
            owave[order] = w
            scales.append(h['CD1_1'])

    for order in range(1,11):
        s = ospex[order]
        for i in range(len(nums)):
            print s[i][s[i]==0].size,numpy.isnan(s[i]).sum()

    scale = 1.7e-5
    w0 = numpy.log10(owave[1][0])
    w1 = numpy.log10(owave[10][-1])
    outwave = numpy.arange(w0,w1,scale)
    outspec = numpy.zeros((outwave.size,10))*numpy.nan
    outvar = outspec.copy()

    corr = numpy.load('orderCorr_BD284211.dat')
    right = None
    rb = None
    rr = None
    for order in range(1,11):
        w = owave[order]
        s = w*0.
        v = w*0.
        for i in range(len(nums)):
            tmp = ndimage.median_filter(ospex[order][i],7)
            s += tmp/ovars[order][i]
            v += 1./ovars[order][i]

        if order==70:
            for i in range(len(nums)):
                pylab.plot(owave[order],ospex[order][i])
                pylab.plot(owave[order],ovars[order][i]**0.5,ls='--')
            pylab.show()
        os = s/v

        ov = 1./v
        r = ov.max()*100

        for j in range(1):
            os = ndimage.median_filter(os,5)
            s = numpy.empty((os.size,len(nums)))
            v = s.copy()
            spec = w*0.
            var = w*0.
            for i in range(len(nums)):
                s[:,i] = ospex[order][i]
                v[:,i] = ovars[order][i]

            S2N = (os-s.T).T/(v.T+ov).T**0.5

            c = abs(S2N)<5.

            s[~c] = numpy.nan
            v[~c] = numpy.nan
            os = numpy.nansum(s/v,1)/numpy.nansum(1./v,1)
            ov = numpy.nansum(1./v,1)**-1
            print c[s==0]
        spec = os
        var = ov

        w0,w1,mod = corr[order]
        mod = sf.genfunc(w,0.,mod)
        spec /= mod
        var /= mod**2

        c = numpy.isnan(spec)
        spec[c] = 0.
        var[c] = 1e9

        c = (w>w0)&(w<w1)

        w = w[c]
        spec = spec[c]
        var = var[c]
        if right is not None:
            left = numpy.median(spec[(w>rb)&(w<rr)])
            spec *= right/left
            var *= (right/left)**2
        try:
            rb = owave[order+1][0]
            rr = w[-1]
            right = numpy.median(spec[(w>rb)&(w<rr)])
        except:
            pass

        lw = numpy.log10(w)
        c = (outwave>=lw[0])&(outwave<=lw[-1])
        mod = interpolate.splrep(lw,spec,k=1)
        outspec[c,order-1] = interpolate.splev(outwave[c],mod)
        mod = interpolate.splrep(lw,var,k=1)
        outvar[c,order-1] = interpolate.splev(outwave[c],mod)
        #pylab.plot(10**outwave,outspec[:,order-1])
        #pylab.plot(10**outwave,outvar[:,order-1]**0.5,ls='--')
    #pylab.show()
    #print outspec.shape
    #df
    spec = numpy.nansum(outspec/outvar,1)/numpy.nansum(1./outvar,1)

    var = numpy.nansum(1./outvar,1)**-1

    return outwave,spec,var
    pylab.plot(10**(outwave),spec/var**0.5)
    pylab.show()

from scipy import ndimage
"""
ow,s,v = extract('STR_2938832465',[28,29],loc=[54,88],wid=3)
s = ndimage.gaussian_filter(s,5.5)
pylab.plot(10**ow,s)
ow2,s2,v2 = extract('STR_2938832465',[28,29],loc=[30,65],wid=3)
s2 = ndimage.gaussian_filter(s2,5.5)
pylab.plot(10**ow2,s2)
pylab.show()
df
"""
import indexTricks as iT
ow,s,v = extract('EEL_J1248+4711',[33,34,35],wid=3)
ow,s,v=ow[:-33],s[:-33],v[:-33]


pylab.plot(10**ow,s)
pylab.xlim(4000.,9000.)
pylab.ylim(0.,1.6)
pylab.xlabel('Observed Wavelength')

print ow.size

pylab.show()
