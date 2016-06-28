import numpy,pylab,cPickle
from scipy import ndimage
import special_functions as sf
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

#d = pyfits.open('EEL1248+4711_0033_2d.fits')
#v = pyfits.open('EEL1248+4711_0033_var.fits')

blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]

left = [1000,1000,700,500,300,0,0,0,0,0]
right = [3800,4000,4000,-1,-1,-1,-1,-1,-1,2500]

def clip(arr,nsig=3.):
    a = arr.flatten()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s

pref = 'BD284211'
nums = [56,57,58]

wid = 5.

ospex = {}
ovars = {}
owave = {}
for order in range(1,11):
    ospex[order] = []
    ovars[order] = []

for num in nums:
    d = pyfits.open('%s_%04d_bgsub.fits'%(pref,num))
    v = pyfits.open('%s_%04d_var.fits'%(pref,num))

    omean = 0
    for order in range(1,11):
        B = blue[order-1]
        R = red[order-1]
        slit = d[order].data.copy()
        vslit = v[order].data.copy()
        vslit[vslit<0.] = 1e9
        h = d[order].header
        x = numpy.arange(slit.shape[1])*1.
        w = 10**(h['CRVAL1']+x*h['CD1_1'])
    
        slice = numpy.median(slit[:,B:R],1)
        m,s = clip(slice)
    
        smooth = ndimage.gaussian_filter(slice,1)
        x = numpy.arange(slice.size)*1.
        fit = numpy.array([0.,smooth.max(),smooth.argmax(),1.])
        fit = sf.ngaussfit(slice,fit)[0]
        ap = numpy.where(abs(x-fit[2])<wid*2.355*fit[3],1.,0.)
        #ap /= ap.sum()
        spec = (slit.T*ap).sum(1)
        vspec = (vslit.T*ap**2).sum(1)
        vspec /= numpy.median(spec)**2
        spec /= numpy.median(spec)
        omean += spec.mean()        
        ospex[order].append(spec)
        ovars[order].append(vspec)
        owave[order] = w
    omean /= 10.
#    for order in range(1,11):
#        ospex[order][-1] /= omean
#        ovars[order][-1] /= omean**2

#for order in range(1,11):
#    for i in range(len(nums)):
#        pylab.plot(ospex[order][i])
#    pylab.show()

#wmodel,fmodel = numpy.loadtxt('fluxCal.dat').T
wmodel,fmodel,jmodel = numpy.loadtxt('../Raw/standard/BD284211.dat').T
fmodel /= fmodel.mean()

from scipy import interpolate,ndimage
model = interpolate.splrep(wmodel,fmodel)

corrections = {}
for order in range(1,11):
    w = owave[order]
    s = w*0.
    v = w*0.
    for i in range(len(nums)):
        s += ospex[order][i]/ovars[order][i]
        v += 1./ovars[order][i]
    os = s/v
    ov = 1./v
    r = ov.max()*100

    spec = w*0.
    var = w*0.
    for i in range(len(nums)):
        s = ospex[order][i]
        v = ovars[order][i]
        c = abs(os-s)<215*(v+ov)**0.5
        spec[c] += (s/v)[c]
        var[c] += 1./v[c]
    spec = spec/var
    var = 1./var
    mspec = interpolate.splev(w,model)
    ratio = spec/mspec
    smo = ndimage.gaussian_filter(ratio,5.)
    #smoModel = interpolate.splrep(w,ratio,s=2*w.size*res)
    knots = numpy.linspace(w[1],w[-2],w.size/100)
    smoModel = interpolate.splrep(w,ratio,t=knots)
    smo = interpolate.splev(w,smoModel)

    W = w[left[order-1]:right[order-1]]
    R = ratio[left[order-1]:right[order-1]]
    c = (abs(W-6890)>30.)&(abs(W-7655)>65)
    c = c&(abs(W-4101)>10)&(abs(W-4340)>10)&(abs(W-4861)>10)&(abs(W-6560)>10)
    c = c&(abs(W-4792)>15)
    c = c&numpy.isfinite(R)
    fitData = numpy.array([W[c],R[c]]).T    
    for i in range(5):
        mod2 = sf.lsqfit(fitData,'chebyshev',7)
        smo2 = sf.genfunc(fitData[:,0],0.,mod2)
        
        res = fitData[:,1]-smo2
        m,s = clip(res)
        pylab.plot(res)
        c = abs(res)<2.5*s
        c = c&(res>-2*s)
        fitData = fitData[c]
    smo2 = sf.genfunc(w,0.,mod2)
    corrections[order] = [W[0],W[-1],mod2]
    continue
    pylab.plot(w,ratio)
    pylab.plot(w,spec)
    pylab.plot(w,mspec)
    pylab.plot(w,smo)
    pylab.plot(w,smo2)
    pylab.xlim([W[0],W[-1]])
    pylab.show()
    c = numpy.isfinite(spec)
    

    #pylab.plot(w,spec)
    #pylab.plot(w,var**0.5)
    #pylab.plot(w,spec/var**0.5)

f = open('orderCorr_BD284211.dat','wb')
cPickle.dump(corrections,f,2)
f.close()
