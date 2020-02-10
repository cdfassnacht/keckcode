import esi
from esi.biastrim import make_bias,biastrim

from esi.flat import *
from esi.straighten import startrace,straighten,fullSolution,getOrders
from esi.wavesolve import solve

from spectra import spectools,offset,measure_width
from spectra.extract import extract
import special_functions as sf

from pickle import dump,load
from math import floor,ceil,fabs
import os,sys

import numpy,scipy,cPickle
from scipy import stats,interpolate,ndimage
from scipy import io as sio
try:
    import pyfits
except:
    from astropy.io import fits as pyfits


def message(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def clip(arr,nsig=3.5):
    a = arr.flatten()
    a = a[numpy.isfinite(a)]
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<nsig*s]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def clip2(arr,nsig=3.,edge=0.01):
    a = arr.flatten()
    a.sort()
    a = a[a.size*edge:a.size*(1.-edge)]
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<nsig*s]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def crFind(img,var,nsig=10.,sigfrac=0.3):
    simg = img/var**0.5
    deriv = ndimage.sobel(simg)
    deriv = ndimage.sobel(deriv)

    mean,std = clip(deriv[deriv!=0.])
    crmask = numpy.where(abs(deriv)>15*std,1,0)
    return crmask
    thresh = nsig*var**0.5
    n = 5
    blkimg = img.repeat(n,0).repeat(n,1)
    deriv = ndimage.laplace(blkimg)*-1.
    deriv[deriv<0] = 0.
    d = iT.resamp(deriv,n)
    m = numpy.where((d>thresh)&(img>thresh),1,0)
    bmap = ndimage.maximum_filter(m,5)
    cond = (bmap==1)&(d>thresh*sigfrac)&(img>thresh*sigfrac)
    m[cond] = 1
    return m


def bgsub(dir,inname,out_prefix,cal_prefix):
    # Where things begin and end....
    blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
    red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]

    readvar = 7.3

    bias = pyfits.open(cal_prefix+"_bias.fits")[0].data.astype(scipy.float32)
    bpm = pyfits.open(cal_prefix+"_bpm.fits")[0].data.astype(scipy.float32)
    flat = pyfits.open(cal_prefix+"_norm.fits")[0].data.astype(scipy.float32)

    orders,y_soln,wideorders = numpy.load(cal_prefix+"_ycor.dat")
    fullsoln = numpy.load(cal_prefix+"_full.dat")

    try:
        masks = numpy.loadtxt(cal_prefix+"_bpmOrders.dat")
        print('BPM opened')
    except:
        print("Making BPM")
        mask = numpy.where(bpm==1.,1.,0.)
        mask = getOrders(mask,y_soln,orders,wideorders)
        for i in range(len(mask)):
            mask[i][0] = numpy.where(mask[i][0]<0.7,numpy.nan,1.)
        f = open(cal_prefix+"_bpmOrders.dat",'wb')
        cPickle.dump(mask,f,2)
        f.close()

    hdu = pyfits.open(dir+inname)
    data = hdu[0].data.copy()
    data = biastrim(data,bias,bpm)/flat


    try:
        back = pyfits.open(out_prefix+"_bg.fits")[0].data.copy()
    except:
        strt = straighten(data,y_soln,orders,wideorders)
        back = scipy.zeros(strt.shape)
        for indx in range(len(orders)):
            i,j = orders[indx]
            bg = stats.stats.nanmedian(strt[i:j],axis=0)
            back[i:j] += bg
        back = curve(back,y_soln,orders,wideorders)
        pyfits.PrimaryHDU(back).writeto(out_prefix+"_bg.fits")

    try:
        strt = pyfits.open(out_prefix+"_strt.fits")[0].data.copy()
    except:
        bgsub = data-back
        sub = bgsub.copy()
        omap = sub*0.
        for loop in range(3):
            message('CR iteration %d ...'%(loop+1))
            map = crFind(sub,data+readvar)
            if map.sum()==0:
                message(' no new CRs found.\n')
                break
            else:
                message(' %d pixels flagged.\n'%(map.sum()))
            inmask = sub.flatten()
            cond = numpy.where(map.ravel()==1)[0]
            inmask[cond] = numpy.where(cond%2==0,1e5,-1e5).astype(inmask.dtype)
            med5 = ndimage.median_filter(inmask.reshape(sub.shape),7)
            med5 *= map
            sub = (1.-map)*sub + med5
            oCRS = bgsub-sub
            omap += map
            pyfits.PrimaryHDU(oCRS).writeto('blah.fits',clobber=True)
        omap = omap>0

        sub = numpy.where(numpy.isfinite(bgsub),bgsub,0.)
        ccdBG = numpy.where(numpy.isfinite(back),back,0.)
        orig = numpy.where(numpy.isfinite(data),data,0.)
        readvar = 2.7**2
        data -= crmask
        strt = straighten(data,y_soln,orders)
        pyfits.PrimaryHDU(strt).writeto(out_prefix+"_strt.fits")


    hdulist = pyfits.HDUList([pyfits.PrimaryHDU()])
    hdulist[0].header = hdu[0].header.copy()
    varhdu = pyfits.HDUList([pyfits.PrimaryHDU()])
    varhdu[0].header = hdu[0].header.copy()
    for indx in range(len(orders)):
        i,j = orders[indx]
        cut = strt[i:j].copy()
        cut[scipy.isinf(cut)] = scipy.nan
        bg = stats.stats.nanmedian(cut,axis=0)

        B = mask[i:j]
        tmp = (cut-bg)*B

        T = tmp[3:-3,blue[indx]:red[indx]]
        T = T[numpy.isfinite(T)]
        avg,std = clip2(T)

        slice = stats.stats.nanmedian(tmp[:,blue[indx]:red[indx]],axis=1)
        slice[:3] = 0.
        slice[-3:] = 0.
        avg,std = clip2(slice[3:-3],3.5,0.05)
        good = numpy.where(slice<2.*std,1.,0.)
        good = ndimage.maximum_filter(good,5)
        good = ndimage.minimum_filter(good,15)

        good = good==1

        xvals = scipy.arange(good.size).astype(scipy.float32)
        fitdata = scipy.empty((slice[good].size,2))
        fitdata[:,0] = xvals[good].copy()
        bgsub = cut.copy()
        for k in range(cut.shape[1]):
            fitdata[:,1] = cut[:,k][good].copy()
            if fitdata[:,1][scipy.isfinite(fitdata[:,1])].size<4:
                continue
            avg,std = clip(fitdata[:,1])
            keep = scipy.where((abs(fitdata[:,1]-avg)<4.*std)&scipy.isfinite(fitdata[:,1]))[0]
            fit = fitdata[keep]
            if fit[:,1][scipy.isfinite(fit[:,1])].size<4:
                continue
            fit = sf.lsqfit(fit,'chebyshev',1)
            bg = sf.genfunc(xvals,0,fit)
            bgsub[:,k] -= bg

        solution = wave_soln[indx][0]
        invsoln = wave_soln[indx][1]
        vals = scipy.arange(cut.shape[1]).astype(scipy.float64)
        wave = sf.genfunc(vals,0.,solution)

        delta = (wave[wave.size-1]-wave[0])/(wave.size-1)
        outc = scipy.arange(wave[0],wave[0]+2,delta)
        outc = outc[:cut.shape[1]].copy()
        pix = sf.genfunc(outc,0,invsoln)
        coords = spectools.array_coords(cut.shape).astype(numpy.float64)
        coords[1] *= 0.
        coords[1] += pix
        bgsub[scipy.isnan(bgsub)] = 0.
        out = ndimage.map_coordinates(bgsub,coords,output=scipy.float64,order=5)
        hdu = pyfits.ImageHDU(out.copy())
        hdu.header.update('CTYPE1','WAVE-LOG')
        hdu.header.update('CRVAL1',wave[0])
        hdu.header.update('CRPIX1',1)
        hdu.header.update('CDELT1',delta)
        hdu.header.update('CTYPE2','LINEAR')
        hdu.header.update('CRVAL2',1)
        hdu.header.update('CRPIX2',1)
        hdu.header.update('CDELT2',1)
        hdu.header.update('CD1_1',delta)
        hdu.header.update('CD2_2',1)
        hdu.header.update('CD1_2',0)
        hdu.header.update('CD2_1',0)
        #hdu.header.update('PC1_1',1)
        #hdu.header.update('PC1_2',0)
        #hdu.header.update('PC2_1',0)
        #hdu.header.update('PC2_2',1)
        hdu.header.update('DC-FLAG',1)
        hdu.header.update('WFITTYPE','LOG_LINEAR')
        hdu.header.update('DISPAXIS',1)
        hdulist.append(hdu)

        cut[scipy.isnan(cut)] = 0.
        out = ndimage.map_coordinates(cut,coords,output=scipy.float64,order=5)

        hdu = pyfits.ImageHDU(out.copy()+5.3)
        hdu.header.update('CTYPE1','WAVE-LOG')
        hdu.header.update('CRVAL1',wave[0])
        hdu.header.update('CRPIX1',1)
        hdu.header.update('CD1_1',delta)
        hdu.header.update('CRVAL2',1)
        hdu.header.update('CRPIX2',1)
        hdu.header.update('CD2_2',1)
        hdu.header.update('DC-FLAG',1)
        hdu.header.update('WFITTYPE','LOG_LINEAR')
        hdu.header.update('DISPAXIS',1)
        varhdu.append(hdu)

    hdulist.verify('fix')
    hdulist.writeto(out_prefix+"_2d.fits")
    varhdu.verify('fix')
    varhdu.writeto(out_prefix+"_var.fits")
