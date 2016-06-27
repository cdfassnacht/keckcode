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

import numpy,scipy
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


def prepare(dir,prefix,bias,stars,hgne,cuar,xe,flat,out_prefix,redoWave=False):
    biasnums = bias.split(",")
    starnums = stars.split(",")
    flatnums = flat.split(",")

    for i in range(len(biasnums)):
        biasnums[i] = dir+prefix+biasnums[i]+".fits"
    for i in range(len(starnums)):
        starnums[i] = dir+prefix+starnums[i]+".fits"
    for i in range(len(flatnums)):
        flatnums[i] = dir+prefix+flatnums[i]+".fits"

    hgne = dir+prefix+hgne+".fits"
    cuar = dir+prefix+cuar+".fits"
    xe = dir+prefix+xe+".fits"

    message("Bias and bad pixel mask ... ")
    try:
        bias = pyfits.open(out_prefix+"_bias.fits")[0].data.copy()
        bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.copy()
        message("opened.\n")
    except:
        bias,bpm = make_bias(biasnums)
        pyfits.PrimaryHDU(bias).writeto(out_prefix+"_bias.fits")
        pyfits.PrimaryHDU(bpm).writeto(out_prefix+"_bpm.fits")
        message("created.\n")

    message("Flat field ... ")
    try:
        flat = pyfits.open(out_prefix+"_flat.fits")[0].data.copy()
        message("opened.\n")
    except:
        flat = make_flat(flatnums,out_prefix)
        pyfits.PrimaryHDU(flat).writeto(out_prefix+"_flat.fits")
        message("created.\n")

    message("Star traces ... ")
    try:
        star = pyfits.open(out_prefix+"_trace.fits")[0].data.astype(scipy.float32)
        message("opened.\n")
    except:
        from esi.straighten import starstack
        star = starstack(starnums,out_prefix)
        pyfits.PrimaryHDU(star).writeto(out_prefix+"_trace.fits")
        message("created.\n")

    message("Order solutions ... ")
    try:
        orders,solutions,wideorders = numpy.load(out_prefix+"_ycor.dat")
        message("opened.\n")
    except:
        orders = find_orders(flat)
        solutions,wideorders = startrace(star,orders)
        file = open(out_prefix+"_ycor.dat","wb")
        dump([orders,solutions,wideorders],file,2)
        file.close()
        message("created.\n")

    message("Normalized flat ... ")
    try:
        respflat = pyfits.open(out_prefix+"_norm.fits")
        message("opened.\n")
    except:
        cflat = response(flat,orders,solutions,wideorders)
        pyfits.PrimaryHDU(cflat).writeto(out_prefix+"_norm.fits")
        message("created.\n")

    try:
        hgne = pyfits.open(out_prefix+"_hgne.fits")[0].data.astype(scipy.float32)
    except:
        print "Straightening HgNe arc"
        hdu = pyfits.open(hgne)
        hgne = hdu[0].data.copy()
        hgne = biastrim(hgne,bias,bpm)
        hgne = straighten(hgne,solutions,orders,wideorders)
        hdu[0].data = hgne
        hdu[0].writeto(out_prefix+"_hgne.fits")

    try:
        cuar = pyfits.open(out_prefix+"_cuar.fits")[0].data.astype(scipy.float32)
    except:
        print "Straightening CuAr arc"
        hdu = pyfits.open(cuar)
        cuar = hdu[0].data.copy()
        cuar = biastrim(cuar,bias,bpm)
        cuar = straighten(cuar,solutions,orders,wideorders)
        hdu[0].data = cuar
        hdu[0].writeto(out_prefix+"_cuar.fits")
    try:
        xe = pyfits.open(out_prefix+"_xe.fits")[0].data.astype(scipy.float32)
    except:
        print "Straightening Xe arc"
        hdu = pyfits.open(xe)
        xe = hdu[0].data.copy()
        xe = biastrim(xe,bias,bpm)
        xe = straighten(xe,solutions,orders,wideorders)
        hdu[0].data = xe 
        hdu[0].writeto(out_prefix+"_xe.fits")

    arcs = {}
    arcs['cuar'] = cuar.copy()
    arcs['hgne'] = hgne.copy()
    arcs['xe'] = xe.copy()

    if redoWave==False:
        try:
            wave_solution = numpy.load(out_prefix+"_wave.dat")
        except:
            redoWave = True
    if redoWave==True:
        print "Determining wavelength solutions"
        wave_solution = solve(arcs,orders)

        file = open(out_prefix+"_wave.dat","wb")
        dump(wave_solution,file,2)
        file.close()


    print "Creating full slit solutions"
    soln = fullSolution(bias.shape,solutions,orders,wideorders,wave_solution)
    file = open(out_prefix+"_full.dat",'wb')
    dump(soln,file,2)
    file.close()

    print ""
    print "**********************"
    print "Preparations complete!"
    print "**********************"
    print ""


def bgsub(dir,inname,diffim,out_prefix,cal_prefix):
    bias = pyfits.open(cal_prefix+"_bias.fits")[0].data.astype(scipy.float32)
    bpm = pyfits.open(cal_prefix+"_bpm.fits")[0].data.astype(scipy.float32)
    flat = pyfits.open(cal_prefix+"_norm.fits")[0].data.astype(scipy.float32)

    file = open(cal_prefix+"_ycor.dat")
    orders = load(file)
    y_soln = load(file)
    file.close()

    hdu = pyfits.open(dir+inname)
    data = hdu[0].data.copy()
    data = biastrim(data,bias,bpm)/flat

    diff = biastrim(pyfits.open(dir+diffim)[0].data.copy(),bias,bpm)/flat
    del bias,bpm,flat

    back = scipy.empty(data.shape)*scipy.nan
    xvals = scipy.arange(back.shape[1])
    for indx in range(len(orders)):
        i,j = orders[indx]
        bottom = i+3
        top = j-3
        top = sf.genfunc(xvals,top,y_soln[indx][0])
        bottom = sf.genfunc(xvals,bottom,y_soln[indx][0])
        for k in range(xvals.size):
            b = floor(bottom[k])
            t = ceil(top[k])
            if t<0:
                back[b:t,k] = scipy.nan
                continue

            if b<0:
                b = 0
            if t>data.shape[0]:
                t = data.shape[0]
            y = scipy.arange(b,t,1.)
            d = data[b:t,k].copy()-diff[b:t,k].copy()
            a,s = clip(d)
            good = abs(d-a)<2.5*s

            d = data[b:t,k].copy()
            b -= 3
            t += 3
            if b<0:
                b = 0
            if t>data.shape[0]:
                t = data.shape[0]
            yout = scipy.arange(b,t,1.)

            if d[good].size<4:
                back[b:t,k] = scipy.nan
                continue
            fitd = scipy.empty((y[good].size,2))
            fitd[:,0] = y[good].copy()
            fitd[:,1] = d[good].copy()
            fit = sf.lsqfit(fitd,'chebyshev',3)

            m = d - sf.genfunc(y,0.,fit)
            a,s = clip(d)
            good = abs(d-a)<2.5*s

            if d[good].size<4:
                back[b:t,k] = sf.genfunc(yout,0.,fit)
                continue
            fitd = scipy.empty((y[good].size,2))
            fitd[:,0] = y[good].copy()
            fitd[:,1] = d[good].copy()
            fit = sf.lsqfit(fitd,'chebyshev',3)

            back[b:t,k] = sf.genfunc(yout,0.,fit)
    back[scipy.isinf(back)] = scipy.nan
    pyfits.PrimaryHDU(back).writeto(out_prefix+"_bg.fits")
    pyfits.PrimaryHDU(data-back).writeto(out_prefix+"_bgsub.fits")

    
def pipeline(dir,inname,out_prefix,cal_prefix):
    bias = pyfits.open(cal_prefix+"_bias.fits")[0].data.astype(scipy.float32)
    bpm = pyfits.open(cal_prefix+"_bpm.fits")[0].data.astype(scipy.float32)
    flat = pyfits.open(cal_prefix+"_norm.fits")[0].data.astype(scipy.float32)

    file = open(cal_prefix+"_ycor.dat")
    orders = load(file)
    y_soln = load(file)
    file.close()

    file = open(cal_prefix+"_wave.dat")
    wave_soln = load(file)
    file.close()

    try:
        mask = pyfits.open(cal_prefix+"_mask.fits")[0].data.copy()
        print 'mask opened'
    except:
        mask = numpy.where(bpm==1.,1.,0.)
        mask = straighten(mask,y_soln,orders)
        mask[mask<0.7] = numpy.nan
        mask[mask>=0.7] = 1.
        pyfits.PrimaryHDU(mask).writeto(cal_prefix+"_mask.fits")

    hdu = pyfits.open(dir+inname)
    data = hdu[0].data.copy()
    data = biastrim(data,bias,bpm)#*flat

    blue = [1500,1400,1300,1200,1100,900,600,200,0,0,0]
    red = [3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]

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

    try:
        back = pyfits.open(out_prefix+"_bg.fits")[0].data.astype(scipy.float32)
    except:
        strt = straighten(data,y_soln,orders)
        back = scipy.zeros(strt.shape)
        for indx in range(len(orders)):
            i,j = orders[indx]
            bg = stats.stats.nanmedian(strt[i:j],axis=0)
            back[i:j] += bg
        back = curve(back,y_soln,orders)
        pyfits.PrimaryHDU(back).writeto(out_prefix+"_bg.fits")

    try:
        strt = pyfits.open(out_prefix+"_strt.fits")[0].data.astype(scipy.float32)
    except:
        bgsub = data-back
        bgsub[scipy.isnan(bgsub)] = 0.
        bgsub[scipy.isinf(bgsub)] = 0.
        model = ndimage.median_filter(bgsub,3)
        diff = bgsub-model
        model = scipy.sqrt(data)
        model[scipy.isnan(model)] = 0.
        model[scipy.isinf(model)] = 0.
        diff[scipy.isnan(diff)] = 0.
        diff[scipy.isinf(diff)] = 0.
        crmask = scipy.where(diff>4.*model,diff,0.)

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
