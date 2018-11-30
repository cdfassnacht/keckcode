from spectra import ycorrect,id_slits,spectools
from esi.biastrim import biastrim

from special_functions import lsqfit,genfunc
from linear import lsqfit as llsqfit#,genfunc

import numpy,scipy,pyfits,pickle
from scipy import ndimage,stats

def starstack(starfiles,out_prefix):
    # Read data from files, biastrim and coadd
    bias = pyfits.open(out_prefix+"_bias.fits")[0].data.astype(scipy.float32)
    bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.astype(scipy.float32)
    stars = []
    for name in starfiles:
        startmp = pyfits.open(name)[0].data.astype(scipy.float32)
        stardata = biastrim(startmp,bias,bpm) 
        stars.append(stardata)
    return numpy.median(stars,0)
#        if first:
#            first = False
#            star = stardata
#        else:
#            star += stardata
#        del stardata,startmp
#    return star

def findpeaks(stars,peak,table,start,end,step,thresh):
    WIDTH = 4
    peak = peak
    standard = peak

    nrows = stars.shape[0]
    xvals = scipy.arange(nrows).astype(scipy.float32)

    solution = {'type':'chebyshev','coeff':scipy.array([[peak]])}
    found = []
    order = 0
    for col in range(start,end,step):
        peak = genfunc(col,0,solution)[0]
        if peak-WIDTH<0 or peak+WIDTH+1>nrows:
            break
        peak = stars[peak-WIDTH:peak+WIDTH+1,col].argmax()+peak-WIDTH
        x = xvals[peak-WIDTH:peak+WIDTH+1]-peak
        flux = stars[peak-WIDTH:peak+WIDTH+1,col]

        crj = ndimage.maximum_filter(flux,5)
        tmppeak = scipy.where(crj==flux)[0]

        """
        if tmppeak.size==0:
            continue
        if abs(tmppeak-flux.size/2.).min()<1.5:
            center = abs(tmppeak-flux.size/2.).argmin()
            if flux.max()!=flux[tmppeak[center]]:
                flux[crj==crj.max()] = scipy.nan
        """

        good = scipy.isfinite(flux)
        if good.size==0:
            continue
        x = x[good]
        flux = flux[good]
        f = flux.sum()
        center = (flux*x).sum()/f
        if abs(center)>1.5 or f<thresh:
            continue
        peak += center
        found.append([col,peak])
        table.append([col,standard,peak])
        if len(found)<8 and len(found)%2**order!=0:
            continue
        if len(found)>8 and abs((col-start)/step)%2**order!=0:
            continue
        if order<3:
            order += 1
        if len(found)>8:
            AA = numpy.array(found)
            solution = lsqfit(AA,'chebyshev',3)
            #bah = llsqfit(AA,'chebyshev',3)
            #import pylab
            #pylab.plot(AA[:,0],genfunc(AA[:,0],0.,solution))
            #pylab.plot(AA[:,0],genfunc(AA[:,0],0.,bah))
            #pylab.plot(AA[:,0],AA[:,1])
            #pylab.show()
        elif len(found)==8:
            solution = lsqfit(scipy.asarray(found),'chebyshev',2)
        elif len(found)==4:
            solution = lsqfit(scipy.asarray(found),'chebyshev',1)
        else:
            solution = {'type':'chebyshev','coeff':scipy.array([[peak]])}
    import pylab
    pylab.plot(AA[:,0],genfunc(AA[:,0],0.,solution))
    pylab.plot(AA[:,0],AA[:,1])
    pylab.show()

    return table



def clip(arr,sig):
    a = arr.copy()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<sig*s]
        if a.size==l:
            return m,s


def startrace(stars,orders):
    from scipy import ndimage,signal
    step = 5
    WIDTH = 6

    nrows = stars.shape[0]
    ncols = stars.shape[1]

    background = ndimage.percentile_filter(stars,10.,(37,1))
    background[numpy.isnan(background)] = 0.
    background = ndimage.median_filter(background,7)
    #import pyfits
    #pyfits.PrimaryHDU(stars-background).writeto('check.fits',clobber=True)
    #df
    smooth = ndimage.gaussian_filter1d(stars,1.5,0)

    stars -= background

    start = stars[:,ncols/2].copy()
    #background = ndimage.percentile_filter(start,10.,37)
    #start -= background
    peaks = ndimage.maximum_filter(start,9)
    peaks2 = ndimage.maximum_filter(start,3)
    
    avg,sig = clip(start,3.)
    thresh = sig*5.
    all = scipy.where((peaks==start)&(peaks>thresh)&(peaks==peaks2))[0]
    all = all.tolist()

    xvals = numpy.arange(nrows)*1.
    yvals = numpy.arange(ncols)*1.
    solutions = []
    for i in range(len(orders)):
        start,end = orders[i]
        while all[0]<start:
            del all[0]
        ordermatches = numpy.empty((0,3))
        while all[0]<end:
            peak = all[0]
            col = ncols/2
            amp = smooth[peak,col]
            line = peak
            matches = []
            while col>0:
                line = int(line)
                if line<WIDTH or line+WIDTH>=nrows or stars[line-1:line+2,col].max()<amp/20.:
                    break
                s1 = slice(line-WIDTH,line+WIDTH+1)
                line = (stars[s1,col]*xvals[s1]).sum()/stars[s1,col].sum()
                matches.append([col,line,peak])
                col -= 1
            col = ncols/2+1
            line = peak
            while col<ncols:
                line = int(line)
                if line<WIDTH or line+WIDTH>=nrows or stars[line-1:line+2,col].max()<amp/20.:
                    break
                s1 = slice(line-WIDTH,line+WIDTH+1)
                if numpy.isnan(stars[s1,col]).any():
                    fitData = numpy.array(matches)[:,:2]
                    fit = lsqfit(fitData,'chebyshev',3)
                    linePos = genfunc(yvals,0.,fit).astype(numpy.int32)
                    while col<nrows-WIDTH and numpy.isnan(stars[s1,col]).any():
                        col += 1
                        s1 = slice(linePos[col]-WIDTH,linePos[col]+WIDTH+1)
                    if col>=nrows-WIDTH:
                        break
                line = (stars[s1,col]*xvals[s1]).sum()/stars[s1,col].sum()
                matches.append([col,line,peak])
                col += 1
            matches = numpy.array(matches)
            soln = lsqfit(matches[:,:2],'cheybshev',3)
            resid = matches[:,1]-genfunc(matches[:,0],0.,soln)
            good = abs(resid)<2.
            ordermatches = numpy.concatenate((ordermatches,matches[good]))
            del all[0]
        matches = ordermatches
        ord1,ord2 = 3,3
        ymap = lsqfit(matches,"chebyshev",ord1,ord2)
        resid = matches[:,2]-genfunc(matches[:,0],matches[:,1],ymap)
        m,s = clip(resid,3.5)
        good = abs(resid-m)<2.5
        matches = matches[good]
        ymap = lsqfit(matches,"chebyshev",ord1,ord2)
        tmp = matches[:,1].copy()
        matches[:,1] = matches[:,2].copy()
        matches[:,2] = tmp.copy()
        ytrue = lsqfit(matches,"chebyshev",ord1,ord2)
        solutions.append([ytrue,ymap])

    return solutions



    peaksByOrder = []
    for i in range(len(orders)):
        start,end = orders[i]
        peaksByOrder.append([])
        for p in all:
            if p>start-5 and p<end+5:
                peaksByOrder[i].append(p)

    solutions = []

    """ Mask cosmic rays """
    stars = stars.copy()
    mask = ndimage.maximum_filter(stars,5)
    stars[mask>63000.] = scipy.nan

    xvals = numpy.arange(nrows)
    for peaks in peaksByOrder:
        matches = []
        
    for p in peaks:
        table = []
        xvals = scipy.arange(nrows).astype(scipy.float32)
        for peak in p:
            if peak-WIDTH<0 or peak+WIDTH+1>nrows:
                continue
            x = xvals[peak-WIDTH:peak+WIDTH+1]-peak
            flux = stars[peak-WIDTH:peak+WIDTH+1,ncols/2]
            good = scipy.isfinite(flux)
            x = x[good]
            flux = flux[good]
            peak += (flux*x).sum()/flux.sum()

            table = findpeaks(stars,peak,table,ncols/2,ncols,step,thresh)
            table = findpeaks(stars,peak,table,ncols/2-step,0,-1*step,thresh)
        transform_table = scipy.asarray(table)


        ord1,ord2 = 3,3
        ytrue = lsqfit(transform_table,"chebyshev",ord1,ord2)
        tmp = transform_table[:,1].copy()
        transform_table[:,1] = transform_table[:,2].copy()
        transform_table[:,2] = tmp.copy()
        del tmp
#        transform_table[:,1],transform_table[:,2] = transform_table[:,2],transform_table[:,1]

        ymap = lsqfit(transform_table,"chebyshev",ord1,ord2)
        solutions.append([ytrue,ymap])

    return solutions

def straighten2(data,solutions,orders):
    output = data*scipy.nan

    data[scipy.isnan(data)] = 0.
    data[scipy.isinf(data)] = 0.
    coords = spectools.array_coords(data.shape)
    x = coords[1].flatten()
    y = coords[0].flatten()
    for i in range(len(solutions)):
        low,high = orders[i]
        ytrue,ymap = solutions[i]
        yforw = genfunc(x,y,ytrue).reshape(coords[0].shape)
        coords[0] = yforw.copy()

        tmpcoords = coords[:,low:high,:].copy()
        tmp = ndimage.map_coordinates(data,tmpcoords,output=scipy.float64,cval=-2**15,order=5)
        output[low:high] = tmp.copy()
    output[output==-2**15] = scipy.nan
    return output


def straighten(data,solutions,orders):
    d = data.copy()
    output = d*numpy.nan
    bad = ~numpy.isfinite(d)
    d[bad] = 0.
    coords = spectools.array_coords(d.shape)
    x = coords[1].flatten()
    y = coords[0].flatten()
    for i in range(len(solutions)):
        low,high = orders[i]
        ytrue,ymap = solutions[i]
        yforw = genfunc(x,y,ytrue).reshape(coords[0].shape)
        cmin = int(yforw[low:high].min())-5
        cmax = int(yforw[low:high].max())+6
        if cmin<0:
            cmin = 0
        if cmax>d.shape[0]:
            cmax = d.shape[0]

        coords[0] = yforw-cmin

        tmpcoords = coords[:,low:high,:]
        tmp = ndimage.map_coordinates(d[cmin:cmax],tmpcoords,output=scipy.float64,cval=-2**15,order=5)
        output[low:high] = tmp.copy()
    output[output==-2**15] = scipy.nan
    return output



def curve(indata,solutions,orders):
    d = indata.astype(scipy.float64)

    mid = scipy.median(d[scipy.isfinite(indata)])
    d[scipy.isnan(indata)] = mid
    output = scipy.zeros(d.shape) 
    coords = spectools.array_coords(d.shape)
    x = coords[1].flatten()
    y = coords[0].flatten()
    yorig = coords[0].copy()
    for i in range(len(solutions)):
        low,high = orders[i]
        ytrue,ymap = solutions[i]
        yforw = genfunc(x,y,ytrue).reshape(coords[0].shape)
        yback = genfunc(x,y,ymap).reshape(coords[0].shape)
        coords[0] = yback.copy()

        cond = (yorig<high+1)&(yorig>low-1)
        ymin = yforw[cond].min()
        ymax = yforw[cond].max()
        if ymin<0:
            ymin = 0
        if ymax>output.shape[0]:
            ymax = output.shape[0]
        tmpcoords = coords[:,ymin:ymax,:].copy()
        tmp = ndimage.map_coordinates(d,tmpcoords,output=scipy.float64,cval=0.,order=5)
        tmp[tmpcoords[0]>high+2] = 0.
        tmp[tmpcoords[0]<low-2] = 0.
        output[ymin:ymax] += tmp.copy()
    return output
