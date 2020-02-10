from keckcode.spectra import ycorrect,id_slits,spectools
from .biastrim import biastrim

from special_functions import lsqfit,genfunc

import numpy,scipy,pickle
from scipy import ndimage,stats
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

def starstack(starfiles,out_prefix):
    # Read data from files, biastrim and coadd
    bias = pyfits.open(out_prefix+"_bias.fits")[0].data.copy()
    bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.copy()
    stars = []
    for name in starfiles:
        startmp = pyfits.open(name)[0].data.copy()
        stardata = biastrim(startmp,bias,bpm) 
        stars.append(stardata)
    return numpy.median(stars,0)


def clip(arr,sig):
    a = arr.copy()
    while 1:
        m,s,l = a.mean(),a.std(),a.size
        a = a[abs(a-m)<sig*s]
        if a.size==l:
            return m,s


def startrace(stars,orders):
    from scipy import ndimage
    if stars.shape[0]>1500:
        WIDTH = 6
        fw = 9
    else:
        WIDTH = 4
        fw = 7

    nrows = stars.shape[0]
    ncols = stars.shape[1]

    # This probably isn't necessary....
    background = ndimage.percentile_filter(stars,10.,(37,1))
    background[numpy.isnan(background)] = 0.
    background = ndimage.median_filter(background,7)

    smooth = ndimage.gaussian_filter1d(stars-background,1.5,0)

    stars -= background
    start = stars[:,int(ncols/2)].copy()

    # Find the traces down the middle of the image
    peaks = ndimage.maximum_filter(start,fw)
    peaks2 = ndimage.maximum_filter(start,3)
    
    avg,sig = clip(start,3.)
    thresh = sig*5.
    all = numpy.where((peaks==start)&(peaks>thresh)&(peaks==peaks2))[0]
    all = all.tolist()

    yvals = numpy.arange(nrows)*1.
    xvals = numpy.arange(ncols)*1.
    wideorders = []
    solutions = []
    # Loop through the orders, following the traces in each order
    for i in range(len(orders)):
        start,end = orders[i]
        while all[0]<start:
            del all[0]
        ordermatches = numpy.empty((0,3))
        while all[0]<end:
            peak = all[0]
            col = int(ncols/2)
            amp = smooth[peak,col]
            line = peak
            matches = []
            # Go towards the blue side of the image
            while col>0:
                line = int(line)
                if line<WIDTH or line+WIDTH>=nrows or stars[line-1:line+2,col].max()<amp/20.:
                    break
                s1 = slice(line-WIDTH,line+WIDTH+1)
                line = (stars[s1,col]*yvals[s1]).sum()/stars[s1,col].sum()
                matches.append([col,line,peak])
                col -= 1
            col = int(ncols/2)+1
            line = peak
            # Go towards the red side of the image
            while col<ncols:
                line = int(line)
                if line<WIDTH or line+WIDTH>=nrows:
                    break
                s1 = slice(line-WIDTH,line+WIDTH+1)
                # The red side has some defects....
                if numpy.isnan(stars[s1,col]).any():
                    if i==0:
                        break
                    fitData = numpy.array(matches)[:,:2]
                    fit = lsqfit(fitData,'chebyshev',3)
                    linePos = genfunc(xvals,0.,fit).astype(numpy.int32)
                    while col<ncols-WIDTH and numpy.isnan(stars[s1,col]).any():
                        col += 1
                        s1 = slice(linePos[col]-WIDTH,linePos[col]+WIDTH+1)
                    if col>=ncols-WIDTH:
                        break
                line = int((stars[s1,col]*yvals[s1]).sum()/stars[s1,col].sum())
                if stars[line-1:line+2,col].max()<amp/20.:
                    break
                matches.append([col,line,peak])
                col += 1

            # Filter garbage; trace should be within ~2pix of the model
            matches = numpy.array(matches)
            soln = lsqfit(matches[:,:2],'cheybshev',3)
            resid = matches[:,1]-genfunc(matches[:,0],0.,soln)
            good = abs(resid)<2.
            ordermatches = numpy.concatenate((ordermatches,matches[good]))
            del all[0]
        matches = ordermatches

        ord1,ord2 = 3,3
        ymap = lsqfit(matches,"chebyshev",ord1,ord2)

        # Clip bad points
        resid = matches[:,2]-genfunc(matches[:,0],matches[:,1],ymap)
        m,s = clip(resid,3.5)
        good = abs(resid-m)<3.5*s
        matches = matches[good]

        ymap = lsqfit(matches,"chebyshev",ord1,ord2)
        tmp = matches[:,1].copy()
        matches[:,1] = matches[:,2].copy()
        matches[:,2] = tmp.copy()
        ytrue = lsqfit(matches,"chebyshev",ord1,ord2)
        solutions.append([ytrue,ymap])

        lower = int(genfunc(xvals,xvals*0.+start-3,ytrue).min())
        upper = int(genfunc(xvals,xvals*0+end+3,ytrue).max())+1
        if lower<0:
            lower = 0
        if upper>stars.shape[1]:
            upper = stars.shape[1]
        wideorders.append([lower,upper])

    return solutions,wideorders


def straighten(data,solutions,orders,wideorders,interp=5):
    d = data.copy()
    shape = d.shape
    output = d*numpy.nan
    bad = ~numpy.isfinite(d)
    d[bad] = 0.
    coords = spectools.array_coords(shape)
    x = coords[1]#.flatten()
    y = coords[0]#.flatten()
    for i in range(len(solutions)):
        low,high = orders[i]
        wlow,whigh = wideorders[i]
        ytrue,ymap = solutions[i]

        yforw = genfunc(x[low:high].ravel(),y[low:high].ravel(),ytrue)
        yforw = yforw.reshape((high-low,shape[1]))

        tmpcoords = coords[:,low:high,:].copy()
        tmpcoords[0] = yforw-wlow
        tmp = ndimage.map_coordinates(d[wlow:whigh],tmpcoords,output=numpy.float64,cval=-2**15,order=interp)
        output[low:high] = tmp.copy()
    output[output==-2**15] = scipy.nan
    return output


def curve(indata,solutions,orders,wideorders,interp=5):
    d = indata.astype(numpy.float64)

    mid = scipy.median(d[scipy.isfinite(indata)])
    d[scipy.isnan(indata)] = mid
    output = scipy.zeros(d.shape) 
    coords = spectools.array_coords(d.shape)
    x = coords[1]#.flatten()
    y = coords[0]#.flatten()
    yorig = coords[0].copy()
    for i in range(len(solutions)):
        low,high = orders[i]
        wlow,whigh = wideorders[i]
        ytrue,ymap = solutions[i]

        yback = genfunc(x[wlow:whigh].ravel(),y[wlow:whigh].ravel(),ymap)
        yback = yback.reshape((whigh-wlow,coords[0].shape[1]))

        tmpcoords = coords[:,wlow:whigh,:].copy()
        tmpcoords[0] = yback-low
        tmp = ndimage.map_coordinates(d[low:high],tmpcoords,output=numpy.float64,cval=0.,order=interp)
        output[wlow:whigh] += tmp.copy()
    return output


def fullSolution(shape,ysoln,orders,wideorders,wsoln):
    import indexTricks as iT
    coords = iT.coords(shape)
    y = coords[0].copy()
    x = coords[1].copy()
    soln = []
    if shape[1]>3000:
        disp = 1.65e-5
    else:
        disp = 2*1.65e-5
    owave = numpy.arange(3.585,4.038,disp)
    for i in range(len(ysoln)):
        low,high = orders[i]
        wlow,whigh = wideorders[i]
        ytrue,ymap = ysoln[i]
        win,wout = wsoln[i]
        tmpw = genfunc(x[0],0.,win)
        ow = owave[(owave>=tmpw[0])&(owave<=tmpw[-1])].copy()
        if ow.size>x.shape[1]:
            diff = ow.size-x.shape[1]
            if diff%2==0:
                ow = ow[int(diff/2):int(diff/-2)]
            else:
                ow = ow[int(diff/2):int((diff+1)/-2)]
        xc = genfunc(ow,0.,wout)
        xc = xc.repeat(high-low).reshape((xc.size,high-low)).T
        yc = genfunc(xc.ravel(),y[low:high,:xc.shape[1]].ravel(),ytrue)
        yc = yc.reshape((high-low,xc.shape[1]))-wlow

        corr = numpy.linspace(-0.2,0.2,high-low)
        xc = (xc.T+corr).T
        #xc = xc.repeat(high-low).reshape(yc.shape[::-1]).T
        soln.append([numpy.array([yc,xc]),ow[0],ow[-1],disp])
    return soln
        

def getOrders(indata,orders,wideorders,fullsoln,interp=5):
    d = indata.astype(numpy.float64)
    d[numpy.isnan(d)] = 0.
    shape = d.shape

    outspex = []
    size = -5
    for i in range(len(orders)):
        low,high = orders[i]
        wlow,whigh = wideorders[i]
        soln,wlo,whi,disp = fullsoln[i]

        slit = ndimage.map_coordinates(d[wlow:whigh],soln,output=numpy.float64,order=interp)
        outspex.append([slit,wlo,whi,disp])
        size += high-low+5
    return outspex
