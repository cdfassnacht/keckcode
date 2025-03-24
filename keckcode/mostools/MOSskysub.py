"""
Background subtract, resample, and coadd 2d spectra.
"""
import time
from mostools import spectools,skysub#correct_telluric,skysub
from special_functions import genfunc,lsqfit

import scipy,numpy
from scipy import optimize,interpolate,ndimage,signal,stats


RESAMPLE = 0        # 1 to only resample
magicnum = -2**15



def quickNanMed(a,ax=0):
    m = numpy.nanmin(a)-1.
    M = numpy.nanmax(a)+1.
    arr = numpy.rollaxis(a,ax,0).repeat(2,0).astype(numpy.float64)
    s = arr[0].size
    n = arr.shape[0]
    x = arr.reshape((n,s)).T.flatten()
    N = numpy.where(numpy.isnan(x))[0]
    new = numpy.where(numpy.arange(N.size)%2==1,m,M)
    x[N] = new
    o = numpy.median(x.reshape((s,n)).T.reshape(arr.shape),0)
    o[o==(m+M)/2] = numpy.nan
    return o


def cr_reject(sub,sky,nsig=4.,contrast=1.,sigfrac=0.3):
    import indexTricks as iT
    sub = sub.copy()
    med5 = ndimage.median_filter(sub,5)
    med5 += sky

    n = 2
    blksub = sub.repeat(n,0).repeat(n,1)
    blksub = ndimage.laplace(blksub)*-1.
    blksub[blksub<0] = 0.
    sub2 = iT.resamp(blksub,n)

    med5[med5<=0.] = 0.01

    noise = med5**0.5
    sigmap = sub2/noise
    med5 = ndimage.median_filter(sigmap,5)
    sigmap -= med5

    map = sigmap.copy()
    map[map<nsig] = 0
    map[map>0] = 1

    med3 = ndimage.median_filter(sub,3)
    med7 = ndimage.median_filter(med3,7)

    med3 -= med7
    med3 /= noise

    med3[med3<0.01] = 0.01

    stars = (map*sigmap)/med3

    stars[stars<contrast] = 0
    stars[stars>0.] = 1

    map *= stars

    gmap = ndimage.maximum_filter(map,5)*sigmap
    gmap[gmap<nsig] = 0
    gmap[gmap>0] = 1
    map = gmap.copy()

    fmap = ndimage.maximum_filter(gmap,5)*sigmap
    fmap[fmap<(nsig*sigfrac)] = 0
    fmap[fmap>0] = 1
    map = fmap.copy()
    return map

def crFind(img,var,nsig=10.,sigfrac=0.3):
    import indexTricks as iT
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


"""
doskysub()
"""
def doskysub(straight,ylen,xlen,sci,yback,sky2x,sky2y,ccd2wave,disp,mswave,offsets,cutoff,readvar,crflag):
    sci = sci.copy()

    # If cutoff is not a float, we are using the blueside
    if type(cutoff)==type([]):
        locutoff,hicutoff = cutoff
    else:
        locutoff = cutoff
        hicutoff = 10400.#7

    nsci = sci.shape[0]
    width = sci.shape[2]

    # Perform telluric correction
    coords = numpy.indices(sci[0].shape)
    x = coords[1].flatten()
    y = coords[0].flatten()


    # Create arrays for output images
    outcoords = numpy.indices((ylen,xlen)).astype(scipy.float64)
    outcoords[1] *= disp
    outcoords[1] += mswave - disp*xlen/2.
    xout = outcoords[1].flatten()
    yout = outcoords[0].flatten()

    out = scipy.zeros((nsci,ylen,xlen))

    if RESAMPLE==1:
        varimage = out.copy()
        bgimage = out.copy()
    else:
        fudge = scipy.ceil(abs(offsets).max())
        bgimage = scipy.zeros((nsci,ylen+fudge,xlen))
        varimage = bgimage.copy()

        bgcoords = spectools.array_coords((ylen+fudge,xlen))
        bgcoords[1] *= disp
        bgcoords[1] += mswave - disp*xlen/2.

    #
    # Cosmic Ray Rejection and Background Subtraction
    #
    yfit = yback.flatten()
    ycond = (yfit>straight-0.4)&(yfit<straight+ylen-0.6)

    coords = numpy.indices(yback.shape).astype(numpy.float32)
    xvals = coords[1].flatten()
    yvals = coords[0].flatten()

    ap_y = scipy.zeros(0)
    aper = scipy.zeros(0)

    n = 5
    re_xout = outcoords[1].flatten()
    re_yout = outcoords[0].flatten()
    re_coords = outcoords.copy()


    for k in range(nsci):
        xfit = genfunc(xvals,yfit-straight,ccd2wave[k])
        zfit = sci[k].flatten()

        x = xfit[ycond]
        y = yfit[ycond]
        z = zfit[ycond]

        # The plus/minus 20 provides a better solution for the edges
        wavecond = (x>locutoff-20.)&(x<hicutoff+20.)
        x = x[wavecond]
        y = y[wavecond]
        z = z[wavecond]

        samp_x = genfunc(re_xout,re_yout,sky2x[k])
        samp_y = genfunc(re_xout,re_yout,sky2y[k])
        re_coords[0] = samp_y.reshape(re_coords[0].shape)
        re_coords[1] = samp_x.reshape(re_coords[1].shape)
        if RESAMPLE==1:
            tmpsci = sci[k].copy()
            tmpOut = ndimage.map_coordinates(tmpsci,re_coords,output=scipy.float64,order=1,cval=-32768)

            background = numpy.sort(tmpOut,0)[5]
            tmpShape = tmpOut.shape
            for i in range(3):
                bg = numpy.tile(background,tmpShape[0]).reshape(tmpShape)
                tmp = numpy.where(tmpOut-bg>4.*(bg+readvar)**0.5,numpy.nan,tmpOut)
                background = quickNanMed(tmp)
                cond = numpy.isnan(background)
                background[cond] = bg[0][cond]
            bgout = numpy.tile(background,tmpShape[0]).reshape(tmpShape)

            tmpCoords = numpy.indices(sci[k].shape).astype(numpy.float32)
            tmpCoords[1] = (genfunc(tmpCoords[1].ravel(),yfit-straight,ccd2wave[k]).reshape(tmpCoords[1].shape)-(mswave - disp*xlen/2.))/disp
            mod = interpolate.splrep(numpy.arange(background.size),background,s=0,k=3)
            ccdBG = interpolate.splev(tmpCoords[1].ravel(),mod).reshape(tmpCoords[1].shape)
        else:
            print "Determining `optimal' sky for image %d"%(k+1)
            bgfit = skysub.skysub(x,y,z,disp)
            print "Subtracting sky"
            background = zfit.copy()
            a = time.time()
            for indx in range(background.size):
                x0 = xfit[indx]
                y0 = yfit[indx]
                if x0<locutoff-10 or x0>hicutoff+10:
                    background[indx] = scipy.nan
                else:
                    background[indx] = interpolate.bisplev(x0,y0,bgfit)
#            bgout = re_coords[0,0]*0.
#            for i in range(bgout.size):
#                bgout[i] = interpolate.bisplev(re_xout[],re_coords[0,0,i],bgfit)
#            bgout = numpy.tile(bgout,re_coords[0].shape[0]).reshape(re_coords[0].shape)
            bgout = re_coords[0]*0.
            for i in range(bgout.shape[0]):
                for j in range(bgout.shape[1]):
 #                   bgout[i,j] = interpolate.bisplev(re_coords[1,i,j],re_coords[0,i,j],bgfit)
                    bgout[i,j] = interpolate.bisplev(re_xout[j],re_coords[0,i,j]+straight,bgfit)
            sub = zfit-background
            sub[scipy.isnan(sub)] = 0.
            sky = sub*0.
            sky[ycond] = sub[ycond]
            sky = sky.reshape(sci[k].shape)
            #sub = sky.copy()
            background[numpy.isnan(background)] = 0.
            ccdBG = sub*0.
            ccdBG[ycond] = background[ycond]
            ccdBG = ccdBG.reshape(sci[k].shape)

        sub = sci[k]-ccdBG
        orig = sub.copy()

        if crflag==True:
            print "Finding and removing cosmic rays"
            for loop in range(5):
                map = crFind(sub,ccdBG+readvar)
                if map.sum()==0:
                    break
                inmask = sub.flatten()
                cond = numpy.where(map.ravel()==1)[0]
                inmask[cond] = numpy.where(cond%2==0,1e5,-1e5).astype(inmask.dtype)
                med5 = ndimage.median_filter(inmask.reshape(sub.shape),7)
                med5 *= map
                sub = (1.-map)*sub + med5
            oCRS = orig-sub
        else:
            oCRS = orig*0.

        re_coords = re_coords*n + n/2
        crs = ndimage.map_coordinates(oCRS.repeat(n,0).repeat(n,1),re_coords,output=scipy.float64,order=1,cval=-32768)

        re_coords[0] = samp_y.reshape(re_coords[0].shape)
        re_coords[1] = samp_x.reshape(re_coords[1].shape)
        tmpOut = ndimage.map_coordinates(sci[k]-oCRS,re_coords,output=scipy.float64,order=5,cval=-32768)
#            tmpBGsub = ndimage.map_coordinates(sci[k]-oCRS-ccdBG,re_coords,output=scipy.float64,order=5,cval=-32768)#-bgout
        tmpBGsub = tmpOut-bgout
        out[k] = tmpOut.copy()

        out[k][xout.reshape(outcoords[1].shape)<locutoff] = scipy.nan
        out[k][xout.reshape(outcoords[1].shape)>hicutoff] = scipy.nan
        out[k][out[k]==-32768] = scipy.nan

        vartmp = out[k].copy()
        vartmp[crs>0.] = numpy.nan
        vartmp[vartmp<0.] = numpy.nan
        varimage[k] = vartmp.copy()

        if RESAMPLE==1:
            background = tmpBGsub[0]*0.
            for i in range(3):
                bg = numpy.tile(background,tmpShape[0]).reshape(tmpShape)
                tmp = numpy.where(tmpBGsub-bg>2.5*(bgout+bg+readvar)**0.5,numpy.nan,tmpBGsub)
                background = quickNanMed(tmp)
                cond = numpy.isnan(background)
                background[cond] = bg[0][cond]
            bgout = numpy.tile(background,tmpShape[0]).reshape(tmpShape)
            bgimage[k] = tmpBGsub.copy()-bgout
        else:
            bgimage[k] = ndimage.map_coordinates(sci[k]-oCRS-ccdBG,re_coords,output=scipy.float64,order=5,cval=-32768)

        bgimage[k][numpy.isnan(out[k])] = numpy.nan


    if RESAMPLE==1:
        return out,bgimage,varimage,outcoords[1,0]
    if bgimage.shape[0]>1:    
        bgimage = quickNanMed(bgimage)
        varimage = quickNanMed(varimage)/nsci
    elif bgimage.ndim==3:
        bgimage = bgimage[0].copy()
        varimage = varimage[0].copy()


    return out,bgimage,varimage,outcoords[1,0]
