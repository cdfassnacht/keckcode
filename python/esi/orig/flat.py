from spectra import ycorrect,id_slits,spectools
from esi.biastrim import biastrim
from esi.straighten import straighten,curve

from special_functions import lsqfit,genfunc

import scipy,pyfits,pickle
from scipy import stats


def make_flat(flatfiles,out_prefix):
    # Read data from files, biastrim and coadd
    weight = 0
    bias = pyfits.open(out_prefix+"_bias.fits")[0].data.astype(numpy.float32)
    bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.astype(numpy.float32)
    for name in flatfiles:
        flattmp = pyfits.open(name)
        flatdata = biastrim(flattmp[0].data.astype(numpy.float32),bias,bpm) 
        exptime = flattmp[0].header['elaptime']

        if weight==0:
            flatavg = flatdata
        else:
            flatavg += flatdata
        weight += exptime
        del flatdata
    del flattmp
    flatavg /= weight

    return flatavg


def response(flat,orders,solutions):
    wide_orders = []
    for i,j in orders:
        wide_orders.append([i-3,j+3])
    norm = straighten(flat,solutions,wide_orders)
    for k in range(len(orders)):
        i,j = orders[k]
        response = stats.stats.nanmedian(norm[i:j],axis=0)
        response[response==0] = 1.
        response[numpy.isnan(response)] = 1.
        i,j = wide_orders[k]
        norm[i:j] /= response
    norm = curve(norm,solutions,wide_orders)
    return norm

def find_orders(flat):
    import pylab
    ncols = flat.shape[1]

    mid = flat[:,ncols/2].copy()
    mid[mid<1.] = mid[mid>1.].min()
    midtrace = numpy.log(mid)
    from scipy import signal,ndimage
    laplace = numpy.array([-1,-1,-1,1,1,1])

    deriv = signal.convolve(midtrace,laplace,mode='same')#/midtrace
    deriv = ndimage.gaussian_filter1d(deriv,1)
    deriv[:5] = 0.
    deriv[-5:] = 0.
    std = deriv.std()
    peaks = ndimage.maximum_filter(deriv,9)
    right = numpy.where((peaks==deriv)&(peaks>std))[0]
    peaks = ndimage.minimum_filter(deriv,9)
    left = numpy.where((peaks==deriv)&(peaks<-1.*std))[0]


    orders = []
    for i in range(len(left)):
        avg = midtrace[left[i]+5:right[i]-5].mean()
        std = midtrace[left[i]+5:right[i]-5].std()

        good = numpy.where(midtrace[left[i]:right[i]]>avg-1.*std)[0]+left[i]
        orders.append([good[0],good[-1]])
    return orders
