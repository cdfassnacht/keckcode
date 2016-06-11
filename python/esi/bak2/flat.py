from spectra import ycorrect,id_slits,spectools
from esi.biastrim import biastrim
from esi.straighten import straighten,curve

from special_functions import lsqfit,genfunc

import numpy,scipy,pyfits,pickle
from scipy import stats


def make_flat(flatfiles,out_prefix):
    # Read data from files, biastrim and coadd
    weight = 0
    bias = pyfits.open(out_prefix+"_bias.fits")[0].data.copy()
    bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.copy()
    imgs = []
    for name in flatfiles:
        flattmp = pyfits.open(name)
        flatdata = biastrim(flattmp[0].data.copy(),bias,bpm) 
        exptime = flattmp[0].header['elaptime']

        imgs.append(flatdata/exptime)
    return numpy.median(imgs,0)


def response(flat,orders,solutions,wideorders):
    from scipy import ndimage,interpolate
    if flat.shape[0]>1500:
        edge = 3
    else:
        edge = 1
    if flat.shape[1]>3000:
        kwidth = 5.
    else:
        kwidth = 3.

    tmpOrders = []
    for low,high in orders:
        tmpOrders.append([low-edge,high+edge])
    norm = straighten(flat,solutions,tmpOrders,wideorders)
    xvals = numpy.arange(norm.shape[1])*1.
    for k in range(len(orders)):
        i,j = orders[k]
        response = stats.stats.nanmedian(norm[i:j],axis=0)
        response[response==0] = 1.
        response[numpy.isnan(response)] = 1.
        smooth = ndimage.gaussian_filter(response,kwidth)
        i,j = tmpOrders[k]
        norm[i:j] /= smooth
    norm = curve(norm,solutions,tmpOrders,wideorders)
    return norm


def find_orders(flat):
    if flat.shape[0]>1500:
        fw = 9
        border = 5
    else:
        fw = 7
        border = 3
    ncols = flat.shape[1]

    mid = flat[:,ncols/2].copy()
    mid[mid<1.] = mid[mid>1.].min()
    midtrace = numpy.log(mid)
    from scipy import signal,ndimage
    laplace = numpy.array([-1,-1,-1,1,1,1])

    deriv = signal.convolve(midtrace,laplace,mode='same')#/midtrace
    deriv = ndimage.gaussian_filter1d(deriv,1)
    deriv[:border] = 0.
    deriv[-border:] = 0.
    std = deriv.std()
    peaks = ndimage.maximum_filter(deriv,fw)
    right = numpy.where((peaks==deriv)&(peaks>std))[0]
    peaks = ndimage.minimum_filter(deriv,fw)
    left = numpy.where((peaks==deriv)&(peaks<-1.*std))[0]

    if right[0]<left[0]:
        right = right[1:]

    orders = []
    for i in range(len(left)):
        avg = midtrace[left[i]+border:right[i]-border].mean()
        std = midtrace[left[i]+border:right[i]-border].std()

        good = numpy.where(midtrace[left[i]:right[i]]>avg-1.*std)[0]+left[i]
        orders.append([good[0],good[-1]])
    return orders
