"""
Helper routines to perform bias subtraction and overscan trimming of LRIS data.
"""

import numpy,pyfits
SIDE = None
def biastrim(hdu):
    """
    Subtracts bias from data and returns the overscan region-subtracted
	image.
    """
    if type(hdu)!=type(pyfits.HDUList([])):
        return oldredside(hdu)
    namps = (len(hdu)-1)/2
    if hdu[0].header['INSTRUME'].strip()=='LRIS':
        order = [2,1,4,3]
        #gain = [1.022,0.955,0.877,0.916]
        gain = [1.255,1.180,1.191,1.162] # From web on 01 Sept 2011
        blue = False
        if namps==1:
            order = [1,2]
            gain = [0.955,0.916]
        if len(hdu)==1:
            order = [0]
            gain = [2.]
    else:
        order = [1,2,3,4]
        gain = [1.55,1.56,1.63,1.70]
        blue = True
        if namps==1:
            order = [1,2]
            gain = [1.56,1.70]

#    xbin,ybin = [int(i) for i in hdu[0].header['CCDSUM'].split()]
    xbin,ybin = eval(hdu[0].header['BINNING'])
    xshape = 4096/xbin
    yshape = 4096/ybin

    oimg = numpy.empty((yshape,xshape))
#    xi,xf = [int(i)+1 for i in hdu[1].header['DATASEC'].split(',')[0][1:].split(':')]
#    yi,yf = [int(i) for i in hdu[1].header['DATASEC'].split(',')[1][:-1].split(':')]
#    xi,yi = xi-1,yi-1
    sign = 1
    i = 0
    xmin,ymin,xmax,ymax = 4096,4096,0,0
    for h in order:
        xi,xf,yi,yf = eval(hdu[h].header['DATASEC'].replace(':',','))
        xi -= 1
        yi -= 1
        d = hdu[h].data[yi:yf,xi:xf][:,::sign].astype(numpy.float32)
        if hdu[0].header['INSTRUME'].strip()=='LRIS':
            b = numpy.median(hdu[h].data[yf+2:,xi:xf].astype(numpy.float32))
        else:
            b = numpy.median(hdu[h].data[yi:yf,-10:].astype(numpy.float32).T,0)
        x1,x2,y1,y2 = eval(hdu[h].header['DETSEC'].replace(':',','))
        if y1>y2:
            y1,y2 = y2,y1
        if x1>x2:
            x1,x2 = x2,x1
        x1 -= 1
        y1 -= 1
        x1,x2,y1,y2 = x1/xbin,x2/xbin,y1/ybin,y2/ybin
        oimg[y1:y2,x1:x2] = (d.T-b).T*gain[h-1]
#        oimg[:,x0:x0+xshape/4] = (hdu[h].data[yi:yf,xi:xf][:,::sign].astype(numpy.float32).T-numpy.median(hdu[h].data[:,-10:].astype(numpy.float32).T,0)).T*gain[i]
        if x1<xmin:
            xmin = x1
        if y1<ymin:
            ymin = y1
        if x2>xmax:
            xmax = x2
        if y2>ymax:
            ymax = y2
        if namps==2:
            sign *= -1
    outdata = oimg.T
    ymin,ymax,xmin,xmax = xmin,xmax,ymin,ymax
    yshape = outdata.shape[0]
    if SIDE is None:
         return {'bottom':outdata[ymin:yshape/2,xmin:xmax].copy(),'top':outdata[yshape/2:ymax,xmin:xmax].copy()}
    elif SIDE.lower()=='bottom':
        return outdata[ymin:yshape/2,xmin:xmax]
    elif SIDE.lower()=='top':
        return outdata[yshape/2:ymax,xmin:xmax]
    else:
        return outdata[ymin:ymax,xmin:xmax]#{'bottom':outdata[ymin:yshape/2,xmin:xmax].copy(),'top':outdata[yshape/2:ymax,xmin:xmax].copy()}
    return outdata


def dummyblueside(data):
    return data.astype(numpy.float32)


def oldredside(data):
    if data.shape[1]==2148:
        bias = (data[:,2:21].mean(axis=1)*18+data[:,2069:2148].mean(axis=1)*80)/98.
        out_data = (data[:,21:2069].T-bias).T
        return out_data
    elif data.shape[1]>3000:
        return blueside(data)

    xsize = data.shape[0]
    bias = numpy.empty((2,data.shape[0]))
    bias[0] = (data[:,1:20].mean(axis=1)*19+data[:,2088:2168].mean(axis=1)*80)/99.
    bias[1] = (data[:,20:38].mean(axis=1)*18+data[:,2168:2248].mean(axis=1)*80)/98.
    out_data = numpy.empty((2046,xsize))

    """
    Mask out the bad columns. Note this might not be appropriate for older
        data (or if the CCDs change).
    """
    mask = (data[0:1490,995]+data[0:1490,997])/2.
    data[0:1490,996] = mask.copy()
    mask = (data[0:1493,999]+data[0:1493,1001])/2.
    data[0:1493,1000] = mask.copy()

    data = data.transpose()

    out_data[0:1023,:] = data[41:1064,:] - bias[0]
    out_data[1023:2046,:] = data[1064:2087,:] - bias[1]

    """
    Fix difference in amplifier gains. This does *not* convert from DN
        to electrons.
    """
    out_data[1023:2046,:] *= 1.0765 # Note this differs from the LRIS
                                        #  website that would suggest 1.0960
    out_data = out_data.transpose()
    return out_data


def blueside(data):
    """
    Subtracts bias from data and returns the overscan region-subtracted image.
    """
    data = data.T
    data = data[:,::-1]
    size = data.shape[1]

    bias = numpy.empty((4,size))
    bias[0] = (data[0:50,:].mean(axis=0) + data[4300:4380,:].mean(axis=0))/2.
    bias[1] = (data[52:101,:].mean(axis=0) + data[4380:4460,:].mean(axis=0))/2.
    bias[2] = (data[102:153,:].mean(axis=0) + data[4460:4540,:].mean(axis=0))/2.
    bias[3] = (data[153:202,:].mean(axis=0) + data[4540:4620,:].mean(axis=0))/2.

    """
    Conversion factor from DN to electrons (from LRIS website)
    """
    gain = [1.55,1.56,1.63,1.70]
    outdata = numpy.empty((4096,4096))
    for i in range(4):
        outstart = i*1024
        datastart = i*1024 + 204
        outdata[outstart:outstart+1024,:] = data[datastart:datastart+1024,:] - bias[i]
        outdata[outstart:outstart+1024,:] *= gain[i]
    del bias

    if SIDE is None:
        return {'bottom':outdata[:2048].copy(),'top':outdata[2048:].copy()}
    elif SIDE.lower()=='bottom':
        return outdata[:2048]
    elif SIDE.lower()=='top':
        return outdata[2048:]
    else:
        return {'bottom':outdata[:2048].copy(),'top':outdata[2048:].copy()}
    return outdata

