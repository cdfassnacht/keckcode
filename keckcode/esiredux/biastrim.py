import numpy
from astropy.io import fits as pyfits

def make_bias(filelist):
    count = 0
    n = 1
    for file in filelist:
        data = pyfits.open(file)[0]
        namps = data.header['NUMAMPS']
        data = data.data.astype(numpy.float32).T
        if data.size<1500*3000:
            n = 2
        if count==0:
            if namps==2:
                bias = numpy.empty((len(filelist),2048/n,4096/n))
            else:
                bias = numpy.empty((len(filelist),2046/n,4096/n))
        if namps==2:
            bias[count] = data[24/n:2072/n].copy()
        else:
            bias[count] = data[14/n:2060/n].copy()
        count += 1
    bias = numpy.median(bias,axis=0)
    diff = bias.copy()
    diff[:1024/n] -= numpy.median(bias[:1024/n])
    diff[1024/n:] -= numpy.median(bias[1024/n:])
    pix = diff.flatten()
    pix.sort()
    std = pix[10000/n:-10000/n].std()
    bpm = numpy.where(abs(diff)>10.*std,numpy.nan,1)
    return bias,bpm


def biastrim(data,bias,bpm):
    # The bias floats, so the overscan region is better than bias frames
    if data.size>1500*3000:
        data = data.T.astype(float)
        b1 = numpy.median(data[2074:2150],0)
        b2 = numpy.median(data[2152:2230],0)
        data[24:1048] -= b1
        data[1048:2072] -= b2

        # There is a small gain difference between the amplifiers
        off = numpy.median((data[1048]/data[1047])[data[1047]!=0])
        off = 1.08
        data[24:1048] *= off
        return data[24:2072]*bpm*1.29
    else:
        data = data.T.astype(float)
        b1 = numpy.median(data[1037:1075],0)
        b2 = numpy.median(data[1076:1115],0)
        data[12:524] -= b1
        data[524:1036] -= b2
        off = numpy.median(data[524]/data[523])
        off = 1.08
        data[12:524] *= off
        return data[12:1036]*bpm*1.29
