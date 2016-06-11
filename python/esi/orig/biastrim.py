import pyfits,numpy

def make_bias(filelist):
    count = 0
    for file in filelist:
        data = pyfits.open(file)[0]
        namps = data.header['NUMAMPS']
        data = data.data.astype(numpy.float32).T
        if count==0:
            if namps==2:
                bias = numpy.empty((len(filelist),2048,4096))
            else:
                bias = numpy.empty((len(filelist),2046,4096))
        if namps==2:
            bias[count] = data[24:2072].copy()
        else:
            bias[count] = data[14:2060].copy()
        count += 1
    bias = numpy.median(bias,axis=0)
    diff = bias.copy()
    diff[:1024] -= numpy.median(bias[:1024])
    diff[1024:] -= numpy.median(bias[1024:])
    pix = diff.flatten()
    pix.sort()
    std = pix[10000:-10000].std()
    bpm = numpy.where(abs(diff)>10.*std,numpy.nan,1)
    return bias,bpm


def biastrim(data,bias,bpm):
    data = data.T
    if bias.shape[0]==2048:
        data = data[24:2072].copy()-bias
        offset = numpy.median(data[1024]/data[1023])
        data[0:1024] *= offset
    else:
        data = data[14:2060].copy()-bias
    return data*bpm
