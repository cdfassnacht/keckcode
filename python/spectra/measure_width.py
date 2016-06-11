import scipy,pyfits
from scipy import ndimage

import special_functions

def clipped(data,clip=3.5):
        d = data.copy()
        mean = d.mean()
        std = d.std()

        delta = 1
        while delta:
                size = d.size
                d = d[abs(d-mean)<clip*std]
                delta = size-d.size
                if delta==size:
                        return 0.
                mean = d.mean()
                std = d.std()
        return mean,std


def measure(spectrum,num=25):
	data = spectrum.astype(scipy.float64)
	size = data.size

#	sorted = scipy.sort(data)
#	zero = sorted[size/20:size/10].mean()
#	sigma = sorted[size/20:size/10].std()

	pixels = scipy.arange(size)*1.

#	thresh = zero+100*sigma

	avg,std = clipped(data,3.)
	thresh = avg + 30.*std

	mask = ndimage.maximum_filter(data,13)
	mask = scipy.where(mask==data,data,0)
	lines = scipy.where((mask>thresh)&(mask<50000.))[0]
	count = 0

	vals = []
	for pos in lines:
		if pos<5 or pos+6>size:
			continue
		fitdata = scipy.zeros((11,2))
		fitdata[:,0] = pixels[pos-5:pos+6]
		fitdata[:,1] = data[pos-5:pos+6]

		par = scipy.zeros(4)
		par[1] = data[pos]
		par[2] = pixels[pos]
		par[3] = 1.

		fit,chi2 = special_functions.ngaussfit(fitdata,par)
		fit[2] -= fitdata[0,0]
		if fit[3]>4. or fit[3]<0.8 or fit[2]<2 or fit[2]>9:
			continue

		vals.append(fit[3])

		count += 1
	vals = scipy.asarray(vals)
	return scipy.median(vals)

