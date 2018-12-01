import spectools,special_functions,scipy
try:
	import pyfits
except:
	from astropy.io import fits as pyfits

def findoffset(scidata,coord,offset,slice=None):

	straight = spectools.resampley(scidata,coord,offset,slice=slice)
	slice = straight.mean(axis=1)

	p = scipy.zeros(4)
	p[0] = slice.min()
	p[1] = slice.max()
	p[2] = slice.argmax()
	p[3] = 1.
	fit,val = special_functions.ngaussfit(slice,p)
	chi2 = val/(slice.size-4.)
	if fit[2]>0 and fit[2]<slice.size:
		return fit[2]
	else:
		return scipy.nan
