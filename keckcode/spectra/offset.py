import scipy
from scipy import stats
try:
	import pyfits
except:
	from astropy.io import fits as pyfits

import special_functions
from . import spectools

def findoffset(scidata,coord,offset,slice=None):

	
	straight = spectools.resampley(scidata,coord,offset,slice=slice)
	slice = stats.stats.median(straight,axis=1)#straight.mean(axis=1)
	xvals = scipy.arange(slice.size)*1.

	slice = slice**4

	sum = slice.sum()
	slice *= xvals

	center = slice.sum()/sum
	return center
