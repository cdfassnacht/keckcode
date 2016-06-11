import spectools,special_functions,scipy
import pyfits

from scipy import stats

def findoffset(scidata,coord,offset,slice=None):

	
	straight = spectools.resampley(scidata,coord,offset,slice=slice)
	slice = stats.stats.median(straight,axis=1)#straight.mean(axis=1)
	xvals = scipy.arange(slice.size)*1.

	slice = slice**4

	sum = slice.sum()
	slice *= xvals

	center = slice.sum()/sum
	return center
