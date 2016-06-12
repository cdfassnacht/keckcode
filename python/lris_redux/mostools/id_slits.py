import scipy
from scipy import ndimage,signal,stats

"""
New version; not well-tested.
"""
def new_id_slits(flat):
        ncols = flat.shape[1]
	midtrace = flat[:,ncols/2].astype(scipy.float64)
	midtrace[midtrace<=0.] = 1.

        laplace = scipy.array([-1.,-1.,-1.,1.,1.,1.])
        deriv = signal.convolve(midtrace,laplace,mode='same')/midtrace
	d2 = abs(signal.convolve(midtrace,[-1.,1.],mode='same'))
        deriv = ndimage.gaussian_filter1d(deriv,1)
        deriv[:5] = 0.
        deriv[-5:] = 0.
	tmp = deriv.copy()
	tmp.sort()
	std = tmp[tmp.size*0.01:tmp.size*0.99].std()
        peaks = ndimage.maximum_filter(deriv,9)
        right = scipy.where((peaks==deriv)&(peaks>std))[0]
        peaks = ndimage.minimum_filter(deriv,9)
        left = scipy.where((peaks==deriv)&(peaks<-1.*std))[0]

        orders = []
	stars = []
        for i in range(len(left)):
                thresh = stats.stats.std(midtrace[left[i]+3:right[i]-3])
		good = scipy.where(d2[left[i]:right[i]]<thresh)[0]+left[i]
                orders.append([good[0],good[-1]])
		if deriv[right[i]]>3.*std:
			stars.append([good[0],good[-1]])
	print orders
	print stars
        return orders


"""
Old version; works well, though sometimes fails if ther is a substantial dip
  in the middle of a slit.
"""
def clip(arr,thresh=3.5):
	"""
	clip(arr,thresh=3.5)

	Simple sigma-clipping algorithm. Returns avg,std of clipped array.
	"""
	a = arr.copy()

	avg,std = a.mean(),a.std()
	while 1:
		size = a.size
		a = a[abs(a-avg)<thresh*std]
		avg,std = a.mean(),a.std()
		if size==a.size:
			break
	return avg,std


def id_slits(flat_data,findstars=True):
	"""
	id_slits(flat_data,small)

	Identifies slits and starboxes. There are issues with correctly
	  classifying short edge slits (which mimic star boxes) and dealing
	  with slits that have substantial dips.

	Inputs:
	  flat_data - 2d array containing the flat data (usually the brighter
	                pixels in a given row)
	  small     - parameter to describe how small star boxes are expected
	                to be (currently unused). Default = 20

	Outputs:
	  slit - list containing pairs describing the locations of the edges
	           of the slits
	  star - list containing pairs describing the locations of the edges
	          of the star boxes
	"""

	y_axis = flat_data.shape[0]

	data = flat_data.mean(axis=1)
	d = data.copy()

	"""
	The slits tend to be demarcated by when the sorted data begins to
	  grow at an accelerating rate; the first derivative tends to be an
	  acceptable proxy, though. The edges are masked for bad pixels/star
	  boxes.
	"""
	srt = scipy.sort(d)
	brk = signal.convolve(srt,[-1.,1.],mode='same')
	pix = brk[brk.size/10:brk.size*9/10].argmin()+brk.size/10

	lowvals = srt[pix]

	d[d<lowvals] = 0.
	d[d>0.] = 1.


	"""
	This needs to be tweaked to properly account for slits at the top and
	  bottom of the mask.
	"""
	edges = signal.convolve(d,[-1.,1.],mode='same')
	left = scipy.where(edges<0)[0]
	right = scipy.where(edges>0)[0]

	slits = []
	for i in range(left.size):
		slits.append([left[i],right[i]-1])

	if findstars is False:
		return slits

	"""
	The star boxes are identified by locating where the slit amplitudes
	  begin to spike. The current criterion is that a slit amplitude is
	  more than one sigma greater than the previous slit.
	"""
	amps = []
	for l,r in slits:
		amps.append(scipy.median(data[l:r]))
	amps = scipy.asarray(amps)
	args = amps.argsort()
	amps.sort()

	indx = amps.size-1
	for i in range(amps.size/2,amps.size):
		std = amps[:i].std()
		if amps[i]>amps[i-1]+std:
			indx = i
			break
	starindx = args[indx:]
	starindx.sort()
	stars = []
	for i in starindx:
		stars.append(slits[i])
	for i in starindx[::-1]:
		del slits[i]

	return slits,stars
