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


def id_slits(flat_data,small=20.):
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

	right_edges = signal.convolve(flat_data,[[-1.],[0.],[0.],[1.]],mode='same')
	right_edges = right_edges.mean(axis=1)
	right_edges[0] = 0.
	right_edges[-1] = 0.

	left_edges = right_edges*-1.
	deriv = abs(right_edges)


	"""
	Find sigma-clipped avg and std, then filter to determine edges
	"""
	avg,std = clip(right_edges[right_edges<1000],5)
	max = ndimage.maximum_filter(right_edges,7)==right_edges
	coords = scipy.where(max&(right_edges>15.*std))[0]
#	mid = right_edges[coords]
#	avg,std = clip(mid,1.5)
#	avg = scipy.median(mid)
#	std = avg/2.
#	right_edges = coords[right_edges[coords]>avg+std].tolist()
	right_edges = coords.tolist()
	print coords

	avg,std = clip(left_edges[left_edges<1000],5)
	max = ndimage.maximum_filter(left_edges,7)==left_edges
	coords = scipy.where(max&(left_edges>15.*std))[0]
#	mid = left_edges[coords]
#	avg,std = clip(mid,1.5)
#	avg = scipy.median(mid)
#	std = avg/2.
#	left_edges = coords[left_edges[coords]>avg+std].tolist()
	left_edges = coords.tolist()
	print coords

	"""
	The edge finding may fail at the start/end of the CCD; this performs a
	  simple test/correction by asserting that the first left edge must
	  come before the first right edge (and likewise, reversed, for the
	  final edge).
	"""
	if right_edges[0]<left_edges[0]:
		left_edges.insert(0,0)
	if left_edges[-1]>right_edges[-1]:
		right_edges.append(y_axis)

	"""
	Check for 'doubled' peaks; these can be caused by reflections, for
	  example.
	"""
	mask = []
	coords = []
	for i in left_edges:
		coords.append(i)
		mask.append(0)
	for i in right_edges:
		coords.append(i)
		mask.append(1)
	coords = scipy.asarray(coords)
	mask = scipy.asarray(mask)
	args = coords.argsort()
	coords.sort()
	mask = mask[args]
	for i in range(mask.size-1):
		if mask[i]==mask[i+1]:
			if mask[i]==1:
				right_edges.remove(coords[i+1])
			else:
				left_edges.remove(coords[i])

	slice = flat_data.mean(axis=1)

	store = []
	separate = []
	widths = []
	slit = []
	star = []
	index = 0
	for i in range(len(right_edges)):
		start = left_edges[i]
		end = right_edges[i]

		if end-start<9:
			continue
		sigma = stats.stats.std(slice[start+3:end-3])
		mean = slice[start+3:end-3].mean()
		if scipy.sqrt(mean)>sigma:
			sigma = scipy.sqrt(mean)
		
		good = scipy.where(deriv[start:end]<sigma)[0]+start
		start,end = good[0],good[-1]+1

		if end-start<4:
			continue
		store.append([start,end,slice[start:end].mean()])
		width = float(end-start)
#		separate.append(slice[start:end].mean())#/width**0.5)
		separate.append(slice[start:end].max())
		widths.append(end-start)


	"""
	Remove the most obvious star box (ie with the greatest maximum
	  flux) then clip to find the 'true' mean/std.
	"""
	separate = scipy.asarray(separate)
	args = separate.argsort()
	separate.sort()
	separate = separate[:-3]
	mean,std = clip(separate)
	cutoff = mean+5*std
	print cutoff
	for start,end,flux in store:
		print flux
		x = flux
		if x>cutoff:
			start -= 3
			end += 3
			if start<0:
				start = 0
			if end>y_axis:
				end = y_axis
			star.append([start,end])
		# Include short slits, but not short slits on the edges that
		#   might be due to clipping on the red side.
		elif end-start<10:
			count = len(star)+len(slit)
			if count==0 or count==len(store)-1:
				continue
			start -= 3
                        end += 3
                        if start<0:
                                start = 0
                        if end>y_axis:
                                end = y_axis
                        star.append([start,end])
		else:
			slit.append([start,end])

	return slit,star		
