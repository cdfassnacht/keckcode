import scipy
try:
	import pyfits
except:
	from astropy.io import fits as pyfits
from scipy import ndimage

def id_slits(flat_data):
	#
	# Identifies slits and starboxes.
	#

	y_axis = flat_data.shape[0]

	# Squash the data to one dimension
	if flat_data.ndim==2:
		slice = flat_data.mean(axis=1)
	else:
		slice = flat_data.copy()

	# We assume that the starboxes are at least 2.5x greater than
	#  the mean of the array.
	star_threshold = slice.mean()*2.5
	slit_threshold = slice.mean()/2.5

	# Pixels with noise levels that greatly differ from the "model" noise
	#  level (ie Poisson shot noise) are set to zero, as are pixels below
	#  the slit threshold. This finds edges, even when two slits are close
	#  together and the edge floor is above the slit threshold.
	original = slice.copy()
	for i in range(y_axis):
		start = i-2
		end = i+3
		if start<0:
			start = 0
		if end>y_axis:
			end = y_axis
		tmp = original[start:end].std()
		sigma = scipy.sqrt(original[i])
		if tmp>5*sigma:
			slice[i] = 0.
		if slice[i]<slit_threshold:
			slice[i] = 0.

	# Narrow the slits by two pixels
	mask = ndimage.minimum_filter1d(slice,5)
	mask[mask>0] = 1
	slice *= mask

	star = []
	slit = []
	start = 0
	i = 0
	while i<y_axis:
		# Check for a rising edge
		if slice[i]<slit_threshold:
			while i<y_axis and slice[i]<slit_threshold:
				i += 1
			start = i

		# Checking for the falling edge of a slit
		while i<y_axis and slice[i]>slit_threshold:
			i += 1
		end = i

		if start==end:
			i += 1
			continue

		sigma = slice[start:end].std()
		if scipy.sqrt(slice[start:end].mean())>sigma:
			sigma = scipy.sqrt(slice[start:end].mean())
		while start>0 and abs(original[start-1]-original[start])<sigma:
			start -= 1
		while end<y_axis and abs(original[end]-original[end-1])<sigma:
			end += 1
		slice[start:end] = original[start:end].copy()
		i = end

		if end-start<6:
			continue

		if slice[start:end].mean()>star_threshold:
			limits = scipy.where(slice[start:end]>star_threshold)
#			end = start + limits[0][-1] - limits[0][0] + 1 
			start += limits[0][0]
			end = start + limits[0].size
			star.append([start,end])
		elif end-start>10:
			slit.append([start,end])

	return slit,star

	slice = flat_data.copy()
	mask = ndimage.maximum_filter(slice,5)
	std = scipy.sqrt(mask)
	diff = abs(slice-mask)
	mask = scipy.where(diff<std,slice,0)
	mask = scipy.where(mask>slit_threshold,1,0)

	slit = []
	star = []
	start = 0
	i = 0
	while i<slice.size:
		# In a valley
		if mask[i]==0:
			i += 1
			continue
		# Not in a valley!
		start = i
		while i<slice.size and mask[i]!=0:
			i += 1
		end = i
		if end-start<6:
			continue
		if slice[start:end].mean()>star_threshold:
			limits = scipy.where(slice[start:end]>star_threshold)
			start += limits[0][0]
			end = start + limits[0].size
			star.append([start,end])
		elif end-start>10:
			slit.append([start,end])

	return slit,star
