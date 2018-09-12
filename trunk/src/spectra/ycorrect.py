import scipy,special_functions
from scipy import ndimage
#import pylab as p

def ycorrect(data):

	# Parameters
	SUMWIDTH = 41	# Width of summing over columns

	y_axis = data.shape[0]
	x_axis = data.shape[1]


	# Determine the location of the central column
	central = x_axis/2

	# Determine the centers of the holes in the center column to use as the
	#  reference for all other columns
	x_min_orig = central - SUMWIDTH/2
	x_max_orig = central + SUMWIDTH/2

	midcol = data[:,x_min_orig:x_max_orig].mean(axis=1)
	central_edges,threshold,star_cutoff = find_holes(midcol)
	transform_table = scipy.zeros((1,3),'f4')

	index = 0
	for peak in central_edges:
		if index:
			transform_table.resize((index+1,3))
		transform_table[index,0] = central
		transform_table[index,1] = peak
		transform_table[index,2] = peak

		index += 1

	offset = scipy.zeros(len(central_edges))

	x_min = x_min_orig
	x_max = x_max_orig
	current_column = central

	while current_column>SUMWIDTH + 20:
		current_column = current_column - SUMWIDTH - 10
		x_min = x_min - SUMWIDTH - 10
		x_max = x_max - SUMWIDTH - 10

		comp_array = data[:,x_min:x_max].mean(axis=1)
		comp_array.clip(min=-1000.,max=star_cutoff)
		derivative = deriv_1d(comp_array)
		derivative = ndimage.gaussian_filter1d(derivative,3)
		derivative = abs(derivative)

		for i in range(offset.size):
			if scipy.isnan(offset[i]):
				continue
			ref = central_edges[i] + offset[i]

			start = int(ref) - 6
			end = start + 13

			if derivative[start:end].max()<threshold:
				offset[i] = scipy.nan
				continue

			fit = find_peak(derivative[start:end])

			if(fit[2]<0 or fit[2]>13 or fit[3]<1 or fit[3]>6):
				offset[i] = scipy.nan
				continue

			peak = fit[2]+float(start)			
			offset[i] = peak - central_edges[i]

			transform_table.resize((index+1,3))

			transform_table[index,0] = current_column
			transform_table[index,1] = central_edges[i]
			transform_table[index,2] = peak

			index += 1

	offset = scipy.zeros(offset.size)

	x_min = x_min_orig
	x_max = x_max_orig
	current_column = central
	while current_column<x_axis - SUMWIDTH - 19:
		current_column = current_column + SUMWIDTH + 10
		x_min = x_min + SUMWIDTH + 10
		x_max = x_max + SUMWIDTH + 10

		comp_array = data[:,x_min:x_max].mean(axis=1)
		comp_array.clip(min=-1000.,max=star_cutoff)
		derivative = deriv_1d(comp_array)
		derivative = ndimage.gaussian_filter1d(derivative,3)
		derivative = abs(derivative)

		for i in range(offset.size):
			if scipy.isnan(offset[i]):
				continue
			ref = central_edges[i] + offset[i]
 
			start = int(round(ref)) - 6
			end = start + 13
 
			if derivative[start:end].max()<threshold:
				offset[i] = scipy.nan
				continue
 
			fit = find_peak(derivative[start:end])

			if(fit[2]<0 or fit[2]>13 or fit[3]<1 or fit[3]>6):
				offset[i] = scipy.nan
				continue
 
			peak = fit[2]+float(start)
			offset[i] = peak - central_edges[i]

			transform_table.resize((index+1,3))

			transform_table[index,0] = current_column
			transform_table[index,1] = central_edges[i]
			transform_table[index,2] = peak

			index += 1

	true_coeffs = special_functions.lsqfit(transform_table,"chebyshev",4,4)

	temp = transform_table[:,1].copy()
	transform_table[:,1] = transform_table[:,2].copy()
	transform_table[:,2] = temp.copy()

	map_coeffs = special_functions.lsqfit(transform_table,"chebyshev",4,4)

	# The xtrue_coeffs give the "true" y value for that pixel. xmap_coeffs
	#  describe where the coordinate should be mapped to in an output grid
	#  (in other words the true pixel value for that y).

	return true_coeffs,map_coeffs


########################
#
# Supporting Functions
#
########################

def find_holes(data):
	sample = data.copy()
	size = sample.size

	# Here's a little hack to "flatten" star boxes
	tmp = scipy.sort(sample)
	star_cutoff = scipy.median(tmp[-30:-10])*0.6
	sample = scipy.where(sample>star_cutoff,star_cutoff,sample)

	derivative = deriv_1d(sample)
	derivative = ndimage.gaussian_filter1d(derivative,3)
	derivative = abs(derivative)


	tmp = scipy.sort(derivative)
	avg = scipy.median(tmp[size/8:size*3/8])
	sigma = tmp[size/8:size*3/8].std()


	threshold = avg + sigma*100.

	edge = []

	count = 0
	while derivative.max()>threshold:
		start = derivative.argmax()-7
		end = derivative.argmax()+8

		if start<0:
			start = 0
		if end>derivative.size:
			end = derivative.size

		fit = find_peak(derivative[start:end])

		if start>7 and end<derivative.size-7:
			edge.append(float(start) + fit[2])

		start -= 3
		end += 3

		if start<0:
			start = 0
		if end>derivative.size:
			end = derivative.size

		derivative[start:end] = 0.

	edge.sort()
	return edge,threshold,star_cutoff



def find_peak(data):
	size = data.shape[0]
	fit_data = scipy.reshape(scipy.arange(0,size,0.5),(size,2))
	fit_data[:,1] = data.copy()

	fit = scipy.zeros(4)
	fit[0] = 0.
	fit[1] = data.max()
	fit[2] = data.argmax()
	fit[3] = 2.
	fit,val = special_functions.ngaussfit(fit_data,fit)

	return fit

def deriv_1d(data):
	size = data.shape[0]
	out = scipy.zeros(size,'f8')
	out[0] = data[1] - data[0]
	for i in range(1,size-1):
		out[i] = (data[i+1]-data[i-1])/2.
	out[size-1] = data[size-1]-data[size-2]

	return out

