from spectra import ycorrect,id_slits,spectools
from esi.biastrim import biastrim

#from special_functions import lsqfit,genfunc
from linear import lsqfit,genfunc

import scipy,pyfits,pickle


def starstack(starfiles,out_prefix):
	# Read data from files, biastrim and coadd
	bias = pyfits.open(out_prefix+"_bias.fits")[0].data.astype(scipy.float32)
	bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.astype(scipy.float32)
	first = True
	for name in starfiles:
		startmp = pyfits.open(name)[0].data.astype(scipy.float32)
		stardata = biastrim(startmp,bias,bpm) 

		if first:
			first = False
			star = stardata
		else:
			star += stardata
		del stardata,startmp
	return star

def findpeaks(stars,peak,table,start,end,step,thresh):
	WIDTH = 4
	peak = peak
	standard = peak

	nrows = stars.shape[0]
	xvals = scipy.arange(nrows).astype(scipy.float32)

	solution = {'type':'chebyshev','coeff':scipy.array([[peak]])}
	found = []
	order = 0
	for col in range(start,end,step):
		peak = genfunc(col,0,solution)[0]
		if peak-WIDTH<0 or peak+WIDTH+1>nrows:
			break
		peak = stars[peak-WIDTH:peak+WIDTH+1,col].argmax()+peak-WIDTH
		x = xvals[peak-WIDTH:peak+WIDTH+1]-peak
		flux = stars[peak-WIDTH:peak+WIDTH+1,col]
		good = scipy.isfinite(flux)
		if good.size==0:
			continue
		x = x[good]
		flux = flux[good]
		f = flux.sum()
		center = (flux*x).sum()/f
		if abs(center)>2 or f<thresh:
			continue
		peak += center
		found.append([col,peak])
		table.append([col,standard,peak])
		if len(found)<8 and len(found)%2**order!=0:
			continue
		if len(found)>8 and abs((col-start)/step)%2**order!=0:
			continue
		if order<3:
			order += 1
		if len(found)>8:
			solution = lsqfit(scipy.asarray(found),'chebyshev',3)
		elif len(found)==8:
			solution = lsqfit(scipy.asarray(found),'chebyshev',2)
		elif len(found)==4:
			solution = lsqfit(scipy.asarray(found),'chebyshev',1)
		else:
			solution = {'type':'chebyshev','coeff':scipy.array([[peak]])}
	return table


def startrace(stars,orders):
	from scipy import ndimage,signal
	step = 5
	WIDTH = 4

	nrows = stars.shape[0]
	ncols = stars.shape[1]


	start = stars[:,ncols/2].copy()
	background = ndimage.percentile_filter(start,10.,37)
	start -= background
	peaks = ndimage.maximum_filter(start,9)
	peaks2 = ndimage.maximum_filter(start,3)
	
	thresh = stars.flatten()
	thresh[scipy.isnan(thresh)] = 0.
	thresh.sort()
	thresh = thresh[20000:-1000000].std()
	all = scipy.where((peaks==start)&(peaks>thresh)&(peaks==peaks2))[0]

	peaks = []
	for i in range(len(orders)):
		start,end = orders[i]
		peaks.append([])
		for p in all:
			if p>start-5 and p<end+5:
				peaks[i].append(p)

	solutions = []
	for p in peaks:
		table = []
		xvals = scipy.arange(nrows).astype(scipy.float32)
		for peak in p:
			if peak-WIDTH<0 or peak+WIDTH+1>nrows:
				continue
			x = xvals[peak-WIDTH:peak+WIDTH+1]-peak
			flux = stars[peak-WIDTH:peak+WIDTH+1,ncols/2]
			good = scipy.isfinite(flux)
			x = x[good]
			flux = flux[good]
			peak += (flux*x).sum()/flux.sum()

			table = findpeaks(stars,peak,table,ncols/2,ncols,step,thresh)
			table = findpeaks(stars,peak,table,ncols/2-step,0,-1*step,thresh)
		transform_table = scipy.asarray(table)

		ytrue = lsqfit(transform_table,"chebyshev",3,3)
		transform_table[:,1],transform_table[:,2] = transform_table[:,2],transform_table[:,1]

		ymap = lsqfit(transform_table,"chebyshev",3,3)
		solutions.append([ytrue,ymap])

	return solutions

def straighten(data,solutions,orders):
	output = data*scipy.nan
	from scipy import ndimage

	data[scipy.isnan(data)] = 0.
	coords = spectools.array_coords(data.shape)
	x = coords[1].flatten()
	y = coords[0].flatten()
	prefilter = True
	for i in range(len(solutions)):
		low,high = orders[i]
		ytrue,ymap = solutions[i]
		yforw = genfunc(x,y,ytrue).reshape(coords[0].shape)
		yback = genfunc(x,y,ymap).reshape(coords[0].shape)
		coords[0] = yforw.copy()

		ymax = scipy.where(yback<high+4)[0].max()
		ymin = scipy.where(yback>low-4)[0].min()
		tmpcoords = coords[:,ymin:ymax,:].copy()
		tmp = ndimage.map_coordinates(data,tmpcoords,output=scipy.float64,cval=0.,order=5,prefilter=prefilter)
		prefilter = False

		output[low:high] = tmp[low-ymin:high-ymin].copy()

	return output
