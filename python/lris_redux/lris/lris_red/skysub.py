"""
Background subtract, resample, and coadd 2d spectra.
"""

from mostools import spectools,correct_telluric,skysub
from special_functions import genfunc

import scipy
from scipy import optimize,interpolate,ndimage,signal,stats
from astropy.io import fits as pyfits

RESAMPLE = 0		# 1 to only resample
magicnum = -2**15



def offsety(data,outshape,offset):
	coords = spectools.array_coords(outshape)
	coords[0] += offset
	return ndimage.map_coordinates(data,coords,output=scipy.float64,cval=-999)


"""
Many of the pixels in the slit stack will be unused border pixels. This avoids
  taking the median of all of those (and therefore saves a bit of time).
"""
def fastmed(arr):
	tmp = arr.copy()
	tmp[scipy.isnan(tmp)] = 0.
	tmp = tmp.sum(0)
	tmp = tmp.sum(0)
	ends = scipy.where(tmp!=0)[0]
	left = ends.min()
	right = ends.max()+1
	tmp = stats.stats.nanmedian(arr[:,:,left:right],0)
	out = arr[0]*scipy.nan
	out[:,left:right] = tmp.copy()
	return out

"""
doskysub()
"""

def doskysub(straight,ylen,xlen,sci,yback,sky2x,sky2y,ccd2wave,disp,mswave,offsets,cutoff,airmass):
	sci = sci.copy()

	# If cutoff is not a float, we are using the blueside
	locutoff = cutoff
	hicutoff = 10400.

	nsci = sci.shape[0]
	width = sci.shape[2]

	# Perform telluric correction
	coords = spectools.array_coords(sci[0].shape)
	x = coords[1].flatten()
	y = coords[0].flatten()

	for k in range(nsci):
		w = genfunc(x,y,ccd2wave[k])
		telluric = correct_telluric.correct(w,airmass[k],disp)
		sci[k] *= telluric.reshape(sci[k].shape)
	del coords,x,y,telluric

	# Create arrays for output images
	outcoords = spectools.array_coords((ylen,xlen))
	outcoords[1] *= disp
	outcoords[1] += mswave - disp*xlen/2.
	xout = outcoords[1].flatten()
	yout = outcoords[0].flatten()

	out = scipy.zeros((nsci,ylen,xlen))

	fudge = scipy.ceil(abs(offsets).max())
	bgimage = scipy.zeros((nsci,ylen+fudge,xlen))
	varimage = bgimage.copy()

	bgcoords = spectools.array_coords((ylen+fudge,xlen))
	bgcoords[1] *= disp
	bgcoords[1] += mswave - disp*xlen/2.

	#
	# Cosmic Ray Rejection and Background Subtraction
	#
	yfit = yback.flatten()
	ycond = (yfit>straight-0.4)&(yfit<straight+ylen-0.6)

	coords = spectools.array_coords(yback.shape)
	xvals = coords[1].flatten()
	yvals = coords[0].flatten()

	ap_y = scipy.zeros(0)
	aper = scipy.zeros(0)
	for k in range(nsci):
		xfit = genfunc(xvals,yfit-straight,ccd2wave[k])
		zfit = sci[k].flatten()

		x = xfit[ycond]
		y = yfit[ycond]
		z = zfit[ycond]

		# The plus/minus 20 provides a better solution for the edges
		wavecond = (x>locutoff-20.)&(x<hicutoff+20.)
		x = x[wavecond]
		y = y[wavecond]
		z = z[wavecond]

		# If only resampling...
		if RESAMPLE==1:
			coords = outcoords.copy()
			samp_x = genfunc(xout,yout,sky2x[k])
			samp_y = genfunc(xout,yout,sky2y[k])
			coords[0] = samp_y.reshape(coords[0].shape)
			coords[1] = samp_x.reshape(coords[1].shape)
			out[k] = scipy.ndimage.map_coordinates(sci[k],coords,output=scipy.float64,order=5,cval=-32768,prefilter=False)
			out[k][xout.reshape(coords[1].shape)<locutoff] = scipy.nan
			out[k][xout.reshape(coords[1].shape)>hicutoff] = scipy.nan
			out[k][out[k]==-32768] = scipy.nan
			continue


		bgfit = skysub.skysub(x,y,z,disp)

		background = zfit.copy()
		for indx in range(background.size):
			x0 = xfit[indx]
			y0 = yfit[indx]
			if x0<locutoff-10 or x0>hicutoff+10:
				background[indx] = scipy.nan
			else:
				background[indx] = interpolate.bisplev(x0,y0,bgfit)
		sub = zfit-background
		sub[scipy.isnan(sub)] = 0.
		sky = sub*0.
		sky[ycond] = sub[ycond]
		sky = sky.reshape(sci[k].shape)
		sub = sky.copy()

		background[scipy.isnan(background)] = 0.

		# Note that 2d filtering may flag very sharp source traces!
		sub = sub.reshape(sci[k].shape)
		sky = ndimage.median_filter(sky,5)
		diff = sub-sky
		model = scipy.sqrt(background.reshape(sci[k].shape)+sky)
		crmask = scipy.where(diff>4.*model,diff,0.)
		sub -= crmask
		sci[k] -= crmask

		# Create straightened slit
		coords = outcoords.copy()
		samp_x = genfunc(xout,yout,sky2x[k])
		samp_y = genfunc(xout,yout,sky2y[k])
		coords[0] = samp_y.reshape(coords[0].shape)
		coords[1] = samp_x.reshape(coords[1].shape)
		out[k] = scipy.ndimage.map_coordinates(sci[k],coords,output=scipy.float64,order=5,cval=magicnum,prefilter=False)
		out[k][xout.reshape(coords[1].shape)<locutoff] = scipy.nan
		out[k][xout.reshape(coords[1].shape)>hicutoff] = scipy.nan
		out[k][out[k]==magicnum] = scipy.nan

		# Output bgsub image
		coords = bgcoords.copy()
		bgy = bgcoords[0].flatten()+offsets[k]
		bgx = bgcoords[1].flatten()
		samp_x = genfunc(bgx,bgy,sky2x[k])
		samp_y = genfunc(bgx,bgy,sky2y[k])
		coords[0] = samp_y.reshape(coords[0].shape)
		coords[1] = samp_x.reshape(coords[1].shape)

		varimage[k] = scipy.ndimage.map_coordinates(sci[k],coords,output=scipy.float64,order=5,cval=magicnum,prefilter=False)

		# Only include good data (ie positive variance, wavelength
		#   greater than dichroic cutoff)
		cond = (bgcoords[0]+offsets[k]<0.)|(bgcoords[0]+offsets[k]>ylen)
		cond = (varimage[k]<=0)|cond
		cond = (bgcoords[1]<locutoff)|(bgcoords[1]>hicutoff)|cond
		varimage[k][cond] = scipy.nan

		bgimage[k] = scipy.ndimage.map_coordinates(sub,coords,output=scipy.float64,order=5,cval=magicnum,prefilter=False)
		bgimage[k][cond] = scipy.nan
		bgimage[k][bgimage[k]==magicnum] = scipy.nan # Shouldn't be
							     #   necessary...

	if RESAMPLE==1:
		return out,bgimage,varimage
	
	bgimage = fastmed(bgimage)
	varimage = fastmed(varimage)/nsci

	return out,bgimage,varimage
