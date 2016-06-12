"""
Helper functions to remove telluric absorption based on a model from Keck.
"""

import scipy
from scipy import io as sio,ndimage,interpolate
import mostools

def correct(inwave,airmass=1.,scale=0.85):
	"""
	correct(inwave,airmass=1.,scale=0.85)

	Computes telluric correction for the A-band and B-band.

	Inputs:
	  inwave  - wavelengths for which corrections should be determined
	  airmass - airmass of spectrum
	  scale   - approximate resolution (in sigma) of science data

	Outputs:
	  *multiplicative* telluric correction for A-band and B-band
	"""
	a = aband(inwave,airmass,scale)
	b = bband(inwave,airmass,scale)

	return a+b


def get_correction(inwave,airmass,scale,wave,data):
	"""
	get_correction(inwave,airmass,scale,wave,data)

	Determines a telluric correction from a model.

	Inputs:
	  inwave  - wavelengths for which corrections should be determined
	  airmass - airmass of science data
	  scale   - approximate resolution (in sigma) of science data
	  wave    - wavelengths of telluric model
	  data    - telluric model

	Outputs:
	  *multiplicative* telluric correction
	"""
	data = data**(airmass**0.55)

	if scale>0.85:
		kernel = scipy.sqrt(scale**2-0.85**2)
		data = ndimage.gaussian_filter1d(data,kernel)
	spline = interpolate.splrep(wave,data,s=0)

	cond = (inwave>wave[0])&(inwave<wave[-1])
	good = inwave[cond]

	correction = interpolate.splev(good,spline)
	output = scipy.ones(inwave.size)

	output[cond] = correction
	return output


def bband(inwave,airmass=1.,scale=0.85):
	"""
	Telluric correction for the B-band
	"""
	path = mostools.__path__[0]
	file = path+"/data/bband.dat"

	bband = sio.read_array(file)

	wave = scipy.power(10.,bband[:,0])
	data = bband[:,1].astype(scipy.float32)

	return get_correction(inwave,airmass,scale,wave,data)


def aband(inwave,airmass=1.,scale=0.85):
	"""
	Telluric correction for the A-band
	"""
	path = mostools.__path__[0]
	file = path+"/data/aband.dat"

	aband = sio.read_array(file)

	wave = aband[:,0]
	data = aband[:,1].astype(scipy.float32)

	return get_correction(inwave,airmass,scale,wave,data)
