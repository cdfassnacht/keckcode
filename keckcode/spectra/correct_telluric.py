import scipy
from scipy import io as sio,ndimage,interpolate
import spectra

def get_correction(inwave,airmass,scale,wave,data):
	data = data**(airmass**0.55)

	if scale>0.85:
		kernel = scipy.sqrt(scale**2-0.85**2)
		data = ndimage.gaussian_filter1d(data,kernel)
	spline = interpolate.splrep(wave,data)

	cond = (inwave>wave[0])&(inwave<wave[-1])
	good = inwave[cond]

	correction = interpolate.splev(good,spline)
	output = scipy.ones(inwave.size)

	output[cond] = correction
	return output


def correct(inwave,airmass=1.,scale=0.85):
	a = aband(inwave,airmass,scale)
	b = bband(inwave,airmass,scale)

	return a+b


def bband(inwave,airmass=1.,scale=0.85):
	path = spectra.__path__[0]
	file = path+"/data/bband.dat"

	bband = sio.read_array(file)

	wave = scipy.power(10.,bband[:,0])
	data = bband[:,1].astype(scipy.float32)

	return get_correction(inwave,airmass,scale,wave,data)


def aband(inwave,airmass=1.,scale=0.85):
	path = spectra.__path__[0]
	file = path+"/data/aband.dat"

	aband = sio.read_array(file)

	wave = aband[:,0]
	data = aband[:,1].astype(scipy.float32)

	return get_correction(inwave,airmass,scale,wave,data)
