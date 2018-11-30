import scipy
try:
	import pyfits
except:
	from astropy.io import fits as pyfits

# resample1d(data,coeffs,axis,const)
# array_coords(shape)
# resampley(data,ycoords,yoffset,cval,mode)
# get_slit(2ddata,slit#)
# logrebin(spec,crval,crpix,cd1)
# fits_logrebin(in,out)



# Resample an image along the axis specified
def resample1d(data,coeffs,axis,const=0.0,offset=0.):
	from special_functions import genfunc
	coords = array_coords(data.shape)
	y = coords[0].flatten()
	x = coords[1].flatten()
	if axis=="Y" or axis=="y" or axis==1:
		y += offset
		coords[0] = genfunc(x,y,coeffs).reshape(data.shape)-offset
	else:
		x += offset
		coords[1] = genfunc(x,y,coeffs).reshape(data.shape)-offset
	del x,y
	return scipy.ndimage.map_coordinates(data,coords,cval=const,output=scipy.float64,order=5)

def array_coords(shape):
	y = shape[0]
	x = shape[1]
	out = scipy.empty((2,y,x))
	t = scipy.arange(y,dtype='f8')
	out[0] = scipy.tile(t,(x,1)).T
	t = scipy.arange(x,dtype='f8')
	out[1] = scipy.tile(t,(y,1))
	return out

def resampley(data,ycoords,yoffset=0.,cval=0.,mode="constant",slice=None):
	from scipy import ndimage
	try:
		coords = array_coords(ycoords.shape)
	except:
		from special_functions import genfunc
		coords = array_coords(data.shape)
		x = coords[1].flatten()
		y = coords[0].flatten()+yoffset
		ycoords = genfunc(x,y,ycoords).reshape(coords[0].shape)
		if slice is not None:
			ycoords = ycoords[slice].copy()
		coords = array_coords(ycoords.shape)
		del x,y
	coords[0] = ycoords.astype(scipy.float64)-yoffset
	return ndimage.map_coordinates(data,coords,output=scipy.float64,mode=mode,cval=cval,order=5)

def cutout(file,outname,slit_number,plane=0):
	f = pyfits.open(file)
	data = f[0].data.copy()
	hdr = f[0].header.copy()
	data[scipy.isnan(data)] = 0.

	if len(data.shape)==3:
		data = data[plane]
	slit,start,bottom = cut_slit(data,slit_number)
	crval = hdr['CRPIX1']-start

	hdu = pyfits.PrimaryHDU(slit)
	hdu.header = hdr
	hdu.header.update('CRPIX1',crval)
	hdu.verify('fix')

	hdu.writeto(outname)
	

def get_slit(data,slit_number):
	a,b,c = cut_slit(data,slit_number)
	return a

def cut_slit(data,slit_number):
	slit_number -= 1

	data[scipy.isnan(data)] = 0.

	slice = data.sum(axis=1)

	# Find bottom of good data
	bottom = 0
	while slice[bottom]==0:
		bottom += 1

	slice = scipy.trim_zeros(slice)

	zeros = scipy.where(slice==0)
	zeros = zeros[0]+bottom

	# Special case for images with only one slit
	if zeros.size==0:
		# Deal with upper zeros
		slice = data.sum(axis=1)
		top = bottom
		while top<slice.size and slice[top]!=0:
			top += 1

		slice = data.sum(axis=0)
		indx = scipy.where(slice!=0)
		start = indx[0][0]
		end = indx[0][-1]+1
		return data[bottom:top,start:end],start,bottom

	borders = scipy.array_split(zeros,zeros.size/5)
	if slit_number>len(borders):
		return scipy.asarray(0),0,0

	if slit_number==0:
		start = 0
		end = borders[0][0]
	elif slit_number==len(borders):
		start = borders[slit_number-1][4]+1
		end = slice.size
	else:
		start = borders[slit_number-1][4]+1
		end = borders[slit_number][0]
	data = data[start:end]
	bottom = start
	slice = data.sum(axis=0)
	indx = scipy.where(slice!=0)
	start = indx[0][0]
	end = indx[0][-1]+1
	return data[:,start:end],start,bottom


def logrebin(spectrum,crval,crpix,cd1):
	from scipy import interpolate
	# First determine the "best" pixel scale
	start = crval+cd1*(1-crpix)
	inwave = scipy.arange(start,start+spectrum.size*cd1,cd1)
	sampwave = scipy.log10(inwave)
	pixscale = (sampwave[-1]-sampwave[0])/sampwave.size
	outwave = scipy.arange(sampwave[0],sampwave[-1],pixscale)
	outwave = scipy.power(10,outwave)

	spline = interpolate.splrep(inwave,spectrum)
	newspec = interpolate.splev(outwave,spline)
	return outwave,newspec

def fits_logrebin(infile,outfile):
	try:
		import pyfits
	except:
		from astropy.io import fits as pyfits
	f = pyfits.open(infile)
	header = f[0].header
	data = f[0].data.copy()
	crval1 = header['crval1']
	crpix = header['crpix1']
	cd = header['cd1_1']
	outwave,outspec = rebin_log(data,crval1,crpix,cd)
	out = pyfits.PrimaryHDU(outspec)
	outwv = scipy.log10(outwave)
	start = outwv[0]
	delt = outwv[1]-outwv[0]
	out.header.update('CRVAL1',start)
	out.header.update('CRPIX1',1)
	out.header.update('CD1_1',delt)
	out.header.update('DC-FLAG',1)
	out.writeto(outfile)

def optimal_smooth(data,var,width=1.5):
	from scipy import signal
	window = 8*width # Go out to 4sigma in the tails
	if window%2==0:
		window += 1
	w = signal.gaussian(window,width)
	tmp = data/var
	smooth = signal.convolve(tmp,w,'same')
	norm = signal.convolve(1./var,w,'same')
	return smooth/norm

def parse_hdr(header):
	crval = header['crval1']
	try:
		crpix = header['crpix1']
	except:
		crpix = 1
	log = 0
	try:
		cd = header['cd1_1']
	except:
		cd = header['cdelt1']
	try:
		log = header['dc-flag']
	except:
		try:
			tmp = header['WFITTYPE']
			if tmp=='LOG_LINEAR':
				log = 1
		except:
			pass 
	return [crval,crpix,cd,log]

def wavelength(filename,ext=0):
	f = pyfits.open(filename)
	hdr = f[ext].header

	hdr_info = parse_hdr(hdr)
	crval = hdr_info[0]
	crpix = hdr_info[1]
	cd = hdr_info[2]
	islog = hdr_info[3]
	npix = hdr['NAXIS1']

	start = crval+(1.-crpix)*cd
	wave = scipy.arange(start,start+npix*cd,cd)
	if islog:
		wave = scipy.power(10.,wave)
	return wave[0:npix]
