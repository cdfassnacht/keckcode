from spectra import spectools
from esi import skysub
import special_functions

import pyfits,scipy,pickle
from scipy import optimize,interpolate,ndimage,signal,stats


# Define saturation level for arclines
SATURATED = 57000.
SCALE1 = 1000.
SCALE2 = 1.

LINES = [5460.735,5577.345,6300.32,6363.81,6533.04,6553.61,6912.62,6923.21,6939.52,7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,7821.51,7841.266,7993.332,8310.719,8399.16,8415.231,8430.17,8791.186,8885.83,8943.395,9337.854]


def myoptimize(p,x,z,scale,model,model2=None):
#	vals = special_functions.genfunc(x,0,p).astype(scipy.float64)
#	mod = interpolate.splev(vals,model)
#	if model2 is not None:
#		mod += interpolate.splev(model2)
#	ratio = scipy.median(z)/scipy.median(mod)
#	del mod

	par = p['coeff'].copy()
	scale = float(scale)
	for i in range(2,par.size):
		par[i] *= scale**i
	# The starting wavelength is scaled *down* by 100...
	par[0] /= SCALE1
	par[1] *= SCALE2
	p['coeff'] = par.copy()
	par = special_functions.unpack_coeff(p)

#	par = scipy.zeros(par0.size+1)
#	par[0] = ratio 
#	par[1:] = par0.copy()
	coeff,ier = optimize.leastsq(skyfitfunc,par,(x,z,scale,model,p,model2),maxfev=100000)
#	coeff = coeff[1:].copy()

	par = special_functions.build_coeff(coeff,p)
	for i in range(2,par['coeff'].size):
		par['coeff'][i] /= scale**i
	# The starting wavelength is scaled *down* by 100...
	par['coeff'][0] *= SCALE1
	par['coeff'][1] /= SCALE2

	return par

def skyfitfunc(p,x,data,scale,model,tmp,model2=None):
#	ratio = p[0]
#	offset = p[1]
#	p = p[1:].copy()

	par = special_functions.build_coeff(p,tmp)
	for i in range(2,par['coeff'].size):
		par['coeff'][i] /= scale**i
	# The starting wavelength is scaled *down* by 100...
	par['coeff'][0] *= SCALE1
	par['coeff'][1] /= SCALE2

	data = data.astype(scipy.float64)
	z = special_functions.genfunc(x,0,par).astype(scipy.float64)

	mod = interpolate.splev(z,model).astype(scipy.float64)
	if model2 is not None:
		mod += interpolate.splev(z,model2).astype(scipy.float64)

#	ratio = data.max()/(mod+scipy.median(data)).max()
#	mod = mod*ratio + scipy.median(data)
#	mod = mod*data.max()/mod.max()# + scipy.median(data)
#	mod = mod*scipy.median(data)/scipy.median(mod)
	ratio = scipy.median(data)/scipy.median(mod)
	mod = mod*ratio #+ offset

	diff = (data-mod)/scipy.sqrt(abs(mod))
	return diff

def arcfitfunc(p,x,y,z,model):
	par = special_functions.unpack_coeff(p)
	diag = rescale_pars(p,x.max(),y.max())
	coeff,ier = optimize.leastsq(doarcfitfunc,par,(x,y,z,model,p),maxfev=100000,diag=diag)
	return special_functions.build_coeff(coeff,p)

# Fitting function for the arcs
def doarcfitfunc(p,xdata,ydata,scidata,model,coeff):
	par = special_functions.build_coeff(p,coeff)
	scidata = scidata.astype(scipy.float64)
	z = special_functions.genfunc(xdata,ydata,par).astype(scipy.float64)
	z = z.reshape((1,scidata.size))

	resample = ndimage.map_coordinates(model,z,output=scipy.float64,cval=-1)

	# Do not fit to bad pixels (the left/right pixels on top/bottom slits)
	z = z[0].round().astype(scipy.int16)
	z[z>=model.size] = model.size-1
	z[z<0] = 0.
	cond = (model[z]>0)&(resample>0)
	return (scidata[cond] - resample[cond])/scipy.sqrt(abs(resample[cond]))

def rescale_pars(p,x,y):
	p = p.copy()
	p['coeff'] *= 0.
	p['coeff'] += 1.
	for i in range(p['coeff'].shape[0]):
		for j in range(p['coeff'].shape[1]):
			if i==1 and j==0:
				continue
			p['coeff'][i,j] /= scipy.power(x,i)
			p['coeff'][i,j] /= scipy.power(y,j)
	return special_functions.unpack_coeff(p).tolist()

def clipped_std(data,clip):
	d = data.copy()
	mean = d.mean()
	std = d.std()

	delta = 1
	while delta:
		size = d.size
		d = d[abs(d-mean)<clip*std]
		delta = size-d.size
		if delta==size:
			return 0.
		mean = d.mean()
		std = d.std()
	return std


def getpeaks(x,z):
	# Peak finding
	tmp = ndimage.maximum_filter(z,9)
	peaks = scipy.where(tmp==z,z,0)

	# Remove non-peaks
	tmp = ndimage.maximum_filter(z,3)
	tmp = scipy.where(tmp==z,1.,0.)
	peaks *= tmp

	# Only choose 10 sigma lines
	avg = ndimage.percentile_filter(z,35.,40)
	std = scipy.sqrt(avg)
	peaks = scipy.where(peaks>avg+5.*std)[0]

	# Centroid the lines, fixing the zeropoint of the gaussians to avg
	fit = scipy.ones(8)
	fit[4] = 0.
	datalines = []
	for i in peaks:
		if i-4<0 or i+5>x.size:
			continue
		fit[0] = avg[i]
		fit[1] = z[i]
		fit[2] = x[i]
		fit[3] = 1.

		fitdata = scipy.empty((9,2))
		fitdata[:,0] = x[i-4:i+5].copy()
		fitdata[:,1] = z[i-4:i+5].copy()

		gfit,chi2 = special_functions.ngaussfit(fitdata,fit)
		datalines.append(gfit[2])

	return scipy.asarray(datalines)

def matchlines(peaks,w,tol,order):
	lines = scipy.asarray(LINES)
	lines.sort()

	wave = special_functions.genfunc(peaks,0.,w)

	gooddata = []
	goodlines = []
	delts = []
	for i in range(wave.size):
		delta = 1e9
		match = None
		for j in lines:
			if abs(j-wave[i])<delta:
				delta = abs(j-wave[i])
				match = j
			else:
				break
		if delta<tol:
			delts.append(delta)
			gooddata.append(peaks[i])
			goodlines.append(match)

	fitdata = scipy.empty((len(gooddata),2))
	fitdata[:,0] = scipy.asarray(gooddata)
	fitdata[:,1] = scipy.asarray(goodlines)

	w = special_functions.lsqfit(fitdata,'chebyshev',order)
	lines = special_functions.genfunc(fitdata[:,0],0.,w)
	error = lines-fitdata[:,1]

	return error,w

# Chi-square of skymodel and data
def push(data,model):
	start = 0
	min = scipy.empty(model.size-data.size)
	div = abs(model)
	while start+data.size<model.size:
		end = start+data.size
		diff = data-model[start:end]
		chi = (diff*diff/div[start:end]).sum()
		min[start] = chi
		start += 1
	return min.min(),min.argmin()


#
# curve = first pix of non-ycor
# sci = biastrimmed, non-ycor science data array for slit
# arc = y-corrected arc for this slit
# yforw = yforw definition for slit
# arcfit = wavelength model
# disp = scale
# mswave = central wavelength
#

def skymatch(curve,sci,arc,yforw,widemodel,finemodel,goodmodel,disp,mswave,cutoff):
	# Do not change input
	sci = sci.copy()
	nsci = sci.shape[0]

	if cutoff==0:
		cutoff = 5000


	# Height and width of slit
	height = arc.shape[0]
	width = arc.shape[1]

	# Straighten science data
	scimodel = scipy.zeros((nsci,height,width))
	for k in range(nsci):
		scimodel[k] = spectools.resampley(sci[k],yforw,curve)

	# Arcline straightening fit
	if height>200:
		xdim = 2
		ydim = 2
	else:
		xdim = 2
		ydim = 2

	#
	# Find approximate dichroic edge
	#
	smooth = ndimage.percentile_filter(scimodel[nsci/2],45.,3)
	tmp = scipy.where(smooth<=0,scipy.nan,smooth)
	slice = stats.stats.nanmean(tmp,1)
	indx = slice.argsort()[slice.size/4]
	skymodel = tmp[indx]
	slice = stats.stats.nanmedian(tmp,0)
	skymodel[scipy.isnan(skymodel)] = slice[scipy.isnan(skymodel)]
	skymodel[scipy.isnan(skymodel)] = 0.
	buffer = scipy.where(skymodel!=0)[0].min()
	gooddata = scipy.trim_zeros(skymodel)
	slice = ndimage.percentile_filter(gooddata,25.,35)
	med = scipy.median(slice)
	indx = scipy.where((gooddata>med/2)&(slice>med/4))[0]
	first = buffer + indx.min()-30
	if first<0:
		first = 0
	if first>0:
		last = first + (7600-cutoff)/disp
		lastwave = 8000
		if last-first>width*2/3:
			last = first+width*2/3
			if last>width:
				last = width
			lastwave = 10400
	else:
		last = width*2/3
		lastwave = 10400


	# Skip the bottom and top of the slit!
	i = 2
	j = height-2

	coords = spectools.array_coords((j-i,width))
	coords[0] += 2.

	# Mask cosmic rays
	cond = abs(scimodel[nsci/2]-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc = scipy.where(cond,smooth,scimodel[nsci/2])
	smooth = ndimage.percentile_filter(newarc,45.,5)
	cond = abs(newarc-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc[cond] = smooth[cond]
	smooth = ndimage.percentile_filter(newarc,45.,5)
	cond = abs(newarc-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc[cond] = smooth[cond]

	# Choose a row without source flux as model row for straightening
	# First mask bad pixels (from grating smile correction, for example)
	newarc[newarc<=0.] = scipy.nan
	slice = stats.stats.nanmean(newarc,1)
	indx = slice.argsort()[slice.size/4]
	arcmid = newarc[indx]
	fitdata = newarc[i:j,10:-10].flatten()
	xdata = coords[1,:,10:-10].flatten()
	ydata = coords[0,:,10:-10].flatten()

	# A hack to remove rows with EXTREMELY strong sources
	badrow = scipy.where(slice>slice[indx]+100.*scipy.sqrt(slice[indx]))[0]
	ycond = ydata==ydata
	for row in badrow:
		ycond = ycond&(ydata!=row)
	cond = (scipy.isfinite(fitdata))&(xdata>first)&ycond
	xdata = xdata[cond]
	ydata = ydata[cond]
	fitdata = fitdata[cond]
	arcmid[scipy.isnan(arcmid)] = 0.

	# Perform the fit
	p = scipy.zeros((xdim+1,ydim+1))
	p[1,0] = 1.
	p = {'coeff':p,'type':"chebyshev"}
	xback = arcfitfunc(p,xdata,ydata,fitdata,arcmid)

	# Do the reverse fit
	xdata = coords[1].flatten()
	ydata = coords[0].flatten()
	yorig = yforw[i:j].flatten()-curve
	newxdata = special_functions.genfunc(xdata,ydata,xback)

	tmp = scipy.zeros((xdata.size,3))
	tmp[:,0] = newxdata.copy()
	tmp[:,1] = ydata.copy()
	tmp[:,2] = xdata.copy()
	xforw = special_functions.lsqfit(tmp,"chebyshev",xdim,ydim)

	# Produce a model of the sky to use for the wavelength calibration
	bgmodel = scipy.zeros((nsci,sci.shape[2]))
	for k in range(nsci):
		tmp = spectools.resample1d(scimodel[k],xforw,"x",-999)
		scimodel[k] = tmp.copy()
		tmp.sort(axis=0)
		tmp[tmp==-999] = scipy.nan
		bg = stats.stats.nanmedian(tmp,axis=0)
		bg = scipy.tile(bg,(height,1))
		tmp[scipy.isnan(tmp)] = bg[scipy.isnan(tmp)]
		tmp[scipy.isnan(tmp)] = 0.
		bgmodel[k] = tmp[tmp.shape[0]/4,:]
		del tmp,bg
	#skymodel = stats.stats.median(bgmodel,0)
	skymodel = bgmodel[nsci/2]

	#
	# Obtain the wavelength solution
	#
	max = 1.0e19
	minwave = mswave-disp*skymodel.size*1.1
	maxwave = mswave+disp*skymodel.size*1.1
	if minwave<3400:
		minwave = 3400
	if maxwave>lastwave:
		maxwave = lastwave


	# Due a quick fit to define initial guesses...first the starting
	#   wavelength
	fudge = disp/800.
	nsteps = 15

	xmodel = skymodel[first:last]
	#xmodel = ndimage.gaussian_filter1d(xmodel,5.)
	p0 = 0.
	p1 = disp
	for l in range(nsteps):
		try_disp = disp + ((l-nsteps/2)*fudge)
		skyfit_x = scipy.arange(minwave,maxwave,try_disp)
		fitmodel = interpolate.splev(skyfit_x,finemodel)
		tratio = scipy.median(xmodel)/scipy.median(fitmodel)
		fitmodel *= tratio
		chi2,off = push(xmodel,fitmodel)
		if chi2<max:
			p1 = try_disp
			p0 = (off-first)*try_disp+minwave
			max = chi2
	scale = skymodel.size
	xdim = 3
	ydim = 2
	p = scipy.zeros((xdim+1,1))
	p[0,0] = p0
	p[1,0] = p1
	sky2x = []
	sky2y = []
	ccd2wave = []
	# Solve for the wavelength solution of each mask.
	for k in range(nsci):
		# Do pixel-by-pixel fitting for first order solution
		fitx = scipy.arange(first,skymodel.size,1.)
		fitdata = bgmodel[k][first:]

		peaks = getpeaks(fitx,fitdata)

		fitdata1 = ndimage.gaussian_filter1d(fitdata,5./disp)
		start = p.copy()
		skycoeff = {'coeff':start,'type':"chebyshev"}
		skycoeff = myoptimize(skycoeff,fitx,fitdata1,scale,widemodel)
		error,skycoeff = matchlines(peaks,skycoeff,3.*disp,3)
		skycoeff = myoptimize(skycoeff,fitx,fitdata,scale,finemodel)
		error,skycoeff = matchlines(peaks,skycoeff,2.*disp,3)
		error,skycoeff = matchlines(peaks,skycoeff,disp,3)
		wlen = special_functions.genfunc(newxdata,0,skycoeff)

		print "  Image %d wavelength error: %5.3f angstroms from %d lines" % (k+1,error.std(),error.size)

		xdim = 5
		ydim = 2

		revmodel = scipy.zeros((wlen.size,3))
		revmodel[:,0] = wlen.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = xdata.copy()
		if k==0:
			sky2x.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		else:
			sky2x.append(special_functions.lsqfitter(revmodel,sky2x[k-1]))

		revmodel[:,2] = yorig.copy()
		if k==0:
			sky2y.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		else:
			sky2y.append(special_functions.lsqfitter(revmodel,sky2y[k-1]))

		revmodel[:,0] = xdata.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = wlen.copy()
		if k==0:
			ccd2wave.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		else:
			ccd2wave.append(special_functions.lsqfitter(revmodel,ccd2wave[k-1]))

	return sky2x,sky2y,ccd2wave
