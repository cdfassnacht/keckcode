"""
Module to determine the wavelength solution of the lris redside using skylines.

TODO: add logging
"""


from mostools import spectools
from lris.lris_red import skysub
import special_functions

import pyfits,scipy
from scipy import optimize,interpolate,ndimage,signal,stats


""" Define saturation level for arclines """
SATURATED = 57000.
SCALE1 = 1000.
SCALE2 = 1.


""" List of skylines to use """
LINES = [5460.735,5577.345,5915.308,5932.864,6257.970,6300.32,6363.81,6533.04,6553.61,6863.971,6871.073,6912.62,6923.21,6939.52,7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,7821.51,7841.266,7993.332,8310.719,8399.16,8415.231,8430.17,8791.186,8885.83,8943.395,8988.384,9038.059,9337.854,9375.977,9419.746,9439.67,9458.524]


"""
Wrapper to optimize.leastsq(). Rescales the parameters to make all of them
  approximately the same order -- this might not be necessary any more.
"""
def myoptimize(p,x,z,scale,model,model2=None):
	par = p['coeff'].copy()
	scale = float(scale)
	for i in range(2,par.size):
		par[i] *= scale**i

	par[0] /= SCALE1
	par[1] *= SCALE2
	p['coeff'] = par.copy()
	par = special_functions.unpack_coeff(p)

	coeff,ier = optimize.leastsq(skyfitfunc,par,(x,z,scale,model,p,model2),maxfev=100000)

	par = special_functions.build_coeff(coeff,p)
	for i in range(2,par['coeff'].size):
		par['coeff'][i] /= scale**i
	par['coeff'][0] *= SCALE1
	par['coeff'][1] /= SCALE2

	return par

"""
The model of the sky/arc spectrum is evaluated at the wavelengths of the data
  pixels, and the difference between the model and the data is returned to the
  fitting routine.
"""
def skyfitfunc(p,x,data,scale,model,tmp,model2=None):
	par = special_functions.build_coeff(p,tmp)
	for i in range(2,par['coeff'].size):
		par['coeff'][i] /= scale**i
	par['coeff'][0] *= SCALE1
	par['coeff'][1] /= SCALE2

	data = data.astype(scipy.float64)
	z = special_functions.genfunc(x,0,par).astype(scipy.float64)

	mod = interpolate.splev(z,model).astype(scipy.float64)
	if model2 is not None:
		mod += interpolate.splev(z,model2).astype(scipy.float64)

	ratio = scipy.median(data)/scipy.median(mod)
	mod = mod*ratio

	diff = (data-mod)/scipy.sqrt(abs(mod))
	return diff


"""
Functions for fitting the arclines for line straigtening.
"""
def arcfitfunc(p,x,y,z,model):
	par = special_functions.unpack_coeff(p)
	diag = rescale_pars(p,x.max(),y.max())
	coeff,ier = optimize.leastsq(doarcfitfunc,par,(x,y,z,model,p),maxfev=100000,diag=diag)
	return special_functions.build_coeff(coeff,p)

def doarcfitfunc(p,xdata,ydata,scidata,model,coeff):
	par = special_functions.build_coeff(p,coeff)
	scidata = scidata.astype(scipy.float64)
	z = special_functions.genfunc(xdata,ydata,par).astype(scipy.float64)
	z = z.reshape((1,scidata.size))

	resample = ndimage.map_coordinates(model,z,output=scipy.float64,cval=-1)

	diff = (scidata - resample)/scipy.sqrt(abs(resample))
	return diff

"""
Rescales parameters to approx. the same order. Probably not necessary....
"""
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


"""
Calculates a clipped std.
"""
def clipped_std(data,clip):
	d = data.copy()

	while 1:
		mean,std = d.mean(),d.std()
		size = d.size
		d = d[abs(d-mean)<clip*std]
		if size==d.size or d.size==0:
			return mean,std
	return mean,std


"""
Finds and centroids peaks in the spectrum.
"""
def findlines(x,z,nsigma=5.):
	""" Start by estimating the sky continuum background. """
	min = ndimage.minimum_filter(z,9)
	minx = scipy.where(min==z)[0]
	min = z[minx]
	fit = scipy.empty((minx.size,2))
	fit[:,0] = minx.astype(scipy.float32)
	fit[:,1] = min.copy()
	bgfit = special_functions.lsqfit(fit,'chebyshev',7)

	""" Trim outliers (usually from blended lines) """
	while 1:
		bg = special_functions.genfunc(minx,0,bgfit)
		oldsize = minx.size
		std = stats.stats.std(min-bg)
		cond = ((min-bg)<3.*std)&((bg-min)<5.*std)
		min = min[cond]
		minx = minx[cond]
		if minx.size==oldsize or minx.size==0:
			break
		fit = scipy.empty((minx.size,2))
		fit[:,0] = minx.astype(scipy.float32)
		fit[:,1] = min.copy()
		bgfit = special_functions.lsqfit(fit,'chebyshev',7)
	bg = special_functions.genfunc(scipy.arange(z.size),0,bgfit)

	"""
	Assuming poisson statistics, a line must be 'nsigma' times the noise
	  level above the background.
	"""
	tmp = ndimage.maximum_filter(z,9)
	threshold = bg+nsigma*(bg**0.5)
	peaks = scipy.where((tmp==z)&(z>threshold))[0]


	""" Centroid the lines, fixing the bias level of the gaussians """
	fit = scipy.ones(8)
	fit[4] = 0.
	datalines = []
	for i in peaks:
		if i-14<0 or i+15>x.size:
			continue
		fit[0] = bg[i]
		fit[1] = z[i]
		fit[2] = x[i]
		fit[3] = 1.

		"""
		Deal with saturated lines appropriately (ie skipping the
		  non-linear pixels)
		"""
		if z[i]<SATURATED:
			fitdata = scipy.empty((9,2))
			fitdata[:,0] = x[i-4:i+5].copy()
			fitdata[:,1] = z[i-4:i+5].copy()
		else:
			fitdata = scipy.empty((25,2))
			fitdata[:,0] = x[i-12:i+13].copy()
			fitdata[:,1] = z[i-12:i+13].copy()
			fitdata = fitdata[fitdata[:,1]<SATURATED*0.8]

		""" A simple weighted sum would probably be robust.... """
		gfit,chi2 = special_functions.ngaussfit(fitdata,fit)
		datalines.append(gfit[2])

	return scipy.asarray(datalines)



def matchlines(peaks,w,tol,order,clean=False):
	"""
	Line-matching routine to associate features in the data with known
	  features.
	"""
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
	diff = lines-fitdata[:,1]

	"""
	Clip the outliers from the fit. Usually the outliers are 'real' matches
	  but do not have good centroids or are blended. Cleaning is probably
	  best utilized at the final iteration. One negative side effect of
	  cleaning is that isolated lines (ie 5577) may get removed and the
	  solution becomes unconstrained in that region.
	"""
	while clean:
		size = diff.size
		cond = abs(diff)<2.*stats.stats.std(diff)
		data = fitdata[:,0][cond]
		if data.size==diff.size:
			break
		lines = fitdata[:,1][cond]
		fitdata = scipy.empty((data.size,2))
		fitdata[:,0] = data.copy()
		fitdata[:,1] = lines.copy()
		w = special_functions.lsqfit(fitdata,'chebyshev',order)
		lines = special_functions.genfunc(fitdata[:,0],0.,w)
		diff = lines-fitdata[:,1]

	return diff,w


"""
Simple chi-square routine for calculation of initial parameters of a linear
model. This can be used instead of a correlation.
"""
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


"""
Debugging function for viewing wavelength solutions.
"""
def debug(wave,fitdata,finemodel,skycoeff=None,offset=0.):
        if skycoeff is not None:
                wave = special_functions.genfunc(wave,0.,skycoeff)
        mod = interpolate.splev(wave,finemodel)
        import pylab
        pylab.plot(wave,fitdata)
        pylab.plot(wave+offset,mod*fitdata.max()/mod.max())
        pylab.show()


"""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
"""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""


"""
skymatch() is the main function for the wavelength solution.

 curve     - first pix of non-ycor slit (ie bottom of curved slit)
 sci       - biastrimmed, non-ycor science data array for slit
 arc       - y-corrected arc for this slit (not used in current implementation)
 yforw     - the y-forward model for this slit
 widemodel - broad wavelength model
 finemodel - matched-resolution wavelength model
 goodmodel - ??
 disp      - output pixel scale
 mswave    - central wavelength
 cutoff    - bluest allowed wavelength
"""
def skymatch(curve,sci,arc,yforw,widemodel,finemodel,goodmodel,disp,mswave,cutoff):
	""" Do not alter input data! """
	sci = sci.copy()
	nsci = sci.shape[0]

	if cutoff==0:
		cutoff = 5000

	height = arc.shape[0]
	width = arc.shape[1]

	
	""" Straighten science data """
	scimodel = scipy.zeros((nsci,height,width))
	for k in range(nsci):
		scimodel[k] = spectools.resampley(sci[k],yforw,curve)


	""" Fit order may need to be changed for very long slits """
	xord = 2
	yord = 2


	"""
	Find approximate dichroic edge. A lot of this is tweaked for LRIS, but
	  it is also somewhat general in that it tries to use the area of
	  dense sky features *without using much of the repeating self-similar
	  patterns of OH lines* for the correlation. The repeating 'forests'
	  can lead to 'jumps' in the wavelength solution.
	"""
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
	if cutoff>6000:
		last = first+width*2/3
		if last>width:
			last = width
		lastwave = 10400


	""" Skip the bottom and top of the slit, which may have artifacts """
	i = 2
	j = height-2

	coords = spectools.array_coords((j-i,width))
	coords[0] += 2.

	
	""" Cosmic-ray removal """
	cond = abs(scimodel[nsci/2]-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc = scipy.where(cond,smooth,scimodel[nsci/2])
	smooth = ndimage.percentile_filter(newarc,45.,5)
	cond = abs(newarc-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc[cond] = smooth[cond]
	smooth = ndimage.percentile_filter(newarc,45.,5)
	cond = abs(newarc-smooth)>3.*scipy.sqrt(abs(smooth))
	newarc[cond] = smooth[cond]

	"""
	Choose a row without source flux as model row for straightening. First
	  mask bad pixels (from grating smile correction, for example).
	"""
	newarc[newarc<=0.] = scipy.nan
	slice = stats.stats.nanmean(newarc,1)
	indx = slice.argsort()[slice.size/4]
	arcmid = newarc[indx]
	fitdata = newarc[i:j,height:-1*height].flatten()
	xdata = coords[1,:,height:-1*height].flatten()
	ydata = coords[0,:,height:-1*height].flatten()

	""" A hack to remove rows with EXTREMELY strong sources. """
	badrow = scipy.where(slice>slice[indx]+100.*scipy.sqrt(slice[indx]))[0]
	ycond = ydata==ydata
	for row in badrow:
		ycond = ycond&(ydata!=row)
	cond = (scipy.isfinite(fitdata))&(xdata>first)&ycond
	xdata = xdata[cond]
	ydata = ydata[cond]
	fitdata = fitdata[cond]
	arcmid[scipy.isnan(arcmid)] = 0.

	""" Perform the line-straightening fit. """
	p = scipy.zeros((xord+1,yord+1))
	p[1,0] = 1.
	p = {'coeff':p,'type':"chebyshev"}
	xback = arcfitfunc(p,xdata,ydata,fitdata,arcmid)

	""" The inverse fit. """
	xdata = coords[1].flatten()
	ydata = coords[0].flatten()
	yorig = yforw[i:j].flatten()-curve
	newxdata = special_functions.genfunc(xdata,ydata,xback)
	tmp = scipy.zeros((xdata.size,3))
	tmp[:,0] = newxdata.copy()
	tmp[:,1] = ydata.copy()
	tmp[:,2] = xdata.copy()
	xforw = special_functions.lsqfit(tmp,"chebyshev",xord,yord)

	"""
	Produce a model of the sky to use for the initial wavelength calibration
	"""
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
	skymodel = bgmodel[nsci/2]  ## Use the middle exposure for the model


	"""
	******************************
	Obtain the wavelength solution
	******************************

	The skylines on the red side tend to be much more stable than using
	  the arcs on the blue side, so the code is somewhat simpler.
	"""

	"""
	We perform a correlation for the initial wavelength solution. The
	  correlation can also be replaced with a chi-square routine to solve
	  for the best (linear) offset. We loop through a range of guesses at
	  the dispersion to find the best. We also only use the part of the
	  sky from first:last to get an initial solution that is optimized for
	  the blue end of the sky spectrum, where there are fewer sky lines.
	"""
	max = 1.0e19
	minwave = mswave-disp*skymodel.size*1.1
	maxwave = mswave+disp*skymodel.size*1.1
	if minwave<3400:
		minwave = 3400
	if maxwave>lastwave:
		maxwave = lastwave


	fudge = disp/800.
	nsteps = 15

	xmodel = skymodel[first:last]  ## Or a smoothed model can be used
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


	p = scipy.zeros((xord+1,1))
	p[0,0] = p0
	p[1,0] = p1
	sky2x = []	# Convert wavelength,Y in the output image to input X
	sky2y = []	# Convert wavelength,Y in the output image to input Y
	ccd2wave = []	# Convert input X,Y to wavelength

	""" Solve for the wavelength solution of each mask independently. """
	for k in range(nsci):
		"""
		We first do pixel-by-pixel fitting for the initial solution
		  then do a line-matching when the solution is pretty good.
		"""
		fitx = scipy.arange(first,skymodel.size,1.)
		fitdata = bgmodel[k][first:]

		peaks = findlines(fitx,fitdata,5.)

		fitdata1 = ndimage.gaussian_filter1d(fitdata,5./disp)
		start = p.copy()
		skycoeff = {'coeff':start,'type':"chebyshev"}
		skycoeff = myoptimize(skycoeff,fitx,fitdata1,scale,widemodel)
		error,skycoeff = matchlines(peaks,skycoeff,3.*disp,3)
		skycoeff = myoptimize(skycoeff,fitx,fitdata,scale,finemodel)
		error,skycoeff = matchlines(peaks,skycoeff,2.*disp,3)

		"""
		Set the 'clean' parameter and determine the final wavelength
		  solution.
		"""
		error,skycoeff = matchlines(peaks,skycoeff,disp,3,True)
		wlen = special_functions.genfunc(newxdata,0,skycoeff)

		print "  Image %d wavelength error: %5.3f angstroms from %d lines" % (k+1,stats.stats.std(error),error.size)


		"""
		Create the full 2d solution. The 'if' statements (should) speed
		  up the fitting by using the previous mask's solution as the
		  starting guess of the solution for the current mask.
		"""
		xord = 5
		yord = 2

		revmodel = scipy.zeros((wlen.size,3))
		revmodel[:,0] = wlen.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = xdata.copy()
		if k==0:
			sky2x.append(special_functions.lsqfit(revmodel,"chebyshev",xord,yord))
		else:
			sky2x.append(special_functions.lsqfitter(revmodel,sky2x[k-1]))

		revmodel[:,2] = yorig.copy()
		if k==0:
			sky2y.append(special_functions.lsqfit(revmodel,"chebyshev",xord,yord))
		else:
			sky2y.append(special_functions.lsqfitter(revmodel,sky2y[k-1]))

		revmodel[:,0] = xdata.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = wlen.copy()
		if k==0:
			ccd2wave.append(special_functions.lsqfit(revmodel,"chebyshev",xord,yord))
		else:
			ccd2wave.append(special_functions.lsqfitter(revmodel,ccd2wave[k-1]))

	return sky2x,sky2y,ccd2wave
