from mostools import spectools
import special_functions

import scipy
from scipy import optimize,interpolate,ndimage,signal,stats,random



""" Define saturation level for arclines """
SATURATED = 57000.



"""
Fitting functions for the wavelength solution. One thing that has been looked
  at only briefly is parametrizing a gain offset between the model and the data
  so that the chi-square was something like:
    (data - gain*model)
"""

def arcfitfunc2(p,x,y,z,model):
	par = special_functions.unpack_coeff(p)
	coeff,ier = optimize.leastsq(doarcfitfunc,par,(x,y,z,model,p),maxfev=100000)
	return special_functions.build_coeff(coeff,p)

# Fitting function for the arcs
def doarcfitfunc(p,xdata,ydata,scidata,model,coeff):
	par = special_functions.build_coeff(p,coeff)
	scidata = scidata.astype(scipy.float64)
	z = special_functions.genfunc(xdata,ydata,par).astype(scipy.float64)
	z = z.reshape((1,scidata.size))
	resample = ndimage.map_coordinates(model,z,output=scipy.float64,cval=-1)
	diff = (scidata - resample)/scipy.sqrt(abs(resample))
	return diff

"""
A wrapper to the optimize.leastsq function.
"""
def myoptimize(p,x,z,scale,model,model2=None):
	par = special_functions.unpack_coeff(p)
	coeff,ier = optimize.leastsq(arcfitfunc,par,(x,z,scale,model,p,model2),maxfev=100000)
	par = special_functions.build_coeff(coeff,p)
	return par

"""
The chi-square function for myoptimize. Compares a model of the arcs with
the arc data.
"""
def arcfitfunc(p,x,data,scale,model,tmp,model2=None):
	par = special_functions.build_coeff(p,tmp)


	""" Test for increasing function... """
	tmp = par['coeff'][1:].copy()
	for i in range(tmp.size):
		tmp[i] *= (i+1.)
	tmp = {'coeff':tmp.copy(),'type':par['type']}
	t = special_functions.genfunc(x,0,par).astype(scipy.float32)
	if t[t<0].size>0:
		return scipy.ones(x.size)*1.e9

	data = data.astype(scipy.float64)
	z = special_functions.genfunc(x,0,par).astype(scipy.float64)

	mod = interpolate.splev(z,model).astype(scipy.float64)
	if model2 is not None:
		mod += interpolate.splev(z,model2).astype(scipy.float64)

	mod = mod*data.max()/mod.max() + scipy.median(data)
	diff = (data-mod)/scipy.sqrt(abs(mod))

	return diff


"""
Other helper functions.
"""
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
Calculates a clipped std.
"""
def clipped_std(data,clip):
	d = data.copy()
	mean,std = d.mean(),d.std()

	while 1:
		size = d.size
		d = d[abs(d-mean)<clip*std]
		if size-d.size==0 or d.size==0:
			return mean,std
		mean,std = d.mean(),d.std()
	return mean,std


"""
Finds and centroids peaks in the spectrum.
"""
def findlines(x,z,sigma=5.):
	""" Peak finding """
	tmp = ndimage.maximum_filter(z,9)
	peaks = scipy.where(tmp==z,z,0)

	""" Only choose >sigma lines """
	avg = ndimage.percentile_filter(z,30.,55)
	tmp = scipy.where(peaks>0.)[0]
	peaks = []
	for i in tmp:
		start = i-150
		end = i+150
		if start<0:
			start = 0
		if end>z.size:
			end = z.size
		mean,std = clipped_std(z[start:end],3.5)
		if z[i]>avg[i]+sigma*std:
			peaks.append(i)

	""" Centroid the lines, fixing the zeropoint of the gaussians to avg """
	fit = scipy.ones(8)
	fit[4] = 0.
	datalines = []
	for i in peaks:
		if i-14<0 or i+15>x.size:
			continue
		fit[0] = avg[i]
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

"""
Linefile should contain the arclines expected to be present in the spectrum.
"""
def getlines(linefile):
	from scipy import io as sio
	file = open(linefile)

	lines = sio.read_array(file)
	file.close()
	lines = lines[:,0]
	lines.sort()

	return lines


def matchlines(peaks,solution,linefile,tol=30.,order=3):
	wave = special_functions.genfunc(peaks,0.,solution)
	lines = getlines(linefile)

	gooddata = []
	goodlines = []
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
			gooddata.append(peaks[i])
			goodlines.append(match)

	fitdata = scipy.empty((len(gooddata),2))
	fitdata[:,0] = scipy.asarray(gooddata)
	fitdata[:,1] = scipy.asarray(goodlines)

	return special_functions.lsqfit(fitdata,'chebyshev',order)

"""
Debugging function for viewing wavelength solutions.
"""
def debug(wave,fitdata,finemodel,skycoeff=None):
	if skycoeff is not None:
		wave = special_functions.genfunc(wave,0.,skycoeff)
	mod = interpolate.splev(wave,finemodel)
	import pylab
	pylab.plot(wave,fitdata)
	pylab.plot(wave,mod*fitdata.max()/mod.max())
	pylab.show()

"""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
"""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""




"""
Arcmatch is the main function of the arcmatch script. 

 curve = first pixel of non-ycor slit (ie bottom of curved slit)
 sci = biastrimmed, non-ycor science data array for slit
 arc = y-corrected arc for this slit
 yforw = yforw definition for slit
 widemodel = broad wavelength model
 finemodel = matched-resolution wavelength model
 skymodel = model of the sky
 disp = scale
 mswave = central wavelength
 extra = includes name of linefile, starting/ending wavelengths
"""
def arcmatch(curve,sci,arc,yforw,widemodel,finemodel,goodmodel,linemodel,disp,mswave,extra,logfile):
	""" Do not alter input data! """
	sci = sci.copy()
	nsci = sci.shape[0]

	""" Parse extra arguments """
	linefile = extra[0]
	cutoff = extra[1]
	if len(extra)==3:
		blue = extra[2]
	else:
		blue = 3000.

	height = arc.shape[0]
	width = arc.shape[1]

	""" Straighten science data """
	scimodel = scipy.zeros((nsci,height,width))
	for k in range(nsci):
		scimodel[k] = spectools.resampley(sci[k],yforw,curve)

	""" Skip the bottom and top of the slit, which may have artifacts """
	i = 4
	j = height-4

	xdim = 2
	ydim = 2

	avg = scipy.median(scipy.median(arc))
	m,std = clipped_std(arc,4.)

	""" The 'fiducial' arc spectrum is taken from the middle of the slit """
	arcmid = stats.stats.nanmedian(arc[height/2-3:height/2+4],0)

	""" First we straighten the lines out. """
	coords = spectools.array_coords((j-i,width))
	coords[0] += 4.
	fitdata = arc[i:j,height:-1*height].flatten()
	xdata = coords[1,:,height:-1*height].flatten()
	ydata = coords[0,:,height:-1*height].flatten()

	"""
	Only use every third row (but at least 8 rows) for the solution to save
	  time
	"""
	nskip = (j-i)/8
	if nskip>3:
		nskip = 3
	cond = ydata%nskip==0
	fitdata = fitdata[cond]
	ydata = ydata[cond]
	xdata = xdata[cond]

	""" Only use pixels around the peaks (another time saver) """
	thresh = scipy.median(arcmid)
	thresh += 5.*thresh**0.5
	mask = ndimage.maximum_filter(arcmid,15)
	mask = fitdata>thresh
	fitdata = fitdata[mask]
	ydata = ydata[mask]
	xdata = xdata[mask]

	p = scipy.zeros((xdim+1,ydim+1))
	p[1,0] = 1.

	p = {'coeff':p,'type':"chebyshev"}
	xback = arcfitfunc2(p,xdata,ydata,fitdata,arcmid)

	"""
	xdata, ydata, yorig, and newxdata hold original/transformed CCD
	  coordinates. We don't need to do a fit to *all* points, so we skip
	  ~every third.
	"""
	xdata = coords[1].flatten()
	ydata = coords[0].flatten()
	yorig = yforw[i:j].flatten()-curve
	cond = (ydata+xdata)%2==0
	xdata = xdata[cond]
	ydata = ydata[cond]
	yorig = yorig[cond]
	newxdata = special_functions.genfunc(xdata,ydata,xback)

	tmp = scipy.zeros((xdata.size,3))
	tmp[:,0] = newxdata.copy()
	tmp[:,1] = ydata.copy()
	tmp[:,2] = xdata.copy()
	xforw = special_functions.lsqfit(tmp,"chebyshev",xdim,ydim)

	"""
	Sky background model isn't needed if there aren't expected to be any
	  sky lines.
	"""
	if cutoff>5200:
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

	newarc = spectools.resample1d(arc,xforw,"x",-999)
	newarc[newarc==-999] = scipy.nan
	arcmid = stats.stats.nanmedian(newarc,axis=0)
	arcmid[scipy.isnan(arcmid)] = 0.

	"""
	We don't want to do the cross-correlation with all pixels; there
	  are reflections and other badness (like higher order lines). We
	  set the beginning of the 'useable' correlation spectrum to 150
	  pixels before the part of the arc with substantial flux. We set
	  the final pixel to be 100 pixels after the last substantial peak.
	  Here, 'substantial' means within 20% of the maximum; this in part
	  assumes that the final peak with 20% of the maximum will be the 5460
	  line, so we are in affect hard-coding that the correlation be carried
	  out to ~ 5500 angstroms.
	"""
	mask = scipy.where(arcmid>avg+4.*std)[0]
	first = mask[0]-150
	mask = ndimage.maximum_filter(arcmid,9)
	mask = scipy.where(mask==arcmid,mask,0)
	last = scipy.where(mask>0.20*arcmid.max())[0][-1]+100

	if first<0.:
		first = 0.
	if last>width:
		last = width

	"""
	Set reasonable starting and ending wavelengths for the correlation
	  routines.
	"""
	minwave = mswave-disp*width*1.1
	maxwave = mswave+disp*width*1.1
	if minwave<2000:
		minwave = 2000.
	if maxwave>8000:
		maxwave = 8000.

	"""
	Due a quick fit to define initial guesses...first the starting
	  wavelength. I've used both a correlation and a chi-square-like
	  function to determine the start; It's still not clear which is
	  most robust! The correlation also uses a broadened model of the
	  arcs to make it more robust wrt deviations from a linear solution.
	"""
	broad = 10.
	fudge = disp/200.
	nsteps = 19
	xmodel = arcmid[first:last].copy()
	xmodel = ndimage.gaussian_filter(xmodel,broad/disp)

	p0 = 0.
	p1 = disp
	offset = 0.
	max = 1e29
	for l in range(nsteps):
		try_disp = disp + ((l-nsteps/2)*fudge)
		skyfit_x = scipy.arange(minwave,maxwave,try_disp)
		fitmodel = interpolate.splev(skyfit_x,widemodel)
		tratio = xmodel.max()/fitmodel.max()
		fitmodel *= tratio
		fitmodel += scipy.median(xmodel)
		chi2,off = push(xmodel,fitmodel)
		if chi2<max:
			p1 = try_disp
			p0 = (off-first)*try_disp+minwave
			offset = off
			max = chi2
	firstwave = minwave+offset*p1

	"""
	We start with a second order fit: lambda = p0 + p1*x + p2*x*x
	   x = 0 is *not* neccessarily the first pix; it is the first *good*
	   pixel. p2 = 0. should be a reasonable first guess.
	"""
	p0 = blue
	first += (p0-firstwave)/p1

	print first, p0, p1
	if first<0:
		p0 = p0-first*p1
		first = 0
	last = first+(5650.-p0)/p1
	if last>width:
		last = width
	first = int(first)
	last = int(last)
	scale = last

	"""
	Create a normalized model of the arc data for the first refinement.
	"""
	fitdata = arcmid[first:last].copy()
	fitx = scipy.arange(fitdata.size).astype(scipy.float64)
	arclines = findlines(fitx,fitdata,20)

	fitmodel = scipy.zeros(3*len(arclines)+1)
	index = 1
	for iter in range(len(arclines)):
		fitmodel[index] = 1.
		fitmodel[index+1] = arclines[iter]
		fitmodel[index+2] = broad/disp
		index += 3
	arcmodel = special_functions.ngauss(fitx,fitmodel)


	""" We begin with a broadened arc-spectrum for the initial fit. """
	xdim = 2
	p = scipy.zeros((xdim+1,1))
	p[0,0] = p0
	p[1,0] = p1
	p = {'coeff':p,'type':"chebyshev"}
	fitdata = arcmid[first:last].copy()
	skycoeff = myoptimize(p,fitx,arcmodel,scale,linemodel)

	debug(fitx,arcmodel,linemodel,skycoeff)

	skycoeff = matchlines(arclines,p,linefile,order=1,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,order=2,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=10.*disp)

	""" 3rd order fit """
	xdim = 3
	p = scipy.zeros((xdim+1,1))
	p[0:skycoeff['coeff'].size,0] = skycoeff['coeff'][:,0].copy()
	skycoeff = {'coeff':p,'type':"chebyshev"}
	skycoeff = myoptimize(skycoeff,fitx,arcmodel,scale,linemodel)

	skycoeff = matchlines(arclines,skycoeff,linefile,order=2)
	skycoeff = matchlines(arclines,skycoeff,linefile)
	skycoeff = matchlines(arclines,skycoeff,linefile)

	skycoeff = matchlines(arclines,skycoeff,linefile,tol=10.)

	""" debug() is useful here to probe the quality of the line fitting """
#	debug(fitx,arcmodel,linemodel,skycoeff)
#	debug(fitx,fitdata,finemodel,skycoeff)

	xvals = scipy.arange(width)-first
	w = special_functions.genfunc(xvals,0,skycoeff)

	cond = (w>=blue)&(w<=cutoff)  # 'good' pixels
	cond1 = scipy.where(cond)[0]  # index of first good pixel

	fitx = scipy.arange(width).astype(scipy.float64)
	fitx = fitx[cond]-first
	fitdata = arcmid[cond]
	skycoeff = myoptimize(skycoeff,fitx,fitdata,scale,finemodel)

	arclines = findlines(fitx,fitdata,10)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=5.*disp)

	xdim = 5
	ydim = 2

	wave = special_functions.genfunc(newxdata-first,0.,skycoeff)
	sky2x = []
	sky2y = []
	ccd2wave = []

	"""
	We can only use the arcs if the cutoff is too blue.
	"""
	if cutoff<5200:
		revmodel = scipy.zeros((wave.size,3))
		revmodel[:,0] = wave.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = xdata.copy()
		sky2x.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))

		revmodel[:,2] = yorig.copy()
		sky2y.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		revmodel[:,0] = xdata.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = wave.copy()
		ccd2wave.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		for k in range(1,nsci):
			sky2x.append(sky2x[0])
			sky2y.append(sky2y[0])
			ccd2wave.append(ccd2wave[0])
		return sky2x,sky2y,ccd2wave


	"""
	Try to use the 5577 line to refine the wavelength solution.
	"""
	wlen = wave.copy()
	xvals = scipy.arange(width)-first
	for k in range(nsci):
		peaks = findlines(xvals,bgmodel[k],5.)
		w = special_functions.genfunc(peaks,0.,skycoeff)
		delta = 5577.345-w
		offset = delta[abs(delta).argmin()]
		if abs(offset)>3.*disp:
			print "Warning: Substantial offset found!"

		wave = wlen+offset
		revmodel = scipy.zeros((wave.size,3))
		revmodel[:,0] = wave.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = xdata.copy()
		sky2x.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		revmodel[:,2] = yorig.copy()
		sky2y.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))
		revmodel[:,0] = xdata.copy()
		revmodel[:,1] = ydata.copy()
		revmodel[:,2] = wave.copy()
		ccd2wave.append(special_functions.lsqfit(revmodel,"chebyshev",xdim,ydim))

	return sky2x,sky2y,ccd2wave
