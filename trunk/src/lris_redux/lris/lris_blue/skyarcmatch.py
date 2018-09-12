"""
skyarcmatch

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

TODO:
  Right now the pixel-by-pixel fitting is disabled because the arcs
    aren't characterized well enough (the relative strengths for a given
    lamp might be, but the ratios between lamps change and this is not
    accounted for yet). However, the line matching is also kind of
    crappy because: (1) the sky has sparse lines with varying intensities
    that shift the centroid (because low level lines that are 'ignored'
    become high and merge with real lines, altering their center--this
    happens for the 5200. line and the 6363 line); (2) the bluest arc
    lines are merged and difficult to centroid; and (3) the 4800 and 4810
    arc lines merge.

"""

from mostools import spectools
import special_functions

import scipy
from scipy import optimize,interpolate,ndimage,signal,stats,random



""" Define saturation level for arclines """
SATURATED = 57000.


""" List of lines used for the line matching """
#SKYLINES = [5577.338,6300.304,6363.78,6553.617,6912.62]
#SKYLINES = [5197.93,5460.735,5577.345,5891.9,6300.32,6363.81,6533.04,6553.61,6912.62,6923.21,6939.52]
#ARCLINES = [4046.563,4077.831,4358.327,4678.149,4680.14,4722.15,4799.912,4810.53,5085.822,5460.735]


# Strong sky lines. Blends of 5198/5200 and NaD.
SKYLINES = [5198.794,5200.286,5460.735,5577.345,5891.9,6300.32,6363.81,6533.04,6553.61,6912.62,6923.21,6939.52]

# Common arc lines. Blend of 4678.149 and 4680.14.
ARCLINES = [3650.158,3654.840,3663.28,4046.563,4077.831,4358.327,4679.149,4722.15,4799.912,4810.53,5085.822,5460.735]


""" Helper function to simply update the processing log """
def keeplog(logfile,entry):
	logfile = open(logfile.name,'a')
	logfile.write(entry)
	logfile.close()

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
A wrapper to the optimize.leastsq function. It scales the wavelength
data to make the non-linear solver more robust.
"""
def myoptimize(p,x,z,model,model2=None):
	par = special_functions.unpack_coeff(p)
	mask = ndimage.maximum_filter(z,3)
	coeff,ier = optimize.leastsq(arcfitfunc,par,(x,z,model,p,model2,mask),maxfev=100000)

	par = special_functions.build_coeff(coeff,p)
	return par

"""
The chi-square function for myoptimize. Compares a model of the arcs with
the arc data.
"""
def arcfitfunc(p,x,data,model,tmp,model2=None,mask=None):
	par = special_functions.build_coeff(p,tmp)

	""" Test for increasing function... """
	tmp = special_functions.genfunc(x,0.,par).astype(scipy.float32)
        diff = signal.convolve(tmp,[1.,-1.],mode='same')[10:-10].copy()
        if diff[diff<0].size>0:
                return scipy.ones(x.size)*1.e9

	data = data.astype(scipy.float64)
	z = special_functions.genfunc(x,0,par).astype(scipy.float64)

	mod = interpolate.splev(z,model).astype(scipy.float64)
	if model2 is not None:
		mod += interpolate.splev(z,model2).astype(scipy.float64)

	mod = mod*data.max()/mod.max() + scipy.median(data)
	diff = (data-mod)/scipy.sqrt(abs(mod))
	return diff

def myopt2(p,x,z,model,skyx,sky,skymodel,r,o):
	par = special_functions.unpack_coeff(p)
	coeff,ier = optimize.leastsq(skyarcfit,par,(x,z,model,skyx,sky,skymodel,p,r,o),maxfev=100000)
	par = special_functions.build_coeff(coeff,p)
	return par

def skyarcfit(p,x,data,model,skyx,sky,skymodel,tmp,r,o):
	par = special_functions.build_coeff(p,tmp)
	""" Test for increasing function... """
	tmp = special_functions.genfunc(x,0.,par).astype(scipy.float32)
	diff = signal.convolve(tmp,[1.,-1.],mode='same')[10:-10].copy()
	if diff[diff<0].size>0:
		return scipy.ones(x.size)*1.e9
	data = data.astype(scipy.float64)

	z = special_functions.genfunc(skyx,0,par).astype(scipy.float64)
	skymod = interpolate.splev(z,skymodel).astype(scipy.float64)
	skymod = (skymod)*r
#	skymod = skymod*scipy.median(sky)/scipy.median(skymod)

	z = special_functions.genfunc(x,0,par).astype(scipy.float64)
	mod = interpolate.splev(z,model).astype(scipy.float64)
	mod = mod*data.max()/mod.max() + scipy.median(data)

	diff1 = (data-mod)/scipy.sqrt(abs(mod))
#	diff2 = (sky-skymod)/scipy.sqrt(abs(sky))**0.25
	diff2 = (sky-skymod)

	return scipy.concatenate((diff1,diff2))

def skyopt(p,x,data,model):
	par = special_functions.unpack_coeff(p).tolist()

	wave = special_functions.genfunc(x,0,p).astype(scipy.float64)
	sky = interpolate.splev(wave,model).astype(scipy.float64)
	ratio = scipy.median(data)/scipy.median(sky)
	offset = 0.
	par.append(ratio)
	par.append(offset)

	coeff,ier = optimize.leastsq(skyfit,par,(x,data,model,p),maxfev=100000)
	return special_functions.build_coeff(coeff,p),coeff[-2],coeff[-1]

def skyfit(p,x,data,model,tmp):
	par = special_functions.build_coeff(p[:-2],tmp)
	""" Test for increasing function... """
	tmp = special_functions.genfunc(x,0.,par).astype(scipy.float32)
        diff = signal.convolve(tmp,[1.,-1.],mode='same')[10:-10].copy()
        if diff[diff<0].size>0:
                return scipy.ones(x.size)*1.e9

	data = data.astype(scipy.float64)

	z = special_functions.genfunc(x,0,par).astype(scipy.float64)
	sky = interpolate.splev(z,model).astype(scipy.float64)
	sky *= p[-2]

	return (sky-data)/scipy.sqrt(abs(sky))

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
			return std
		mean = d.mean()
		std = d.std()
	return std

"""
Finds and centroids peaks in the spectrum.
"""
def findlines(x,z,sigma=3.):
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
		std = clipped_std(z[start:end],3.5)
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

"""
Match arclines with peaks.
"""
def matchlines(peaks,solution,linefile,tol=30.,order=3,offset=False):
	wave = special_functions.genfunc(peaks,0.,solution)
	try:
		lines = getlines(linefile)
	except:
		lines = linefile

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

	if offset:
		a = special_functions.genfunc(fitdata[:,0],0.,solution)
		return fitdata[:,1]-a

	return special_functions.lsqfit(fitdata,'chebyshev',order)

"""
Match arclines and skylines jointly.
"""
def doublematch(peaks,w,linefile,skyin,tol=30.,order=3,shift=0.,logfile=None):
	lines = scipy.asarray(skyin)
	lines.sort()

	p = peaks[0]
	wave = special_functions.genfunc(p,0.,w)
	goodmatches = []
	for i in range(wave.size):
		delta = 1e9
		match = None
		for j in lines:
			if abs(j-wave[i])<delta:
				delta = abs(j-wave[i]+shift)
				match = j
			else:
				break
		if delta<tol:
			goodmatches.append([p[i],match])

	p = peaks[1]
	wave = special_functions.genfunc(p,0.,w)
	lines = scipy.asarray(linefile)
#	lines = getlines(linefile)
	lines.sort()

	for i in range(wave.size):
		delta = 1e9
		match = None
		for arc in lines:
			if abs(arc-wave[i])<delta:
				delta = abs(arc-wave[i]+shift)
				match = arc
			else:
				break
		if delta<tol:
			goodmatches.append([p[i],match])

	fitdata = scipy.asarray(goodmatches)
	fit = special_functions.lsqfit(fitdata,'chebyshev',order)
	match = special_functions.genfunc(fitdata[:,0],0.,fit)
	diff = match-fitdata[:,1]
	error = stats.stats.std(diff)
	if logfile is not None:
		logentry = "\tWavelength solution error: %5.3f angstroms from %d lines\n" % (error,diff.size)
		for i in range(diff.size):
			logentry += "\t\t%7.2f   %7.2f   %5.2f\n" % (fitdata[i,1],match[i],diff[i])
		keeplog(logfile,logentry)
	return fit,error


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
	"""
	arcmatch(curve,sci,arc,yforw,widemodel,finemodel,goodmodel,linemodel
	          ,disp,mswave,extra,logfile)

	Creates the full 2d distortion solution for a slit.

	Inputs:
	  curve     - first pixel (ie bottom) of curved slit
	  sci       - biastrimmed, non-ycorrected science data
	  arc       - y-corrected arc for this slit
	  yforw     - yforw definition for slit
	  widemodel - broad wavelength model
	  finemodel - matched-resolution wavelength model
	  goodmodel - model of the sky
	  linemodel - broad wavelength model (why are there two...?)
	  disp      - scale
	  mswave    - central wavelength
	  extra     - includes name of linefile, starting/ending wavelengths
	  logfile   - a file object

	Outputs:
	  
	"""

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
	cutoff = 5650.

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
	std = clipped_std(arc,4.)

	""" The 'fiducial' arc spectrum is taken from the middle of the slit """
	arcmid = stats.stats.nanmedian(arc[height/2-2:height/2+3],0)

	""" First we straighten the lines out. """
	coords = spectools.array_coords((j-i,width))
	coords[0] += 4.
	fitdata = arc[i:j,height:-1*height].flatten()
	xdata = coords[1,:,height:-1*height].flatten()
	ydata = coords[0,:,height:-1*height].flatten()

	""" Only use every third row for the solution (to save time) """
	cond = ydata%3==0
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

	p = {'coeff':p,'type':"polynomial"}
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
	mask = scipy.where(mask==arcmid,mask,0)[:arcmid.argmax()+460/disp]
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
	  function to determine the start; it's still not clear which is
	  most robust! The correlation also uses a broadened model of the
	  arcs to make it more robust wrt deviations from a 1st order solution.
	"""
	broad = 9.
	fudge = disp/1000.
	nsteps = 35
	xmodel = arcmid[first:last].copy()
#	xmodel = ndimage.gaussian_filter(xmodel,broad/disp)
	
	p0 = 0.
	p1 = disp
	offset = 0.
	max = 1e29
	c = []
	b = []
	fitdata = arcmid[first:last].copy()
	fitx = scipy.arange(fitdata.size).astype(scipy.float64)
	arclines = findlines(fitx,fitdata,10)

	fitmodel = scipy.zeros(3*len(arclines)+1)
	index = 1
	for iter in range(len(arclines)):
		fitmodel[index] = 1.
		fitmodel[index+1] = arclines[iter]
		fitmodel[index+2] = broad/disp
		index += 3
	xmodel = special_functions.ngauss(fitx,fitmodel)
	for l in range(nsteps):
		try_disp = disp + ((l-nsteps/2)*fudge)
		skyfit_x = scipy.arange(minwave,maxwave,try_disp)
		fitmodel = interpolate.splev(skyfit_x,linemodel)
#		tratio = xmodel.max()/fitmodel.max()
#		fitmodel *= tratio
#		fitmodel += scipy.median(xmodel)
		corr = signal.correlate(xmodel,fitmodel,mode='valid')
		chi2,off = 1./corr.max(),corr.argmax()
#		chi2,off = push(xmodel,fitmodel)
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

	if first<0:
		p0 = p0-first*p1
		first = 0
	last = first+(5650.-p0)/p1
	if last>width:
		last = width
	first = int(first)
	last = int(last)
	keeplog(logfile,'\tInitial Solution: %f, %f\n' % (p0,p1))

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

	""" We begin with a fit to the normalized spectra. """
	xdim = 2
	p = scipy.zeros((xdim+1,1))
	p[0,0] = p0
	p[1,0] = p1
	p = {'coeff':p,'type':"polynomial"}
	fitdata = arcmid[first:last].copy()
	skycoeff = myoptimize(p,fitx,arcmodel,linemodel)
#	debug(fitx,fitdata,finemodel,skycoeff)

	skycoeff = matchlines(arclines,p,linefile,order=1,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,order=2,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=10.*disp)
	""" Refinement of fit, still with normalized spectra. """
	xdim = 3
	p = scipy.zeros((xdim+1,1))
	p[0:skycoeff['coeff'].size,0] = skycoeff['coeff'][:,0].copy()
	skycoeff = {'coeff':p,'type':"chebyshev"}
#	skycoeff = myoptimize(skycoeff,fitx,arcmodel,linemodel)
	skycoeff = myoptimize(skycoeff,fitx,fitdata,finemodel)

#	debug(fitx,fitdata,finemodel,skycoeff)
	""" Match lines, progressively constraining the fit """
	skycoeff = matchlines(arclines,skycoeff,linefile,order=2,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=10.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=5.*disp)
	skycoeff = matchlines(arclines,skycoeff,linefile,tol=5.*disp)
#	debug(fitx,fitdata,finemodel,skycoeff)

	xvals = scipy.arange(width)-first
	w = special_functions.genfunc(xvals,0,skycoeff)

	"""
	The wavelength solution might wrap around on itself (ie it might not be
	  a strictly increasing function). This ensures that only the
	  increasing part of the current solution is used. It might be more
	  reasonable to put something in the chi-square to enforce this....
	  I've added it to arcfitfunc for now....
	val = w.argmax()
	if cutoff>w.max():
		cutoff = w.max()
	w[val:] = cutoff + 1.
	"""

	cond = (w>=blue)&(w<=cutoff)  # 'good' pixels
	cond1 = scipy.where(cond)[0]  # index of first good pixel

	fitx = scipy.arange(width).astype(scipy.float64)
	fitx = fitx[cond]-first
	fitdata = arcmid[cond]
#	debug(fitx,fitdata,finemodel,skycoeff)
#	skycoeff = myoptimize(skycoeff,fitx,fitdata,finemodel)

#	xdim = 5
#	p = scipy.zeros((xdim+1,1))
#	p[0:skycoeff['coeff'].size,0] = skycoeff['coeff'][:,0].copy()
#	skycoeff = {'coeff':p,'type':"polynomial"}
#	skycoeff = myoptimize(skycoeff,fitx,fitdata,finemodel)

	

#	debug(fitx,fitdata,finemodel,skycoeff)

	coeff = skycoeff['coeff'].copy()
	arclamppeaks = findlines(fitx,fitdata)
	skycoeff = matchlines(arclamppeaks,skycoeff,linefile,tol=5.*disp)

	startx = fitx.copy()
	startdata = fitdata.copy()
	wave = special_functions.genfunc(startx,0.,skycoeff)
	model = interpolate.splrep(wave,startdata,s=0)

	sky2x = []
	sky2y = []
	ccd2wave = []

	def findoffset(p,x,data,model,coeff):
		xvals = x+p[0]
		w = special_functions.genfunc(xvals,0.,coeff)
		mod = interpolate.splev(w,model)
		mod = mod*data.max()/mod.max()
		diff = (mod-data)/scipy.sqrt(abs(mod)) 
		return diff

	"""
	If the highest line is greater than 3x the next highest, there
	  was probably some sky activity. This could throw off the
	  pixel-based fitting, so we should choose line-matching
	  instead. A better solution would be to 'dampen' the active
	  line (we can't just alter the residuals because the sky bias
	  and gain offsets were calculated using the line), but that is
	  non-trivial to accomplish without shifting the central
	  wavelength. TURNED OFF--LINE MATCHING ALWAYS USED!
	"""
	LINEMATCH = False
	for k in range(nsci):
		sky = bgmodel[k]
		x_vals = scipy.arange(width)
		skypeaks = findlines(x_vals,sky,3)
		amps = ndimage.map_coordinates(sky,[skypeaks],output=scipy.float64,order=5)
		amps.sort()
		if amps[-1]>2.*amps[-2]:
			LINEMATCH = True
	# Use LINEMATCH until arcs are better dealt with
	LINEMATCH = True

	"""
	We attempt to determine an offset between the arclamp frame and the sky
	  image due to flexure (or if the arcs were taken in the afternoon...).
	  First we use the arc solution to determine the offset of the bluest
	  sky line. We then shift the arc lines to correct for the offset and
	  find a solution for just the sky. We extrapolate the sky solution to
	  the reddest arc line, determine the arc line's offset and use this as
	  the final offset between sky and arc. We then solve for the full
	  solution using the sky and corrected arc lines.
	"""
	for k in range(nsci):
		sky = bgmodel[k]

		skycoeff['coeff'] = coeff.copy()
		x_vals = scipy.arange(width)-first

		skypeaks = findlines(x_vals,sky,5)
		"""
		Start by trying to match to the 5200 line; if it is too faint
		  (and therefore was not extracted as a skypeak), use the
		  5577 line. 5197.93,5198.794 -- the line is a blend and the
		  left value is the non-blended wavelength. One issue is that
		  the line exhibits auroral brightening--which appears to
		  change the effective wavelength because the fainter line at
		  5200. is the one that brightens! BAH! We'll just use 5577.
		"""
#		sline = 5198.794
		sline = 5577.345
		highpeaks = findlines(x_vals,sky,7)
		w = special_functions.genfunc(highpeaks,0.,skycoeff)
		delta = sline-w
		if abs(delta).min()>10.*disp:
			sline = 5577.345
			delta = sline-w
		w = special_functions.genfunc(x_vals,0.,skycoeff)
		skycond = (w>sline-40.)&(w<sline+40.)
		off = delta[abs(delta).argmin()]/disp

		tmp = special_functions.genfunc(skypeaks,0.,skycoeff)
		entry = '\tSky lines found for image %d:' % (k+1)
		for t in tmp:
			entry += ' %7.2f' % t
		entry += "\n"
		keeplog(logfile,entry)

		min = 10.
		step = None
		for i in scipy.arange(-3.,3.,0.1):
			arcpeaks = arclamppeaks+i
			tmpcoeff,err = doublematch([skypeaks,arcpeaks],skycoeff,ARCLINES,[sline],tol=10.*disp)
			tmpcoeff,err = doublematch([skypeaks,arcpeaks],tmpcoeff,ARCLINES,SKYLINES,tol=10.*disp)
			tmpcoeff,err = doublematch([skypeaks,arcpeaks],tmpcoeff,ARCLINES,SKYLINES,tol=5.*disp)
			if err<min:
				min = err
				outcoeff = tmpcoeff.copy()
				step = i
		min = 10.
		for i in scipy.arange(step-0.5,step+0.5,0.01):
			arcpeaks = arclamppeaks+i
			tmpcoeff,err = doublematch([skypeaks,arcpeaks],outcoeff,ARCLINES,SKYLINES,tol=5.*disp)
			if err<min:
				min = err
				skycoeff = tmpcoeff.copy()
				arcoffset = i
		keeplog(logfile,"\tArc offset for image %d: %6.3f\n" % (k+1,arcoffset))
		skypeaks = findlines(x_vals,sky,2.)
		arcpeaks = arclamppeaks+arcoffset
		skycoeff,e = doublematch([skypeaks,arcpeaks],skycoeff,ARCLINES,SKYLINES,tol=5.*disp,logfile=logfile)

		w = special_functions.genfunc(x_vals,0.,skycoeff)
		skycond = (w>5500.)&(w<6700.)
		skyx = x_vals[skycond]
		skyz = sky[skycond]
		tmp = skyz.copy()
		tmp.sort()
		med = tmp[tmp.size*0.3]
		skyx = skyx[skyz>med]
		skyz = skyz[skyz>med]
		tmpcoeff,ratio,offset = skyopt(skycoeff,skyx,skyz,goodmodel)
		skycoeff = myopt2(skycoeff,fitx+arcoffset,fitdata,finemodel,skyx,skyz,goodmodel,ratio,offset)


		""" Evaluate wavelength of all pixels in straightened frame """
		wave = special_functions.genfunc(newxdata-first,0,skycoeff)

		"""
		Now that we have the wavelength solution for the rectified slit
		  we can create the full 2d solutions for:
			x(lambda,pos)
			y(lambda,pos)
			lamba(x,y)
		  where x and y are CCD coordinates and lambda and pos are the
		  'true' coordinates (wavelength and spatial position along
		  the slit).
		"""
		xdim = 5
		ydim = 2

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

#	return sky2x,sky2y,[xforw,xback] # This is for non-2d subtraction
	return sky2x,sky2y,ccd2wave      # Thsi is for 2d subtraction
