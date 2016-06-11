import scipy,special_functions
from scipy import signal,stats,ndimage

# From an input 2d array, find and extract spectral traces.
#   Returns a list with entries that include fit parameters, the extracted
#   spectrum, a smoothed spectrum, and a simple model of the noise.

# (TUNEABLE) CONSTANTS
WIDTH = 10.     # Number of pixels to fit profile to = 2*WIDTH + 1
NSIG = 1.5      # Width (in sigmas) to extract
FILTSIZE = 7    # Wiener filter smoothing size
NOISE = 5.	# Signal to noise limit

def extract(data,varimg,width=WIDTH,nsig=NSIG,noise=NOISE):
	WIDTH = width
	NSIG = nsig
	NOISE = noise

	data = data.copy()
	spectra = []

	# Replace nan with zero
	data[scipy.isnan(data)] = 0.
	varimg[scipy.isnan(varimg)] = 0.

	# Create model of real flux. We ignore the slit ends, which may have
	#  artifacts from the resampling.
	slit = data[:,8:-8].astype(scipy.float32)
	var = varimg[:,8:-8]

	# OK...so negative-variance also isn't good; set these pixels to zero
	var[var<0] = 0

	# Create noise models
	sigmaimg = slit/scipy.sqrt(var)
	highpix = scipy.where(sigmaimg>1.5,sigmaimg,0.)
	source_columns = highpix.sum(axis=0)

	# MASKING DISABLED (this would take only columns with lotsa flux...)
#	mask = scipy.where(source_columns>4.,1.,scipy.nan)
	mask = source_columns*0.

	# Condition 1, dealing with bad pixels
	if (var==0).any():
		cond = var==0
		var[cond] = scipy.nan
		slit[cond] = scipy.nan
		mask = scipy.where(cond,0,1)
		flux = scipy.nansum(slit/var,axis=1)/scipy.nansum(1./var,axis=1)
		noise = scipy.sqrt(scipy.nansum(var,axis=1))/mask.sum(axis=1)
	# Condition 2, no masking
	elif scipy.nansum(mask)==0:
		flux = (slit/var).sum(axis=1)/(1./var).sum(axis=1)
		noise = scipy.sqrt(var.sum(axis=1))/mask.size
	# Condition 3, masking
	else:
		fluxmodel = slit*mask
		noisemodel = var*mask

		noise = scipy.sqrt(scipy.nansum(noisemodel,axis=1))/scipy.nansum(mask)
		flux = stats.stats.nanmean(fluxmodel,axis=1)

	# A smooth S/N estimate for the slit
#	sig2noise = ndimage.gaussian_filter1d(flux,1)/noise

	row = scipy.arange(flux.size)
	model = flux.copy()
	nspec = 10	# Maximum number of attempts
	while nspec:
		nspec -= 1

		# Fit a gaussian around the peak of the S/N model
		start = model.argmax()-WIDTH
		end = model.argmax()+WIDTH+1
		if start<0:
			start = 0.
		if end>model.size:
			end = model.size

		fitarr = model[start:end]
		p = scipy.zeros(4)
		p[1] = fitarr.max()
		p[2] = fitarr.argmax()
		p[3] = 2.

		fit,val = special_functions.ngaussfit(fitarr,p)
		chi2 = val/(fitarr.size-3)
		fit[2] += start

		# If the centroid doesn't lie on the slit, get use the edge pix
		midcol = fit[2].round()
		if midcol>=flux.size:
			midcol = flux.size-1
		elif midcol<0:
			midcol = 0
		# Require a reasonable S/N and width
		if fit[3]>fitarr.size/2. or fit[3]<0.85:
			break
		elif fit[0]>0 and fit[1]<NOISE*noise[midcol]:
			break
		elif fit[0]<0 and fit[1]-fit[0]<NOISE*noise[midcol]:
			break
		else:
			fit[1] += fit[0]
			fit[0] = 0.
			# Subtract away a model of the source
			source = special_functions.ngauss(row,fit)
			model -= scipy.where(source>noise,source,0.)

			# Skip residuals!
			if fit[2]<flux.size and fit[1]<scipy.sqrt(flux[fit[2]]):
				continue
			fit[1] = 1.
			weight = special_functions.ngauss(row,fit)
			cond = (row>fit[2]-fit[3]*NSIG)&(row<fit[2]+fit[3]*NSIG)
			weight = scipy.where(cond,weight,0)
			weight /= weight.sum()
			spec = weight*data.T
			spec = spec.sum(axis=1)
			varspec = weight*varimg.T
			varspec = varspec.sum(axis=1)
			spec[varspec==0] = 0.
			smooth = signal.wiener(spec,FILTSIZE,varspec)
			smooth[scipy.isnan(smooth)] = 0.
			spectra.append([fit,spec,smooth,varspec])
	return spectra
