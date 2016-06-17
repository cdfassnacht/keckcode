"""
Pipeline to reduce LRIS redside spectra.
 Inputs:
   prefix      - Filename prefix, usually 'lred'
   dir         - Directory with input files, ie "raw/" (backslash is important!)
   science     - Numbers of the science files, ie "0023,0024,0029,0030"
   arc         - Number of arc file, ie "0025"
   flats       - Numbers of flat files, ie "0026,0027,0028"
   outprefix   - Prefix for output files, ie "mask2_11oct04"
   useflat     - 1 to use flat data from a previous run, otherwise 0
   usearc      - 1 to use arc data from a previous run, otherwise 0
   cache       - 1 to cache data to disk (useful for blueside with RAM<2GB)
   offsets     - a list/array of relative offsets between masks (in pixels); this will be unnecessary if the stars remained in the starboxes for all masks.

"""


import lris
from lris.lris_biastrim import redside as biastrim

from lris.lris_red.flat import *
from lris.lris_red.skymatch import skymatch as wavematch
from lris.lris_red.skysub import doskysub

from mostools import spectools,offset,measure_width
from mostools.extract import extract
import special_functions

from math import ceil,fabs
import pickle

import numpy as np
import scipy
from scipy import stats,interpolate,ndimage
from scipy import io as sio
from astropy.io import fits as pyfits


""" A control routine to encapsulate the pipeline. """
def lris_pipeline(prefix,dir,scinames,arcname,flatnames,out_prefix,useflat=0,usearc=0,cache=0,offsets=None):
	print "Processing mask",out_prefix


	nsci = len(scinames)

	print "Preparing flatfields"
	if useflat==1:
		yforw,yback,slits,starboxes,flatnorm = flatload(out_prefix)
	else:
		yforw,yback,slits,starboxes,flatnorm = flatpipe(flatnames,out_prefix)
	axis1 = flatnorm.shape[0]
	axis2 = flatnorm.shape[1]

	"""
	Read lamps data from the arclamp file; this is unnecssary for the red
	  side unless the line fitting is altered to better calibrate the blue
	  end (ie for 460 dichroic data).
	"""
	print "Preparing arcs for line identification"
	if usearc==1:
		arcdata = biastrim(pyfits.open(arcname)[0].data)
		arcname = out_prefix+"_arc.fits"
		arc_tmp = pyfits.open(arcname)
		arc_ycor = arc_tmp[0].data.astype(scipy.float32)
		lamps = arc_tmp[0].header['LAMPS']
		del arc_tmp
	else:
		arc_tmp = pyfits.open(arcname)
		arcdata = arc_tmp[0].data.copy()
		lamps = arc_tmp[0].header['LAMPS']
		del arc_tmp
		arcdata = biastrim(arcdata)
		arc_ycor = spectools.resampley(arcdata,yforw).astype(scipy.float32)
		arcname = out_prefix+"_arc.fits"
		arc_hdu = pyfits.PrimaryHDU(arc_ycor)
		arc_hdu.header.update('LAMPS',lamps)
		arc_hdu.writeto(arcname)
		del arc_hdu

	"""
	Skysubtraction, centering, &c. may work better if there is some slop on
	  the sides of the data (the slit definitions have been created to only
	  include good data, so 'bad' edges are rejected).
	"""
	wide_stars = []
	for i,j in starboxes:
		mod = scipy.where((yback<j)&(yback>i))
		a = mod[0].min()-3
		b = mod[0].max()+3
		if a<0:
			a = 0
		if b>axis1:
			b = axis1
		wide_stars.append([a,b])


	print "Bias trimming and flatfielding science data"
	scidata = scipy.zeros((nsci,axis1,axis2),'f4')
	center = scipy.zeros((nsci,len(starboxes)),'f4')
	flux = scipy.zeros((nsci),'f4')
	airmass = []
	for i in range(nsci):
		filename = scinames[i]
		scitmp = pyfits.open(filename)

		scidatatmp = scitmp[0].data.copy()
		scidatatmp = biastrim(scidatatmp).astype(scipy.float32)

		"""
		The biastrim routine should take care of bad columns, but this
		  is just in case; we do a simple linear interpolation over
		  bad columns.
		"""
		bad = scipy.where(scidatatmp>56000.)
		nbad = bad[0].size
		for k in range(nbad):
			y = bad[0][k]
			x = bad[1][k]
			if x==0:
				x1 = x+1
				x2 = x+1
			elif x == scidatatmp.shape[1]-1:
				x1 = x-1
				x2 = x-1
			else:
				x1 = x-1
				x2 = x+1
			scidatatmp[y,x] = \
			    (scidatatmp[y,x1]+scidatatmp[y,x2])/2.

		"""
		Apply the flatfield and copy the data into the working array.
		"""		  
		scidatatmp = scidatatmp/flatnorm
		scidata[i,:,:] = scidatatmp.copy()

		"""
		Copy key header keywords; note that old data might not have
		  MSWAVE or DICHNAME keywords.
		"""
		try:
			mswave = scitmp[0].header['MSWAVE']
		except:
			mswave = 6500.
		disperser = scitmp[0].header['GRANAME']
		airmass.append(scitmp[0].header['AIRMASS'])
		try:
			dichroic = scitmp[0].header['DICHNAME']
		except:
			dichroic = None

		"""
		This should give a reasonable estimate of the sky level; the
		  program does a dumb scaling (to the level of exposure with the
		  highest sky level)
		"""
		flux[i] = scipy.sort(scipy.ravel(scidatatmp))[scidatatmp.size/4]

		"""
		Centroid stars in starboxes to find shifts between mask
		  exposures.
		"""
		for j in range(len(starboxes)):
			a,b = starboxes[j]
			m,n = wide_stars[j]
			a -= 4
			b += 4
			m -= 2
			n += 2
			if m<0:
				m = 0
			if n>scidatatmp.shape[0]:
				n = scidatatmp.shape[0]
			if a<0:
				a = 0
			if b>yforw.shape[0]:
				b = yforw.shape[0]
			center[i,j] = offset.findoffset(scidatatmp[m:n],yforw[a:b],m)

		del scitmp
		del scidatatmp
	del flatnorm

	"""
	This implements the mechanism for manually entering offsets (if for
	  example we dithered the stars out of the starboxes).
	"""
	if offsets is not None:
		center = scipy.asarray(offsets)
	else:
		center = stats.stats.nanmean(center,axis=1)

	"""
	Perform the flux scaling and set the offsets relative to each other.
	"""
	print "Normalizing Fluxes"
	cmax = center.max()
	fmax = flux.max()
	for i in range(center.size):
		center[i] -= cmax
		ratio = fmax/flux[i]
		scidata[i] *= ratio
	cmax = ceil(fabs(center.min()))


	"""
	Set the output scale (and approximate input scale), as well as blue
	  cutoff limits.
	"""
	if disperser=="150/7500":
		scale = 4.8
	elif disperser=="300/5000":
		scale = 2.45
	elif disperser=="400/8500":
		scale = 1.85
	elif disperser=="600/5000":
		scale = 1.25
	elif disperser=="600/7500":
		scale = 1.25
	elif disperser=="600/10000":
		scale = 1.25
	elif disperser=="831/8200":
		scale = 0.915
	elif disperser=="900/5500":
		scale = 0.85
	elif disperser=="1200/7500":
		scale = 0.64

	if dichroic=='mirror':
		redcutoff = 4000.
		dich_file = ''
	elif dichroic=='460':
		redcutoff = 4600.  # I haven't checked this...
		dich_file = '460'
	elif dichroic=='500':
		redcutoff = 5000.  # I haven't checked this...
		dich_file = '500'
	elif dichroic=='560':
		redcutoff = 5500.
		dich_file = '560'
	elif dichroic=='680':
		redcutoff = 6700.
		dich_file = '680'
	else:
		redcutoff = 3500.
		dich_file = ''

	"""
	Determine the y-size of the output arrays. We also find an estimate of
	  the mask resolution while looping through. Only every seventh slit
	  is examined to expedite the process.
	"""
	nsize = 0
	csize = 0
	wide_slits = []
	linewidth = []
	for i,j in slits:
		csize += int(j-i+cmax) + 5
		nsize += j-i+5
		mod = scipy.where((yback>i)&(yback<j))
		a = mod[0].min()-4
		b = mod[0].max()+4
		if a<0:
			a = 0
		if b>axis1:
			b = axis1
		wide_slits.append([a,b])
		if len(wide_slits)%7==0:
			linewidth.append(measure_width.measure(arc_ycor[(i+j)/2,:],15))
	csize -= 5
	nsize -= 5

	linewidth = scipy.median(scipy.asarray(linewidth))

	print "Loading wavelength model"
	lris_path = lris.__path__[0]

	filename = lris_path+"/data/uves_sky.model"
	infile = open(filename,"r")
	wavecalmodel = pickle.load(infile)
	infile.close()
	wave = scipy.arange(3400.,10400.,0.1)

	"""
	We make the sky spectrum slightly more realistic by taking into account
	  the dichroic cutoff. This mainly helps with matching the 5577 line
	  for the 560 dichroic. It would be nice if the response of the
	  instrument was somewhat well characterized for all grating/dichroic
	  combinations....
	"""
	if dich_file!='':
		filename = lris_path+"/data/dichroics/dichroic_"+dich_file+"_t.dat"
		infile = open(filename,"r")
		#input = sio.read_array(infile)
		input = np.loadtxt(infile)
		infile.close()
		spline = interpolate.splrep(input[:,0],input[:,1],s=0)
		dich = interpolate.splev(wave,spline)
		dich[wave<4500.] = 1.
		dich[wave>8800.] = 1.
		del input,spline
	else:
		dich = scipy.ones(wave.size)

	"""
	Create two sky spectrum spline models. One is a 'fine' model matched
	  to the resolution of the instrumental setup. The other is a widened
	  model for coarse wavelength matching.
	"""
	wavemodel = interpolate.splev(wave,wavecalmodel)
	finemodel = ndimage.gaussian_filter1d(wavemodel,linewidth*scale/0.1)
	wavemodel = ndimage.gaussian_filter1d(finemodel,5./0.1)
	finemodel *= dich
	finemodel = interpolate.splrep(wave,finemodel,s=0)
	wavemodel *= dich
	widemodel = interpolate.splrep(wave,wavemodel,s=0)
	goodmodel = finemodel
	del dich,wave,wavemodel

	""" See extract.py; sets default extraction width. """
	extractwidth = 10

	print "Creating output arrays"

	"""
	We choose an output array size that *should* be large enough to contain
	  all of the valid data (given reasonable assumptions about how far
	  the slits are placed from the center of the mask). We could also
	  decrease the needed size by enforcing the blue limit....
	"""
	outlength = int(axis2*1.6)
	out = scipy.zeros((nsci,nsize,outlength))*scipy.nan
	out2 = scipy.zeros((2,csize,outlength))*scipy.nan

	"""
	For systems with limited RAM, it might make sense to cache the output
	  arrays to disk. This increases the time it takes to run but may be
	  necessary and also allows the progress of the reduction to be
	  monitored.
	"""
	if cache:
		import os
		print "Caching..."
		strtfile = out_prefix+"_TMPSTRT.fits"
		bgfile = out_prefix+"_TMPBSUB.fits"
		try:
			os.remove(strtfile)
		except:
			pass

		outfile = pyfits.PrimaryHDU(out)
		outfile.header.update('CTYPE1','LINEAR')
		outfile.header.update('CRPIX1',1)
		outfile.header.update('CRVAL1',mswave-(0.5*out2.shape[2])*scale)
		outfile.header.update('CD1_1',scale)
		outfile.header.update('CTYPE2','LINEAR')
		outfile.header.update('CRPIX2',1)
		outfile.header.update('CRVAL2',1)
		outfile.header.update('CD2_2',1)
		if nsci>1:
			outfile.header.update('CTYPE3','LINEAR')
			outfile.header.update('CRPIX3',1)
			outfile.header.update('CRVAL3',1)
			outfile.header.update('CD3_3',1)
		outfile.writeto(strtfile)
		del outfile,out

		try:
			os.remove(bgfile)
		except:
			pass

		outfile = pyfits.PrimaryHDU(out2)
		outfile.header.update('CTYPE1','LINEAR')
		outfile.header.update('CRPIX1',1)
		outfile.header.update('CRVAL1',mswave-(0.5*out2.shape[2])*scale)
		outfile.header.update('CD1_1',scale)
		outfile.header.update('CTYPE2','LINEAR')
		outfile.header.update('CRPIX2',1)
		outfile.header.update('CRVAL2',1)
		outfile.header.update('CD2_2',1)
		outfile.header.update('CTYPE3','LINEAR')
		outfile.header.update('CRPIX3',1)
		outfile.header.update('CRVAL3',1)
		outfile.header.update('CD3_3',1)
		outfile.writeto(bgfile)
		del outfile,out2

	"""
	Loop through all of the slits, determining the wavelength solution and
	  performing the background subtraction. It might be more robust to
	  determine all wavelength solutions, then jointly determine a 'master'
	  solution.... posc stores the current (starting) position of the
	  coadded array, and posn stores the current position of the straight
	  array.
	"""
	posc = 0
	posn = 0
	count = 1

	""" Debugging feature; set to 1 to skip background subtraction """
	lris.lris_red.skysub.RESAMPLE = 0
	""" Extract 1d spectra? """
	do_extract = False

	for k in range(len(slits)):
		i,j = slits[k]
		a,b = wide_slits[k]

		""" Debugging feature; change number to skip initial slits """
		if count<1:
			count += 1
			continue

		print "Working on slit %d (%d to %d)" % (count,i,j)
		# Determine the wavelength solution
		sky2x,sky2y,ccd2wave = wavematch(a,scidata[:,a:b],arc_ycor[i:j],yforw[i:j],widemodel,finemodel,goodmodel,scale,mswave,redcutoff)
		# Resample and background subtract
		print 'Doing background subtraction'
		#scidata[0,a:b] = arcdata[a:b] # This line may be a debugging step that MWA put in.  See what happens with it missing.
		strt,bgsub,varimg = doskysub(i,j-i,outlength,scidata[:,a:b],yback[a:b],sky2x,sky2y,ccd2wave,scale,mswave,center,redcutoff,airmass)

		# Store the resampled 2d spectra
		h = strt.shape[1]
		if cache:
			file = pyfits.open(strtfile,mode="update")
			out = file[0].data
		out[:,posn:posn+h] = strt.copy()
		if cache:
			file.close()
			del file,out
		posn += h+5

		if lris.lris_red.skysub.RESAMPLE:
			count += 1
			continue

		# Store the resampled, background subtracted 2d spectra
		h = bgsub.shape[0]
		if cache:
			file = pyfits.open(bgfile,mode="update")
			out2 = file[0].data
		out2[0,posc:posc+h] = bgsub.copy()
		out2[1,posc:posc+h] = varimg.copy()
		if cache:
			file.close()
			del file,out2
		posc += h+5


		# Find and extract object traces
		if do_extract:
			print '  Extracting object spectra'
			tmp = scipy.where(scipy.isnan(bgsub),0.,bgsub)
			filter = tmp.sum(axis=0)
			mod = scipy.where(filter!=0)
			start = mod[0][0]
			end = mod[0][-1]+1
			del tmp
			slit = bgsub[:,start:end]
			spectra = extract(slit,varimg[:,start:end],extractwidth)
			num = 1
			crval = mswave-(0.5*bgsub.shape[1]-start)*scale
			for spec in spectra:
				for item in spec:
					if item.size==4:
						hdu = pyfits.PrimaryHDU()
						hdu.header.update('CENTER',item[2])
						hdu.header.update('WIDTH',item[3])
						hdulist = pyfits.HDUList([hdu])
					else:
						thdu = pyfits.ImageHDU(item)
						thdu.header.update('CRVAL1',crval)
						thdu.header.update('CD1_1',scale)
						thdu.header.update('CRPIX1',1)
						thdu.header.update('CRVAL2',1)
						thdu.header.update('CD2_2',1)
						thdu.header.update('CRPIX2',1)
						thdu.header.update('CTYPE1','LINEAR')
						hdulist.append(thdu)
					outname = out_prefix+"_spec_%02d_%02d.fits" % (count,num)
					hdulist.writeto(outname)
					num += 1

		count += 1


	""" Output 2d spectra """
	if cache:
		file = pyfits.open(bgfile)
		out2 = file[0].data.copy()
		del file
	tmp = out2[0].copy()
	tmp = scipy.where(scipy.isnan(tmp),0,1)
	mod = scipy.where(tmp.sum(axis=0)!=0)
	start = mod[0][0]
	end = mod[0][-1]+1
	del tmp

	outname = out_prefix+"_bgsub.fits"
	outfile = pyfits.PrimaryHDU(out2[0,:,start:end])
	outfile.header.update('CTYPE1','LINEAR')
	outfile.header.update('CRPIX1',1)
	outfile.header.update('CRVAL1',mswave-(0.5*out2.shape[2]-start)*scale)
	outfile.header.update('CD1_1',scale)
	outfile.header.update('CRPIX2',1)
	outfile.header.update('CRVAL2',1)
	outfile.header.update('CD2_2',1)
	outfile.writeto(outname)
	hdr = outfile.header.copy()

	outname = out_prefix+"_var.fits"
	outfile = pyfits.PrimaryHDU(out2[1,:,start:end])
	outfile.header=hdr
	outfile.writeto(outname)
	del out2,hdr

	if cache:
		file = pyfits.open(strtfile)
		out = file[0].data.copy()
		del file
	for i in range(nsci):
		outname = out_prefix+"_straight_%d.fits" % (i+1)
		outfile = pyfits.PrimaryHDU(out[i,:,start:end])
		outfile.header.update('CTYPE1','LINEAR')
		outfile.header.update('CRPIX1',1)
		outfile.header.update('CRVAL1',mswave-(0.5*out.shape[2]-start)*scale)
		outfile.header.update('CD1_1',scale)
		outfile.header.update('CRPIX2',1)
		outfile.header.update('CRVAL2',1)
		outfile.header.update('CD2_2',1)
		#if nsci>1:
		#	outfile.header.update('CRPIX3',1)
		#	outfile.header.update('CRVAL3',1)
		#	outfile.header.update('CD3_3',1)	
		outfile.writeto(outname)
		del outfile

	del out
