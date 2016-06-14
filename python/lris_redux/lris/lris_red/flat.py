"""
Module to process the LRIS redside flat files. flatpipe() processes new files
  and flatload() loads pre-existing files.
"""

import scipy,pickle
from astropy.io import fits as pyfits

def flatpipe(flatfiles,out_prefix):
	"""
	flatpipe(flatfiles,out_prefix)

	Processes the spectroscopic flats by:
	  coadding
	  determining y-distortion
	  finding science slits and star boxes
	  creating normalized response flat

	Input:
	  flatfiles  - a list of the files to be processed
	  out_prefix - the output prefix of the distortion maps, normalized
	                 flat, and slit descriptors

	Output:
	  pickled description of the slits, star boxes, and distortion
	    coefficients
	  distortion maps
	  normalized flats
	"""

	from mostools import ycorrect,id_slits,spectools
	from lris.lris_biastrim import redside as biastrim
	from special_functions import lsqfit,genfunc
	from scipy import ndimage,stats

	# Read data from files, biastrim and coadd
	weight = 0
	for name in flatfiles:
		flattmp = pyfits.open(name)
		flatdata = biastrim(flattmp[0].data.astype(scipy.float32))
		exptime = flattmp[0].header['elaptime']

		if weight==0:
			flatavg = flatdata
		else:
			flatavg += flatdata
		weight += exptime
		del flatdata
	del flattmp
	flatavg /= weight

	# Find the y-correction
	ytrue,ymap = ycorrect.ycorrect(flatavg)

	# Generate the distortion maps
	coords = spectools.array_coords(flatavg.shape)
	x = coords[1].flatten()
	y = coords[0].flatten()
	yforw = genfunc(x,y,ytrue).reshape(coords[0].shape)
	yback = genfunc(x,y,ymap).reshape(coords[0].shape)
	del coords,x,y

	# Create straight image to find slits/normalization
	flat_ycor = spectools.resampley(flatavg,yforw,cval=0.)
	del flatavg

	# Identify slits and star boxes using brighter columns for each
	#   (straightened) row
	search = scipy.sort(flat_ycor,1)
	width = search.shape[1]
	slits,starboxes = id_slits.id_slits(search[:,width*2/3:width*9/10])

	print "Starbox Locations..."
	for i,j in starboxes:
		print "[:,%d:%d]" % (i,j)

	print "Slit boundaries defined by..."
	for i,j in slits:
		print "[:,%d:%d]" % (i,j)

	# Normalize the straightened flatfield
	flatmodel = scipy.ones(flat_ycor.shape,scipy.float32)
	for i,j in slits:
		data = flat_ycor[i:j,:].copy()
		data[data==0] = scipy.nan
		slice = stats.stats.nanmedian(data,axis=0)
		slice[scipy.isnan(slice)] = stats.stats.nanmedian(slice) 
		norm = ndimage.gaussian_filter1d(slice,23.)
		norm[norm<=0] = 1.
		top = i-3
		bottom = j+3
		if top<0:
			top = 0
		if bottom>flatmodel.shape[0]:
			bottom = flatmodel.shape[0]
		flatmodel[top:bottom,:] = flat_ycor[top:bottom,:]/norm
	del flat_ycor

	# Create normalized flatfield by resampling to curved frame
	flatnorm = spectools.resampley(flatmodel,yback,cval=1.)
	del flatmodel

	# Output normalized flat
	flatname = out_prefix+"_flat.fits"
	pyfits.PrimaryHDU(flatnorm.astype(scipy.float32)).writeto(flatname)

	# Output coefficient/slit definition arrays
	outname = out_prefix+"_yforw.fits"
	pyfits.PrimaryHDU(yforw.astype(scipy.float32)).writeto(outname)
	outname = out_prefix+"_yback.fits"
	pyfits.PrimaryHDU(yback.astype(scipy.float32)).writeto(outname)
	outname = out_prefix+"_ygeom.dat"
	outfile = open(outname,"w")
	pickle.dump(slits,outfile)
	pickle.dump(starboxes,outfile)
	pickle.dump(ytrue,outfile)
	pickle.dump(ymap,outfile)
	outfile.close()

	return yforw.astype(scipy.float32),yback.astype(scipy.float32),slits,starboxes,flatnorm.astype(scipy.float32)



def flatload(out_prefix):
	"""
	Loads the outputs of flatpipe() from a previous run.
	"""
	# Open FITS flatfile
	flatname = out_prefix+"_flat.fits"
	flattmp = pyfits.open(flatname)
	flatnorm = flattmp[0].data.astype(scipy.float32)
	del flattmp

	# Get distortion arrays
	inname = out_prefix+"_yforw.fits"
	tmp = pyfits.open(inname)
	yforw = tmp[0].data.astype(scipy.float32)
	inname = out_prefix+"_yback.fits"
	tmp = pyfits.open(inname)
	yback = tmp[0].data.astype(scipy.float32)
	del tmp
	# Open pickle data
	inname = out_prefix+"_ygeom.dat"
	infile = open(inname,"r")
	slits = pickle.load(infile)
	starboxes = pickle.load(infile)
	infile.close()

	return yforw,yback,slits,starboxes,flatnorm
