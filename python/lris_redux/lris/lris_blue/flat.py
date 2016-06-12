"""
Module to process the LRIS blue side flat files. flatpipe() processes new files
  and flatload() loads pre-existing files.

NOTES:	The pipeline treats the two CCDs as independent (since they are...).
	The flat fields are broken (there are too many reflections on the blue
	  side to make them worthwhile), so no normalization is performed.
"""

import scipy,pyfits,pickle

def flatpipe(flatfiles,out_prefix):
	"""
	flatpipe(flatfiles,out_prefix)

	Processes the spectroscopic flats by:
	  coadding
	  determining y-distortion
	  finding science slits and star boxes

	Input:
	  flatfiles  - a list of the files to be processed
	  out_prefix - the output prefix of the distortion maps, normalized
			 flat, and slit descriptors

	Output:
	  pickled description of the slits, star boxes, and distortion
	    coefficients
	  distortion maps
	"""
	from lris.lris_biastrim import blueside as biastrim
	from mostools import ycorrect,id_slits,spectools
	from special_functions import lsqfit,genfunc


	# Read data from files, biastrim and coadd
	weight = 0
	for name in flatfiles:
		flattmp = pyfits.open(name)
		flatdata = biastrim(flattmp[0].data.copy())
		exptime = flattmp[0].header['elaptime']

		if weight==0:
			flatavg = flatdata
		else:
			flatavg += flatdata
		weight += exptime
		del flatdata
	del flattmp
	flatavg /= weight

	# The center of the mask...
	YMID = 2048


	# Find the y-correction
	coords = spectools.array_coords(flatavg[:YMID].shape)
	x = coords[1].flatten()
	y = coords[0].flatten()
	ytrue = {}
	ymap = {}
	yforw = {}
	yback = {}

	# Bottom
	a,b = ycorrect.ycorrect(flatavg[:YMID])
	yforw['bottom'] = genfunc(x,y,a).reshape(coords[0].shape)
	yback['bottom'] = genfunc(x,y,b).reshape(coords[0].shape)
	ytrue['bottom'] = a
	ymap['bottom'] = b
	# Top
	a,b = ycorrect.ycorrect(flatavg[YMID:])
	yforw['top'] = genfunc(x,y,a).reshape(coords[0].shape)
	yback['top'] = genfunc(x,y,b).reshape(coords[0].shape)
	ytrue['top'] = a
	ymap['top'] = b
	del coords,x,y

	# Create straight image to find slits
	flat_ycor = {}
	flat_ycor['bottom'] = spectools.resampley(flatavg[:YMID],yforw['bottom'],cval=0.)
	flat_ycor['top'] = spectools.resampley(flatavg[YMID:],yforw['top'],cval=0.)
	del flatavg

	slits = {}
	starboxes = {}
	for i in ['bottom','top']:
		add = 0
		if i=='top':
			add = YMID
		# Identify slits and star boxes using brighter columns for each
		#   (straightened) row
		search = scipy.sort(flat_ycor[i],1)
		del flat_ycor[i]
		width = search.shape[1]
		slits[i],starboxes[i] = id_slits.id_slits(search[:,width*7/10:width*9/10])
		print "Starbox Locations for %s:" % i
		for a,b in starboxes[i]:
			print "[:,%d:%d]" % (a+add,b+add)
		print ""
		print "Slit Locations for %s:" % i
		for a,b in slits[i]:
			print "[:,%d:%d]" % (a+add,b+add)
		print ""

		# Output coefficient arrays
		outname = out_prefix+"_yforw_%s.fits" % i
		pyfits.PrimaryHDU(yforw[i]).writeto(outname)
		outname = out_prefix+"_yback_%s.fits" % i
		pyfits.PrimaryHDU(yback[i]).writeto(outname)

	# Output distortion solutions and slit definitions
	outname = out_prefix+"_ygeom.dat"
	outfile = open(outname,"w")
	pickle.dump(slits,outfile)
	pickle.dump(starboxes,outfile)
	pickle.dump(ytrue,outfile)
	pickle.dump(ymap,outfile)
	outfile.close()

	return yforw,yback,slits,starboxes

# For flatfield reductions that have already been performed.
def flatload(out_prefix):
	"""
	Loads the outputs of flatpipe() from a previous run.
	"""

	# Get distortion arrays
	yforw = {}
	yback = {}
	for i in ['bottom','top']:
		inname = out_prefix+"_yforw_%s.fits" % i
		tmp = pyfits.open(inname)
		yforw[i] = tmp[0].data.astype(scipy.float32)
		inname = out_prefix+"_yback_%s.fits" % i
		tmp = pyfits.open(inname)
		yback[i] = tmp[0].data.astype(scipy.float32)
	# Open pickle data
	inname = out_prefix+"_ygeom.dat"
	infile = open(inname,"r")
	slits = pickle.load(infile)
	starboxes = pickle.load(infile)
	infile.close()

	return yforw,yback,slits,starboxes
