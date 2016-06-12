"""
lris_pipeline(prefix,dir,science,arc,flats,out_prefix,useflat,usearc,cache,offsets)

Pipeline to reduce LRIS red or blueside spectra. Automatically performs almost
  *all* operations, including: removing the instrumental signature (bias,
  overscan), determining the distortion solution, finding the slit locations,
  determining the wavelength solution, 'optimal' background subtraction,
  resampling, and source extraction.

Inputs:
  prefix    - Filename prefix, ie "lred" or "lblue"
  dir       - Directory with input files, ie "raw/"
  science   - Numbers of the science files, ie "0023,0024,0029,0030"
  arc       - Number of arc file, ie "0025"
  flats     - Numbers of flat files, ie "0026,0027,0028"
  outprefix - Prefix for output files, ie "mask2_11oct04"
  useflat   - 1 to use flat data from a previous run
  usearc    - 1 to use arc data from a previous run
  cache     - 1 to cache data to disk (useful for blueside with RAM<2GB)
  offsets   - a list/array of relative y-offsets between masks (in pixels)

Outputs:
  straightened, wavelength calibrated, cosmic-ray cleaned 2d spectra
  straightened, wavelength calibrated, cosmic-ray cleaned, background
    subtracted, coadded 2d spectra
  2d variance spectra
  1d extracted spectra, including a smoothed spectrum and variance spectrum
"""

import pyfits

def lris_pipeline(prefix,dir,science,arc,flats,out_prefix,useflat=0,usearc=0,cache=0,offsets=None):
	""" Batch files will have a prefix """
	if prefix is not None:
		arcname = dir+prefix+arc+".fits"

		scinames = []
		flatnames = []

		for flatn in flats.split(","):
			flatnames.append(dir + prefix + flatn + ".fits")
		for scin in science.split(","):
			scinames.append(dir + prefix + scin + ".fits")
	else:
		arcname = arc
		scinames = science
		flatnames = flats

	"""
	We use the arclamp file to determine the name of the slitmask and which
	  instrument (redside or blueside) we are using.
	"""
	try:
		instrument = pyfits.getval(arcname,'INSTRUME')
	except:
		instrument = "LRISOLD"
	maskname = pyfits.getval(arcname,'SLITNAME')


	"""
	Ensure that all files are for the same slitmask and instrument. If they
	  do not match the arc slitmask/instrument name, the program prints a
	  warning and exits.
	"""
	for name in flatnames:
		if instrument!="LRISOLD":
			inst = pyfits.getval(name,'INSTRUME')
			if inst!=instrument:
				print "Warning: Instrument mismatch in file %s" % name
				return None
		mask = pyfits.getval(name,'SLITNAME')
		if mask!=maskname:
			print "Warning: Mask mismatch in file %s" % name
			return None

	for name in scinames:
		if instrument!="LRISOLD":
			inst = pyfits.getval(name,'INSTRUME')
			if inst!=instrument:
				print "Warning: Instrument mismatch in file %s" % name
				return
			mask = pyfits.getval(name,'SLITNAME')
			if mask!=maskname:
				print "Warning: Mask mismatch in file %s" % name
				return


	"""
	Choose the correct pipeline to import, depending on whether we have red
	  or blue side data.
	"""
	if instrument=="LRISBLUE":
		from lris.lris_blue.lris_blue_pipeline import lris_pipeline as pipeline
	else:
		from lris.lris_red.lris_red_pipeline import lris_pipeline as pipeline

	pipeline(prefix,dir,scinames,arcname,flatnames,out_prefix,useflat,usearc,cache,offsets)
