"""
Tools for creating arc models.
"""

import special_functions as sf

import scipy
from scipy import io as sio
from math import sqrt,pi

def make_arc(filename,dispersion,wave):
	"""
	make_arc(filename,dispersion,wave)

	Create models of the arclamps used.

	Inputs:
	  filename   - name of file containing arc wavelengths and amplitudes
	  dispersion - width of arcs
	  wave       - wavelength scale on which to evaluate arcs

	Outputs:
	  resolution-matched arc spectrum and 3x broadened arc spectrum
	"""
	f = open(filename)
	lines = sio.read_array(f)
	f.close()

	n = lines.shape[0]
	fitmodel = scipy.zeros(3*n+1)
	fitmodel2 = scipy.zeros(3*n+1)
	index = 1
	for i in range(n):
		fitmodel[index] = lines[i,1]
		fitmodel[index+1] = lines[i,0]
		fitmodel[index+2] = dispersion
		fitmodel2[index] = 1.
		fitmodel2[index+1] = lines[i,0]
		fitmodel2[index+2] = dispersion*3.
		index += 3
	return sf.ngauss(wave,fitmodel),sf.ngauss(wave,fitmodel2)


def make_linelist(files,cutoff,fluxlimit,linefile):
	"""
	make_linelist(files,cutoff,fluxlimit,linefile)

	Cleans/combines linelists.

	Inputs:
	  files     - list containing the names of the arc files
	  cutoff    - maximum wavelength
	  fluxlimit - minimum flux
	  linefile  - name of output file

	Outputs:
	  outputs the data for 'good' arc lines to a single file.
	"""
	outfile = open(linefile,'w')
	if fluxlimit is None:
		fluxlimit = 0.
	for name in files:
		f = open(name)
		lines = sio.read_array(f)
		f.close()
		for i,j in lines:
			if i>cutoff or j<fluxlimit:
				continue
			else:
				outfile.write("%8.3f    %8.3f\n" % (i,j))
	outfile.close()


def response(wave,file1,file2):
	""" Not used. """
	from scipy import interpolate
	f = open(file1)
	resp1 = sio.read_array(f)
	f.close()
	x = resp1[:,0]
	z = resp1[:,1]
	spline = interpolate.splrep(x,z,s=0)
	resp1 = interpolate.splev(wave,spline)

	f = open(file2)
	resp2 = sio.read_array(f)
	f.close()
	x = resp2[:,0]
	z = resp2[:,1]
	spline = interpolate.splrep(x,z,s=0)
	resp2 = interpolate.splev(wave,spline)

	return resp1/resp2
