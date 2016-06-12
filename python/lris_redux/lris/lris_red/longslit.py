"""
This code doesn't work (yet).
"""

import lris
from lris.lris_biastrim import redside as biastrim

import lris_red
from lris_red.flat import *
from lris_red.skymatch import skymatch as wavematch
from lris_red.skysub import doskysub

from spectra import spectools,offset,measure_width
from spectra.extract import extract
import special_functions

from pickle import dump,load
from math import ceil,fabs
import os

import scipy,pyfits
from scipy import stats,interpolate,ndimage
from scipy import io as sio

# Pipeline to reduce LRIS red or blueside spectra
# Inputs:
#   prefix	- Filename prefix, ie "lred" or "lblue"
#   dir		- Directory with input files, ie "raw/"
#   science     - Numbers of the science files, ie "0023,0024,0029,0030"
#   arc		- Number of arc file, ie "0025"
#   flats	- Numbers of flat files, ie "0026,0027,0028"
#   outprefix   - Prefix for output files, ie "mask2_11oct04"
#   useflat	- 1 to use flat data from a previous run
#   usearc	- 1 to use arc data from a previous run
#   cache	- 1 to cache data to disk (useful for blueside with RAM<2GB)
#   offsets     - a list/array of relative offsets between masks (in pixels)

def lris_pipeline(prefix,dir,science,arc,flats,out_prefix,useflat=0,usearc=0,cache=0,offsets=None):
	print "Processing mask",out_prefix

	scinums = science.split(",")
	flatnums = flats.split(",")

	for i in range(len(flatnums)):
		flatnums[i] = dir + prefix + flatnums[i] + ".fits"
	scinames = []
	for i in range(len(scinums)):
		name = dir + prefix + scinums[i] + ".fits"
		scinames.append(name)
	arcname = dir + prefix + arc + ".fits"

	nsci = len(scinums)


	print "Preparing flatfields"
	if useflat==1:
		yforw,yback,slits,starboxes,flatnorm = flatload(out_prefix)
	else:
		yforw,yback,slits,starboxes,flatnorm = flatpipe(flatnums,out_prefix)
	axis1 = flatnorm.shape[0]
	axis2 = flatnorm.shape[1]

	print "Preparing arcs for line identification"
	if usearc==1:
		arcname = out_prefix+"_arc.fits"
		arc_tmp = pyfits.open(arcname)
		arc_ycor = arc_tmp[0].data.astype(scipy.float32)
		lamps = arc_tmp[0].header['LAMPS']
		del arc_tmp
	else:
		arcname = dir + prefix + arc + ".fits"
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

		#Remove screwed columns (this should already be done though...)
		bad = scipy.where(scidatatmp>56000.)
		nbad = bad[0].size
		for k in range(nbad):
			y = bad[0][k]
			x = bad[1][k]
			scidatatmp[y,x] = (scidatatmp[y,x-1]+scidatatmp[y,x+1])/2.
		# Don't flatfield blueside data
		scidatatmp = scidatatmp/flatnorm
		scidata[i,:,:] = scidatatmp.copy()

		try:
			mswave = scitmp[0].header['MSWAVE']
		except:
			mswave = 6500.
		if len(slits)==1:
			try:
				mswave = scitmp[0].header['WAVELEN']
			except:
				pass
		disperser = scitmp[0].header['GRANAME']
		airmass.append(scitmp[0].header['AIRMASS'])

		# Old data mightn't have a dichroic keyword!
		try:
			dichroic = scitmp[0].header['DICHNAME']
		except:
			dichroic = None

		flux[i] = scipy.sort(scipy.ravel(scidatatmp))[scidatatmp.size/4]
		for j in range(len(starboxes)):
			a,b = starboxes[j]
			m,n = wide_stars[j]
			a -= 4
			b += 4
			m -= 2
			n += 2
			center[i,j] = offset.findoffset(scidatatmp[m:n],yforw[a:b],m)

		del scitmp
		del scidatatmp
	del flatnorm

	if offsets is not None:
		center = scipy.asarray(offsets)
	else:
		center = stats.stats.nanmean(center,axis=1)

	center[scipy.isnan(center)] = 0.
	print "Normalizing Fluxes"
	cmax = center.max()
	fmax = flux.max()
	for i in range(center.size):
		center[i] -= cmax
		ratio = fmax/flux[i]
		scidata[i] *= ratio
	cmax = ceil(fabs(center.min()))

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
		redcutoff = 4600.
		dich_file = '460'
	elif dichroic=='500':
		redcutoff = 5000.
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

	nsize = 0
	csize = 0
	wide_slits = []
	linewidth = []
	slits = [[1150,1250]]
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
		if len(wide_slits)%7==0 or len(slits)==1:
			linewidth.append(measure_width.measure(arc_ycor[(i+j)/2,:],15))
	csize -= 5
	nsize -= 5

	linewidth = scipy.median(scipy.asarray(linewidth))


	print "Loading wavelength model"
	lris_path = lris_red.__path__[0]

	filename = lris_path+"/uves_sky.model"
	infile = open(filename,"r")
	wavecalmodel = load(infile)
	infile.close()
	wave = scipy.arange(3400.,10400.,0.1)

	if dich_file!='':
		filename = lris_path+"/dichroics/dichroic_"+dich_file+"_t.dat"
		infile = open(filename,"r")
		input = sio.read_array(infile)
		infile.close()
		spline = interpolate.splrep(input[:,0],input[:,1])
		dich = interpolate.splev(wave,spline)
		dich[wave<4500.] = 1.
		dich[wave>8800.] = 1.
		del input,spline
	else:
		dich = scipy.ones(wave.size)
	wavemodel = interpolate.splev(wave,wavecalmodel)
	finemodel = ndimage.gaussian_filter1d(wavemodel,linewidth*scale/0.1)
	wavemodel = ndimage.gaussian_filter1d(finemodel,5./0.1)
	finemodel *= dich
	finemodel = interpolate.splrep(wave,finemodel)
	wavemodel *= dich
	widemodel = interpolate.splrep(wave,wavemodel)
	goodmodel = finemodel
	del dich,wave,wavemodel

	extractwidth = 10

	print "Creating output arrays"
	outlength = int(axis2*1.6)
	out = scipy.zeros((nsci,nsize,outlength),scipy.float32)*scipy.nan
	out2 = scipy.zeros((2,csize,outlength),scipy.float32)*scipy.nan

	if cache:
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

	posc = 0
	posn = 0
	count = 1
	for k in range(len(slits)):
		i,j = slits[k]
		a,b = wide_slits[k]
		##
		if count<1:
			count += 1
			continue
		##
		print "Working on slit %d (%d to %d)" % (count,i,j)
		sky2x,sky2y,ccd2wave = wavematch(a,scidata[:,a:b],arc_ycor[i:j],yforw[i:j],widemodel,finemodel,goodmodel,scale,mswave,redcutoff)

		strt,bgsub,varimg = doskysub(i,j-i,outlength,scidata[:,a:b],yback[a:b],sky2x,sky2y,ccd2wave,scale,mswave,center,redcutoff,airmass)
		
		h = strt.shape[1]
		if cache:
			file = pyfits.open(strtfile,mode="update")
			out = file[0].data
		out[:,posn:posn+h] = strt.copy()
		if cache:
			file.close()
			del file,out
		posn += h+5

##
#		lris_red.skysub.RESAMPLE = 1
#		count += 1
#		continue
##


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
##
#		count += 1
#		continue
##
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

##
#	file = pyfits.open(bgfile)
#	file.writeto(out_prefix+"_save.fits")
#	return
##

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
	outname = out_prefix+"_straight.fits"
	outfile = pyfits.PrimaryHDU(out[:,:,start:end])
	outfile.header.update('CTYPE1','LINEAR')
	outfile.header.update('CRPIX1',1)
	outfile.header.update('CRVAL1',mswave-(0.5*out.shape[2]-start)*scale)
	outfile.header.update('CD1_1',scale)
	outfile.header.update('CRPIX2',1)
	outfile.header.update('CRVAL2',1)
	outfile.header.update('CD2_2',1)
	if nsci>1:
		outfile.header.update('CRPIX3',1)
		outfile.header.update('CRVAL3',1)
		outfile.header.update('CD3_3',1)	
	outfile.writeto(outname)

	del out,outfile
