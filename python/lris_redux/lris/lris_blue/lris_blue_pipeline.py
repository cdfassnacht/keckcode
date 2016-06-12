"""
Pipeline to reduce LRIS blue-side spectra
Inputs:
  prefix      - Filename prefix, ie "lblue"
  dir	 - Directory with input files, ie "raw/"
  science     - Numbers of the science files, ie "0023,0024,0029,0030"
  arc	 - Number of arc file, ie "0025"
  flats       - Numbers of flat files, ie "0026,0027,0028"
  outprefix   - Prefix for output files, ie "mask2_11oct04"
  useflat     - 1 to use flat data from a previous run
  usearc      - 1 to use arc data from a previous run
  cache       - 1 to cache data to disk (useful for blueside with RAM
  offsets     - a list/array of relative offsets between masks (in pixels)
  logfile     - name of the output logfile (out_prefix.log is used otherwise)
"""

import lris
from lris.lris_biastrim import blueside as biastrim

from lris.lris_blue.flat import *
from lris.lris_blue.skysub import doskysub
from lris.lris_blue.arc import make_arc,make_linelist

from mostools import spectools,offset,measure_width
from mostools.extract import extract
import special_functions

from math import ceil,fabs
import time,pickle

import scipy,pyfits
from scipy import stats,interpolate,ndimage
from scipy import io as sio

"""
The top and the bottom of the mask are treated as separate detectors (since
  they *are* separate detectors). The following are helper functions to load
  the arc image or the distortion solution for the top or bottom. This allows
  us to not have to store both top and bottom in memory at the same time.
"""
def get_arc(out_prefix,side='top'):
	name = out_prefix+"_arc_%s.fits" % side
	return pyfits.open(name)[0].data.astype(scipy.float32)

def get_yforw(out_prefix,side='top'):
	name = out_prefix+"_yforw_%s.fits" % side
	return pyfits.open(name)[0].data.astype(scipy.float32)

def get_yback(out_prefix,side='top'):
	name = out_prefix+"_yback_%s.fits" % side
	return pyfits.open(name)[0].data.astype(scipy.float32)

"""
Main pipeline. The blueside currently includes logging.
"""
def lris_pipeline(prefix,dir,scinames,arcname,flatnames,out_prefix,useflat=0,usearc=0,cache=0,offsets=None,logfile=None):
	# Create a logfile for this session
	if logfile is None:
		logfile = open('%s.log' % out_prefix,'w')
	else:
		logfile = open(logfile,'w')
	stime = time.strftime("%d/%m/%y %H:%M:%S")
	logfile.write('%s\n' % stime)

	print "Processing mask",out_prefix
	logfile.write('Processing mask %s\n' % out_prefix)


	""" Prepare image names. """

	nsci = len(scinames)
	YMID = 2048  # offset for the second detector

	print "Preparing flatfields"
	if useflat==1:
		logfile.write('Using pre-made flats\n')
		yforw,yback,slits,starboxes = flatload(out_prefix)
	else:
		logfile.write('Making new flats\n')
		yforw,yback,slits,starboxes = flatpipe(flatnames,out_prefix)


	print "Preparing arcs for line identification"
	if usearc==1:
		logfile.write('Using pre-made arcs\n')
		arc_ycor = {}
		for i in ['bottom','top']:
			arcname = out_prefix+"_arc_%s.fits" % i
			arc_tmp = pyfits.open(arcname)
			arc_ycor[i] = arc_tmp[0].data.astype(scipy.float32)
			lamps = arc_tmp[0].header['LAMPS']
			filter = arc_tmp[0].header['BLUFILT']
			del arc_tmp
	else:
		logfile.write('Making new arcs\n')
		""" Load arc data from fits file """
		arc_tmp = pyfits.open(arcname)
		arcdata = arc_tmp[0].data.copy()

		""" Determine which lamps were used """
		lamps = arc_tmp[0].header['LAMPS']
		try:
			filter = arc_tmp[0].header['BLUFILT']
		except:
			filter = 'clear'
		del arc_tmp

		""" Process arcs for the top and bottom separately """
		arcdata = biastrim(arcdata)
		arc_ycor = {}
		arc_ycor['bottom'] = spectools.resampley(arcdata[:YMID],yforw['bottom']).astype(scipy.float32)
		arcname = out_prefix+"_arc_bottom.fits"
		arc_hdu = pyfits.PrimaryHDU(arc_ycor['bottom'])
		arc_hdu.header.update('LAMPS',lamps)
		arc_hdu.header.update('BLUFILT',filter)
		arc_hdu.writeto(arcname)

		arc_ycor['top'] = spectools.resampley(arcdata[YMID:],yforw['top']).astype(scipy.float32)
		arcname = out_prefix+"_arc_top.fits"
		arc_hdu = pyfits.PrimaryHDU(arc_ycor['top'])
		arc_hdu.header.update('LAMPS',lamps)
		arc_hdu.header.update('BLUFILT',filter)
		arc_hdu.writeto(arcname)

		del arc_hdu,arcdata

	axis1 = 4096
	axis2 = 4096


	"""
	We create 'wide' starboxes that describe the minimum and maximum
	  y-position of the slit in the *unstraightened* frame.
	"""
	wide_stars = {}
	logfile.write('\n')
	for l in ['bottom','top']:
		wide_stars[l] = []
		bmax = arc_ycor[l].shape[0]
		logfile.write('Star boxes for %s\n' % l)
		for i,j in starboxes[l]:
			logfile.write('[:,%d:%d]\n' % (i,j))

			mod = scipy.where((yback[l]<j)&(yback[l]>i))
			a = mod[0].min()-3	# We include a small amount of
			b = mod[0].max()+3	#  padding for resampling.
			if a<0:
				a = 0
			if b>bmax:
				b = bmax
			wide_stars[l].append([a,b])

	print "Bias trimming science data"
	nstars = len(starboxes['bottom'])+len(starboxes['top'])
	scidata = scipy.zeros((nsci,axis1,axis2),'f4')
	center = scipy.zeros((nsci,nstars),'f4')
	flux = scipy.zeros((nsci),'f4')
	for i in range(nsci):
		filename = scinames[i]
		scitmp = pyfits.open(filename)

		scidatatmp = scitmp[0].data.copy()
		scidatatmp = biastrim(scidatatmp).astype(scipy.float32)

		"""
		Remove screwed columns (this should already be done by biastrim
		  though...).
		"""
		bad = scipy.where(scidatatmp>56000.)
		nbad = bad[0].size
		for k in range(nbad):
			y = bad[0][k]
			x = bad[1][k]
			scidatatmp[y,x] = (scidatatmp[y,x-1]+scidatatmp[y,x+1])/2.
		"""
		We don't flatfield blueside data because of ghosts and
		  reflections. Milan Bogosavljevic has data that show that
		  flatfielding is important if using very blue data--the bluest
		  end looks like it has fringes! I should add a flag to allow
		  flatfielding to be turned on....
		"""
		scidata[i,:,:] = scidatatmp.copy()

		"""
		The try/except code is because sometimes the data just don't
		  have these keywords (I think this might be correlated to
		  stopping exposures early, though I'm not sure). Plus old data
		  might not have used the dichroic at all....
		"""
		try:
			disperser = scitmp[0].header['GRISNAME']
		except:
			pass
		try:
			dichroic = scitmp[0].header['DICHNAME']
		except:
			dichroic = None

		"""
		We use the first quartile for the flux normalization.
		"""
		flux[i] = scipy.sort(scipy.ravel(scidatatmp))[scidatatmp.size/4]

		"""
		The starboxes are used to determine relative y-shifts between
		  mask exposures. If the offsets keyword was used, this will
		  be ignored.
		"""
		for l in ['bottom','top']:
			if offsets is not None:
				continue
			for j in range(len(starboxes[l])):
				a,b = starboxes[l][j]
				m,n = wide_stars[l][j]
				a -= 4
				b += 4
				m -= 2
				n += 2
				if a<0:
					a = 0
				if l=='top':
					m += YMID
					n += YMID
				center[i,j] = offset.findoffset(scidatatmp[m:n],yforw[l][a:b],m)

		del scitmp
		del scidatatmp

	if offsets is not None:
		center = scipy.asarray(offsets)
	else:
		center = stats.stats.nanmean(center,axis=1)
	center[scipy.isnan(center)] = 0.

	print "Normalizing Fluxes"
	cmax = center.max()
	fmax = flux.max()
	logfile.write('\nMask pixel and flux offsets\n')
	logfile.write('-----------------------------\n')
	logfile.write('Mask   Pixels   Flux\n')
	for i in range(center.size):
		center[i] -= cmax
		ratio = fmax/flux[i]
		scidata[i] *= ratio
		logfile.write('%4d   %6.2f   %4.2f\n' % (i,center[i],ratio))
	cmax = ceil(fabs(center.min()))
	logfile.write('\n')

	if disperser=="300/5000":
		scale = 1.41
		mswave = 5135.
	elif disperser=="400/3400":
		scale = 1.05
		mswave = 3990.
	elif disperser=="600/4000":
		scale = 0.63
		mswave = 4590.
	elif disperser=="1200/3400":
		scale = 0.24
		mswave = 3505.

	if dichroic=='mirror':
		bluecutoff = 0.
		dich_file = ''
	elif dichroic=='460':
		bluecutoff = 4650.
		dich_file = '460'
	elif dichroic=='500':
		bluecutoff = 5100.
		dich_file = '500'
	elif dichroic=='560':
		bluecutoff = 5650.
		dich_file = '560'
	elif dichroic=='680':
#		bluecutoff = 6800.
		bluecutoff = 5650.
		dich_file = '680'
	else:
		bluecutoff = 8000.
		dich_file = ''


	"""
	We create 'wide' slits that describe the minimum and maximum y-position
	  of the slit in the *unstraightened* frame. We also determine the mask
	  resolution using every seventh slit.
	"""
	nsize = 0  # Size of straightened mask
	csize = 0  # Size of coadded mask
	wide_slits = {}
	linewidth = []
	print "Determining mask resolution"
	for l in ['bottom','top']:
		wide_slits[l] = []
		bmax = arc_ycor[l].shape[0]
		logfile.write('Slits for %s (%d total)\n' % (l,len(slits[l])))
		for i,j in slits[l]:
			logfile.write('[:,%d:%d]\n' % (i,j))
			csize += int(j-i+cmax) + 5
			nsize += j-i+5
			mod = scipy.where((yback[l]>i)&(yback[l]<j))
			a = mod[0].min()-4
			b = mod[0].max()+4
			if a<0:
				a = 0
			if b>bmax:
				b = bmax
			wide_slits[l].append([a,b])
			if len(wide_slits[l])%7==0:
				linewidth.append(measure_width.measure(arc_ycor[l][(i+j)/2,:],15))
	csize -= 5
	nsize -= 5
	logfile.write("\n\n")
	logfile.close()

	linewidth = scipy.median(scipy.asarray(linewidth))


	""" We can temporarily delete the top CCD from memory """
	yforw = yforw['bottom']
	yback = yback['bottom']
	arc_ycor = arc_ycor['bottom']


	"""
	Create the arclamp model by using only the blue lamps that were turned
	  on. Turning off the red lamps (Ne and Ar) reduces reflections on the
	  blue side, and it is difficult to use these lines for the wavelength
	  solution because 2nd order blue lines show up starting at ~5800.
	"""
	print "Loading wavelength model"
	lris_path = lris.__path__[0]
	lamps = lamps.split(',')
	wave = scipy.arange(2000.,8000.,0.1)
	filenames = []
	if lamps[0]=='1':
		filenames.append(lris_path+"/data/bluearcs/hg.dat")
	if lamps[3]=='1':
		filenames.append(lris_path+"/data/bluearcs/cd.dat")
	if lamps[4]=='1':
		filenames.append(lris_path+"/data/bluearcs/zn.dat")
	fluxlimit = None
	if filter=='SP580' and bluecutoff>5650:
		cutoff = 5650. 
	else:
		fluxlimit = 150.
		cutoff = bluecutoff
	linefile = out_prefix+"_lines.dat"
	make_linelist(filenames,cutoff,fluxlimit,linefile)

	"""
	The relative amplitudes of the lines in the hg, cd, and zn.dat files
	  are more appropriate for the bluer grisms. A separate linelist is
	  used for the 300 grism and assumes all three lamps were on. This is
	  one of the problems with (1) not knowing the throughput for each
	  setup and (2) not having stable lamps.
	"""
	if disperser=="300/5000":
		filename = lris_path+"/data/bluearcs/300_lines.dat"
		arc,lines = make_arc(filename,linewidth*scale,wave)
	else:
		arc,lines = make_arc(linefile,linewidth*scale,wave)
	finemodel = interpolate.splrep(wave,arc,s=0)
	smooth = ndimage.gaussian_filter1d(arc,9./0.1)
	widemodel = interpolate.splrep(wave,smooth,s=0)
	linemodel = interpolate.splrep(wave,lines,s=0)

	filename = lris_path+"/data/uves_sky.model"
	infile = open(filename,"r")
	wavecalmodel = pickle.load(infile)
	infile.close()
	wave = scipy.arange(3400.,10400.,0.1)
	""" We attempt to model the dichroic cutoff for the sky mode. """
	if dichroic=='680' and disperser=="300/5000":
		filename = lris_path+"/data/dichroics/dichroic_680_t.dat"
		infile = open(filename,"r")
		input = sio.read_array(infile)
		infile.close()
		input[:,1] = 1. - input[:,1]
		spline = interpolate.splrep(input[:,0],input[:,1],s=0)
		dich = interpolate.splev(wave,spline)
		dich[wave<4500.] = 1.
		dich[wave>8800.] = 1.
		filename = lris_path+"/data/grisms/grism_300.dat"
		infile = open(filename,"r")
		input = sio.read_array(infile)
		infile.close()
		spline = interpolate.splrep(input[:,0],input[:,1],s=0)
		eff = interpolate.splev(wave,spline)
		eff[wave<5100.] = 1.
		eff[wave>7200.] = 1.
		dich *= eff
		del input,spline
	else:
		dich = scipy.ones(wave.size)
	wave = scipy.arange(3400.,10400.,0.1)
	wavemodel = interpolate.splev(wave,wavecalmodel)
	goodmodel = ndimage.gaussian_filter1d(wavemodel,linewidth*scale/0.12)
	goodmodel *= dich
	goodmodel = interpolate.splrep(wave,goodmodel,s=0)

	extra = [linefile,cutoff,3000]
	extractwidth = 15

	del arc,wave,smooth

	"""
	Use the skyarcmatch routine if the 300grism is employed, otherwise
	  just match the arclines.
	"""
	if dichroic=='680' and disperser=="300/5000":
		from lris.lris_blue.skyarcmatch import arcmatch as wavematch
		extra2 = [linefile,6850,3500]
	else:
		from lris.lris_blue.arcmatch import arcmatch as wavematch
		extra2 = extra

	"""
	This could be improved by making the arrays the (pre-determined) size
	  stipulated by the red and blue cutoffs.
	"""
	print "Creating output arrays"
	outlength = int(axis2*1.6)
	out = scipy.zeros((nsci,nsize,outlength),scipy.float32)*scipy.nan
	out2 = scipy.zeros((2,csize,outlength),scipy.float32)*scipy.nan

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


	logfile = open(logfile.name,'a')
        logfile.write('Beginning Wavelength Solution and Resampling\n')
        logfile.write('--------------------------------------------\n')
        logfile.close()


	"""
        Loop through all of the slits, determining the wavelength solution and
          performing the background subtraction. It might be more robust to
          determine all wavelength solutions, then jointly determine a 'master'
          solution.... posc stores the current (starting) position of the
          coadded array, and posn stores the current position of the straight
          array. off is 0 while looping over the bottom and YMID for the top. n
	  is the number of bottom slits, so that the 1st top slit is n+1.
        """
	nbottom = len(slits['bottom'])
	nslits = nbottom+len(slits['top'])
	posc = 0
	posn = 0
	count = 1
	off = 0
	n = 0
	narrow = slits['bottom']
	wide = wide_slits['bottom']

	""" Debugging feature; set to 1 to skip background subtraction """
	lris.lris_blue.skysub.RESAMPLE = 0
	for k in range(nslits):
		"""
		When we have finished all of the bottom slits, switch
		  parameters over to their top values.
		"""
		if k==nbottom:
			arc_ycor = get_arc(out_prefix)
			yforw = get_yforw(out_prefix)
			yback = get_yback(out_prefix)
			n = nbottom
			off = YMID
			narrow = slits['top']
			wide = wide_slits['top']
		i,j = narrow[k-n]
		a,b = wide[k-n]

		""" Debugging feature; change number to skip initial slits """
		if count<1:
			count += 1
			continue
		
		print "Working on slit %d (%d to %d)" % (count,i+off,j+off)
		logfile = open(logfile.name,'a')
		logfile.write("Working on slit %d (%d to %d)\n" % (count,i+off,j+off))
		logfile.close()
		sky2x,sky2y,ccd2wave = wavematch(a,scidata[:,a+off:b+off],arc_ycor[i:j],yforw[i:j],widemodel,finemodel,goodmodel,linemodel,scale,mswave,extra,logfile)
		logfile = open(logfile.name,'a')
		logfile.write("\n")
		logfile.close()
		strt,bgsub,varimg = doskysub(i,j-i,outlength,scidata[:,a+off:b+off],yback[a:b],sky2x,sky2y,ccd2wave,scale,mswave,center,extra2)

		""" Store the resampled 2d spectra """
		h = strt.shape[1]
		if cache:
			file = pyfits.open(strtfile,mode="update")
			out = file[0].data
		out[:,posn:posn+h] = strt.copy()
		if cache:
			file.close()
			del file,out
		posn += h+5

		if lris.lris_blue.skysub.RESAMPLE==1:
			count += 1
			continue

		""" Store the resampled, background subtracted 2d spectra """
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


		""" Find and extract object traces """
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



	""" Output 2d spectra"""
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
