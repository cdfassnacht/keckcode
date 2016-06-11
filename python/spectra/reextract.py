import extract,pyfits,sys,spectools,scipy

def write_slits(spectra,crpix,crval,scale,out_prefix,slitnum):
	num = 1
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
				thdu.header.update('CRPIX1',crpix)
				thdu.header.update('CRVAL2',1)
				thdu.header.update('CD2_2',1)
				thdu.header.update('CRPIX2',1)
				thdu.header.update('CTYPE1','LINEAR')
				hdulist.append(thdu)
		outname = out_prefix+"_spec_%02d_%02d.fits" % (slitnum,num)
		hdulist.writeto(outname)
		num += 1


scifile = sys.argv[1]
varfile = sys.argv[2]
prefix = sys.argv[3]

data = pyfits.open(scifile)[0].data.copy()
hdr = pyfits.open(scifile)[0].header.copy()
vardata = pyfits.open(varfile)[0].data.copy()

input = {}
input['slits'] = None
input['width'] = 10
input['noise'] = 5.

for i in range(4,len(sys.argv)):
	key,val = sys.argv[i].lower().split('=')
	input[key] = val

width = float(input['width'])
noise = float(input['noise'])

crpix = hdr['CRPIX1']
crval = hdr['CRVAL1']
disp = hdr['CD1_1']

if input['slits'] is not None:
	slits = input['slits'].split(',')
	for i in slits:
		slitnum = int(i)
		slit,start,bottom = spectools.cut_slit(data,slitnum)
		if slit.size==1:
			continue
		var = spectools.get_slit(vardata,slitnum)
		spectra = extract.extract(slit,var,width,noise=noise)

		pix = crpix - start

		write_slits(spectra,pix,crval,disp,prefix,slitnum)
else:
	slitnum = 1
	while 1:
		slit,start,bottom = spectools.cut_slit(data,slitnum)
		if slit.size==1:
			break
		var = spectools.get_slit(vardata,slitnum)
		spectra = extract.extract(slit,var,width,noise=noise)

		pix = crpix - start

		write_slits(spectra,pix,crval,disp,prefix,slitnum)
		slitnum += 1
