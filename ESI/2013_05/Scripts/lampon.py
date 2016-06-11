import pyfits

d = pyfits.open('../raw/e130514_0060.fits')
cuar = pyfits.open('../raw/e130514_0002.fits')[0].data.copy()

d[0].data += cuar
d[0].header['LAMPCU1'] = 'on'
d.writeto('../raw/e130514_2001.fits')#,clobber=True)
