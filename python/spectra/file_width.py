import scipy,pyfits,special_functions
from spectra import spectools

def measure(file,num=25):
	f = pyfits.open(file)

	data = f[0].data.astype(scipy.float64)/1000.
	size = data.size
	wave = spectools.wavelength(file)
	sorted = scipy.sort(data)
	zero = sorted[size/20:size/10].mean()
	sigma = sorted[size/20:size/10].std()

	thresh = zero+100*sigma

	count = 0
	search = data.copy()

	vals = scipy.zeros(num)
	place = scipy.zeros(num)
	while count<num:
		max = search.max()
		if max<thresh:
			break

		pos = search.argmax()


		search[pos-5:pos+6] = 0.

		if pos<5 or pos+6>size:
			continue
		fitdata = scipy.zeros((11,2))
		fitdata[:,0] = wave[pos-5:pos+6]
		fitdata[:,1] = data[pos-5:pos+6]

		par = scipy.zeros(4)
		par[1] = max
		par[2] = wave[pos]
		par[3] = wave[pos]-wave[pos-1]

		fit,chi2 = special_functions.ngaussfit(fitdata,par)
		if chi2>4:
			continue
		model = special_functions.ngauss(fitdata[:,0],fit)
		pylab.plot(wave[pos-5:pos+6],model)

		feature = wave[pos]
		width = fit[3]*299800/feature

		vals[count] = width
		place[count] = feature

		count += 1
	args = place.argsort()
	vals = vals[args]
	place = place[args]

	return scipy.median(vals)

