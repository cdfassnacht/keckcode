from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# gaussian function to fit to Ha emission line
def gauss(x, x0, a, b, c, d):
    return a * np.exp(-((x-x0)**2)/(2*(b**2))) + c * x + d

def kin_extract(CUBE,
                sky_ra_min, sky_ra_max, sky_dec_min, sky_dec_max,
                signal_ra_min, signal_ra_max, signal_dec_min, signal_dec_max,
                guess, bounds):
    
    wave = []
    err = []

    # read in data cube
    cube = fits.open(CUBE)

    # define parameters from header
    # number of pixels across
    n1 = cube[0].header['naxis1']
    # starting/lowest wavelength observed
    wl = cube[0].header['crval1']
    # change in wavelength per pixel
    dw = cube[0].header['cdelt1']
    
    # access and get information about primaryHDU -- cube1.shape() = (slice #, y, x)
    cube1 = cube[0].data

    # define intensity spectrum over a region with signal
    val = cube1[signal_ra_min:signal_ra_max, signal_dec_min:signal_dec_max , :]
    data = val.sum(axis=0).sum(axis=0)

    # define clean sky region for error/noise estimates
    sky = cube1[sky_ra_min:sky_ra_max, sky_dec_min:sky_dec_max, :]

    for i in range(n1):
        # create list of all measured wavelengths
        val = dw * i + wl
        wave.append(val)

        # estimate error from a clean sky region
        error = np.std(sky[:,:,i])
        err.append(error)


    popt, pcov = curve_fit(gauss, wave, data, p0=guess, sigma=err, bounds=(0, bounds))

    fit = []
    for i in range(0,len(wave)):
        fit.append(gauss(wave[i],*popt))

    # plot data imported from fits file
    plt.plot(wave, data, color='gray', label = 'Raw data', alpha=0.5)

    # plot gaussian fit
    plt.plot(wave, fit, color='k', label = 'Gaussian fit')


    # add labels and title to plot
    plt.xlabel('$\lambda$ [nm]')
    plt.ylabel('flux')
    plt.axis([popt[0] - 15, popt[0] + 15, 0.90 * min(data), 1.1 * max(data)])
    plt.legend()
    plt.show()