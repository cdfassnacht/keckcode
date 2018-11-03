from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter

from emission_line_fitting import gauss

def integrated_dispersion(CUBE,
                sky_ra_min, sky_ra_max, sky_dec_min, sky_dec_max,
                signal_ra_min, signal_ra_max, signal_dec_min, signal_dec_max,
                guess, bounds, inst_disp, filt, plot=True):
    
    c = 2.998 * 10**5
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
    
    filtereddata = gaussian_filter(cube1, sigma = [filt, filt, 0])

    # define intensity spectrum over a region with signal
    val = filtereddata[signal_ra_min:signal_ra_max, signal_dec_min:signal_dec_max , :]
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
    
    unc = np.sqrt(np.diag(pcov))
    
    wavelength = popt[0]
    sig = popt[2]
    vdisp = c * (sig/wavelength)
    
    corrected_vdisp = np.sqrt(vdisp**2 - inst_disp**2)
    err = unc[2] * c / wavelength
    
    if plot == True:
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
    
    return corrected_vdisp, err