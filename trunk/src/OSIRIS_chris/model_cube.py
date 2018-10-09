import astropy
import numpy as np
from astropy.io import fits
from astropy.modeling import functional_models
from astropy.modeling import fitting
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
import numpy.ma as ma

from emission_line_fitting import gauss

def TT_point_spread(TIP_TILT, TT_ra_min, TT_ra_max, TT_dec_min, TT_dec_max, filt):
    tip_tilt = fits.open(TIP_TILT)
    tip_tilt1 = tip_tilt[0].data
    
    tip_tilt1 = gaussian_filter(tip_tilt1, sigma = [filt, filt, 0])

    val = tip_tilt1[TT_ra_min:TT_ra_max, TT_dec_min:TT_dec_max, :]

    spatial = val.sum(axis=2)

    
    # calculate xgrid + ygrid
    x, y = np.mgrid[:len(spatial), :len(spatial[0])]

    test = astropy.modeling.functional_models.Gaussian2D(amplitude = 1000, x_mean = 0.5 * len(spatial),
                                                         y_mean = 0.5 * len(spatial[0]))
    f1 = fitting.LevMarLSQFitter()
    t = f1(test, x, y, spatial)
    
    avg_fwhm = 0.5 * (t.x_fwhm+t.y_fwhm)

    return (t.x_fwhm, t.y_fwhm, avg_fwhm)

def skyline_spec(INST):

    CUBE = fits.open(INST)

    # define parameters from header

    # number of pixels across
    n1 = CUBE[0].header['naxis1']
    # starting/lowest wavelength observed
    wl = CUBE[0].header['crval1']
    # change in wavelength per pixel
    dw = CUBE[0].header['cdelt1']

    # create list of all measured wavelengths
    wave = []
    for i in range(0, n1):
        val = dw * i + wl
        wave.append(val)
    
    # access and get information about primaryHDU -- (slice #, y, x)
    cube1 = CUBE[0].data

    val = cube1[:, :, :]
    #data = val
    data = val.sum(axis=0).sum(axis=0)

    return wave, data

def inst_disp(INST, wave, data, p0_w, bound_w, p0_a, bound_a):
    CUBE = fits.open(INST)
    cube1 = CUBE[0].data
    n1 = CUBE[0].header['naxis1']
    
    # estimate error from a clean sky region
    sky = cube1[28:36,45:60,:]
    err = []

    for i in range(0,n1):
        error = np.std(sky[:,:,i])
        err.append(error)
        
    # assign and bounds for parameters
    popt, pcov = curve_fit(gauss, wave, data, p0=(p0_w, p0_a, 0.5, 0., 0.006), sigma=err,  
                           bounds=(0,[bound_w, bound_a, 5, 0.5 , .01]))

    # calculate intensity of fit
    fit = []
    for i in range(0,len(wave)):
        fit.append(gauss(wave[i],*popt))

    xmin = min(wave)
    xmax = max(wave)
    ymin = 1.05 * min(data)
    ymax = 1.05 * max(data)

    # plot data imported from fits file
    plt.plot(wave, data, label = 'Loaded from file')

    # plot gaussian fit
    plt.plot(wave, fit, label = 'gaussian fit')


    # add labels and title to plot
    plt.xlabel('wavelength')
    plt.ylabel('intensity')
    plt.title('Raw Data')
    plt.axis([xmin, xmax, ymin, ymax])
    plt.legend()
    plt.show()

    # calculate velocity dispersion
    c = 2.998 * 10**5
    sigma = popt[2]
    wavelength = popt[0]
    vdisp = c * (sigma/wavelength)
    
    return vdisp

def modeled_gauss(x, x0, A, b):
    return A * np.exp(-((x-x0)**2)/(2*(b**2)))

def model_data_cube(params, combined_mask, TT_psf, inst, wave_list, filt, deltara, deltadec):

    comb_mask = []
    for i in combined_mask:
        val = [i]*5
        comb_mask.append(val)
    c = 2.998 * 10**5

    y = ma.array(params, mask = comb_mask)
    modeled_cube = []

    for i in y:
        wavelength = i[0]
        amp = i[1]
        spread = (inst * wavelength)/c
        val = modeled_gauss(wave_list, wavelength, amp, spread)
        modeled_cube.append(val)
        
    modeled_cube = [modeled_cube[deltadec*k:deltadec*(k+1)] for k in range(0, deltara)]


    s = 2
    w = TT_psf
    t = (((w - 1)/2)-0.5)/s
    filtereddata_psf = gaussian_filter(modeled_cube, sigma = s, truncate = t)
    filtereddata = gaussian_filter(filtereddata_psf, sigma = [filt, filt, 0])
    
    return filtereddata

def correct_dispersion(dispersion, deltara, deltadec, wave_list, model, guess, bounds, combined_mask, scale, plot=False):

    vdisp_vals = []
    c = 2.998 * 10**5

    for j in range(0, deltara):
        for k in range(deltadec):
            popt, pcov = curve_fit(gauss, wave_list, model[j][k], p0 = guess, bounds = (0, bounds))
        
            wavelength = popt[0]
            amp = popt[1]
            sig = popt[2]

       
            vdisp = c * (sig/wavelength)
            vdisp_vals.append(vdisp)
        
    corrections = ma.array(vdisp_vals, mask = combined_mask)
    corrections = [corrections[deltadec*k:deltadec*(k+1)] for k in range(0, deltara)]
    corrections = ma.masked_array(corrections, mask=combined_mask)
    
    disp2 = dispersion**2 - corrections**2
    disp2[disp2 < 0] = 0.01
    disp = np.sqrt(disp2)

    if plot == True:
        plt.imshow(corrections, origin='lower', cmap='jet', extent = (0, deltadec * scale, 0, deltara * scale))
        plt.colorbar()
        plt.show()
    
    return disp