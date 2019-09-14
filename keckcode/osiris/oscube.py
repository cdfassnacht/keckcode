"""

oscube.py

"""

import numpy as np
from scipy.ndimage import filters
from astropy.io import fits as pf
from specim import imfuncs as imf
from specim import specfuncs as ss
from specim.imfuncs.wcshdu import WcsHDU

# ===========================================================================


class OsCube(imf.Image):
    """
    A class used to visualize and analyze OSIRIS data
    """

    def __init__(self, indat, verbose=True):
        """
        Loads in an OSIRIS data cube that was processed through the
        standard data reduction pipeline and, possibly, has had some
        additional processing done.

        Inputs:
          indat - either the name of an input fits file or a HDU
        """

        """ Read the data into an Image structure """
        super(OsCube, self).__init__(indat, verbose=verbose)  # Python 2.7
        # super().__init__(self, indat, verbose=verbose) # Python 3 syntax
        print('Number of wavelength slices: %d' % self.header['naxis1'])

        """ Set up the wavelength vector based on the header information """
        hdr = self.header
        nwav = hdr['naxis1']
        self.wav = np.arange(nwav)
        self.wav = hdr['crval1'] + self.wav * hdr['cdelt1']

        """ Convert from nm to angstroms, including in the header """
        self.wav *= 10.
        hdr['crval1'] *= 10.
        hdr['cdelt1'] *= 10.

        """ Save the sizes of the array in easily accessible form """
        self.xsize = self.data.shape[0]
        self.ysize = self.data.shape[1]
        self.wsize = self.data.shape[2]

        """ Set default values """
        self.cube = None
        self.moment0 = None
        self.meanspec = None
        self.varspec = None

    # -----------------------------------------------------------------------

    # Actually, don't use this syntax below unless required.  See the
    #  @property calls in dispparam.py
    #
    # imslice = property(fget=get_imslice, fset=set_imslice)

    # -----------------------------------------------------------------------

    def make_wcs2dhdr(self):
        """

        Takes the input WCS information for the cube (3 dimesions, including
        wavelength) and converts it into 2-dimensional spatial-only WCS
        information for an image slice or the cube compressed along the
        spatial direction.
        """

        """ Get the WCS information """
        outwcs = self.wcsinfo.celestial.swapaxes(0, 1)
        outhdr = outwcs.to_header()

        """ Add other important info """
        hdr = self.header
        kwlist = ['bunit', 'bscale', 'bzero', 'itime', 'coadds', 'object',
                  'sfilter', 'sscale', 'telescop', 'instrume', 'elaptime',
                  'filter']
        for k in kwlist:
            if k.upper() in hdr.keys():
                outhdr[k] = hdr[k]

        return outhdr

    # -----------------------------------------------------------------------

    def set_imslice(self, imslice, hext=0, display=True, mode='xy',
                    **kwargs):
        """
        Sets the 2-dimension slice to use for the display functions.

        Inputs:
          imslice - which image slice to use.
        """

        """ Get the number of dimensions in the input image """
        hdr = self.header
        if 'NAXIS' in hdr.keys():
            ndim = hdr['naxis']
        else:
            raise KeyError('No NAXIS keyword in fits header')

        """
        The OSIRIS data-reduction pipeline produces a 3-dimensional data cube.
        Check this
        """
        if ndim != 3:
            print('')
            print('ERROR: Expected a 3-dimensional data cube but found '
                  '%d dimensions' % ndim)
            print('')
            raise ValueError

        """
        Select the image slice to use.  Right now this assumes the standard
        axis order produced by the pipeline, namely that CTYPE1 is WAVE,
        CTYPE2 is RA, and CTYPE3 is DEC.  This means that the data array
        structure is hard-wired to this assumption
        """

        if imslice >= self.wsize:
            print('')
            print('ERROR: Requested an image slice outside the available '
                  'range')
            print('Maximum available slice value is %d' % (self.wsize-1))
            print('')
            raise IndexError

        """
        Make a wcs header for the slice, i.e., a set of 2d wcs information
        without the spectral information.
        """
        w2dhdr = self.make_wcs2dhdr()

        """
        Actually make the slice.
        The transpose is to get RA and Dec into the order that the WcsHDU
         container expects them to be
        The header information.
        """
        # self['slice'] = WcsHDU(np.transpose(self.data[:, :, imslice]))
        self['slice'] = WcsHDU(np.transpose(self.data[:, :, imslice]), w2dhdr)
        self.prevdext = hext

        """ Display the image slice if requested """
        self.found_rms = False
        if display:
            self.display(dmode='slice', mode=mode,
                         title='Image Slice %d (zero-indexed)' % imslice,
                         **kwargs)

    # -----------------------------------------------------------------------

    def select_cube(self, wlim=None, xlim=None, ylim=None, wmode='slice',
                    hext=0, verbose=False):
        """
        Creates a cube to analyze, defined by ranges in x, y, and wavelength.
        """

        """ Use default values if none are requested """
        if xlim is not None:
            xmin = xlim[0]
            xmax = xlim[1]
        else:
            xmin = 0
            xmax = self.xsize
        if ylim is not None:
            ymin = ylim[0]
            ymax = ylim[1]
        else:
            ymin = 0
            ymax = self.ysize
        if wlim is not None:
            wmin = wlim[0]
            wmax = wlim[1]
        else:
            wmin = 0
            wmax = self.wsize

        if verbose:
            print('')
            print('Creating a cube from original data with ranges:')
            print('  x:      %d - %d' % (xmin, xmax))
            print('  y:      %d - %d' % (ymin, ymax))
            print('  lambda: %d - %d (slice number)' % (wmin, wmax))

        """ Select the cube """
        cube = self.data[xmin:xmax, ymin:ymax, wmin:wmax]
        cubehdr = self.header.copy()
        cubehdr['crpix1'] -= wmin
        cubehdr['crpix2'] -= ymin
        cubehdr['crpix3'] -= xmin

        """ Return the results """
        return cube, cubehdr

    # -----------------------------------------------------------------------

    def compress_spec(self, wlim=None, xlim=None, ylim=None, wmode='slice',
                      combmode='sum', hext=0, display=True, verbose=True,
                      **kwargs):
        """

        Compresses the data cube along the spectral dimension, but only
        for image slices between some minimum and maximum wavelengths.
        These wavelength limits (wlim) can be set either by the slice
        number or the actual wavelength in Angstrom [wavelength mode is NOT
        yet implemented].
        Setting wlim=None (the default) will use the full wavelength range

        The compression can be done either as a sum or as a median.

        The result is a 2-dimensional spatial image, which is stored in the
        data container and, thus, can be easily displayed.
        """

        """ First select the image slices to use """
        if wmode == 'wavelength':
            print('NOTE: wavelength mode has not yet been implemented')
            return
        else:
            if wlim is None:
                wlim = (0, self.wsize-1)
            minslice = wlim[0]
            maxslice = wlim[1]
            wavmin = self.wav[minslice]
            wavmax = self.wav[maxslice]

        if verbose:
            print('')
            print('Data cube will be compressed along the spectral direction')
            print('  Image slice range: %d - %d' % (minslice, maxslice))
            print('  Corresponding wavelength range: %8.2f - %8.2f'
                  % (wavmin, wavmax))

        """ Create a temporary cube container """
        cube, cubehdr = self.select_cube(wlim, xlim, ylim, hext=hext,
                                         verbose=verbose)

        """ Compress the temporary cube along the spectral axis """
        if combmode == 'median':
            self['slice'] = WcsHDU(np.transpose(np.median(cube, axis=2)))
        else:
            self['slice'] = WcsHDU(np.transpose(cube.sum(axis=2)))
        self.prevdext = hext

        """ Display the result if requested """
        if display:
            self.display(dmode='slice', **kwargs)

        """ Clean up """
        del(cube)

    # -----------------------------------------------------------------------

    def whitelight(self, combmode='sum', hext=0, display=True, verbose=True,
                   **kwargs):
        """
        Creates, and displays if requested, the "white light" image.
        This is a specialized version of the compress_spec method, which
        compresses the cube along* the spectral direction, but for the full
        data cube.
        """

        if verbose:
            print('Creating white light image')

        """ Set up for running compress_spec on the full data cube """
        self.compress_spec(wmode='slice', combmode=combmode,
                           hext=hext, display=display, **kwargs)

    # -----------------------------------------------------------------------

    def make_varspec(self, maskfile, verbose=False):
        """
        Steps through the spectral slices and computes the statistics of
        the illuminated region and stores the clipped mean and variance
        for each slice
        """

        """ Get the mask from the mask file """
        mhdu = pf.open(maskfile)
        mask = mhdu[0].data.astype(bool)

        """ Set up the containers to store the info """
        mean = np.zeros(self.wsize)
        var = np.zeros(self.wsize)

        count = 0
        for i in range(self.wsize):
            self.set_imslice(i, display=False)
            self['slice'].sigma_clip(mask=mask, verbose=False)
            mean[i] = self['slice'].mean_clip
            r = self['slice'].rms_clip
            var[i] = r**2
            if verbose:
                print(count, mean[i], r)
            count += 1

        self.meanspec = mean
        self.varspec = var

    # -----------------------------------------------------------------------

    def make_1dspec(self, x, y, hext=0, display=True,
                    skyx=None, skyy=None, debug=False, **kwargs):
        """
        Takes a spaxel or a rectangular region, designated by the (x, y)
         coordinates and extracts the spectral information into a Spec1d
         container.
        Also plots the 1d spectrum if requested.

        Returns: the Spec1d container (instance).

        """

        """ Set the data set to use """
        cube = self.data

        if isinstance(x, int) and isinstance(y, int):
            flux = cube[x, y, :]
            npix = 1
        elif (isinstance(x, list) or isinstance(x, tuple)) and \
                (isinstance(y, list) or isinstance(y, tuple)):
            """ Need to check lengths """
            xmin = int(x[0])
            xmax = int(x[1])
            ymin = int(y[0])
            ymax = int(y[1])
            flux = cube[xmin:xmax, ymin:ymax, :].sum(axis=0)
            flux = flux.sum(axis=0)
            npix = (xmax - xmin) * (ymax - ymin)
        if debug:
            print('npix: %d' % npix)
            print(self.wav.size, flux.size)

        """
        Make the variance spectrum if it has been requested by
        setting skyx and skyy
        """
        if skyx is not None and skyy is not None:
            skycube, hdr = self.select_cube(xlim=skyx, ylim=skyy,
                                            hext=hext)
            var = npix * skycube.var(axis=(0, 1))
        else:
            var = None

        """
        Make, and display if requested, the final spectrum.
        """
        spec = ss.Spec1d(wav=self.wav, flux=flux, var=var)

        if display:
            spec.plot(**kwargs)

        """ Clean up and return """
        del flux
        if var is not None:
            del var
        return spec

    # -----------------------------------------------------------------------

    def click_1dspec(self, hext=0, xysmooth=1, **kwargs):
        """
        An interactive interface to set_1dspec.  Produces the 1D spectrum
        associated with the spaxel that is clicked.
        """

        self.start_interactive()
        self.set_1dspec(int(self.xclick), int(self.yclick), **kwargs)

    # -----------------------------------------------------------------------

    def smooth_xy(self, kwidth, hext=0, outfile=None):
        """
        Smooths the cube over the two spatial dimensions.
        The smoothing is a circular gaussian with sigma=kwidth, and where
         kwidth is given in pixels
        """

        """ Smooth the data """
        data = self.data
        cube = filters.gaussian_filter(data, sigma=[kwidth, kwidth, 0])

        """ Put the smoothed data into a new HDU """
        hdr = self.header.copy()
        hdr['history'] = 'Data have been spatially smoothed (gaussian)'
        hdr['smoothw'] = ('%5.1f' % kwidth,
                          'Gaussian smoothing kernel width')
        phdu = pf.PrimaryHDU(cube, hdr)

        """ Save the smoothed cube in an output file if desired """
        if outfile:
            print('')
            print('Wrote smoothed data cube to %s' % outfile)
            phdu.writeto(outfile, overwrite=True)
            print('')

        """ Return the smoothed cube """
        hdu = pf.HDUList()
        hdu.append(phdu)
        return hdu

    # -----------------------------------------------------------------------

    def make_moment0(self, wlim, xlim, ylim, wmode='slice', combmode='sum',
                     hext=0, display=True, verbose=True, **kwargs):
        """

        Makes the moment 0 (i.e., total flux) map of the science region
        of the data cube.  The science region is defined by a wavelength
        range (wlim), x-pixel range (xlim), and y-pixel range (ylim).

        Essentially this method is just a front-end for the (barely) more
        generic compress_spec method, with a perhaps easier to remember
        name
        """

        """ Call compress_spec to create the 2-d moment0 image """
        self.compress_spec(wlim, xlim, ylim, wmode=wmode, hext=hext,
                           combmode=combmode, display=display,
                           verbose=verbose, **kwargs)

        """ Save the data in a moment0 attribute """
        self.moment0 = self.data.copy()

    # -----------------------------------------------------------------------

    def save_xyl(self, outfits, outtext=None, **kwargs):
        """

        Saves the data cube as a fits file but with the spectral axis
         as axis 3 (in fits parlance) rather than axis 1.
         The OSIRIS data reduction pipeline produces a data cube with
         FITS data axis order of (lambda, y, x).  This method produces
         an output cube with a FITS data axis order of (x, y, lambda)
        If the optional outtext parameter is set to a filename, then an
         output text file will be created in the format expected by
         Francesca Rizzo's modeling code

        Required inputs:
          outfits - name of output fits file

        Optional inputs:
          outtext - if set to a string value, then save an output text file
                    in the format expected by Francesca Rizzo's modeling code.
                    The default value (None) means that no output file is
                    created.
        """

        """
        Create a copy of the data cube and swap the axes
        """
        cube, hdr0 = self.select_cube(verbose=True, **kwargs)
        tmp = cube.copy()
        tmp2 = np.swapaxes(tmp, 0, 2)

        """
        Create a new header and write the output file
        """
        tmphdu = pf.PrimaryHDU(tmp2)
        hdr = tmphdu.header
        """ Non-coordinate keywords """
        klist = ['bunit', 'bscale', 'bzero', 'itime', 'coadds', 'sampmode',
                 'numreads', 'saturate', 'instr', 'pscale', 'object',
                 'pa_spec', 'sfilter', 'telescop', 'instrume', 'targwave',
                 'airmass', 'wcsaxes']
        for k in klist:
            try:
                hdr[k] = hdr0[k]
            except KeyError:
                continue
        """
        Coordinate keywords.
        For now we're assuming that the input was (Wavelength, Dec, RA) and
         the output will be (RA, Dec, Wavelength)
        """
        crlist = ['ctype', 'cunit', 'crval', 'crpix', 'cdelt']
        for k in crlist:
            hdr['%s1' % k] = hdr0['%s3' % k]
            hdr['%s2' % k] = hdr0['%s2' % k]
            hdr['%s3' % k] = hdr0['%s1' % k]
        hdr['cunit3'] = 'Ang'
        """ PC matrix """
        hdr['pc1_1'] = hdr0['pc%d_%d' % (self.raaxis, self.raaxis)]
        hdr['pc1_2'] = hdr0['pc%d_%d' % (self.raaxis, self.decaxis)]
        hdr['pc2_1'] = hdr0['pc%d_%d' % (self.decaxis, self.raaxis)]
        hdr['pc2_2'] = hdr0['pc%d_%d' % (self.decaxis, self.decaxis)]

        print('')
        print('Saving to output file %s' % outfits)
        tmphdu.writeto(outfits, overwrite=True)

        """ Create an output text file if requested """
        if outtext is not None:
            print('Creating kinematics data file %s' % outtext)
            f = open(outtext, 'w')
            f.write('#Cube\n')
            f.write('nch %d\n' % tmp2.shape[0])
            f.write('ctype wave\n')
            f.write('cunit angstrom\n')
            f.write('cd3 %f\n' % hdr['cdelt3'])
            f.write('cv3 %f\n' % hdr['crval3'])
            f.write('crpix3 %d\n' % hdr['crpix3'])
            f.write('crota %7.2f\n' % self.impa)
            f.close()

        """ Clean up before exiting """
        del(tmp, tmphdu)
