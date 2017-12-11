"""

osiris_funcs.py

"""

import numpy as np
import imfuncs as imf
import spec_simple as ss

# ===========================================================================

class osCube(imf.Image):
    """
    A class used to visualize and analyze OSIRIS data
    """

    def __init__(self, infile, verbose=True):
        """
        Loads in an OSIRIS data cube that was processed through the
        standard data reduction pipeline and, possibly, has had some
        additional processing done.
        """

        """ Read the data into an Image structure """
        imf.Image.__init__(self, infile, verbose=verbose)

        """ Set a default image plane to plot """
        self.set_imslice(0, display=False)

        """ Set up the wavelength vector based on the header information """
        hdr = self.hdu[0].header
        nwav = hdr['naxis1']
        self.wav = np.arange(nwav)
        self.wav = hdr['crval1'] + self.wav * hdr['cdelt1']
        """ Convert to angstroms from nm """
        self.wav *= 10.

        """ Save the sizes of the array in easily accessible form """
        self.xsize = self.hdu[0].data.shape[0]
        self.ysize = self.hdu[0].data.shape[1]
        self.wsize = self.hdu[0].data.shape[2]

        """ Set default values """
        self.cube = None

    # -----------------------------------------------------------------------

    def set_imslice(self, imslice=0, hext=0, display=True, mode='xy',
                    **kwargs):
        """
        Sets the 2-dimension slice to use for the display functions.

        Inputs:
          imslice - which image slice to use.
        """

        """ Get the number of dimensions in the input image """
        hdr = self.hdu[hext].header
        if 'NAXIS' in hdr.keys():
            ndim = hdr['naxis']
        else:
            raise KeyError

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

        nwave = hdr['naxis1']
        if imslice >= nwave:
            print('')
            print('ERROR: Requested an image slice outside the available '
                  'range')
            print('Maximum available slice value is %d' % (nwave-1))
            print('')
            raise IndexError

        self.data = np.transpose(self.hdu[hext].data[:, :, imslice].copy())

        """ Display the image slice if requested """
        self.found_rms = False
        if display:
            self.display(title='Image Plane %d (zero-indexed)' % imslice,
                         mode=mode, **kwargs)

    # -----------------------------------------------------------------------

    def compress_spec(self, wmin, wmax, wmode='slice', combmode='sum',
                      hext=0, display=True, **kwargs):
        """
        Compresses the data cube along the spectral dimension, but only
        for image slices between some minimum and maximum wavelengths.
        These wavelength limits (wmin, wmax) can be set either by the slice
        number or the actual wavelength in Angstrom [wavelength mode is NOT
        yet implemented].

        The compression can be done either as a sum or as a median.

        The result is a 2-dimensional spatial image, which is stored in the
        data container and, thus, can be easily displayed.
        """

        """ First select the image slices to use """
        if wmode == 'wavelength':
            print('NOTE: wavelength mode has not yet been implemented')
            return
        else:
            minslice = wmin
            maxslice = wmax
            wavmin = self.wav[minslice]
            wavmax = self.wav[maxslice]
        print('')
        print('Data cube will be compressed along the spectral direction')
        print('  Image slice range: %d - %d' % (minslice, maxslice))
        print('  Corresponding wavelength range: %8.2f - %8.2f'
              % (wavmin, wavmax))

        """ Create a temporary cube container """
        cube = self.hdu[hext].data[:, :, minslice:(maxslice+1)].copy()

        """ Compress the temporary cube along the spectral axis """
        if combmode == 'median':
            self.data = np.transpose(np.median(cube, axis=2))
        else:
            self.data = np.transpose(cube.sum(axis=2))

        """ Display the result if requested """
        if display:
            self.display(**kwargs)

        """ Clean up """
        del(cube)

    # -----------------------------------------------------------------------

    def whitelight(self, combmode='sum', hext=0, display=True, verbose=True,
                   **kwargs):
        """
        Creates, and displays if requested, the "white light" image.
        This is a specialized version of the compress_spec method, which
        compresses the cube along the spectral direction, but for the full
        data cube.
        """

        if verbose:
            print('Creating white light image')

        """ Set up for running compress_spec on the full data cube """
        wmax = self.hdu[hext].data.shape[2] - 1
        self.compress_spec(wmin=0, wmax=wmax, wmode='slice', combmode=combmode,
                           hext=hext, display=display, **kwargs)

    # -----------------------------------------------------------------------

    def set_1dspec(self, x, y, hext=0, display=True, **kwargs):
        """
        Takes a spaxel designated by its (x, y) coordinates and extracts the
        spectral information into a Spec1d container.
        Also plot the 1d spectrum if requested.
        """

        flux = self.hdu[hext].data[x, y, :]
        print(self.wav.size, flux.size)
        spec = ss.Spec1d(wav=self.wav, flux=flux)

        if display:
            spec.plot(**kwargs)

        return spec

    # -----------------------------------------------------------------------

    def select_cube(self, wmin=None, wmax=None, xmin=None, xmax=None,
                    ymin=None, ymax=None, wmode='slice', hext=0):
        """
        Creates a cube to analyze, defined by ranges in x, y, and wavelength.
        """

        """ Use default values if none are requested """
        if xmin == None:
            xmin = 0
        if xmax == None:
            xmax = self.xsize
        if ymin == None:
            ymin = 0
        if ymax == None:
            ymax = self.ysize
        if wmin == None:
            wmin = 0
        if wmax == None:
            wmax = self.wsize

        """ Select the cube """
        self.cube = self.hdu[0].data[xmin:xmax, ymin:ymax, wmin:wmax]
