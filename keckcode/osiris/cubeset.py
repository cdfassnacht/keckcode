"""

cubeset - Defines a CubeSet class that contains code to handle operations
          on several IFU data cubes, e.g., coaddition

"""

import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.io import fits as pf
from .oscube import OsCube


class CubeSet(list):
    """

    The Cubeset class is effectively just a list of Spec1d instances that
    includes operations that act on the whole list (e.g., coadds)

    """

    # ------------------------------------------------------------------------

    def __init__(self, inlist, informat=None, verbose=True):
        """

        There are three ways to create a Cubeset instance

        Option 1: List of filenames
        Option 2: List of already existing OsCubeinstances
        Option 3: A filename, where the file contains within it a list
         of input files

        """

        """
        First check to see if inlist is a list, of either filenames (strings)
         or OsCube instances (options 1 and 2).
        If it's not a list, check to see if it is a single string, in which
         case it may be a file that contains a list of filenames (option 3)
        """

        if isinstance(inlist, list):
            """ Option 1 or 2: list of files or of OsCube instances """
            if isinstance(inlist[0], str):
                for f in inlist:
                    try:
                        cube = OsCube(f)
                    except IOError:
                        print('')
                        print('Could not open requested file: %s' % f)
                        return
                    self.append(cube)

            elif isinstance(inlist[0], OsCube):
                for cube in inlist:
                    self.append(cube)

            else:
                raise TypeError('input list must be a list of filenames or'
                                ' OsCube instances')

        elif isinstance(inlist, str):
            """ Option 3: a file containing a list of filenames """
            try:
                if informat is not None:
                    intab = ascii.read(inlist, guess=False, format=informat)
                else:
                    intab = ascii.read(inlist)
            except IOError:
                print('')
                print('Could not open requested file: %s' % inlist)
                print('')
                return

            for f in intab.columns[0]:
                try:
                    cube = OsCube(f)
                except IOError:
                    print('')
                    print('Could not open requested file: %s' % f)
                    return
                self.append(cube)

    # ------------------------------------------------------------------------

    def clean(self, maskfile='default', maskdir='../Clean',
              outdir='../Clean', debug=False, **kwargs):
        """

        Runs the clean algorithm for each cube in the CubeSet

        """

        """ Set up the mask file """
        if maskfile == 'default':
            s0 = self[0]
            mfile = 'mask_%s_%s.fits' % (s0.filt, s0.lenslet)
            maskfile = os.path.join(maskdir, mfile)
            if debug:
                print(maskfile)

        """
        Read the mask file into the first cube, and then copy it to the
        subsequent cubes
        """
        self[0].read_maskfile(maskfile)
        for i in range(1, len(self)):
            self[i].mask = self[0].mask.copy()

        """
        Loop through the files, cleaning each one and saving to the
        designated directory
        """
        print('')
        for f in self:
            basename = os.path.basename(f.infile)
            outfile = os.path.join(outdir, basename)
            varfile = outfile.replace('.fits', '.varspec')
            if debug:
                print('Input file:   %s' % f.infile)
                print('Output file:  %s' % outfile)
                print('Varspec file: %s' % varfile)
            f.make_varspec(outfile=varfile)
            f.clean(**kwargs)
            f.save_drp('clean', outfile)

    # ------------------------------------------------------------------------

    def coadd(self, configfile, outfile, testslice='default', testonly=False,
              slroot='tmp', wmin=None, wmax=None, verbose=True, **kwargs):
        """

        Coadd the cubes through calls to swarp
        (NOT YET IMPLEMENTED)

        """

        """ Set the actual range of wavelength slices to use """
        if wmin is None:
            wmin = 0
        if wmax is None:
            wmax = self[0].wsize

        """
        Make a test file, with the swarped coadd of a single slice.
        This is used to set the size of the output array
        """
        if testslice == 'default':
            testslice = int((wmin + wmax) / 2.)

        for i, c in enumerate(self):
            c.set_imslice(testslice, display=False)
            outname = '%s_%02d.fits' % (slroot, i)
            c['slice'].writeto(outname, overwrite=True)
        os.system('swarp %s*fits -c swarp_coseye.config' % slroot)

        """
        Get the relevent information out of the test file, and then delete
        the temporary files.
        The reason for using the WCS call is to get the WCS information in 
         the header of the swarped file into a standard format
        """
        tmphdu = pf.open('coadd.fits')
        hdr2d = (wcs.WCS(tmphdu[0].header)).to_header()
        dim2d = tmphdu[0].data.shape
        if testonly:
            return
        else:
            os.system('rm %s*fits' % slroot)
            os.system('rm coadd*fits')

        """ Split all the cubes into their slices """
        if(verbose):
            print('Splitting the cubes into slices')
        for i, c in enumerate(self):
            outroot = '%s_%i' % (slroot, i)
            c.slice_cube(outroot=outroot)
