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

    def __init__(self, inlist, informat=None, indir=None, maskfile=None,
                 verbose=True):
        """

        There are three ways to create a Cubeset instance

        Option 1: List of filenames
        Option 2: List of already existing OsCubeinstances
        Option 3: A filename, where the file contains within it a list
         of input files
         - Within this option is a specialized format for the coadd, which
           can be chosen by setting informat='coadd'
           For this option, the file must have the following columns:
             File G CRVAL1 CRVAL2 CRPIX1 CRPIX2
           where the column names are fairly self-explanatory except for G,
           which is either 0 or 1 and indicates whether the file is
           good (G=1) or bad (G=0)

        """

        """ Set default values """
        self.info = None

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
                if informat == 'coadd' or informat is None:
                    intab = ascii.read(inlist)
                else:
                    intab = ascii.read(inlist, guess=False, format=informat)
            except IOError:
                print('')
                print('Could not open requested file: %s' % inlist)
                print('')
                return

            """
            Read in the files.
            If informat is 'coadd', then use the G column to only select
             the good files
            """
            if informat == 'coadd':
                goodtab = intab[intab['G'] == 1]
                flist = goodtab['File']
                self.info = goodtab
            else:
                flist = intab.columns[0]
                self.info = intab
            for f in flist:
                if indir is 'fromfile':
                    osirisdir = os.getenv('osiris')
                    obsdir = (f.split('_')[0]).replace('s', '20')
                    infile = os.path.join(osirisdir, obsdir, 'Clean', f)
                elif indir is not None:
                    infile = os.path.join(indir, f)
                else:
                    infile = f
                try:
                    cube = OsCube(infile)
                except IOError:
                    print('')
                    print('Could not open requested file: %s' % infile)
                    return
                self.append(cube)

    # ------------------------------------------------------------------------

    def clean(self, maskfile='default', maskdir='../Clean',
              outdir='../Clean', debug=False, **kwargs):
        """

        Runs the clean algorithm for each cube in the CubeSet

        """

        """
        Read the mask file into the first cube, and then copy it to the
        subsequent cubes
        """
        self[0].read_maskfile(maskfile, maskdir)
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

    def coadd(self, configfile, outfile, wlim=None, testslice='default',
              testonly=False, slroot='slice', verbose=True, **kwargs):
        """

        Coadd the cubes through calls to swarp
        (NOT YET IMPLEMENTED)

        """

        """
        Get the range of wavelength slices to extract from the cube.
        The default is to use the full wavelength range.
        """
        if wlim is not None:
            wmin = wlim[0]
            wmax = wlim[1]
        else:
            wmin = 0
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
            outwht = outname.replace('.fits', '_wht.fits')
            hdr = c['slice'].header
            if 'ELAPTIME' in hdr.keys():
                hdr['exptime'] = hdr['elaptime']
            elif 'ITIME' in hdr.keys():
                hdr['exptime'] = hdr['itime'] / 1000.
            c['slice'].writeto(outname, overwrite=True)
            pf.PrimaryHDU(c.mask.astype(int),
                          hdr).writeto(outwht, overwrite=True)
        os.system('swarp %s*fits -c swarp_coseye.config' % slroot)

        """ If the testonly mode has been requested, quit here """
        if testonly:
            return

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
            c.slice_cube(wlim=wlim, outroot=outroot)
