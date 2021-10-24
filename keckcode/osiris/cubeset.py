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

    def __init__(self, inlist, informat=None, indir=None, maskfile='default',
                 maskdir=None, verbose=True):
        """

        There are three ways to create a Cubeset instance

        Option 1: List of filenames
        Option 2: List of already existing OsCubeinstances
        Option 3: A filename, where the file contains within it a list
         of input files
         - Within this option is a specialized format for the coadd, which
           can be chosen by setting informat='coadd'
           For this option, the file must have the following columns:
             File CRVAL1 CRVAL2 CRPIX1 CRPIX2
           where the column names are fairly self-explanatory

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
                    tmp = ((os.path.basename(f)).split('.')[0]).split('_')
                    frame, filt, lenslet = tmp[1:4]
                    if maskfile == 'default':
                        maskfile = 'mask_%s_%s.fits' % (filt, lenslet)
                    if indir is not None:
                        infile = os.path.join(indir, f)
                    else:
                        infile = f
                    if maskdir is not None:
                        maskpath = os.path.join(maskdir, maskfile)
                    else:
                        maskpath = maskfile
                    try:
                        cube = OsCube(infile, maskfile=maskpath, verbose=False)
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
                goodmask = (intab['CRPIX1'] > 0.) & (intab['CRPIX2'] > 0.)
                goodtab = intab[goodmask]
                flist = goodtab['File']
                self.info = goodtab
            else:
                flist = intab.columns[0]
                self.info = intab
            for f, info in zip(flist, goodtab):
                tmp = ((os.path.basename(f)).split('.')[0]).split('_')
                obsdir = tmp[0].replace('s', '20')
                frame, filt, lenslet = tmp[1:4]
                if indir is 'fromfile':
                    osirisdir = os.getenv('osiris')
                    infile = os.path.join(osirisdir, obsdir, 'Clean', f)
                elif indir is not None:
                    infile = os.path.join(indir, f)
                else:
                    infile = f
                try:
                    cube = OsCube(infile, maskfile=maskfile,verbose=False)
                    crpix = [info['CRPIX1'], info['CRPIX2']]
                    crval = [info['CRVAL1'], info['CRVAL2']]
                    cube.update_wcs_from_2d(crpix, crval)
                except IOError:
                    print('')
                    print('Could not open requested file: %s' % infile)
                    return
                self.append(cube)

    # ------------------------------------------------------------------------

    def clean(self, maskfile='default', maskdir='../Clean', dsdir='../DarkSub',
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
            varfile = outfile.replace('.fits', '_varspec.fits')
            if dsdir is not None:
                dsfile = os.path.join(dsdir, basename)
                varcubefile = varfile.replace('varspec', 'varcube')
            else:
                dsfile = None
            if debug:
                print('Input file:   %s' % f.infile)
                print('Output file:  %s' % outfile)
                print('Varspec file: %s' % varfile)
                if dsfile is not None:
                    print('Darksub file: %s' % dsfile)
                    print('Varcube file: %s' % varcubefile)
            f.make_varspec(outfile=varfile, outformat='fitstab')
            f.clean(**kwargs)
            f.save_drp(outfile, 'clean')
            if dsfile is not None:
                f.make_varcube('darksub', dsfile=dsfile, outfile=varcubefile)

    # ------------------------------------------------------------------------

    def coadd(self, lensroot, configfile, centpos=None, whttype='mask',
              varsuff='varcube', wlim=None,
              testslice='default', swarp='swarp', testonly=False,
              slroot='slice', verbose=True,
              **kwargs):
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

        """ Get variance information if given """
        if whttype == 'varcube':
            varlist = []
            for i, cube in enumerate(self):
                if cube.infile is not None:
                    inname = cube.infile.replace('.fits', '_%s.fits' % varsuff)
                    varlist.append(inname)
                else:
                    raise IOError('No filename associated with cube %d' % i)
            varcube = CubeSet(varlist)
        else:
            varcube = None

        """
        Make a test file, with the swarped coadd of a single slice.
        This is used to set the size of the output array
        """
        if testslice == 'default':
            testslice = int((wmin + wmax) / 2.)

        for i, c in enumerate(self):
            """ Get the data in the test slice """
            c.set_imslice(testslice, display=False)
            outname = '%s_%03d_%02d.fits' % (slroot, testslice, i)

            """ Get information from the fits header """
            hdr = c['slice'].header
            if 'ELAPTIME' in hdr.keys():
                hdr['exptime'] = hdr['elaptime']
            elif 'ITIME' in hdr.keys():
                hdr['exptime'] = hdr['itime'] / 1000.
            else:
                hdr['exptime'] = 1.

            """ Save the science data """
            c['slice'].writeto(outname)

            """ Set the input weight type for swarp """
            mask = c.mask.astype(int)
            if whttype == 'mask':
                wsuff = '_wht.fits'
                wstr = 'MAP_WEIGHT'
                whtdat = mask
            elif whttype == 'varspec':
                wsuff = '_wht.fits'
                wstr = 'MAP_WEIGHT'
                if c.varspec is None:
                    c.make_varspec()
                whtdat = mask * c.varspec[testslice]
            elif whttype == 'varcube':
                wsuff = '_rms.fits'
                wstr = 'MAP_RMS'
                varcube[i].set_imslice(testslice, display=False)
                whtdat = varcube[i]['slice']
            else:
                raise ValueError('whttype must be one of "varcube", "varspec"'
                                 ' or "mask"')

            """ Define and save the weight image"""
            outwht = outname.replace('.fits', wsuff)
            pf.PrimaryHDU(whtdat, hdr).writeto(outwht, overwrite=True)

            """ Save the mask and exposure time images """
            outmask = outname.replace('.fits', '_mask.fits')
            texp = mask * hdr['exptime']
            outtexp = outmask.replace('mask', 'texp')
            pf.PrimaryHDU(mask, hdr).writeto(outmask, overwrite=True)
            pf.PrimaryHDU(texp, hdr).writeto(outtexp, overwrite=True)

        """ Run swarp on the science and ancillary files """
        if whttype == 'mask':
            keyvals = '-WEIGHT_TYPE %s  -WEIGHT_SUFFIX %s' % (wstr, wsuff)
        os.system('%s %s*fits -c %s %s' % (swarp, slroot, configfile, keyvals))

        addkeys = '-WEIGHT_TYPE NONE -COMBINE_TYPE SUM'
        addkeys = '%s -RESAMPLING_TYPE BILINEAR' % addkeys
        mkeys = '-IMAGEOUT_NAME %s_mask.fits' % lensroot
        os.system('%s %s*mask.fits -c %s %s %s' %
                  (swarp, slroot, configfile, addkeys, mkeys))
        # os.system('%s %s*fits -c %s %s' % (swarp, slroot, configfile, keyvals))

        """ If the testonly mode has been requested, quit here """
        if testonly:
            # os.system('rm %s*fits' % slroot)
            return

        """
        Get the relevant information out of the test file, and then delete
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
