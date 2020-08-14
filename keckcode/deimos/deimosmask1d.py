"""

deimosmask1d  - Code to define a DeimosMask1d class, which is used to
                visualize and process the spec1d*fits files that come
                out of the PypeIt data reduction pipeline for DEIMOS data

"""

import numpy as np
from collections import OrderedDict

from astropy.io import fits as pf
from astropy.table import Table, join

from specim.specfuncs import spec1d, specset1d

import sys
pyversion = sys.version_info.major

# ---------------------------------------------------------------------------


class DeimosMask1d(OrderedDict):
    """

    Class to visualize and perhaps process the 1d spectra that are contained
    in the spec1d*fits that come out of the PypeIt data reduction pipeline
    for DEIMOS data

    """

    def __init__(self, indat, specdict=None, verbose=True):
        """

        Instantiate a DeimosMask1d object

        Inputs:
         indat    - One of the following:
                    1. a spec1d*fits file that has been produced by the PypeIt
                       data reduction pipeline
                    2. an astropy Table containing, at a minimum, the following
                       columns: 'slitid', 'objid', 'spatloc', 'fwhm'
                       In this case, the optional specdict parameter must be
                        provided
         specdict - A dictionary containing the spectra
        """

        """ Set up the empty container by calling the superclass """
        if pyversion == 2:
            super(DeimosMask1d, self).__init__()
        else:
            super().__init__()

        """ Set default values """
        self.hdr0 = None
        self.tabnames = ['det', 'slitid', 'objid', 'spatloc', 'fwhm']

        """ Read in the spectra and other information """
        if isinstance(indat, str):
            self.read_fits(indat, verbose=verbose)
        elif isinstance(indat, Table):
            self.read_table(indat, specdict, verbose=verbose)
        else:
            raise TypeError('Input data must either be the name of a fits '
                            'file or an astropy Table')

        """ Sort the slitinfo table """
        self.slitinfo.sort(['det', 'slitid', 'objid', 'spatloc'])

    # -----------------------------------------------------------------------

    def read_fits(self, infile, verbose=True):
        """

        Reads the spectra and other information from an input fits file
         that has been produced by the pypeit pipeline.

        Inputs:
         infile - name of input fits file

        """

        """
        Load the data from the input file and get information from the
        primary header
        """
        hdu = pf.open(infile)
        self.hdr0 = hdu[0].header
        self.nspec = self.hdr0['nspec']
        self.slitinfo = \
            Table(np.zeros((self.nspec, len(self.tabnames))),
                  names=self.tabnames, dtype=[int, int, int, int, float])

        """ Load the spectra into a list of Spec1d objects """
        colnames = ['opt_wave', 'opt_counts', 'opt_counts_ivar',
                    'opt_counts_sky']
        if verbose:
            print('Reading %d spectra from:\n  %s' % (self.nspec, infile))
        for i in range(self.nspec):
            """ Get information about the object """
            tmphdu = hdu[i+1]
            hdr = tmphdu.header
            det = int(hdr['name'].split('-')[2][3:])
            self.slitinfo['det'][i] = det
            self.slitinfo['slitid'][i] = hdr['slitid']
            spatloc = int(hdr['name'].split('-')[0][4:])
            self.slitinfo['spatloc'][i] = spatloc
            self.slitinfo['objid'][i] = hdr['objid']
            self.slitinfo['fwhm'][i] = hdr['fwhm']

            """ Load in the spectrum """
            spec = spec1d.Spec1d(tmphdu, colnames=colnames, verbose=False)
            mask = spec['var'] > 0.
            spec['var'][mask] = 1. / spec['var'][mask]
            spec['var'][~mask] = 25. * spec['var'].max()
            specid = '%d_%d_%d_%d' % (det, hdr['slitid'], hdr['objid'], spatloc)
            self[specid] = spec

    # -----------------------------------------------------------------------

    def read_table(self, specinfo, specdict, verbose=True):
        """

        Reads the spectra and other information from an astropy Table (for
        the auxillary information) and a dictionary containing the spectra

        """

        """ Check to make sure that we have the spectra """
        if specdict is None:
            raise ValueError('You must provide a dictionary containing the '
                             'input spectra')

        """ Make sure that the specinfo table has the correct column """
        for name in self.tabnames:
            if name not in specinfo.colnames:
                raise KeyError('Missing %s column in input table' % name)
        self.slitinfo = specinfo.copy()

        """ Prepare to read in the data """
        self.nspec = len(specinfo)

        """ Import the spectra """
        if verbose:
            print('Reading %d spectra from input table and spectrum '
                  'dictionary' % self.nspec)
        for info in specinfo:
            specid = '%d_%d_%d_%d' % (info['det'], info['slitid'],
                                      info['objid'], info['spatloc'])
            if specid in specdict.keys():
                self[specid] = specdict[specid]

    # -----------------------------------------------------------------------

    def __str__(self):
        """

        Gives a cleaner output when doing a print statement for the object

        """

        return('DeimosMask1d with %d spectra' % self.nspec)

    # -----------------------------------------------------------------------

    def __repr__(self):
        """

        Gives a cleaner representation of the object

        """

        return('DeimosMask1d with %d spectra' % self.nspec)

    # -----------------------------------------------------------------------

    def plot(self, specid, title='default', **kwargs):
        """

        Plots the spectrum for the slit designated by the specid paramters

        Inputs:
           specid    - slit ID, a string containing info on the detector,
                       slit, object, and location of trace, e.g. '1_120_1_98'
        """
        if title == 'default':
            title = 'Extracted spectrum for %s' % specid
        return self[specid].plot(title=title, **kwargs)

    # -----------------------------------------------------------------------

    def smooth(self, specid, filtwidth, title='default', **kwargs):
        """

        Plots a smoothed version of the spectrum designated by the
         specid parameter

        Inputs:
         specid    - slit ID, a string containing info on the detector,
                       slit, object, and location of trace, e.g. '1_120_1_98'
         filtwidth - width of kernel to be used for smoothing

        """

        if title == 'default':
            title = 'Smoothed spectrum for %s (%d-pix kernel)' \
                % (specid, filtwidth)
        self[specid].smooth(filtwidth, title=title, **kwargs)

    # -----------------------------------------------------------------------

    def mark_lines(self, lines, z, index, **kwargs):
        """
        Draws a plot with spectral line marks on a given spectra.

        Inputs:
          index - index of spectrum to plot (between 0 and nspec-1)
          z     - redshift estimate
          lines - dictionary. Key should be the line name, value should be
                  boolean

        """
        if lines['em'] and lines['strongem']:
            lines.pop('strongem')

        for k,v in lines.items():
            if v:
                self[index].mark_lines(k, z, **kwargs)

    # -----------------------------------------------------------------------

    def save(self, outfile):
        """

        Reformat the data back into the format of the input file, and then
        save it to the output fits file.

        """

        """
        Set up the Primary HDU and put relevant information into its header
        """
        phdu = pf.PrimaryHDU()
        keys = ['ra', 'dec', 'target', 'dispname', 'decker', 'binning',
                'mjd', 'airmass', 'exptime', 'instrume', 'telescop',
                'lon-obs', 'lat-obs', 'alt-obs']
        if self.hdr0 is not None:
            for k in keys:
                if k in self.hdr0.keys():
                    phdu.header[k] = self.hdr0[k]
        phdu.header['nspec'] = self.nspec
        hdulist = pf.HDUList([phdu,])

        """ Loop through the spectra """
        for info in self.slitinfo:
            name = '%d_%d_%d_%d' % (info['det'], info['slitid'], info['objid'],
                                    info['spatloc'])
            extname = 'SPAT%04d-SLIT%04d-DET%02d' % \
                (info['spatloc'], info['slitid'], info['det'])
            spec = self[name]

            """
            Save the variance information, since the output wants inverse
            variance
            """
            mask = spec['var'] > 0
            tmpvar = spec['var'].copy()

            """ Temporarily rename columns """
            incols = ['wav', 'flux', 'var', 'sky']
            outcols = ['opt_wave', 'opt_counts', 'opt_counts_ivar',
                       'opt_counts_sky']
            for i, o in zip(incols, outcols):
                spec.rename_column(i, o)
            spec['opt_counts_ivar'][mask] = 1./ tmpvar[mask]

            """ Put the spectrum into a binary table extension """
            exthdu = pf.BinTableHDU(spec.as_array(), name=extname)

            """ Put relevant information into the table header """
            for name in self.slitinfo.colnames:
                exthdu.header[name] = info[name]
            exthdu.header['name'] = extname
            hdulist.append(exthdu)

            """ Put the column names back to their original values """
            for i, o in zip(incols, outcols):
                spec.rename_column(o, i)
            spec['var'] = tmpvar.copy()

        """ Save the fits file """
        hdulist.writeto(outfile, overwrite=True)

    # -----------------------------------------------------------------------

    def coadd(self, other, outfile=None):
        """

        Coadds the spectra in this slitmask with the corresponding ones
        from other exposures using the same slitmask.

        Inputs:
         other   - either a single DeimosMask1d object or a list of
                   DeimosMask1d objects
         outfile - name for output file

        """

        """ Standardize the format of the input data """
        if not isinstance(other, (list, tuple, np.ndarray)):
            other = [other,]

        """ Match the objects in the masks """
        print('')
        print('Extracted spectra in exposure 1: %d' % self.nspec)
        if self.hdr0 is not None:
            exptime = self.hdr0['exptime']
        else:
            exptime = None
        for i, exp_i in enumerate(other):
            print('Extracted spectra in exposure %d: %d' %
                  ((i+2), exp_i.nspec))
            if exptime is not None:
                if exp_i.hdr0 is not None:
                    exptime += exp_i.hdr0['exptime']
            if i == 0:
                matchtab = join(self.slitinfo, exp_i.slitinfo,
                                keys=['det', 'slitid'])
                diff = matchtab['spatloc_1'] - matchtab['spatloc_2']
            else:
                matchtab = join(matchtab, exp_i.slitinfo,
                                keys=['det', 'slitid'])
                diff = matchtab['spatloc_1'] - matchtab['spatloc']
                for name in ['objid', 'spatloc', 'fwhm']:
                    matchtab.rename_column(name, '%s_%d' % (name, (i+2)))
            mask = np.abs(diff) <= 3.
            matchtab = matchtab[mask]
        print('')
        print('Number of matched spectra in all %d exposures: %d'
              % (len(other)+1, len(matchtab)))

        """ Coadd the matched spectra """
        outspec = {}
        for info in matchtab:
            specid = '%d_%d_%d_%d' % \
                (info['det'], info['slitid'], info['objid_1'],
                 info['spatloc_1'])
            speclist = [self[specid],]
            for i, exp_i in enumerate(other):
                j = i + 2
                specid_i = '%d_%d_%d_%d' % \
                    (info['det'], info['slitid'], info['objid_%d' % j],
                     info['spatloc_%d' % j])
                # ADD RESAMPLING SINCE WAVELENGTH SCALES ARE NOT IDENTICAL
                speclist.append(exp_i[specid_i])
            ss = specset1d.SpecSet1d(spec1dlist=speclist)
            outspec[specid] = ss.coadd(doplot=False, verbose=False)
            # print(specid)

        """
        Convert the dictionary of coadded spectra into a DeimosMask1d object
        and write it to an output file if requested
        """
        for name in ['objid', 'spatloc', 'fwhm']:
            matchtab.rename_column('%s_1' % name, name)
        outtab = matchtab['det', 'slitid', 'objid', 'spatloc', 'fwhm']
        outmask = DeimosMask1d(outtab, specdict=outspec)
        if self.hdr0 is not None:
            outmask.hdr0 = self.hdr0
            if exptime is not None:
                outmask.hdr0['exptime'] = exptime
        if outfile is not None:
            outmask.save(outfile)
            print('Saved %d coadded spectra in %s' % (len(matchtab), outfile))

        del(outspec, matchtab, outtab)
        return outmask
