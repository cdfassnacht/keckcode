"""

osimset.py

Defines the OsImSet class, which is just a specialized version of the
CCDset class, with a few add-ons that are specific to the OSIRIS imager

"""
import os
import sys
import numpy as np
from math import sqrt

from cdfutils import coords
from specim.imfuncs.wcshdu import WcsHDU
from ccdredux.ccdset import CCDSet

pyversion = sys.version_info.major

# ---------------------------------------------------------------------------


class OsImSet(CCDSet):
    """
    Specialized version of the generic CCDSet class, with a few specific
    add-ons that are specific to the Keck OSIRIS imager
    """

    def __init__(self, inlist, wcstype='koa', is_sci=True, texpkey='truitime',
                 gainkey='sysgain', indir=None, obsdate=None, **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nOsImSet: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Create the input filelist from the passed parameters """
        if isinstance(inlist, dict):
            filelist = self.make_filelist([inlist], obsdate, indir=indir)
        elif isinstance(inlist[0], dict):
            filelist = self.make_filelist(inlist, obsdate, indir=indir)
        else:
            """
            For any other data types, let CCDSet (called through the "super"
            calls below) do the type checking.
            """
            filelist = inlist

        """ Set up the empty CCDSet container by calling the superclass """
        if pyversion == 2:
            super(OsImSet, self).__init__(filelist, texpkey=texpkey,
                                         gainkey=gainkey, indir=indir, **kwargs)
        else:
            super().__init__(filelist, texpkey=texpkey, gainkey=gainkey,
                             indir=indir, **kwargs)

        """
        Copy some osiris-specific keywords into the more standard versions
        """
        for hdu in self:
            hdr = hdu.header
            if texpkey.upper() in hdr.keys():
                hdr['exptime'] = hdr[texpkey]

        """ Set up some default values """

        """
        Fix the WCS header info if necessary.
        The KOA processing adds WCS keywords based on the old version of the
         detector and so they don't reflect the new detector parameters
         (still valid as of Oct 2022)
        NOTE: This doesn't correct for the flip, since that is taken care of
         in the osim_funcs.py pipeline.
        """
        # print(wcstype)
        if wcstype is None:
            pass
            # for hdu in self:
            #     hdu.pixscale = 0.01
        elif wcstype == 'raw' and is_sci:
            for hdu in self:
                hdr = hdu.header
                hdr['cdelt1'] = -0.01 / 3600.
                hdr['cdelt2'] = 0.01 / 3600.
                hdr['crval1'] = hdr['ra']
                hdr['crval2'] = hdr['dec']
                hdr['crpix1'] = hdr['naxis1'] / 2. + 0.5
                hdr['crpix2'] = hdr['naxis2'] / 2. + 0.5
                hdr['ctype1'] = 'RA---TAN'
                hdr['ctype2'] = 'DEC--TAN'
                hdu.read_wcsinfo(hdr)
                hdu.pixscale = 0.01
                if 'PA_IMAG' in hdr.keys():
                    hdu.pc = coords.rot_to_pcmatrix(-1. * hdr['pa_imag'],
                                                    verbose=False)
                hdu.crpix = (np.array(hdu.data.shape)[::-1]) / 2. + 0.5
        elif wcstype == 'koa' and is_sci:
            for hdu in self:
                hdu.pixscale = 0.01
                if 'PA_IMAG' in hdr.keys():
                    hdu.pc = coords.rot_to_pcmatrix(-1. * hdr['pa_imag'],
                                                    verbose=False)
                # print(type(hdu.pc))
                hdu.crpix = (np.array(hdu.data.shape)[::-1]) / 2. + 0.5
        elif is_sci:
            for hdu in self:
                hdu.pixscale = 0.01
                if 'PA_IMAG' in hdr.keys():
                    hdu.pc = coords.rot_to_pcmatrix(-1. * hdr['pa_imag'],
                                                    verbose=False)

    #  ------------------------------------------------------------------------

    def make_filelist(self, assnlist, obsdate, indir='auto'):
        """

        Makes a list of file names based on an input directory and an OSIRIS
         association number.
        If the frames parameter is None (the default), then the filelist will
         contain all exposures in the association.  If, however, it is set to
         be a list, tuple, or ndarray, then the output filelist will only
         contain the desired frame numbers.

        """

        """ Make the date string from the provided obsdate """
        if obsdate is not None:
            datestr = '%s_%s_%s' % (obsdate[:4], obsdate[4:6], obsdate[6:8])
        else:
            raise TypeError('\nmake_filelist: an obsdate that is not None '
                            'must be provided\n\n')
        if indir is None:
            indir = '.'
        elif indir == 'auto':
            indir = os.path.join(os.getenv('sharpdat'), 'Raw', datestr)
        else:
            pass

        """ Create a filelist from the inputs """
        filelist = []
        print(indir)
        for i in assnlist:
            """ Check the passed parameters """
            for k in ['assn', 'frames']:
                if k not in i.keys():
                    raise KeyError('\nmake_filelist: dict must contain both'
                                   '"assn" and "frames" keys\n\n')
            """
            Loop through association and frame numbers to create file list
            """
            assn = i['assn']
            for j in i['frames']:
                filebase = 'i%s_a%03d%03d' % (obsdate[2:], assn, j)
                filelist.append(os.path.join(indir, '%s.fits' % filebase))

        return filelist

    #  ------------------------------------------------------------------------

    def make_outlist(self, intext, outtext, outdir=None):
        """
        Makes an output filelist from the input filelist (if one exists) by
        replacing the appropriate part of each input filename with the
        desired output file designator.
        For example, for an input list of ff*fits this method could produce
        corresponding bgsub*fits or ff*_wht.fits.
        """

        outfiles = []
        for f in self.datainfo['basename']:
            if outdir is not None:
                outf = os.path.join(outdir, f.replace(intext, outtext))
            else:
                outf = f.replace(intext, outtext)
            outfiles.append(outf)
        return outfiles

    #  ------------------------------------------------------------------------

    def make_kai_log(self, outfile, outdir=None):
        """

        Makes a text file containing observing log information derived from
        the headers of the fits files.  The format of the log is that used
        by the KAI data reduction pipeline.

        Inputs:
          outfile - Name of output text file
          outdir  - Location for output text file.  The default (None) will
                    save the text file in the current working directory

        """

    #  ------------------------------------------------------------------------

    def make_var_and_bpm(self, flatfile, inpref, bpmglobfile=None, caldir=None,
                         bpmpref='bpm', texppref='texp', onespref='ones',
                         nsig=10., outdir=None, verbose=True):
        """
        Creates the variance files and bad pixel masks associated with each
        image in this OsImSet object.
        Additionally, makes exposure-time maps and datasets that are
         identically equal to 1.0 everywhere, in both cases with WCS
         information copied from their associated science images, if
         these extra files are requested.
        """

        """ Load the flat-field file and global bad pixel mask (if any) """
        self.load_calib(bpmglobfile=bpmglobfile, flatfile=flatfile,
                        caldir=caldir, verbose=verbose, headverbose=False)

        """ Set up the file names for the output files """
        varfiles = self.make_outlist(inpref, 'var', outdir=outdir)
        bpmfiles = self.make_outlist(inpref, bpmpref, outdir=outdir)
        if texppref is not None:
            texpfiles = self.make_outlist(inpref, texppref, outdir=outdir)
        if onespref is not None:
            onesfiles = self.make_outlist(inpref, onespref, outdir=outdir)

        if verbose:
            print('')
        for i, f in enumerate(self):
            """
            Make the exposure time map and, if requested, save it
            """
            texp = f.make_texp('exptime')
            if texppref is not None:
                texp.writeto(texpfiles[i])
                if verbose:
                    print('Wrote exposure time image: %s' % texpfiles[i])
            """
            Make the variance frames by applying the flat-field and exposure
            time corrections to the data a second time
            """
            v = f.copy()
            v /= texp
            v /= self.flat
            v.sigma_clip()
            v.writeto(varfiles[i])
            if verbose:
                print('Variance file mean = %5.3f  ==> rms = %5.3f' %
                      (v.mean_clip, sqrt(v.mean_clip)))
                print('Wrote variance image: %s' % varfiles[i])
            """ Make the bad pixel mask for this individual data set """
            tmpbpm = f.make_bpm('sci', var=v, nsig=nsig)
            bpm = WcsHDU(self.bpmglobal * tmpbpm, wcsverb=False)
            bpm.header['object'] = 'Bad Pixel Mask for %s' % f.infile
            bpm.writeto(bpmfiles[i])
            if verbose:
                print('Wrote bad pixel mask: %s' % bpmfiles[i])
            """ Make the ones file if requested """
            if onespref is not None:
                ones = f.make_ones()
                # ones.writeto(onesfiles[i])
                # if verbose:
                #     print('Wrote ones image: %s' % onesfiles[i])
            else:
                ones = 1.
            if verbose:
                print('')
            del texp, v, bpm, tmpbpm, ones
