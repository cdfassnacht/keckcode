import sys
import os

from pyraf import iraf as ir

from kai import instruments
from kai.reduce import data
from kai.reduce import kai_util
from ..ao_img.aoset import AOSet

""" Define global variables for the two possible instruments """
osiris = instruments.OSIRIS()
nirc2 = instruments.NIRC2()

pyversion = sys.version_info.major


class KaiSet(AOSet):
    """

    Class to run KAI functions on related sets of files.
    Several of the early steps are taken care of by methods in the AOSet
     class, which is the parent class

    """

    def __init__(self, inlist, inst, obsdate, indir=None, gzip=False,
                 verbose=True, **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nKaiSet: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Set up the KaiSet container by calling the superclass """
        if pyversion == 2:
            super(KaiSet, self).__init__(inlist, inst, obsdate, indir=indir,
                                         gzip=gzip, verbose=verbose, **kwargs)
        else:
            super().__init__(inlist, inst, obsdate, indir=indir, gzip=gzip,
                             verbose=verbose, **kwargs)

        """ Get the instrument in KAI format """
        self.inst = None
        try:
            self.get_instrument(inst)
        except ValueError:
            print('')
            print('Could not create kaiset object')
            print('')
            return

    #  ------------------------------------------------------------------------

    def get_instrument(self, instrument):
        """

        Makes sure that either NIRC2 or OSIRIS has been selected

        """
        if instrument.lower() == 'osiris' or instrument.lower() == 'osim':
            self.inst = osiris
        elif instrument.lower() == 'nirc2':
            self.inst = nirc2
        else:
            print('')
            raise ValueError('get_instrument: instrument must be '
                             '"osiris" or "nirc2"\n')

    #  ------------------------------------------------------------------------

    def add_def_hdrinfo(self):
        """

        Adds keywords that are needed for later procrssing to the headers of
         the images

        """

        """ Loop through the images """
        for hdu in self:

            """ Get the central wavelength of the filter being used """
            hdr = hdu.header
            defwave = self.inst.get_central_wavelength(hdr)

            """ Add the new wavelength header cards """
            hdr['effwave'] = defwave
            hdr['cenwave'] = defwave
            hdr['camname'] = 'narrow'

            try:
                nonlinSky = hdr['skylev']
            except KeyError:
                nonlinSky = 0.
            coaddkey = self.inst.hdr_keys['coadds']
            try:
                coadds = hdr[coaddkey]
            except KeyError:
                coadds = 1
            satLevel = (coadds * self.inst.get_saturation_level()) - nonlinSky
            hdr['satlevel'] = satLevel

    #  ------------------------------------------------------------------------

    def clean_cosmicrays(self, obsfilt, outlist, filelist=None, verbose=True):
        """

        Cleans cosmic rays from the images by (eventually) calling the iraf
        task noao.imred.crutils.cosmicrays.  The wrapper for that call is
        the clean_cosmicrays function in the KAI data.py code.

        Inputs:
         obsfilt  - Filter used for the observations, e.g., 'Kp'
         outlist  - List of output files to store the cosmic ray masks
         filelist - The default behavior, which is executed if filelist is None,
                     is to call the data.py code with the filenames stored
                     in this object's datainfo['basename'] column.  However,
                     the user can provide a list of filenames to use instead.
                     In that case, filelist should be a list object containing
                     strings that are the input filenames.
        """

        """ Check that outlist has the proper length """
        if len(outlist) != self.nfiles:
            raise IndexError('clean_cosmicrays: outlist length does not match'
                             ' number of input files')

        """ Set the list of input filenames """
        if filelist is not None:
            if len(filelist) != len(outlist):
                raise IndexError(
                    'clean_cosmicrays: input list length does not match'
                    ' number of output files')
            inlist = filelist
        else:
            inlist = self.datainfo['basename']

        """ Loop through the files, creating cosmic ray masks """
        if verbose:
            print('')
            print('Creating cosmic ray masks')
            print('-------------------------')
        for i, j in zip(inlist, outlist):
            if verbose:
                print('%s ---> %s' % (i, j))
            if os.path.isfile(j):
                os.remove(j)
            data.clean_cosmicrays(i, j, obsfilt.lower())

    #  ------------------------------------------------------------------------

    def dewarp(self, xdistmap, ydistmap, inpref, outpref, whtpref, drizlog,
               fixDAR=True, use_koa_weather=False):
        """

        Corrects for the optical distortion of the files, i.e., "dewarps" them,
        by drizzling them using the provided distortion map.

        Inputs:
        xdistmap - x distortions as a function of position
        ydistmap - y distortions as a function of position
        inpref   - prefix in the input filenames to be replaced in the output
                   filenames
        outpref  - prefix to be used for the output files
        whtpref  - prefix to be used for the output weight files
        drizlog  - name of the logfile to store the text generated by the
                   drizzle process
        fixDAR   - Fix differential atmospheric distortion? Default is True
        use_koa_weather - ?

        """

    #  ------------------------------------------------------------------------

    def make_coo(self, refSrc, strSrc, check_loc=True, inpref='ce', outpref='c',
                 outdir=None, update_fits=True, cent_box=12,
                 update_from_AO=True):

        """ Get reference AO stage position from first image """
        hdr0 = self[0].header
        aotsxyRef = [float(hdr0['aotsx']), float(hdr0['aotsy'])]

        """ Set up table columns to store info """
        self.datainfo['scale'] = 0.0
        self.datainfo['inst_angle'] = 0.0
        self.datainfo['phi'] = 0.0
        self.datainfo['aotsx'] = 0.0
        self.datainfo['aotsy'] = 0.0
        self.datainfo['dx'] = 0.0
        self.datainfo['dy'] = 0.0
        self.datainfo['xref'] = 0.0
        self.datainfo['yref'] = 0.0

        """ Set up output files that have registration information in them """
        cfiles = self.make_outlist(inpref, outpref, outdir=outdir)

        for hdu, info, outfile in zip(self, self.datainfo, cfiles):
            hdr = hdu.header

            radec = [float(hdr['RA']), float(hdr['DEC'])]
            aotsxy = [float(hdr['aotsx']), float(hdr['aotsy'])]
            # aotsxy = kai_util.getAotsxy(hdr) # This line gave a weird error
            info['aotsx'] = aotsxy[0]
            info['aotsy'] = aotsxy[1]

            """ Determine the image's PA and plate scale """
            phi = self.inst.get_position_angle(hdr)
            scale = self.inst.get_plate_scale(hdr)
            info['phi'] = phi
            info['scale'] = scale

            """ Determine the instrument angle w.r.t. the AO bench. """
            inst_angle = self.inst.get_instrument_angle(hdr)
            info['inst_angle'] = inst_angle

            """ Calculate the pixel offsets from the reference image """
            # d_xy = kai_util.radec2pix(radec, phi, scale, radecRef)
            if update_from_AO:
                d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef,
                                           inst_angle=inst_angle)
            else:
                d_xy = [0, 0]
            info['dx'] = d_xy[0]
            info['dy'] = d_xy[1]

            """ In the new image, find the REF and STRL coords """
            xref = refSrc[0] + d_xy[0]
            yref = refSrc[1] + d_xy[1]
            xstr = strSrc[0] + d_xy[0]
            ystr = strSrc[1] + d_xy[1]
            info['xref'] = xref
            info['yref'] = yref

            print('makecoo: xref, yref start = %6.2f %6.2f' % (xref, yref))

            """ Re-center stars to get exact coordinates """
            infile = info['infile']
            if check_loc:
                text = ir.imcntr(infile, xref, yref, cbox=cent_box, Stdout=1)
                values = text[0].split()
                xref = float(values[2])
                yref = float(values[4])

                text = ir.imcntr(infile, xstr, ystr, cbox=cent_box, Stdout=1)
                values = text[0].split()
                xstr = float(values[2])
                ystr = float(values[4])
                print('clean_makecoo: xref, yref final = %6.2f %6.2f}'
                      % (xref, yref))

            """
            Write reference star x,y to fits header and save info in a new
            file
            """
            if update_fits:
                hdr.set('XREF', "%.3f" % xref, 'Cross Corr Reference Src x')
                hdr.set('YREF', "%.3f" % yref, 'Cross Corr Reference Src y')
                hdr.set('XSTREHL', "%.3f" % xstr, 'Strehl Reference Src x')
                hdr.set('YSTREHL', "%.3f" % ystr, 'Strehl Reference Src y')
                hdu.save(outfile)

            """ Also put the xref, yref into into a *.coo file """
            outcoo = outfile.replace('.fits', '.coo')
            with open(outcoo, 'w') as f:
                f.write('%7.2f  %7.2f\n' % (xref, yref))

            """
            Make a temporary rotated coo file, in case there are any data sets
            with various PAs; needed for xregister; remove later
            """
            xyRef_rot = kai_util.rotate_coo(xref, yref, phi)
            xref_r = xyRef_rot[0]
            yref_r = xyRef_rot[1]

            # xyStr_rot = kai_util.rotate_coo(xstr, ystr, phi)
            # xstr_r = xyStr_rot[0]
            # ystr_r = xyStr_rot[1]

            outrcoo = outfile.replace('.fits', '.rcoo')
            with open(outrcoo, 'w') as f:
                f.write('%7.2f  %7.2f\n' % (xref_r, yref_r))
