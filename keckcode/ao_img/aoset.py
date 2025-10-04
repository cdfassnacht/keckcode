import os
import sys
import math
import numpy as np

from cdfutils import coords
from ccdredux.ccdset import CCDSet

pyversion = sys.version_info.major


class AOSet(CCDSet):
    """

    Class to run functions reflecting specifics of AO/NIR imaging on related
     sets of files
    A lot of this class is just setting things up to call methods in the
     CCDSet class, which is the parent class

    """

    def __init__(self, inlist, instrument, obsdate=None, indir=None, gzip=False,
                 frameroot='default', wcstype=None, is_sci=True, verbose=True,
                 **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nAOSet: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Make the date string from the provided obsdate """
        if obsdate is not None:
            self.obsdate = obsdate
            datestr = '%s_%s_%s' % (obsdate[:4], obsdate[4:6], obsdate[6:8])
        else:
            self.obsdate = None
            datestr = None

        """ Set up the SHARP-specific input directory if requested """
        if indir is None:
            indir = '.'
        elif indir == 'auto':
            if datestr is not None:
                indir = os.path.join(os.getenv('sharpdat'), 'Raw', datestr)
            else:
                raise ValueError('If choosing "auto" for indir then obdate '
                                 'must be provided')
        else:
            pass
        if verbose:
            print('Reading files from %s' % indir)

        """
        Set instrument-specific parameters that are needed for make_filelist
        and for the superclass
        """
        if instrument == 'osiris' or instrument == 'osim':
            self.instrument = 'osiris'
            texpkey = 'truitime'
            gainkey = 'sysgain'
        elif instrument == 'nirc2':
            self.instrument = 'nirc2'
            texpkey = 'elaptime'
            gainkey = 'gain'
        else:
            raise ValueError('Instrument must be either osiris or nirc2')

        """ Create the input filelist from the passed parameters """
        if isinstance(inlist, dict):
            filelist = self.make_filelist([inlist], obsdate,
                                          frameroot=frameroot, gzip=gzip)
        elif isinstance(inlist[0], dict):
            filelist = self.make_filelist(inlist, obsdate, frameroot=frameroot,
                                          gzip=gzip)
        else:
            """
            For any other data types, let CCDSet (called through the "super"
            calls below) do the type checking.
            """
            filelist = inlist

        """ Set up the AOSet container by calling the superclass """
        if pyversion == 2:
            super(AOSet, self).__init__(filelist, texpkey=texpkey,
                                        gainkey=gainkey, indir=indir, **kwargs)
        else:
            super().__init__(filelist, texpkey=texpkey, gainkey=gainkey,
                             indir=indir, **kwargs)

        """ Set some default values """
        self.caldir = None
        self.darkdir = None
        self.flatdir = None
        self.maskdir = None
        self.reduxdir = None

        """ Set instrument-specific values """
        if instrument == 'osiris' or instrument == 'osim':
            """ For OSIRIS, mask out left-most 300 columns """
            self.instmask = np.ones(self[0].shape)
            self.instmask[:, :300] = 0
        elif instrument == 'nirc2':
            self.instmask = None

        """
        Copy some keywords into the more standard versions
        """
        for hdu in self:
            hdr = hdu.header
            if texpkey.upper() in hdr.keys():
                hdr['exptime'] = hdr[texpkey]

        """
        For OSIRIS data, fix the WCS header info if necessary.
        The KOA processing adds WCS keywords based on the old version of the
         OSIRIS detector and so they don't reflect the new detector parameters
         (still valid as of Oct 2022)
        NOTE: This doesn't correct for the flip, since that is taken care of
         in the osim_funcs.py pipeline.
        """
        # print(wcstype)
        if wcstype is None or self.instrument != 'osiris':
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
                hdr = hdu.header
                hdu.pixscale = 0.01
                if 'PA_IMAG' in hdr.keys():
                    hdu.pc = coords.rot_to_pcmatrix(-1. * hdr['pa_imag'],
                                                    verbose=False)
                # print(type(hdu.pc))
                hdu.crpix = (np.array(hdu.data.shape)[::-1]) / 2. + 0.5
        elif is_sci:
            for hdu in self:
                hdr = hdu.header
                hdu.pixscale = 0.01
                if 'PA_IMAG' in hdr.keys():
                    hdu.pc = coords.rot_to_pcmatrix(-1. * hdr['pa_imag'],
                                                    verbose=False)

    #  ------------------------------------------------------------------------

    @staticmethod
    def make_filelist_osim(assnlist, obsdate, frameroot='default', suff='fits'):
        """

        Makes a list of file names based on an input directory and an OSIRIS
         association number.
        If the frames parameter is None (the default), then the filelist will
         contain all exposures in the association.  If, however, it is set to
         be a list, tuple, or ndarray, then the output filelist will only
         contain the desired frame numbers.

        """

        """ Create a filelist from the inputs """
        filelist = []
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
                if frameroot is not None:
                    if frameroot == 'default':
                        filebase = 'i%s_a%03d%03d' % (obsdate[2:], assn, j)
                    else:
                        filebase = '%s%03d%03d' % (frameroot, assn, j)
                else:
                    filebase = '%03d%03d' % (assn, j)
                filelist.append('%s.%s' % (filebase, suff))
        # print(filelist)
        return filelist

    #  ------------------------------------------------------------------------

    @staticmethod
    def make_filelist_nirc2(inlist, frameroot='default', suff='fits'):

        """ Create a filelist from the inputs """
        filelist = []
        for i in inlist:
            """ Check the passed parameters """
            for k in ['frames']:
                if k not in i.keys():
                    raise KeyError('\nmake_filelist: dict must contain a'
                                   '"frames" key\n\n')
            """
            Loop through frame numbers to create file list
            """
            for j in i['frames']:
                filebase = 'n%04d' % j
                filelist.append('%s.%s' % (filebase, suff))

        return filelist

    #  ------------------------------------------------------------------------

    def make_filelist(self, inlist, obsdate, frameroot='default', gzip=False):
        """

        Makes a list of file names based on an input directory and an OSIRIS
         association number.
        If the frames parameter is None (the default), then the filelist will
         contain all exposures in the association.  If, however, it is set to
         be a list, tuple, or ndarray, then the output filelist will only
         contain the desired frame numbers.

        """

        """ Set the file extension """
        if gzip:
            suff = 'fits.gz'
        else:
            suff = 'fits'

        """ Get the filelist, which depends on the instrument being used """
        if self.instrument == 'osiris':
            filelist = self.make_filelist_osim(inlist, obsdate, suff=suff,
                                               frameroot=frameroot)
        else:
            filelist = self.make_filelist_nirc2(inlist, frameroot=frameroot,
                                                suff=suff)

        return filelist

    #  ------------------------------------------------------------------------

    def set_kai_dirs(self, obsfilt=None):
        """

        Sets names of directories in the KAI data reduction pipeline standard
        directory structure

        """

        """ Set up directory names """
        # pwd = os.getcwd()
        # self.caldir = os.path.join(pwd, 'calib')
        self.caldir = 'calib'
        self.darkdir = os.path.join(self.caldir, 'darks')
        self.flatdir = os.path.join(self.caldir, 'flats')
        self.maskdir = os.path.join(self.caldir, 'masks')
        if obsfilt is not None:
            self.reduxdir = os.path.join(obsfilt, 'sci_%s' % self.obsdate)

        for d in [self.caldir, self.darkdir, self.flatdir, self.maskdir]:
            if not os.path.isdir(d):
                os.makedirs(d)
        if self.reduxdir is not None:
            if not os.path.isdir(obsfilt):
                os.makedirs(obsfilt)
            if not os.path.isdir(self.reduxdir):
                os.makedirs(self.reduxdir)

    #  ------------------------------------------------------------------------

    def set_caldirs(self, caldir=None, obsfilt=None):
        """

        Sets the locations, for either input or output, of the calibration
        and reduction directories.  There are three options for the caldir
        parameter:
           None         - in this case, use the current directory
           'kaidefault' - set things up with the directory structure used for
                          the KAI data reduction pipeline
           [any other string] - All files will be in the directory designated
                                by the string

        """

        if caldir is None:
            self.caldir = '.'
            self.darkdir = '.'
            self.flatdir = '.'
            self.maskdir = '.'
        elif caldir == 'kaidefault':
            self.set_kai_dirs(obsfilt=obsfilt)
        elif isinstance(caldir, str):
            self.caldir = caldir
            self.darkdir = caldir
            self.flatdir = caldir
            self.maskdir = caldir
        else:
            print('')
            raise TypeError('caldir must be either a string or None')

    #  ------------------------------------------------------------------------

    def create_dark(self, outname, caldir=None, reject='sigclip',
                    nlow=1, nhigh=1):
        """

        Creates a dark-frame using the ccdredux make_dark functionality

        """

        """ Set up output directory """
        if self.darkdir is not None:
            darkdir = self.darkdir
        else:
            self.set_caldirs(caldir=caldir)
            darkdir = self.darkdir

        """ Get the output filenames """
        if outname[-4:] != 'fits':
            outfile = os.path.join(darkdir, '%s.fits' % outname)
        else:
            outfile = os.path.join(darkdir, outname)
        outlist = os.path.join(darkdir, 'dark.lis')

        """ """
        self.make_bias(outfile=outfile, reject=reject, nlow=nlow, nhigh=nhigh)

    #  ------------------------------------------------------------------------

    def create_flat(self, outname, lamps_off=None, caldir=None, normalize=None,
                    indark=None, inflat=None, reject='minmax', nlow=1, nhigh=1,
                    **kwargs):
        """

        Creates a flat-field frame following the KAI recipe.  This approach
        uses frames taken with the flat-field lamps on, which are part of this
        AOSet instance, and possibly frames taken with the flat-field lamps
        off, which are provided via the lamps_off parameter.  The lamps_off
        data should also be a AOSet instance if it is not None.

        """

        """ Check validity of input """
        if lamps_off is not None:
            if not isinstance(lamps_off, AOSet):
                raise ValueError('\nThe lamps_off parameter must be an AOSet\n')
            if len(lamps_off) != len(self):
                raise IndexError('\nDifferent number of files in lamps_off and'
                                 ' lamps_on\n')

        """ Set up directory names """
        if self.flatdir is None:
            self.set_caldirs(caldir=caldir)
        if indark is not None:
            dark = os.path.join(self.darkdir, indark)
        else:
            dark = None

        """ Get the output filenames """
        if outname[-4:] != 'fits':
            outfile = os.path.join(self.flatdir, '%s.fits' % outname)
        else:
            outfile = os.path.join(self.flatdir, outname)
        onfits = os.path.join(self.flatdir, 'lampsOn.fits')
        offfits = os.path.join(self.flatdir, 'lampsOff.fits')
        normfits = os.path.join(self.flatdir, 'flatNotNorm.fits')
        onlist = os.path.join(self.flatdir, 'on.lis')
        offlist = os.path.join(self.flatdir, 'off.lis')
        normlist = os.path.join(self.flatdir, 'onNorm.lis')

        """
        Make the flat, using the lamps-off frames if provided, but otherwise
        using just the lamps-on frames
        """
        if lamps_off is not None:
            """ First make the combined lamps-off and lamps-on frames """
            lamps_off.make_flat(outfile=offfits, normalize=normalize,
                                bias=dark, flat=inflat,
                                reject=reject, nlow=nlow, nhigh=nhigh, **kwargs)
            self.make_flat(outfile=onfits, normalize=normalize, bias=dark,
                           flat=inflat, reject=reject, nlow=nlow,
                           nhigh=nhigh, **kwargs)

            """ Now loop through paired on/off exposures, taking differences """
            nfile = open(onlist, 'w')
            ffile = open(offlist, 'w')
            count = 0
            for n, f in zip(self, lamps_off):
                """ Add the filenames to the onlist and offlist"""
                if 'basename' in self.datainfo.keys():
                    nfile.write('%s\n' % self.datainfo['basename'][count])
                if 'basename' in lamps_off.datainfo.keys():
                    ffile.write('%s\n' % lamps_off.datainfo['basename'][count])
                count += 1

                """ Replace the lamps-on data with the difference data """
                n.data = n.data - f.data
            nfile.close()
            ffile.close()

            """ Make the final flat """
            self.make_flat(outfile=outfile, normalize=normalize, reject=reject,
                           nlow=nlow, nhigh=nhigh, **kwargs)

        else:
            """ Put the filenames into the onlist """
            if 'basename' in self.datainfo.keys():
                nfile = open(onlist, 'w')
                for info in self.datainfo:
                    nfile.write('%s\n' % info['basename'])

            """ Make the final flat """
            self.make_flat(outfile=outfile, normalize=normalize,
                           bias=dark, flat=inflat, reject=reject,
                           nlow=nlow, nhigh=nhigh, **kwargs)

    #  ------------------------------------------------------------------------

    def create_sky(self, outname, obsfilt, skyscale=True,
                   inflat=None, reject='sigclip', nlow=1, nhigh=1):
        """

        Makes a sky file, broadly following the KAI algorithm, but with the
         possibility of applying the flat-field file to the sky frames
         before combining them.

         NOT FINISHED YET.  DON'T USE!
        """

        """ Set up directory names """
        pwd = os.getcwd()
        caldir = os.path.join(pwd, 'calib')
        flatdir = os.path.join(caldir, 'flats')
        skydir = os.path.join(caldir, obsfilt, 'sky')
        # skydir = os.path.join(caldir, 'sky_%s' % self.obsdate)
        if not os.path.isdir(caldir):
            os.makedirs(caldir)
        if not os.path.isdir(skydir):
            os.makedirs(skydir)

        """ Get full path to inflat if it is not None """
        if inflat is not None:
            flatfile = os.path.join(flatdir, 'flat_')

    #  ------------------------------------------------------------------------

    @staticmethod
    def get_scale(instrument, hdr):
        """
        Set the plate scale, depending on which instrument is being used
        """
        if instrument == 'osiris':
            scale = 0.00995
        elif instrument == 'nirc2':
            scales = {"narrow": 0.009952,
                      "medium": 0.019829,
                      "wide": 0.039686}

            scale = scales[hdr['CAMNAME']]
        else:
            raise ValueError('ERROR: get_scale. instrument must be "osiris"'
                             ' or "nirc2"')

        return scale

    #  ------------------------------------------------------------------------

    @staticmethod
    def get_pa_from_rotposn(instrument, hdr):
        """

        Gets the position angle of the detector with respect to the north axis
        based on the value of the ROTPOSN keyword

        """

        """ Make sure that required keys are in the header """
        keylist = ['date-obs', 'rotposn', 'instangl']
        for k in keylist:
            if k.upper() not in hdr.keys():
                raise KeyError('ERROR: get_pa_from_rotposn: %s not found in '
                               ' header.' % k)

        """ 
        Set the value of the offset between the imager and the spectrograph,
         if using OSIRIS 
        """
        if instrument == 'osiris':
            year = int(hdr['date-obs'][:4])
            month = int(hdr['date-obs'][5:7])
            if year < 2022 or (year == 2022 and month < 5):
                imager_offset = 42.5
            else:
                imager_offset = 43.4
        else:
            imager_offset = 0.

        """ Determine the instrument PA and imager instrument angle"""
        phi = hdr['rotposn'] - hdr['instangl'] + imager_offset
        return phi, (hdr['instangl'] - imager_offset)

    #  ------------------------------------------------------------------------

    @staticmethod
    def aotsxy2pix(aotsxy, scale, aotsxyRef, inst_angle=0.0):
        """
        Determine pixel shifts from AOTSX and AOTSY positions.
        :param aotsxy:
        :param scale:
        :param aotsxyRef:
        :param inst_angle:
        :return:
        """

        x = aotsxy[0]
        y = aotsxy[1]

        """ AOTSX,Y are in units of mm. Conversion is 0.727 mm/arcsec """
        d_x = (x - aotsxyRef[0]) / 0.727
        d_y = (y - aotsxyRef[1]) / 0.727
        d_x = d_x * (1.0 / scale)
        d_y = d_y * (1.0 / scale)

        """ Rotate to the instrument PA """
        cosa = np.cos(np.radians(-inst_angle))
        sina = np.sin(np.radians(-inst_angle))

        rot_matrix = np.array([[cosa, sina], [-sina, cosa]])
        coo_ao = np.array([d_x, d_y])
        coo_inst = rot_matrix.dot(coo_ao)

        d_x = coo_inst[0]
        d_y = coo_inst[1]

        return [d_x, d_y]

    #  ------------------------------------------------------------------------

    def calc_offsets_aots(self, refpos):
        """

        Calculates offsets between the first images and the others, in pixels,
        using the values in the AOTSX and AOTSY header keywords in the files

        Inputs:
         refpos - (x, y) position in the first image of a relatively bright and
          compact object such as a star or a lensing galaxy


        """

        """ Set up table columns to store info """
        self.datainfo['inst_angle'] = 0.0
        self.datainfo['phi'] = 0.0
        self.datainfo['aotsx'] = 0.0
        self.datainfo['aotsy'] = 0.0
        self.datainfo['dx'] = 0.0
        self.datainfo['dy'] = 0.0
        self.datainfo['xref0'] = 0.0
        self.datainfo['yref0'] = 0.0

        """ Get reference AO stage position from the first image """
        hdr0 = self[0].header
        if 'AOTSX' in hdr0.keys() and 'AOTSY' in hdr0.keys():
            aotsxyRef = [float(hdr0['aotsx']), float(hdr0['aotsy'])]
        else:
            raise KeyError('ERROR: calc_offsets_aots.  Missing AOTSX and/or '
                           'AOTSY keys in header for %s' %
                           self.datainfo['basename'][0])

        """ Loop through the files, finding the positions of the ref source """
        for hdu, info in zip(self, self.datainfo):

            """ Check to make sure that all keywords are in the header """
            hdr = hdu.header
            keylist = ['date-obs', 'aotsx', 'aotsy', 'rotposn', 'instangl']
            for k in keylist:
                if k.upper() not in hdr.keys():
                    raise KeyError('%s not found in %s header.' %
                                   (k, info['basename']))

            """ Get the AO stage positions for the current image """
            aotsxy = [float(hdr['aotsx']), float(hdr['aotsy'])]
            info['aotsx'] = aotsxy[0]
            info['aotsy'] = aotsxy[1]

            """ Determine the imager PA and plate scale """
            phi, inst_angle = self.get_pa_from_rotposn(self.instrument, hdr)
            scale = self.get_scale(self.instrument, hdr)
            info['phi'] = phi

            """ Calculate the pixel offsets from the reference image """
            d_xy = self.aotsxy2pix(aotsxy, scale, aotsxyRef,
                                   inst_angle=inst_angle)
            info['dx'] = d_xy[0]
            info['dy'] = d_xy[1]

            """ In the new image, find the REF coords """
            xref = refpos[0] + d_xy[0]
            yref = refpos[1] + d_xy[1]
            info['xref0'] = xref
            info['yref0'] = yref

    #  ------------------------------------------------------------------------

    def rotate_coo(self, x, y, phi):
        """
        Rotate an object's coordinates in the *.coo files for data sets
        containing images at different PAs.
        """

        """ Rotate around center of image, and keep origin at center """
        xin = self[0].data.shape[1] / 2
        yin = self[0].data.shape[0] / 2.
        xout = self[0].data.shape[1] / 2
        yout = self[0].data.shape[0] / 2.

        cos = math.cos(math.radians(phi))
        sin = math.sin(math.radians(phi))

        xrot = (x - xin) * cos - (y - yin) * sin + xout
        yrot = (x - xin) * sin + (y - yin) * cos + yout

        return [xrot, yrot]

    #  ------------------------------------------------------------------------

    def make_coo(self, refSrc, strSrc, check_loc=True, inpref='ce', outpref='c',
                 outdir=None, det_thresh=4., coo_update='aots', doplot=False,
                 rmax=40, plotfile='make_coo_check.png', verbose=True):

        """ Get reference AO stage position and RA, Dec from first image """
        hdr0 = self[0].header
        radecRef = [float(hdr0['RA']), float(hdr0['DEC'])]

        """
        Set up table columns to store info 
        (most are now set up now in the call to calc_offsets_aots)
        """
        self.datainfo['xref'] = 0.0
        self.datainfo['yref'] = 0.0

        """ Calculate the initial offsets """
        self.calc_offsets_aots(refSrc)

        """ Set up output files that have registration information in them """
        cfiles = self.make_outlist(inpref, outpref, outdir=outdir)

        """
        Loop through the files, refining the positions of the ref source and
        writing the information to the output files
        """
        if verbose:
            print('Making the *.coo and *.rcoo files')
            print('---------------------------------')
        for hdu, info, outfile in zip(self, self.datainfo, cfiles):
            hdr = hdu.header

            radec = [float(hdr['RA']), float(hdr['DEC'])]

            """ Determine the image's PA and plate scale """
            phi, inst_angle = self.get_pa_from_rotposn(self.instrument, hdr)
            scale = self.get_scale(self.instrument, hdr)

            """ Calculate the pixel offsets from the reference image """
            if coo_update == 'aots':
                pass  # Offsets already calculated from calc_offsets_aots
                # d_xy = kai_util.aotsxy2pix(aotsxy, scale, aotsxyRef,
                #                            inst_angle=inst_angle)
            # elif coo_update == 'radec':
            #     d_xy = kai_util.radec2pix(radec, phi, scale, radecRef)
            #     info['dx'] = d_xy[0]
            #     info['dy'] = d_xy[1]
            else:
                d_xy = [0, 0]
                info['dx'] = d_xy[0]
                info['dy'] = d_xy[1]

            """ In the new image, find the REF and STRL coords """
            xref = info['xref0']
            yref = info['yref0']
            xstr = strSrc[0] + info['dx']
            ystr = strSrc[1] + info['dy']

            basename = info['basename']
            if verbose:
                print('make_coo %s: xref, yref start = %6.2f %6.2f' %
                      (basename, info['xref0'], info['yref0']))

            """ Re-center stars to get exact coordinates """
            if check_loc:
                hdu.im_moments(xref, yref, rmax=rmax, detect_thresh=det_thresh)
                xref = hdu.imex_mux
                yref = hdu.imex_muy

                hdu.im_moments(xstr, ystr, rmax=rmax, detect_thresh=det_thresh)
                xstr = hdu.imex_mux
                ystr = hdu.imex_muy

                if verbose:
                    print('make_coo %s: xref, yref final = %6.2f %6.2f'
                          % (basename, xref, yref))
            info['xref'] = xref
            info['yref'] = yref

            """
            Calculate reference star location in rotated frame, in case it's
            needed later, i.e., if there are files with different PAs
            """
            xyRef_rot = self.rotate_coo(xref, yref, phi)
            xref_r = xyRef_rot[0]
            yref_r = xyRef_rot[1]

            xyStr_rot = self.rotate_coo(xstr, ystr, phi)
            xstr_r = xyStr_rot[0]
            ystr_r = xyStr_rot[1]

            """
            Write reference star x,y to fits header and save info in a new
            file
            """
            hdr.set('XREF', "%.3f" % xref, 'Cross Corr Reference Src x')
            hdr.set('YREF', "%.3f" % yref, 'Cross Corr Reference Src y')
            hdr.set('XREF_R', "%.3f" % xref_r,
                    'Reference Src x in rotated frame')
            hdr.set('YREF_R', "%.3f" % yref_r,
                    'Reference Src y in rotated frame')
            hdr.set('XSTREHL', "%.3f" % xstr, 'Strehl Reference Src x')
            hdr.set('YSTREHL', "%.3f" % ystr, 'Strehl Reference Src y')
            hdu.save(outfile)

            """ Put the xref, yref into into a *.coo file """
            outcoo = outfile.replace('.fits', '.coo')
            with open(outcoo, 'w') as f:
                f.write('%7.2f  %7.2f\n' % (xref, yref))

            """
            Put the rotated coordinates into a *.rcoo file
            """
            # outrcoo = outfile.replace('.fits', '.rcoo')
            # with open(outrcoo, 'w') as f:
            #     f.write('%7.2f  %7.2f\n' % (xref_r, yref_r))

        if doplot:
            self.datainfo['x_cent'] = self.datainfo['xref'].copy()
            self.datainfo['y_cent'] = self.datainfo['yref'].copy()
            self.datainfo['imsize'] = 2 * rmax + 1
            self.plot_multipanel(self.datainfo, mode='xy', markpos=True,
                                 outfile=plotfile)
