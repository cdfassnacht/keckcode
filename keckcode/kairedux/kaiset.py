import os
import sys
from ccdredux.ccdset import CCDSet
from kai import instruments

""" Define global variables for the two possible instruments """
osiris = instruments.OSIRIS()
nirc2 = instruments.NIRC2()

pyversion = sys.version_info.major


class KaiSet(CCDSet):
    """

    Class to run KAI functions on related sets of files

    """

    def __init__(self, inlist, inst, obsdate, indir=None, gzip=False,
                 verbose=True, **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nKaiSet: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Get the instrument in KAI format """
        self.instrument = None
        try:
            self.get_instrument(inst)
        except ValueError:
            print('')
            print('Could not create kaiset object')
            print('')
            return

        """ Set instrument-specific parameters."""
        if self.instrument == osiris:
            texpkey = 'truitime'
            gainkey = 'sysgain'
        else:
            texpkey = 'elaptime'
            gainkey = 'gain'

        """ Make the date string from the provided obsdate """
        if obsdate is not None:
            self.obsdate = obsdate
            datestr = '%s_%s_%s' % (obsdate[:4], obsdate[4:6], obsdate[6:8])
        else:
            raise TypeError('\n an obsdate that is not None '
                            'must be provided\n\n')

        """ Create the input filelist from the passed parameters """
        if isinstance(inlist, dict):
            filelist = self.make_filelist([inlist], obsdate, gzip=gzip)
        elif isinstance(inlist[0], dict):
            filelist = self.make_filelist(inlist, obsdate, gzip=gzip)
        else:
            """
            For any other data types, let CCDSet (called through the "super"
            calls below) do the type checking.
            """
            filelist = inlist

        """ Set up the SHARP-specific input directory if requested """
        if indir is None:
            indir = '.'
        elif indir == 'auto':
            indir = os.path.join(os.getenv('sharpdat'), 'Raw', datestr)
        else:
            pass
        if verbose:
            print('Reading files from %s' % indir)

        """ Set up the KaiSet container by calling the superclass """
        if pyversion == 2:
            super(KaiSet, self).__init__(filelist, texpkey=texpkey,
                                         gainkey=gainkey, indir=indir, **kwargs)
        else:
            super().__init__(filelist, texpkey=texpkey, gainkey=gainkey,
                             indir=indir, **kwargs)

    #  ------------------------------------------------------------------------

    def get_instrument(self, instrument):
        """

        Makes sure that either NIRC2 or OSIRIS has been selected

        """
        if instrument.lower() == 'osiris' or instrument.lower() == 'osim':
            self.instrument = osiris
        elif instrument.lower() == 'nirc2':
            self.instrument = nirc2
        else:
            print('')
            raise ValueError('get_instrument: instrument must be '
                             '"osiris" or "nirc2"\n')

    #  ------------------------------------------------------------------------

    @staticmethod
    def make_filelist_osim(assnlist, obsdate, suff='fits'):
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
                filebase = 'i%s_a%03d%03d' % (obsdate[2:], assn, j)
                filelist.append('%s.%s' % (filebase, suff))

        return filelist

    #  ------------------------------------------------------------------------

    @staticmethod
    def make_filelist_nirc2(inlist, suff='fits'):

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

    def make_filelist(self, inlist, obsdate, gzip=False):
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
        if self.instrument == osiris:
            filelist = self.make_filelist_osim(inlist, obsdate, suff)
        else:
            filelist = self.make_filelist_nirc2(inlist, suff)

        return filelist

    #  ------------------------------------------------------------------------

    def create_dark(self, outname, reject='sigclip', nlow=1, nhigh=1):
        """

        Creates a dark-frame using the ccdredux make_dark functionality

        """

        """ Set up directory names """
        pwd = os.getcwd()
        caldir = os.path.join(pwd, 'calib')
        darkdir = os.path.join(caldir, 'darks')
        if not os.path.isdir(caldir):
            os.makedirs(caldir)
        if not os.path.isdir(darkdir):
            os.makedirs(darkdir)

        """ Get the output filenames """
        if outname[-4:] != 'fits':
            outfile = os.path.join(darkdir, '%s.fits' % outname)
        else:
            outfile = os.path.join(darkdir, outname)
        outlist = os.path.join(darkdir, 'dark.lis')

        """ """
        self.make_bias(outfile=outfile, reject=reject, nlow=nlow, nhigh=nhigh)

    #  ------------------------------------------------------------------------

    def create_flat(self, outname, lamps_off=None, normalize=None,
                    inflat=None, reject='sigclip', nlow=1, nhigh=1):
        """

        Creates a flat-field frame following the KAI recipe.  This approach
        uses frames taken with the flat-field lamps on (which are part of this
        KaiSet instance, and possibly frames taken with the flat-field lamps
        off, which are provided via the lamps_off parameter.  The lamps_off
        data should also be a KaiSet instance if it is not None.

        """

        """ Check validity of input """
        if lamps_off is not None:
            if not isinstance(lamps_off, KaiSet):
                raise ValueError('\nThe lamps_off parameter must be a KaiSet\n')
            if len(lamps_off) != len(self):
                raise IndexError('\nDifferent number of files in lamps_off and'
                                 ' lamps_on\n')

        """ Set up directory names """
        pwd = os.getcwd()
        caldir = os.path.join(pwd, 'calib')
        flatdir = os.path.join(caldir, 'flats')
        if not os.path.isdir(caldir):
            os.makedirs(caldir)
        if not os.path.isdir(flatdir):
            os.makedirs(flatdir)

        """ Get the output filenames """
        if outname[-4:] != 'fits':
            outfile = os.path.join(flatdir, '%s.fits' % outname)
        else:
            outfile = os.path.join(flatdir, outname)
        onfits = os.path.join(flatdir, 'lampsOn.fits')
        offfits = os.path.join(flatdir, 'lampsOff.fits')
        normfits = os.path.join(flatdir, 'flatNotNorm.fits')
        onlist = os.path.join(flatdir, 'on.lis')
        offlist = os.path.join(flatdir, 'off.lis')
        normlist = os.path.join(flatdir, 'onNorm.lis')

        """
        Make the flat, using the lamps-off frames if provided, but otherwise
        using just the lamps-on frames
        """
        if lamps_off is not None:
            """ First make the combined lamps-off and lamps-on frames """
            lamps_off.make_flat(outfile=offfits, normalize=normalize,
                                flatfile=inflat,
                                reject=reject, nlow=nlow, nhigh=nhigh)
            self.make_flat(outfile=onfits, normalize=normalize, flatfile=inflat,
                           reject=reject, nlow=nlow, nhigh=nhigh)

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
                           nlow=nlow, nhigh=nhigh)

        else:
            """ Put the filenames into the onlist """
            if 'basename' in self.datainfo.keys():
                nfile = open(onlist, 'w')
                for info in self.datainfo:
                    nfile.write('%s\n' % info['basename'])

            """ Make the final flat """
            self.make_flat(outfile=outfile, normalize=normalize,
                           flatfile=inflat, reject=reject, nlow=nlow,
                           nhigh=nhigh)

    #  ------------------------------------------------------------------------

    def create_sky(self, outname, obsfilt, skyscale=True,
                   inflat=None, reject='sigclip', nlow=1, nhigh=1):
        """

        Makes a sky file, broadly following the KAI algorithm, but with the
         possibility of applying the flat-field file to the sky frames
         before combining them.
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
