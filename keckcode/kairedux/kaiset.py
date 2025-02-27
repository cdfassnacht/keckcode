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
        if outname[:-4] != 'fits':
            outfile = os.path.join(darkdir, '%s.fits' % outname)
        else:
            outfile = os.path.join(darkdir, outname)
        outlist = os.path.join(darkdir, 'dark.lis')

        """ """
        self.make_bias(outfile=outfile, reject=reject, nlow=nlow, nhigh=nhigh)
