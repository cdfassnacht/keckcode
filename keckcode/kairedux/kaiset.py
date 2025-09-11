import sys
import os
from kai import instruments
from kai.reduce import data
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

        """ Get the instrument in KAI format """
        self.instrument = None
        try:
            self.get_instrument(inst)
        except ValueError:
            print('')
            print('Could not create kaiset object')
            print('')
            return

        """ Set up the KaiSet container by calling the superclass """
        if pyversion == 2:
            super(KaiSet, self).__init__(inlist, inst, obsdate, indir=indir,
                                         gzip=gzip, verbose=verbose, **kwargs)
        else:
            super().__init__(inlist, inst, obsdate, indir=indir, gzip=gzip,
                             verbose=verbose, **kwargs)

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

    def add_def_hdrinfo(self):
        """

        Adds keywords that are needed for later procrssing to the headers of
         the images

        """

        """ Loop through the images """
        for hdu in self:

            """ Get the central wavelength of the filter being used """
            hdr = hdu.header
            defwave = self.instrument.get_central_wavelength(hdr)

            """ Add the new wavelength header cards """
            hdr['effwave'] = defwave
            hdr['cenwave'] = defwave
            hdr['camname'] = 'narrow'

            try:
                nonlinSky = hdr['skylev']
            except KeyError:
                nonlinSky = 0.
            coaddkey = self.instrument.hdr_keys['coadds']
            try:
                coadds = hdr[coaddkey]
            except:
                coadds = 1
            satLevel = (coadds * instrument.get_saturation_level()) - nonlinSky
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
