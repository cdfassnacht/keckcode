import os
import sys
from ccdredux.ccdset import CCDSet
from kai import instruments
from keckcode.ao_img.aoset import AOSet

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

