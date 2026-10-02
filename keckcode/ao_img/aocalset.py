"""

Class AOCalSet

This class takes information about calibration exposures and uses that
information to produce the calibration files that eventually will be applied
to the science frames

"""

import os
import sys
import math
import numpy as np

from specim.imfuncs.wcshdu import WcsHDU
from ccdredux.ccdset import CCDSet
from .aoset import AOSet

pyversion = sys.version_info.major

class AOCalSet(dict):
    """

    In order to create an AOCalSet instance you need to pass a dict of dicts.

    """

    def __init__(self, calinfo, **kwargs):
        """

        calinfo must be a dict that contains information about the calibration
        exposures and other information about the observing run

        """

        """ Check input format """
        if not isinstance(calinfo, dict):
            print('')
            print('Input msut be a dict')
            print('')
            raise TypeError('Input msut be a dict')

        """ Make sure that the absolutely required keys are present """
        reqkeys = ['obsdate', 'instrument']
        for k in reqkeys:
            if k not in calinfo.keys():
                print('')
                raise KeyError('Missing required key %s' % k)

        """ Set up default values for the main parameters """
        self.set_default_params()

        """ Load the parameters from the passed dictionary """
        for k in calinfo.keys():
            self[k] = calinfo[k]

        """ Set base keys for darkinfo, flatinfo, etc. """
        self['basekeys'] = ['name', 'frames']
        if calinfo['instrument'] == 'osiris' or caldata['instrument'] == 'osim':
            self['basekeys'].append('assn')

    # ------------------------------------------------------------------------

    def set_default_params(self):
        """
        Set default values for the main parameters
        """

        self['rawdir'] = None
        self['caldir'] = None
        self['suffix'] = None
        self['darkinfo'] = None
        self['flatinfo'] = None
        self['skyinfo'] = None
        self['dark4mask'] = None
        self['flat4mask'] = None
        self['dark4flat'] = None
        self['flat4sky'] = None
        self['root4sky'] = None
        self['bpmsig'] = 5.
        self['forkai'] = True

    # ------------------------------------------------------------------------

    def check_callist(self, listname, dictkeys):
        """

        Checks the type and dictionary keys of the lists of calibratioin files

        """

        """ Make sure a valid listname has been passed """
        try:
            callist = self[listname]
        except KeyError:
            print('')
            raise KeyError('check_callist got unexpect listname: %s' % listname)

        """ Check callist type and modify if necessary """
        if isinstance(callist, dict):
            newlist = [callist]
        elif isinstance(callist, (list, tuple, np.ndarray)):
            newlist = list(callist)
        else:
            raise TypeError('\nCalibration list must be one of the following:\n'
                            'dict, list, numpy array, or tuple\n\n')

        """
        Make sure each element is a dict and that the dict contains the expected
        keys
        """

        for assndict in newlist:
            if isinstance(assndict, dict):
                for k in dictkeys:
                    if k not in assndict.keys():
                        raise KeyError('\nCalibration list is missing expected '
                                       '%s key.\n\n' % k)
            else:
                raise TypeError('\nCalibration list must contain dict objects'
                                '\n\n')

        return newlist

    # ------------------------------------------------------------------------

    def make_dark(self, **kwargs):
        """

        Makes a dark frame given either an input list of integer frame numbers
        (for NIRC2) or a dict or list of dicts containing 'assn' and 'frames'
        keywords (for OSIRIS)

        """

        """ Check the darkinfo format """
        dkeys = self['basekeys']
        darklist = self.check_callist('darkinfo', dkeys)

        """ Create the dark(s) """
        for info in darklist:
            print('')
            """ Get the output file name """
            if info['name'][-4:] != 'fits':
                outfile = '%s.fits' % info['name']
            else:
                outfile = info['name']
            obj = outfile[:-5]
            print('Creating the dark file: %s' % outfile)

            """ Make the dark file """
            darkset = AOSet(info, self['instrument'], self['obsdate'],
                            indir=self['rawdir'], is_sci=False, wcsverb=False)
            darkset.create_dark(outfile, caldir=self['caldir'], outobj=obj)
            print('')

        print('===========================================================')
        print('   Finished creating dark frames')
        print('===========================================================')

    def make_flat(self, bpm=None, **kwargs):
        """
        Makes a flat frame given either an input list of integer frame numbers
        (for NIRC2) or a dict or list of dicts containing 'assn' and 'frames'
        keywords (for OSIRIS)

        """

        """ Check the flatinfo format """
        fkeys = self['basekeys']
        fkeys.append('obsfilt')
        flatlist = self.check_callist('flatinfo', fkeys)

        """ Set default value """
        normalize = 'sigclip'

        """ Create the flat(s) """
        allflats1 = []
        for info in flatlist:
            print('')

            """ Make an AOSet holder for the lamps-on frames """
            print('Reading flat-field frames (lamps on)')
            flats_on = AOSet(info, self['instrument'], self['obsdate'],
                             indir=self['rawdir'], is_sci=False, wcsverb=False)

            """ Make a lamps-off holder if requested """
            if 'offframes' not in info.keys():
                flats_off = None
            else:
                if self['instrument'] == 'osiris':
                    tmpdict = {'assn': info['assn'],
                               'frames': info['offframes']}
                else:
                    tmpdict = {'frames': info['offframes']}
                print('Reading flat-field frames (lamps off)')
                flats_off = AOSet(tmpdict, self['instrument'], self['obsdate'],
                                  indir=self['rawdir'], is_sci=False,
                                  wcsverb=False)

            """ Make object masks if requested """

            """ Make the flat-field file """
            outfile = '%s_%s.fits' % (info['name'], info['obsfilt'])
            flats_on.create_flat(outfile, lamps_off=flats_off,
                                 normalize=normalize, indark=self['dark4flat'],
                                 bpm=bpm, caldir=self['caldir'], **kwargs)

            allflats1.append('%s.fits' % info['name'])

        del fkeys

        print('===========================================================')
        print('   Finished creating (initial) flat-field frames')
        print('===========================================================')


    def make_cals(self):
        """

        Runs the various methods to create the calibration files, for example,
        make_dark, make_flat, etc.

        """

        if self['darkinfo'] is not None:
            self.make_dark()

        bpm0 = None
        if self['flatinfo'] is not None:
            self.make_flat(bpm=bpm0)
# def make_calfiles(obsdate, darkinfo, flatinfo, skyinfo, dark4mask, flat4mask,
#                   instrument, skyflatinfo=None, rawdir=None, caldir=None,
#                   dark4flat=None, dark4sky=None, flat4sky=None,
#                   root4sky=None, bpmsig=5., suffix=None, forkai=True, **kwargs):

