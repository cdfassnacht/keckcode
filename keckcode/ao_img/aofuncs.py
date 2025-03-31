""" Import basic modules

The matplotlib.use('Agg') line prevents any display functionality by
matplotlib, which is necessary when running the code in a Docker container.
It must come before any import statement that will also lead to an
import of matplotlib, which it needs to be before the specim and kai imports
"""
import warnings

import numpy
import numpy as np
import os

import matplotlib
matplotlib.use('Agg')

from specim.imfuncs.wcshdu import WcsHDU

from .aoset import AOSet

""" Turn off header deprecation warnings """
warnings.filterwarnings('ignore', category=UserWarning, append=True)


def dict_to_framelist(inlist, inst, rootname=None, suffix=None):
    """

    Utility function to make a list of frame numbers from a list of
    dictionaries, where the dictionaries must have two entries:
      'assn'   - the association number
      'frames' - the list of good frames in the association
    For example, if the following frames are good: i220605_a007002,
      i220605_a007003, and i220605_a007004, then the dictionary
      for this association would be {'assn': 7, 'frames': np.arange(2,5)}

    Inputs:
     assnlist - A list of dictionaries, one for each association
     rootname - The base name that is in common for all of the frame names.
                In the example above, the rootname would be "i220605_a"

    Output:
     framelist - a single list of all the good frames from the input
                 list of dictionaries.

    """

    """ First check the assnlist data type """
    iserror = False
    dictlist = []
    if inlist is not None:
        if isinstance(inlist, dict):
            dictlist = [inlist]
        elif isinstance(inlist, list):
            if isinstance(inlist[0], dict):
                dictlist = inlist
            else:
                print('ERROR: first element of assnlist is not a dict')
                iserror = True
        else:
            print('ERROR: inlist is not either a dict or a list')
            iserror = True
    else:
        print('ERROR: inlist is None')
        iserror = True
    if iserror:
        raise TypeError('\ndict_to_framelist: inlist must be a dict or a '
                        'list of dicts.\n')

    """ Create the framelist """
    framelist = []
    for i in dictlist:
        """ Check the passed parameters """
        if inst == 'osiris' or inst == 'osim':
            keylist = ['assn', 'frames']
            assn = i['assn']
        elif inst == 'nirc2':
            keylist = ['frames']
            assn = None
        else:
            raise ValueError('Instrument must be osiris or nirc2')
        for k in keylist:
            if k not in i.keys():
                raise KeyError('\ndict_to framelist: dict is missing the '
                               '%s key\n\n' % k)
        """
        Loop through association and frame numbers to create file list
        """
        for j in i['frames']:
            if inst == 'osiris' or inst == 'osim':
                framename = '%s%03d%03d' % (rootname, assn, j)
                if suffix is not None:
                    framename = '%s_%s' % (framename, suffix)
            else:
                framename = j
            framelist.append(framename)

    return framelist


def inlist_to_framelist(inlist, instrument, obsdate, suffix=None):
    """

    Takes an input list that is either a list of integers (for NIRC2) or
    a list of associations (a list of dicts, for OSIRIS) and converts it
    to the proper framelist format needed for the KAI functions.

    Inputs:
       inlist - A list of either integer frame numbers (for NIRC2) or dict
                objects, with a 'assn' and 'frames' keywords (for ORIRIS)
       instrument - either 'nirc2' or 'osiris'

    Output:
       framelist - a list of frames, in the proper format for KAI

    """

    """ Convert the input list into a frame list"""
    is_error = False
    tmplist = None
    framelist = None
    if isinstance(inlist, (dict, int)):
        tmplist = [inlist]
    # NOTE: in python 3, "range" is a type, but in python 2, which this
    #       code runs in, the range function produces a list output
    # elif isinstance(inlist, range):
    #    framelist = inlist
    elif isinstance(inlist, (tuple, numpy.ndarray, list)):
        tmplist = inlist
    else:
        is_error = True

    if is_error is not True:
        if len(tmplist) == 0:
            framelist = range(0, 0)
        else:
            el1 = tmplist[0]
            if isinstance(el1, dict):
                if instrument == 'osiris' or instrument == 'osim':
                    frameroot = 'i%s_a' % obsdate[2:]
                    framelist = dict_to_framelist(tmplist, instrument,
                                                  frameroot, suffix=suffix)
                elif instrument == 'nirc2':
                    framelist = dict_to_framelist(tmplist, instrument)
                else:
                    is_error = True
            elif instrument == 'nirc2' and isinstance(el1, int):
                framelist = tmplist
            else:
                is_error = True
    if is_error:
        print('')
        print('ERROR: Input frame(s) [inlist parameter] must be one of the '
              'following')
        print('   1. An int (NIRC2)')
        print('   2. A dict with "assn" and "frames" keywords (OSIRIS)')
        print('   3. A list of either ints (NIRC2) or dicts (OSIRIS)')
        print('   4. A tuple or numpy array containing integers (NIRC2)')
        print('   5. A range, i.e., range(151,157) - (NIRC2)')
        print('')
        raise TypeError()
    else:
        return framelist


def check_callist(callist, dictkeys):
    """

    Checks the type and dictionary keys of the lists of calibratioin files

    """

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


def make_dark(darklist, obsdate, instrument, rawdir='../raw', suffix=None):
    """

    Makes a dark frame given either an input list of integer frame numbers
    (for NIRC2) or a dict or list of dicts containing 'assn' and 'frames'
    keywords (for OSIRIS)

    """

    """ Get the output file name """
    outfile = '%s.fits' % darklist['name']
    print('Creating the dark file: %s' % outfile)

    """ Make the dark file """
    darkset = AOSet(darklist, instrument, obsdate, indir=rawdir, wcsverb=False)
    darkset.create_dark(outfile)


def make_flat(flatlist, obsdate, instrument, rawdir='../raw', inflat=None,
              suffix=None):
    """

    Makes a flat-field file

    Inputs:
     onlist     - A list of integers (for NIRC2) or dicts (for OSIRIS)
                   designating the frames associated with the exposures taken
                   when the dome lamp was on.  If twilight or sky exposures are
                   being used instead, then they should be designated here.
     offlist    - The list that designates the frames for which the dome lamp
                   was turned off.  If no dome flats were taken and, thus, the
                   onlist parameter contains twilight flats or sky frames, then
                   just set this offlist parameter to None or to range(0, 0)
     outfile    - Output file name
     instrument - Either 'nirc2' or 'osiris'

    """

    """ Make a KaiSet holder for the lamps-on frames """
    flats_on = AOSet(flatlist, instrument, obsdate, indir=rawdir, wcsverb=False)

    """ Make a lamps-off holder if requested """
    if 'offframes' not in flatlist.keys():
        flats_off = None
    else:
        if instrument == 'osiris':
            tmpdict = {'assn': flatlist['assn'],
                       'frames': flatlist['offframes']}
        else:
            tmpdict = {'frames': flatlist['offframes']}
        flats_off = AOSet(tmpdict, instrument, obsdate, indir=rawdir)

    """ Make the flat-field file """
    outfile = '%s_%s.fits' % (flatlist['name'], flatlist['obsfilt'])
    print('')
    print('Creating the flat-field file: %s' % outfile)
    print('--------------------------------------------')
    flats_on.create_flat(outfile, lamps_off=flats_off, normalize='sigclip',
                         inflat=inflat)


def make_calfiles(obsdate, darkinfo, flatinfo, skyinfo, dark4mask, flat4mask,
                  instrument, root4sky=None, suffix=None):
    """
    
    Makes all of the calibration files
    
    Inputs:
     darkinfo   - A single dict or a list of dicts, where each one contains the
                  information needed to make one dark file.  A list of dicts
                  is needed if dark exposures with several different exposure
                  times (e.g., 7 sec and 180 sec) were taken.  If all of the
                  dark exposures had the same exposure time, then darkinfo
                  will just be a single dict rather than a list of dicts.
                  The keys for each dict are:
                   inlist     - a list of integer frame numbers (for NIRC2) or
                                dicts (for OSIRIS) defining the exposures needed
                                for this given dark.  For the OSIRIS data, the
                                dicts should contain 'assn' and 'frames' 
                                keywords
                   outfile    - output filename 
     flatinfo   - A single dict or a list of dicts, where each one contains the
                  information needed to make one flat file.  A list of dicts
                  is needed if flat-field exposures in several different filters
                  (e.g., Kp and H) were taken.  If flats were taken through only
                  one filter, then flatinfo will just be a single dict rather
                  than a list of dicts.
                  The keys for each dict are:
                   band       -
                   onlist     -
                   offlist    - [OPTIONAL]
                   outfile    -
     dark4mask  -
     flat4mask  -
     instrument - either 'nirc2' or 'osiris'

    """

    """ Set up the base keys that should be in all of the input dicts """
    basekeys = ['name', 'frames']
    if instrument == 'osiris' or instrument == 'osim':
        basekeys.append('assn')

    """ Create the dark(s) as long as darkinfo is not None"""
    if darkinfo is not None:
        """ Check the darkinfo format """
        dkeys = list(basekeys)
        darklist = check_callist(darkinfo, dkeys)

        """ Create the dark(s) """
        for info in darklist:
            make_dark(info, obsdate, instrument, suffix=suffix)
        del dkeys

    """ Create the flat(s) as long as flatinfo is not None"""
    allflats1 = []
    if flatinfo is not None:
        """ Check the flatinfo format """
        fkeys = list(basekeys)
        fkeys.append('obsfilt')
        flatlist = check_callist(flatinfo, fkeys)

        """ Create the flat(s) """
        for info in flatlist:
            make_flat(info, obsdate, instrument, suffix=suffix)
            allflats1.append('%s.fits' % info['name'])

    """
    Make a sky frame
    """
    allflats2 = []
    if skyinfo is not None:
        """ Check the skyinfo format """
        skeys = list(basekeys)
        skeys.append('obsfilt')
        skeys.append('type')
        skylist = check_callist(skyinfo, skeys)

        """ Create the sky flat file(s) and then the final combined flat """
        for info in skylist:
            if root4sky is not None:
                inflatfile = '%s_%s.fits' % (root4sky, info['obsfilt'])
                inflat = os.path.join('calib', 'flats', inflatfile)
            else:
                inflat = None
            make_flat(info, obsdate, instrument, inflat=inflat, suffix=suffix)
            allflats2.append('%s.fits' % info['name'])

            """ Create the final flat"""
            if root4sky is not None:
                flat1 = WcsHDU(inflat, wcsverb=False)
                flat2file = '%s_%s.fits' % (info['name'], info['obsfilt'])
                flat2path = os.path.join('calib', 'flats', flat2file)
                flat2 = WcsHDU(flat2path, wcsverb=False)
                finalflat = flat1 * flat2
                outfile = os.path.join('calib', 'flats', 'flat_%s.fits' %
                                       info['obsfilt'])
                finalflat.writeto(outfile)

        """ Create the sky """
        # for info in skylist:
        #     make_sky(info, obsdate, instrument, suffix=suffix)

    """
    Make the bad pixel mask, which KAI calls the 'supermask' from a dark and
     a flat.
    Use a long-exposure dark (>20 sec) for this.
    Note: the code assumes that the input files are found in calib/darks
     and calib/flats
    """
    # print('')
    # print('Making a bad pixel mask (what KAI calls a supermask)')
    # print('---------------------')
    # calib.makemask(dark4mask, flat4mask, 'supermask.fits', instrument=inst)
