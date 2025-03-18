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
import sys
import glob

import matplotlib
matplotlib.use('Agg')

from astropy.io import ascii
from specim.imfuncs.wcshdu import WcsHDU

from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments

from . import kaiset

""" Turn off header deprecation warnings """
warnings.filterwarnings('ignore', category=UserWarning, append=True)

""" Define global variables for the two possible instruments """
osiris = instruments.OSIRIS()
nirc2 = instruments.NIRC2()


def get_instrument(instrument):
    """

    Makes sure that either NIRC2 or OSIRIS has been selected

    """
    if instrument.lower() == 'osiris':
        return osiris
    elif instrument.lower() == 'nirc2':
        return nirc2
    else:
        print('')
        raise ValueError('get_instrument: instrument must be '
                         '"osiris" or "nirc2"\n')


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
        if inst == osiris:
            keylist = ['assn', 'frames']
            assn = i['assn']
        elif inst == nirc2:
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
            if inst == osiris:
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

    """ Get the instrument """
    try:
        inst = get_instrument(instrument)
    except ValueError:
        return

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
                if inst == osiris:
                    frameroot = 'i%s_a' % obsdate[2:]
                    framelist = dict_to_framelist(tmplist, inst, frameroot,
                                                  suffix=suffix)
                elif inst == nirc2:
                    framelist = dict_to_framelist(tmplist, inst)
                else:
                    is_error = True
            elif inst == nirc2 and isinstance(el1, int):
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


def makelog_and_prep_images(year, instrument, rawdir='../raw'):
    """

    Make an electronic log from all the files in the ../raw/ directory.
    The file will be called kai.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda

    """

    """ Set up the instrument and make the log """
    try:
        inst = get_instrument(instrument)
    except ValueError:
        print('')
        print('ERROR in makelog_and_prep_images: '
              'invalid value for instrument parameter')
        print('')
        return

    print('Making instrument log from files')
    kai_util.makelog(rawdir, instrument=inst)

    """
    If you are reducing OSIRIS, you need to flip the images first.
    """
    if inst == osiris:
        print('Flipping data (needed for OSIRIS images)')
        raw_files = glob.glob('%s/i*.fits' % rawdir)
        osiris.flip_images(raw_files)

    """ Download weather data we will need. """
    print('Downloading weather data for year %s' % year)
    dar.get_atm_conditions(year)

    return


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

    # """ Create the framelist in the proper format """
    # darkframes = inlist_to_framelist(darklist, instrument, obsdate,
    #                                suffix=suffix)
    # print('Frames for dark: %s' % darklist['name'])
    # print(darkframes)

    """ Make the dark file """
    # calib.makedark(darkframes, outfile, instrument=inst)
    darkset = kaiset.KaiSet(darklist, instrument, obsdate, indir=rawdir,
                            wcsverb=False)
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
    flats_on = kaiset.KaiSet(flatlist, instrument, obsdate, indir=rawdir,
                             wcsverb=False)

    """ Make a lamps-off holder if requested """
    if 'offframes' not in flatlist.keys():
        flats_off = None
    else:
        if flats_on.instrument == osiris:
            tmpdict = {'assn': flatlist['assn'],
                       'frames': flatlist['offframes']}
        else:
            tmpdict = {'frames': flatlist['offframes']}
        flats_off = kaiset.KaiSet(tmpdict, instrument, obsdate, indir=rawdir)

    """ Create the framelist in the proper format """
    # onframes = inlist_to_framelist(flatlist, instrument, obsdate, suffix=suffix)
    # if 'offframes' not in flatlist.keys():
    #     offframes = range(0, 0)
    # else:
    #     if inst == osiris:
    #         tmpdict = {'assn': flatlist['assn'],
    #                    'frames': flatlist['offframes']}
    #     else:
    #         tmpdict = {'frames': flatlist['offframes']}
    #     offframes = inlist_to_framelist(tmpdict, instrument, obsdate,
    #                                     suffix=suffix)

    """ Make the flat-field file """
    outfile = '%s_%s.fits' % (flatlist['name'], flatlist['obsfilt'])
    print('')
    print('Creating the flat-field file: %s' % outfile)
    print('--------------------------------------------')
    # calib.makeflat(onframes, offframes, outfile, instrument=inst)
    flats_on.create_flat(outfile, lamps_off=flats_off, normalize='sigclip',
                         inflat=inflat)


def make_sky(skylist, obsdate, instrument, suffix=None):
    """

    Makes a dark frame given either an input list of integer frame numbers
    (for NIRC2) or a dict or list of dicts containing 'assn' and 'frames'
    keywords (for OSIRIS)

    """

    try:
        inst = get_instrument(instrument)
    except ValueError:
        return

    """ Create the framelist in the proper format """
    skyframes = inlist_to_framelist(skylist, instrument, obsdate,
                                    suffix=suffix)
    print(skyframes)

    """ Make the sky file """
    outfile = '%s.fits' % skylist['name']
    print('Creating the sky file: %s' % outfile)
    sky.makesky(skyframes, obsdate, skylist['obsfilt'], instrument=inst)


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

    """ Get the instrument """
    try:
        inst = get_instrument(instrument)
    except ValueError:
        return

    """ Set up the base keys that should be in all of the input dicts """
    basekeys = ['name', 'frames']
    if inst == osiris:
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
                finalflat.writeto('flat_%s.fits' % info['obsfilt'])

        """ Create the sky """
        for info in skylist:
            make_sky(info, obsdate, instrument, suffix=suffix)

    """
    Make the bad pixel mask, which KAI calls the 'supermask' from a dark and
     a flat.
    Use a long-exposure dark (>20 sec) for this.
    Note: the code assumes that the input files are found in calib/darks
     and calib/flats
    """
    print('')
    print('Making supermask.fits')
    print('---------------------')
    calib.makemask(dark4mask, flat4mask, 'supermask.fits', instrument=inst)


def name_checker(a, b):
    length = len(a) + len(b)
    if length > 12+8:
        if length != 25:
            print("Check your ABCs, length is " + str(length))
            print("(Length should be 25 or below 21)")
            sys.exit()


def reduce(target, obsdate, inlist, obsfilt, refSrc, instrument, suffix=None,
           skyscale=False, usestrehl=False, dockerun=False):
    """
    Do the full data reduction.

    Inputs:
      target     - root name of the target object
      obsdate    - 8-digit observation date in yyyymmdd format (e.g., 20201101)
      inlist     -
      refSrc     - make sure you use the position in the _flipped_ image.
      instrument -
      suffix     - any suffix that should be added to the frame names that
                    are generated from the assnlist, e.g., 'flip'.
                   The default (None) means do not add a suffix
      skyscale   -
      usestrehl  -
      dockerrun  - is this function being called within a Docker run in which
                   the weather files have not yet been downloaded.
                   Set to True to download the weather files.
                   Default is False
    """	

    #    -- If you have more than one position angle, make sure to
    #       clean them seperatly.
    #    -- Strehl and Ref src should be the pixel coordinates of a bright
    #       (but non saturated) source in the first exposure of sci_files.
    #    -- If you use the OSIRIS image, you must include the full filename
    #    in the list.

    """ Get the input framelist """
    sci_frames = inlist_to_framelist(inlist, instrument, obsdate, suffix=suffix)

    """ 
    Download weather data that will be used in later steps.
    NOTE: This step is only needed if running the code within docker, and not 
     even then if the "makelog_and_prep_images" function has been run within 
     the same docker session that is being used to run this "go" function  
    """
    obsyear = obsdate[:4]
    if dockerun:
        print('Downloading weather data for year %s' % obsyear)
        dar.get_atm_conditions(obsyear)

    """ Make the list of science frames from the input assn list"""
    #  name_checker(epoch,target) - don't need this anymore
    print('')
    print('Science frames to be cleaned and combined')
    print('-----------------------------------------')
    for frame in sci_frames:
        print('%s' % frame)
    print('')
    print(os.getcwd())
    print('')

    """ Check the instrument """
    try:
        inst = get_instrument(instrument)
    except ValueError:
        print('')
        print('ERROR: Invalid choice of instrument parameter')
        print('')
        return

    """ For this target, use the sky created for 2022_06_05 """
    # sky.makesky(sky_frames, obsdate, obsfilt, instrument=osiris)
    print('Calibrating and cleaning the input files')
    print('----------------------------------------')
    data.clean(sci_frames, obsdate, obsfilt, refSrc, refSrc, field=target,
               instrument=inst, skyscale=skyscale)
    if usestrehl:
        data.calcStrehl(sci_frames, obsfilt, field=target, instrument=inst)
        combwht = 'strehl'
        submaps = 0  # UCB group sets this to 3 for their images (have stars)
    else:
        combwht = None
        submaps = 0
    print('')
    print('Combining the calibrated files')
    print('------------------------------')
    data.combine(sci_frames, obsfilt, obsdate, field=target,
                 trim=0, weight=combwht, submaps=submaps, instrument=inst)


def finalize(target, obsdate, inlist, obsfilt, refradec, instrument,
             combdir='default', suffix=None):
    """
    Get the combined drizzled image into its final format.
    Inputs:
      target     - root name of the target object
      obsdate    - 8-digit observation date in yyyymmdd format (e.g., 20201101)
      inlist     -
      refSrc     - make sure you use the position in the _flipped_ image.
      instrument -
      suffix     - any suffix that should be added to the frame names that
                    are generated from the assnlist, e.g., 'flip'.
                   The default (None) means do not add a suffix

    """

    """ Start by setting up the relevant filenames """
    if combdir is None:
        combdir = '.'
    elif combdir == 'default':
        combdir = '%s/combo' % os.getcwd()
    combroot = os.path.join(combdir, 'mag%s_%s_%s' % (obsdate, target, obsfilt))
    scifile = '%s.fits' % combroot
    sigfile = '%s_sig.fits' % combroot
    outroot = os.path.join(combdir, '%s_kai_%s_%s' % (target, obsdate, obsfilt))
    outsci = '%s.fits' % outroot
    outwht = '%s_wht.fits' % outroot

    """ Make the list of science frames from the input file information """
    sci_frames = inlist_to_framelist(inlist, instrument, obsdate, suffix=suffix)

    """ Read in the science and "sig" files """
    sciin = WcsHDU(scifile)
    whtin = WcsHDU(sigfile)

    """ Set up the instrument, for getting plate scale, etc. """
    try:
        inst = get_instrument(instrument)
    except ValueError:
        print('')
        print('ERROR: Invalid choice of instrument parameter')
        print('')
        return

    """
    Get some information from the drizzle output and put it into SHARP
    standardized form
    """
    hdr = sciin.header
    if 'ITIME' in hdr.keys():
        if inst == osiris:
            tottime = hdr['itime'] / 1000.
        else:
            tottime = hdr['itime']
        hdr['elaptime'] = tottime
        hdr['exptime'] = tottime
        avgtime = tottime / len(sci_frames)
    else:
        avgtime = 1.
    if 'FILTER' not in hdr.keys():
        hdr['filter'] = inst.get_filter_name(hdr)
    if 'NDRIZIM' in hdr.keys():
        hdr['ncombine'] = hdr['ndrizim']
    for i, f in enumerate(sci_frames):
        hdr.set('orig%03d' % (i + 1), f, 'Original input file %d' % (i + 1))

    """
    Use the WcsHDR functionality to:
      1. Set the pixel scale to the appropriate value for the instrument an mode
      2. Convert the pixel units to e-/sec
    """
    pixscale = inst.get_plate_scale(hdr)
    gain = inst.get_gain(hdr)
    if avgtime > 1.:
        texp = avgtime
    else:
        texp = -1.
    newsci = sciin.process_data(gain=gain, texp=texp, pixscale=pixscale)
    newwht = whtin.process_data(pixscale=pixscale)

    """ Put the correct WCS info into the file, if requested """
    if refradec is not None:
        """
        In the final combined science image, assign the refradec values to the
        correct ref pixel in the drizzed image.  The appropriate pixel location
        is stored in the mag[obsdate]_[target]_[obsfilt].coo file
        """
        coofile = '%s.coo' % combroot

        """ Read the reference pixel information from the *.coo file """
        print('')
        print('Reading reference pixel information from %s' % coofile)
        try:
            cootab = ascii.read(coofile, names=['x', 'y'])
            refcoo = [cootab['x'][0], cootab['y'][0]]
            print('Reference pixel: %.2f %.2f' % (refcoo[0], refcoo[1]))
        except IOError:
            print('')
            print('Could not find %s' % coofile)
            print('')
            print('Current directory is: %s' % os.getcwd())
            print('')
            raise IOError

        """ Set the WCS values in the science image """
        newsci.crpix = refcoo
        newsci.crval = refradec
        newwht.crpix = refcoo
        newwht.crval = refradec

    """ Save the updated file with a new name """
    keeplist = ['object', 'telescop', 'instrume', 'filter', 'date-obs',
                'bunit', 'binfo_1', 'binfo_2', 'gain', 'ogain', 'exptime',
                'semester', 'progpi', 'progid', 'elaptime', 'ncombine']
    for i in range(len(sci_frames)):
        keeplist += ['orig%03d' % (i + 1)]
    newsci.writeto(outfile=outsci, keeplist=keeplist)

    """
    Convert the drizzle "sig" file to a weight file in the style of
    astrodrizzle, where the weight in each pixel is the effective exposure time
    in that pixel.
    """
    newwht.data *= avgtime
    """
    Also fix the DATASEC header key for the wht image so that ds9 will display
    it properly.
    """
    whdr = newwht.header
    whdr['datasec'] = '[1:%d,1:%d]' % (hdr['naxis1'], hdr['naxis2'])

    newwht.writeto(outfile=outwht)
