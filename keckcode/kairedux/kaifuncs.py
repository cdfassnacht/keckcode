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
from datetime import datetime

import matplotlib
matplotlib.use('Agg')

from astropy.io import ascii
from astropy.io import fits
from pyraf import iraf as ir
from specim.imfuncs.wcshdu import WcsHDU

from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import bfixpix
from kai.reduce import data
from kai.reduce import dar
from kai.reduce import kai_util
from kai.reduce import util
from kai import instruments

from ..ao_img.aoset import AOSet
from ..ao_img import aofuncs as aofn

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
        raw_files.sort()
        print(raw_files)
        osiris.flip_images(raw_files)

    """ Download weather data we will need. """
    print('Downloading weather data for year %s' % year)
    dar.get_atm_conditions(year)

    return


def check_callist(callist, dictkeys):
    """

    Checks the type and dictionary keys of the lists of calibration files

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
                  instrument, dark4flat=None, root4sky=None, suffix=None):
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

    """ Make darks and flats via a call to ao_funcs """
    aofn.make_calfiles(obsdate, darkinfo, flatinfo, skyinfo, dark4mask,
                       flat4mask, instrument, dark4flat=dark4flat,
                       root4sky=root4sky, suffix=suffix)

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

        """ NOTE: the skyflat is created in the call to ao_funcs above """

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

def apply_calib():
    print('foo')


def kaiclean2(files, nite, wave, refSrc, strSrc, badColumns=None, field=None,
              angOff=0.0, cent_box=12,
              fixDAR=True, use_koa_weather=False, clean_dir=None,
              instrument=instruments.default_inst, check_ref_loc=True,
              update_from_AO=True):
    """
    This is the second part of the original KAI clean function.  The
    functionality of the first part has been replaced by various methods in
    the AOSet class.

    This program should be run from the reduce/ directory.
    Example directory structure is:
    calib/
        flats/
        flat_kp.fits
        flat.fits (optional)
        masks/
        supermask.fits
    kp/
        sci_nite1/
        sky_nite1/
        sky.fits

    All output files will be put into clean_dir (if specified, otherwise
    ../clean/) in the following structure:
    kp/
        c*.fits
        distort/
        cd*.fits
        weight/
        wgt*.fits

    The clean directory may be optionally modified to be named
    <field_><wave> instead of just <wave>. So for instance, for Arches
    field #1 data reduction, you might call clean with: field='arch_f1'.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "kp"), used as
        a wavelength suffix
    field : str, default=None
        Optional prefix for clean directory and final
        combining. All clean files will be put into <field_><wave>. You
        should also pass the same into combine(). If set to None (default)
        then only wavelength is used.
    skyscale : bool, default=False
        Whether or not to scale the sky files to the common median.
        Turn on for scaling skies before subtraction.
    skyfile : str, default=''
        An optional file containing image/sky matches.
    angOff : float, default = 0
        An optional absolute offset in the rotator
        mirror angle for cases (wave='lp') when sky subtraction is done with
        skies taken at matching rotator mirror angles.
    cent_box : int (def = 12)
        the box to use for better centroiding the reference star
    badColumns : int array, default = None
        An array specifying the bad columns (zero-based).
        Assumes a repeating pattern every 8 columns.
    fixDAR : boolean, default = True
        Whether or not to calculate DAR correction coefficients.
    use_koa_weather : boolean, default = False
        If calculating DAR correction, this keyword specifies if the atmosphere
        conditions should be downloaded from the KOA weather data. If False,
        atmosphere conditions are downloaded from the MKWC CFHT data.
    clean_dir : str, optional
        Directory where clean files will be stored. By default,
        assumes that clean files will be stored in '../clean'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """

    # Make sure directory for current passband exists and switch into it
    util.mkdir(wave)
    os.chdir(wave)

    # Determine directory locations
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')

    sciDir = waveDir + '/sci_' + nite + '/'
    util.mkdir(sciDir)
    ir.cd(sciDir)

    # Setup the clean directory
    cleanRoot = rootDir + 'clean/'

    # Check if user has specified a specific clean directory
    if clean_dir is not None:
        cleanRoot = util.trimdir(os.path.abspath(clean_dir) + '/')

    if field is not None:
        clean = cleanRoot + field + '_' + wave + '/'
    else:
        clean = cleanRoot + wave + '/'

    distort = clean + 'distort/'
    weight = clean + 'weight/'
    masks = clean + 'masks/'

    util.mkdir(cleanRoot)
    util.mkdir(clean)
    util.mkdir(distort)
    util.mkdir(weight)
    util.mkdir(masks)

    # Open a text file to document sources of data files
    data_sources_file = open(clean + 'data_sources.txt', 'a')

    try:
        # Bad pixel mask
        _supermask = redDir + 'calib/masks/supermask.fits'

        # Determine the reference coordinates for the first image.
        # This is the image for which refSrc is relevant.
        firstFile = instrument.make_filenames([files[0]], prefix='bp')[0]
        hdr1 = fits.getheader(firstFile, ignore_missing_end=True)
        radecRef = [float(hdr1['RA']), float(hdr1['DEC'])]
        aotsxyRef = kai_util.getAotsxy(hdr1)

        # Prep drizzle stuff
        # Get image size from header - this is just in case the image
        # isn't 1024x1024 (e.g., NIRC2 sub-arrays). Also, if it's
        # rectangular, choose the larger dimension and make it square
        imgsizeX = int(hdr1['NAXIS1'])
        imgsizeY = int(hdr1['NAXIS2'])

        distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr1)
        if imgsizeX >= imgsizeY:
            imgsize = imgsizeX
        else:
            imgsize = imgsizeY
        data.setup_drizzle(imgsize)

        ##########
        # Loop through the list of images
        ##########
        for f in files:
            # Define filenames
            _bp = instrument.make_filenames([f], prefix='bp')[0]
            _cd = instrument.make_filenames([f], prefix='cd')[0]
            _ce = instrument.make_filenames([f], prefix='ce')[0]
            _cc = instrument.make_filenames([f], prefix='c')[0]
            _wgt = instrument.make_filenames([f], prefix='wgt')[0]
            _statmask = instrument.make_filenames([f], prefix='stat_mask')[0]
            _crmask = instrument.make_filenames([f], prefix='crmask')[0]
            _mask = instrument.make_filenames([f], prefix='mask')[0]
            _pers = instrument.make_filenames([f], prefix='pers')[0]
            _max = _cc.replace('.fits', '.max')
            _coo = _cc.replace('.fits', '.coo')
            _rcoo = _cc.replace('.fits', '.rcoo')
            _dlog_tmp = instrument.make_filenames([f], prefix='driz')[0]
            _dlog = _dlog_tmp.replace('.fits', '.log')

            out_line = '{0} from {1} ({2})\n'.format(_cc, _bp,
                                                     datetime.now())
            data_sources_file.write(out_line)

            # Clean up if these files previously existed
            util.rmall([_cd, _ce, _cc,
                        _wgt, _statmask, _crmask, _mask, _pers, _max, _coo,
                        _rcoo, _dlog])

            # Make a static bad pixel mask ###
            # _statmask = supermask + bad columns
            data.clean_get_supermask(_statmask, _supermask, badColumns)

            # Fix cosmic rays and make cosmic ray mask. ###
            data.clean_cosmicrays(_bp, _crmask, wave)

            # Combine static and cosmic ray mask ###
            # This will be used in combine later on.
            # Results are stored in _mask, _mask_static is deleted.
            data.clean_makemask(_mask, _crmask, _statmask, wave,
                                instrument=instrument)

            # Drizzle individual file ###
            data.clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                               fixDAR=fixDAR, instrument=instrument,
                               use_koa_weather=use_koa_weather)

            # Make .max file ###
            # Determine the non-linearity level. Raw data level of
            # non-linearity is 12,000 but we subtracted
            # off a sky which changed this level. The sky is
            # scaled, so the level will be slightly different
            # for every frame.
            hdr = fits.getheader(_bp, ignore_missing_end=True)
            try:
                nonlinSky = hdr['skylev']
            except KeyError:
                nonlinSky = 0.
            coadds = fits.getval(_bp, instrument.hdr_keys['coadds'])
            satLevel = (coadds * instrument.get_saturation_level()) - nonlinSky
            file(_max, 'w').write(str(satLevel))

            # Rename and clean up files ###
            ir.imrename(_bp, _cd)

            # Make the *.coo file and update headers ###
            # First check if PA is not zero
            phi = instrument.get_position_angle(hdr)

            data.clean_makecoo(_ce, _cc, refSrc, strSrc, aotsxyRef, radecRef,
                               instrument=instrument, check_loc=check_ref_loc,
                               cent_box=cent_box, update_from_AO=update_from_AO)

            # Move to the clean directory ###
            util.rmall([clean + _cc, clean + _coo, clean + _rcoo,
                        distort + _cd, weight + _wgt,
                        clean + _ce, clean + _max,
                        masks + _mask, _ce])

            os.rename(_cc, clean + _cc)
            os.rename(_cd, distort + _cd)
            os.rename(_wgt, weight + _wgt)
            os.rename(_mask, masks + _mask)
            os.rename(_max, clean + _max)
            os.rename(_coo, clean + _coo)
            os.rename(_rcoo, clean + _rcoo)

        data_sources_file.close()
    finally:
        # Move back up to the original directory
        ir.cd('../')

    # Change back to original directory
    os.chdir('../')


def kaiclean(files, nite, wave, refSrc, strSrc, badColumns=None, field=None,
             skyscale=False, skyfile=None, angOff=0.0, cent_box=12,
             fixDAR=True, use_koa_weather=False,
             raw_dir=None, clean_dir=None,
             instrument=instruments.default_inst, check_ref_loc=True,
             update_from_AO=True):
    """
    Clean near infrared NIRC2 or OSIRIS images.

    This program should be run from the reduce/ directory.
    Example directory structure is:
    calib/
        flats/
        flat_kp.fits
        flat.fits (optional)
        masks/
        supermask.fits
    kp/
        sci_nite1/
        sky_nite1/
        sky.fits

    All output files will be put into clean_dir (if specified, otherwise
    ../clean/) in the following structure:
    kp/
        c*.fits
        distort/
        cd*.fits
        weight/
        wgt*.fits

    The clean directory may be optionally modified to be named
    <field_><wave> instead of just <wave>. So for instance, for Arches
    field #1 data reduction, you might call clean with: field='arch_f1'.

    Parameters
    ----------
    files : list of int
        Integer list of the files. Does not require padded zeros.
    nite : str
        Name for night of observation (e.g.: "nite1"), used as suffix
        inside the reduce sub-directories.
    wave : str
        Name for the observation passband (e.g.: "kp"), used as
        a wavelength suffix
    field : str, default=None
        Optional prefix for clean directory and final
        combining. All clean files will be put into <field_><wave>. You
        should also pass the same into combine(). If set to None (default)
        then only wavelength is used.
    skyscale : bool, default=False
        Whether or not to scale the sky files to the common median.
        Turn on for scaling skies before subtraction.
    skyfile : str, default=''
        An optional file containing image/sky matches.
    angOff : float, default = 0
        An optional absolute offset in the rotator
        mirror angle for cases (wave='lp') when sky subtraction is done with
        skies taken at matching rotator mirror angles.
    cent_box : int (def = 12)
        the box to use for better centroiding the reference star
    badColumns : int array, default = None
        An array specifying the bad columns (zero-based).
        Assumes a repeating pattern every 8 columns.
    fixDAR : boolean, default = True
        Whether or not to calculate DAR correction coefficients.
    use_koa_weather : boolean, default = False
        If calculating DAR correction, this keyword specifies if the atmosphere
        conditions should be downloaded from the KOA weather data. If False,
        atmosphere conditions are downloaded from the MKWC CFHT data.
    raw_dir : str, optional
        Directory where raw files are stored. By default,
        assumes that raw files are stored in '../raw'
    clean_dir : str, optional
        Directory where clean files will be stored. By default,
        assumes that clean files will be stored in '../clean'
    instrument : instruments object, optional
        Instrument of data. Default is `instruments.default_inst`
    """

    # Make sure directory for current passband exists and switch into it
    util.mkdir(wave)
    os.chdir(wave)

    # Determine directory locatons
    waveDir = os.getcwd() + '/'
    redDir = util.trimdir(os.path.abspath(waveDir + '../') + '/')
    rootDir = util.trimdir(os.path.abspath(redDir + '../') + '/')

    sciDir = waveDir + '/sci_' + nite + '/'
    util.mkdir(sciDir)
    ir.cd(sciDir)

    # Set location of raw data
    rawDir = rootDir + 'raw/'

    # Check if user has specified a specific raw directory
    if raw_dir is not None:
        rawDir = util.trimdir(os.path.abspath(raw_dir) + '/')

    # Setup the clean directory
    cleanRoot = rootDir + 'clean/'

    # Check if user has specified a specific clean directory
    if clean_dir is not None:
        cleanRoot = util.trimdir(os.path.abspath(clean_dir) + '/')

    if field is not None:
        clean = cleanRoot + field + '_' + wave + '/'
    else:
        clean = cleanRoot + wave + '/'

    distort = clean + 'distort/'
    weight = clean + 'weight/'
    masks = clean + 'masks/'

    util.mkdir(cleanRoot)
    util.mkdir(clean)
    util.mkdir(distort)
    util.mkdir(weight)
    util.mkdir(masks)

    # Open a text file to document sources of data files
    data_sources_file = open(clean + 'data_sources.txt', 'a')

    try:
        # Setup flat. Try wavelength specific, but if it doesn't
        # exist, then use a global one.
        flatDir = redDir + 'calib/flats/'
        flat = flatDir + 'flat_' + wave + '.fits'
        if not os.access(flat, os.F_OK):
            flat = flatDir + 'flat.fits'

        # Bad pixel mask
        _supermask = redDir + 'calib/masks/supermask.fits'

        # Determine the reference coordinates for the first image.
        # This is the image for which refSrc is relevant.
        firstFile = instrument.make_filenames([files[0]], rootDir=rawDir)[0]
        hdr1 = fits.getheader(firstFile, ignore_missing_end=True)
        radecRef = [float(hdr1['RA']), float(hdr1['DEC'])]
        aotsxyRef = kai_util.getAotsxy(hdr1)

        # Setup a Sky object that will figure out the sky subtraction
        skyDir = waveDir + 'sky_' + nite + '/'
        skyObj = data.Sky(sciDir, skyDir, wave, scale=skyscale,
                          skyfile=skyfile, angleOffset=angOff,
                          instrument=instrument)

        # Prep drizzle stuff
        # Get image size from header - this is just in case the image
        # isn't 1024x1024 (e.g., NIRC2 sub-arrays). Also, if it's
        # rectangular, choose the larger dimension and make it square
        imgsizeX = int(hdr1['NAXIS1'])
        imgsizeY = int(hdr1['NAXIS2'])

        distXgeoim, distYgeoim = instrument.get_distortion_maps(hdr1)
        if imgsizeX >= imgsizeY:
            imgsize = imgsizeX
        else:
            imgsize = imgsizeY
        data.setup_drizzle(imgsize)

        ##########
        # Loop through the list of images
        ##########
        for f in files:
            # Define filenames
            _raw = instrument.make_filenames([f], rootDir=rawDir)[0]
            _cp = instrument.make_filenames([f])[0]
            _ss = instrument.make_filenames([f], prefix='ss')[0]
            _ff = instrument.make_filenames([f], prefix='ff')[0]
            _ff_f = _ff.replace('.fits', '_f.fits')
            _ff_s = _ff.replace('.fits', '_s.fits')
            _bp = instrument.make_filenames([f], prefix='bp')[0]
            _cd = instrument.make_filenames([f], prefix='cd')[0]
            _ce = instrument.make_filenames([f], prefix='ce')[0]
            _cc = instrument.make_filenames([f], prefix='c')[0]
            _wgt = instrument.make_filenames([f], prefix='wgt')[0]
            _statmask = instrument.make_filenames([f], prefix='stat_mask')[0]
            _crmask = instrument.make_filenames([f], prefix='crmask')[0]
            _mask = instrument.make_filenames([f], prefix='mask')[0]
            _pers = instrument.make_filenames([f], prefix='pers')[0]
            _max = _cc.replace('.fits', '.max')
            _coo = _cc.replace('.fits', '.coo')
            _rcoo = _cc.replace('.fits', '.rcoo')
            _dlog_tmp = instrument.make_filenames([f], prefix='driz')[0]
            _dlog = _dlog_tmp.replace('.fits', '.log')

            out_line = '{0} from {1} ({2})\n'.format(_cc, _raw,
                                                     datetime.now())
            data_sources_file.write(out_line)

            # Clean up if these files previously existed
            util.rmall([_cp, _ss, _ff, _ff_f, _ff_s, _bp, _cd, _ce, _cc,
                        _wgt, _statmask, _crmask, _mask, _pers, _max, _coo,
                        _rcoo, _dlog])

            # Copy the raw file to local directory ###
            ir.imcopy(_raw, _cp, verbose='no')

            # Make persistance mask ###
            # - Checked images, this doesn't appear to be a large effect.
            # clean_persistance(_cp, _pers, instrument=instrument)

            # Sky subtract ###
            # Get the proper sky for this science frame.
            # It might be scaled or there might be a specific one for L'.
            sky = skyObj.getSky(_cp)

            ir.imarith(_cp, '-', sky, _ss)

            # Flat field ###
            ir.imarith(_ss, '/', flat, _ff)

            # Make a static bad pixel mask ###
            # _statmask = supermask + bad columns
            data.clean_get_supermask(_statmask, _supermask, badColumns)

            # Fix bad pixels ###
            # Produces _ff_f file
            bfixpix.bfixpix(_ff, _statmask)
            util.rmall([_ff_s])

            # Fix cosmic rays and make cosmic ray mask. ###
            data.clean_cosmicrays(_ff_f, _crmask, wave)

            # Combine static and cosmic ray mask ###
            # This will be used in combine later on.
            # Results are stored in _mask, _mask_static is deleted.
            data.clean_makemask(_mask, _crmask, _statmask, wave,
                                instrument=instrument)

            # Background Subtraction ###
            bkg = data.clean_bkgsubtract(_ff_f, _bp)

            # Drizzle individual file ###
            data.clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                               fixDAR=fixDAR, instrument=instrument,
                               use_koa_weather=use_koa_weather)

            # Make .max file ###
            # Determine the non-linearity level. Raw data level of
            # non-linearity is 12,000 but we subtracted
            # off a sky which changed this level. The sky is
            # scaled, so the level will be slightly different
            # for every frame.
            nonlinSky = skyObj.getNonlinearCorrection(sky)

            coadds = fits.getval(_ss, instrument.hdr_keys['coadds'])
            satLevel = (coadds * instrument.get_saturation_level()) - nonlinSky - bkg
            file(_max, 'w').write(str(satLevel))

            # Rename and clean up files ###
            ir.imrename(_bp, _cd)
            # util.rmall([_cp, _ss, _ff, _ff_f])

            # Make the *.coo file and update headers ###
            # First check if PA is not zero
            hdr = fits.getheader(_raw, ignore_missing_end=True)
            phi = instrument.get_position_angle(hdr)

            data.clean_makecoo(_ce, _cc, refSrc, strSrc, aotsxyRef, radecRef,
                               instrument=instrument, check_loc=check_ref_loc,
                               cent_box=cent_box, update_from_AO=update_from_AO)

            # Move to the clean directory ###
            util.rmall([clean + _cc, clean + _coo, clean + _rcoo,
                        distort + _cd, weight + _wgt,
                        clean + _ce, clean + _max,
                        masks + _mask, _ce])

            os.rename(_cc, clean + _cc)
            os.rename(_cd, distort + _cd)
            os.rename(_wgt, weight + _wgt)
            os.rename(_mask, masks + _mask)
            os.rename(_max, clean + _max)
            os.rename(_coo, clean + _coo)
            os.rename(_rcoo, clean + _rcoo)

            # This just closes out any sky logging files.
            # skyObj.close()
        data_sources_file.close()
    finally:
        # Move back up to the original directory
        # skyObj.close()
        ir.cd('../')

    # Change back to original directory
    os.chdir('../')


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
    kaiclean(sci_frames, obsdate, obsfilt, refSrc, refSrc, field=target,
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
