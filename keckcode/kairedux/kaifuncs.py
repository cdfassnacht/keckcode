# #
# General Notes:
# -- If you need help on the individual function calls,
#    then in the pyraf prompt, import the module and
#    then print the documentation for that function:
#    --> print kai.kailog.__doc__
#    --> print range.__doc__
#
##################################################

""" Import basic modules """
import warnings
import numpy as np
import os
import sys
import glob

"""
The matplotlib.use('Agg') line prevents any display functionality by 
matplotlib, which is necessary when running the code in a Docker container.
It must come before any import statement that will also lead to an
import of matplotlib, which it needs to be before the specim and kai imports
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm

""" Import additional helper modules """
from astropy.io import ascii
from specim.imfuncs.wcshdu import WcsHDU

""" Import KAI modules """
from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import util
from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments

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


def assn_to_framelist(assnlist, rootname, suffix=None):
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
    inlist = []
    if assnlist is not None:
        if isinstance(assnlist, dict):
            inlist = [assnlist]
        elif isinstance(assnlist, list):
            if isinstance(assnlist[0], dict):
                inlist = assnlist
            else:
                print('ERROR: first element of assnlist is not a dict')
                iserror = True
        else:
            print('ERROR: assnlist is not either a dict or a list')
            iserror = True
    else:
        print('ERROR: assnlist is None')
        iserror = True
    if iserror:
        raise TypeError('\nassn_to_framelist: assnlist must be a dict or a '
                        'list of dicts.\n')

    """ Create the framelist """
    framelist = []
    for i in inlist:
        """ Check the passed parameters """
        for k in ['assn', 'frames']:
            if k not in i.keys():
                raise KeyError('\nassn_to framelist: dict must contain both'
                               '"assn" and "frames" keys\n\n')
        """
        Loop through association and frame numbers to create file list
        """
        assn = i['assn']
        for j in i['frames']:
            framename = '%s%03d%03d' % (rootname, assn, j)
            if suffix is not None:
                framename = '%s_%s' % (framename, suffix)
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
    if isinstance(inlist, (dict, int)):
        framelist = [inlist]
    elif isinstance(inlist, list):
        el1 = inlist[0]
        if inst == osiris and isinstance(el1, dict):
            frameroot = 'i%s_a' % obsdate[2:]
            framelist = assn_to_framelist(inlist, frameroot, suffix=suffix)
        elif inst == nirc2 and isinstance(el1, int):
            framelist = inlist
        else:
            print('')
            print('ERROR: wrong datatype for inlist parameter.')
            if inst == osiris:
                print('For OSIRIS data type must be a dict or a list of dicts')
            else:
                print('For NIRC2 data type must be an int or a list of ints')
            raise TypeError()
    else:
        print('')
        print('ERROR: Input frame(s) [inlist parameter] must be one of the '
              'following')
        print('   1. An int (NIRC2)')
        print('   2. A dict with "assn" and "frames" keywords')
        print('   3. A list of either ints (NIRC2) or dicts (OSIRIS)')
        print('')
        raise TypeError()

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


def make_dark(darklist, obsdate, outfile, instrument, suffix=None):
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
    darkframes = inlist_to_framelist(darklist, instrument, obsdate,
                                     suffix=suffix)

    """ Make the dark file """
    print('Creating the dark file' % outfile)
    calib.makedark(darkframes, outfile, instrument=inst)


def make_flat(onlist, offlist, obsdate, outfile, instrument, suffix=None):
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

    try:
        inst = get_instrument(instrument)
    except ValueError:
        return

    """ Create the framelist in the proper format """
    onframes = inlist_to_framelist(onlist, instrument, obsdate, suffix=suffix)
    offframes = inlist_to_framelist(offlist, instrument, obsdate, suffix=suffix)

    """ Make the dark file """
    print('Creating the flat-field file' % outfile)
    calib.makeflat(onframes, offframes, outfile, instrument=inst)


def make_calfiles(obsdate, darkinfo, flatinfo, dark4mask, flat4mask,
                  instrument, suffix=None):
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

    """ Create the dark(s) as long as darkinfo is not None"""
    if darkinfo is not None:
        """ Check the darkinfo format """
        if isinstance(darkinfo, dict):
            darklist = [darkinfo]
        elif isinstance(darkinfo, list):
            if isinstance(darkinfo[0], dict):
                darklist = darkinfo
            else:
                raise TypeError
        else:
            raise TypeError

        """ Create the dark(s) """
        darkkeys = ['inlist', 'outfile']
        for info in darklist:
            for key in darkkeys:
                if key not in info.keys():
                    raise KeyError
            make_dark(info['inlist'], obsdate, info['outfile'], instrument,
                      suffix=suffix)

    """ Create the flat(s) as long as flatinfo is not None"""
    if flatinfo is not None:
        """ Check the flatinfo format """
        if isinstance(flatinfo, dict):
            flatlist = [flatinfo]
        elif isinstance(flatinfo, list):
            if isinstance(flatinfo[0], dict):
                flatlist = flatinfo
            else:
                raise TypeError
        else:
            raise TypeError

        """ Create the flat(s) """
        flatkeys = ['onlist', 'outfile']
        for info in flatlist:
            for key in flatkeys:
                if key not in info.keys():
                    raise KeyError
            if 'offlist' not in info.keys():
                offlist = range(0, 0)
            else:
                offlist = info['offlist']
            make_flat(info['onlist'], offlist, obsdate, info['outfile'],
                      instrument, suffix=suffix)

    """
    Make the 'supermask' from a dark and a flat.
    Use a long-exposure dark (>20 sec) for this.
    Note: the code assumes that the input files are found in calib/darks
     and calib/flats
    """
    print('Making supermask.fits')
    calib.makemask(dark4mask, flat4mask, 'supermask.fits', instrument=inst)


def plot_image(imagePath, flip=False):

    # Initializing the Image
    img = fits.getdata(imagePath)
   
    # Get image dimensions and make relative to reference
    x_axis = np.arange(img.shape[0], dtype=float)
    y_axis = np.arange(img.shape[1], dtype=float)

    # Extent of image to be plotted in imshow
    extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]
    
    # Flips image in case it's backwards
    if flip:
        x_axis *= -1
        img = np.flip(img, axis=1)
        extent = [x_axis[-1], x_axis[0], y_axis[0], y_axis[-1]]
    
    # Plotting prerequisites
    vmin = 10
    vmax = 1e5
    norm = LogNorm(vmin, vmax)
       
    # Plot the image
    plt.figure(figsize=(10, 8))

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent,
               origin="lower")
    
    # Plot titles, etc.
    plt.colorbar(label='Starlist Magnitude (mag)')
    plt.xlabel('Pixel Coordinates (pixel)')
    plt.ylabel('Pixel Coordinates (pixel)')
    plt.axis('equal')
    plt.title(imagePath.split("/")[-1])

    return


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
             suffix=None):
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
    combdir = '%s/combo' % os.getcwd()
    combroot = os.path.join(combdir, 'mag%s_%s_%s' % (obsdate, target, obsfilt))
    scifile = '%s.fits' % combroot
    sigfile = '%s_sig.fits' % combroot
    outroot = os.path.join(combdir, '%s_kai_%s_%s' % (target, obsdate, obsfilt))
    outsci = '%s.fits' % outroot
    outwht = '%s_wht.fits' % outroot

    """ Make the list of science frames from the input assn list"""
    frameroot = 'i%s_a' % obsdate[2:]
    sci_frames = assn_to_framelist(inlist, frameroot, suffix=suffix)

    """ Read in the science and "sig" files """
    sciin = WcsHDU(scifile)
    whtin = WcsHDU(sigfile)

    """
    Get some information from the drizzle output and put it into SHARP
    standardized form
    """
    hdr = sciin.header
    if 'ITIME' in hdr.keys():
        tottime = hdr['itime'] / 1000.
        hdr['elaptime'] = tottime
        hdr['exptime'] = tottime
        avgtime = tottime / len(sci_frames)
    else:
        avgtime = 1.
    if 'FILTER' not in hdr.keys() and 'IFILTER' in hdr.keys():
        hdr['filter'] = hdr['ifilter']
    if 'NDRIZIM' in hdr.keys():
        hdr['ncombine'] = hdr['ndrizim']
    for i, f in enumerate(sci_frames):
        hdr.set('orig%03d' % (i + 1), f, 'Original input file %d' % (i + 1))

    """
    Use the WcsHDR functionality to:
      1. Set the pixel scale to 0.01 arcsec
      2. Convert the pixel units to e-/sec
    """
    pixscale = 0.01
    if 'SYSGAIN' in hdr.keys():
        gain = hdr['sysgain']
    else:
        gain = -1.
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
    keeplist = ['object', 'telescope', 'instrume', 'filter', 'date-obs',
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
