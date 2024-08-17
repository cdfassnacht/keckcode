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

def inlist_to_framelist(inlist):
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

##########
# Make electronic logs
#    - run this first thing for a new observing run.
##########
def makelog_and_prep_images(year, instrument, rawdir='../raw'):
    """

    Make an electronic log from all the files in the ../raw/ directory.
    The file will be called kai.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda

    """

    """ Set up the instrument and make the log """
    if instrument.lower() == 'osiris':
        inst = osiris
    elif instrument.lower() == 'nirc2'
        inst = nirc2
    else:
        print('')
        raise ValueError('makelog_and_prep_images: instrument must be '
                         '"osiris" or "nirc2"\n')

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

def make_dark(darklist, outfile, instrument, texp=None):
    """

    Makes a dark frame given either an input list of integer frame numbers
    (for NIRC2) or a dict or list of dicts containing 'assn' and 'frames'
    keywords (for OSIRIS)

    """

    """ Create the framelist in the proper format """
    darkframes = inlist_to_framelist(darklist)

    """ Make the dark file """
    if texp is not None:
        darktime = float(texp)
        print('Making the dark for the %.1f sec exposures' % darktime)
    calib.makedark(darkframes, outfile, instrument=instrument)

#def go_calib(darkdict, flatondict, darkflatoffdict=None):
def make_calfiles():
    """Do the calibration reduction.

    @author Jessica Lu
    @author Sylvana Yelda
    """

    ####################
    #
    # Calibration files:
    #     everything created under calib/
    #
    ####################
    # Darks - created in subdir darks/
    #  - darks needed to make bad pixel mask
    #  - store the resulting dark in the file name that indicates the
    #    integration time
    #  - If you use the OSIRIS image, you must include the full filename
    #   in the list. 

    # Use the following format for NIRC2
    # darkFiles = range(23, 27 + 1)
    # Use the following format for OSIRIS
    darkFiles = ['i201105_a020{0:03d}_flip'.format(ii) for ii in range(11, 16)]
    
    # remove instrument = osiris if nirc2.
    print('Making dark for 6sec exposures')
    calib.makedark(darkFiles, 'Dark_006sec.fits', instrument=osiris)

    darkFiles = ['i201105_a020{0:03d}_flip'.format(ii) for ii in range(16, 21)]
    print('Making dark for 180sec exposures')
    calib.makedark(darkFiles, 'dark_180sec.fits', instrument=osiris)


    # Flats - created in subdir flats/
    offFiles = ['i201105_a020{0:03d}_flip'.format(ii) for ii in range(21, 41, 2)]
    onFiles  = ['i201105_a020{0:03d}_flip'.format(ii) for ii in range(22, 42, 2)]
    calib.makeflat(onFiles, offFiles, 'flat_Kp.fits', instrument=osiris)
    # The above's name scheme is 'flat_<filter>_<dichroic>.fits'
    
    # Masks (assumes files were created under calib/darks/ and calib/flats/)
    # Use a long exposure mask (>20 sec) for this makemask.
    print('Making supermask.fits')
    calib.makemask('Dark_180sec.fits', 'flat_Kp.fits', 'supermask.fits',
                   instrument=osiris)


def plot_image(imagePath, flip = False):

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
        img = np.flip(img, axis = 1)
        extent = [x_axis[-1], x_axis[0], y_axis[0], y_axis[-1]]
    
    # Plotting prerequisites
    vmin = 10
    vmax = 1e5
    norm = LogNorm(vmin, vmax)
       
    
    # Plot the image
    plt.figure(figsize=(10,8))

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent,
               origin = "lower")
    
    # Plot titles, etc.
    plt.colorbar(label='Starlist Magnitude (mag)')
    plt.xlabel('Pixel Coordinates (pixel)')
    plt.ylabel('Pixel Coordinates (pixel)')
    plt.axis('equal')
    plt.title(imagePath.split("/")[-1])

    return


def name_checker(a,b):
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

    """ Set the instrument """
    if instrument.lower() == 'osiris':
        inst = osiris
        listtype = dict
    elif instrument.lower() == 'nirc2'
        inst = nirc2
        listtype = int
    else:
        print('')
        raise ValueError('makelog_and_prep_images: instrument must be '
                         '"osiris" or "nirc2"\n')

    """ Get the input framelist """
    sci_frames = inlist_to_framelist(inlist)

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
    sci_frames = assn_to_framelist(assnlist, frameroot, suffix=suffix)

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
