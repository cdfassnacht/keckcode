# #
# General Notes:
# -- If you need help on the individual function calls,
#    then in the pyraf prompt, import the module and
#    then print the documentation for that function:
#    --> print kai.kailog.__doc__
#    --> print range.__doc__
#
##################################################

# Turn off header deprecation warnings
import warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)

"""
These two lines prevent any display functionality by matplotlib, which
is necessary when running the code in a Docker container
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LogNorm

import numpy as np
import os, sys
import glob

# Import KAI modules
from kai.reduce import calib
from kai.reduce import sky
from kai.reduce import data
from kai.reduce import util
from kai.reduce import dar
from kai.reduce import kai_util
from kai import instruments

""" Define a global variable """
osiris = instruments.OSIRIS()


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


##########
# Make electronic logs
#    - run this first thing for a new observing run.
##########
def makelog_and_prep_images(year, instrument=osiris, rawdir='../raw'):
    """

    Make an electronic log from all the files in the ../raw/ directory.
    The file will be called kai.log and stored in the same directory.

    @author Jessica Lu
    @author Sylvana Yelda

    """
    print('Making instrument log from files')
    kai_util.makelog(rawdir, instrument=instrument)

    """
    If you are reducing OSIRIS, you need to flip the images first.
    """
    if instrument == osiris:
        print('Flipping data (needed for OSIRIS images)')
        raw_files = glob.glob('%s/i*.fits' % rawdir)
        osiris.flip_images(raw_files)

    """ Download weather data we will need. """
    print('Downloading weather data for year %s' % year)
    dar.get_atm_conditions(year)

    return

##########
# Reduce
##########
#def go_calib(darkdict, flatondict, darkflatoffdict=None):
def go_calib():
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

    plt.imshow(img, cmap='gist_heat_r', norm=norm, extent=extent, origin = "lower")
    
    # Plot titles, etc.
    plt.colorbar(label='Starlist Magnitude (mag)')
    plt.xlabel('Pixel Coordinates (pixel)')
    plt.ylabel('Pixel Coordinates (pixel)')
    plt.axis('equal')
    plt.title(imagePath.split("/")[-1])
    
    
    return

def name_checker(a,b):
    length = len(a) + len(b)
    if(length > 12+8):
        if(length != 25):
            print("Check your ABCs, length is " + str(length))
            print("(Length should be 25 or below 21)")
            sys.exit()


def go(target, obsdate, assnlist, obsfilt, refSrc, usestrehl=False,
       dockerun=False):
    """
    Do the full data reduction.

    Inputs:
      target   -
      obsdate  -
      assnlist -
      refSrc   - make sure you use the position in the _flipped_ image.
    """	

    ##########
    # Kp-band reduction
    ##########

    #    -- If you have more than one position angle, make sure to
    #       clean them seperatly.
    #    -- Strehl and Ref src should be the pixel coordinates of a bright
    #       (but non saturated) source in the first exposure of sci_files.
    #    -- If you use the OSIRIS image, you must include the full filename in the list. 

    """ 
    Download weather data we will need.
    This step is only needed if running the code within docker, and not even
     then if the "makelog_and_prep_images" function has been run within the same 
     docker session that is being used to run this "go" function  
    """
    obsyear = obsdate[:4]
    if dockerun:
        print('Downloading weather data for year %s' % obsyear)
        dar.get_atm_conditions(obsyear)

    """ Make the list of science frames from the input assn list"""
    #  name_checker(epoch,target) - don't need this any more
    frameroot = 'i%s' % obsdate[2:]
    sci_frames = assn_to_framelist(assnlist, frameroot)
    # sci_files = ['i201105_a007{0:03d}_flip'.format(ii) for ii in range(2, 11)]

    """ For this target, use the sky created for 2022_06_05 """
    # sky.makesky(sky_files, target, 'Kp', instrument=osiris)
    data.clean(sci_frames, obsdate, obsfilt, refSrc, refSrc, field=target,
               instrument=osiris)
    if usestrehl:
        data.calcStrehl(sci_frames, obsfilt, field=target, instrument=osiris)
        combwht = 'strehl'
    else:
        combwht = None
    data.combine(sci_frames, obsfilt, obsdate, field=target,
                 trim=0, weight=combwht, submaps=3, instrument=osiris)
                     
                      