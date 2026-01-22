import sys
import os

import numpy as np
from astropy.io import fits as pf

from kai import instruments
from kai.reduce import data
# from kai.reduce import kai_util
from kai.reduce import util
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

    def __init__(self, inlist, instrument, obsdate, indir=None, gzip=False,
                 verbose=True, **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nKaiSet: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Set up the KaiSet container by calling the superclass """
        if pyversion == 2:
            super(KaiSet, self).__init__(inlist, instrument, obsdate,
                                         indir=indir,
                                         gzip=gzip, verbose=verbose, **kwargs)
        else:
            super().__init__(inlist, instrument, obsdate, indir=indir,
                             gzip=gzip, verbose=verbose, **kwargs)

        """ Get the instrument in KAI format """
        self.inst = None
        try:
            self.get_instrument(instrument)
        except ValueError:
            print('')
            print('Could not create kaiset object')
            print('')
            return

    #  ------------------------------------------------------------------------

    def get_instrument(self, instrument):
        """

        Makes sure that either NIRC2 or OSIRIS has been selected

        """
        if instrument.lower() == 'osiris' or instrument.lower() == 'osim':
            self.inst = osiris
        elif instrument.lower()[:5] == 'nirc2':
            self.inst = nirc2
        else:
            print('')
            raise ValueError('get_instrument: instrument must be '
                             '"osiris" or "nirc2"\n')

    #  ------------------------------------------------------------------------

    def add_def_hdrinfo(self, inpref='bp', maxpref='c', verbose=True):
        """

        Adds keywords that are needed for later procrssing to the headers of
         the images.
        Also write the non-linearity levels out to a set of *.max files

        """

        """ Create a list of filenames for the output *.max files """
        maxlist = self.make_outlist([inpref, '.fits'], [maxpref, '.max'])

        """ Loop through the images """
        if verbose:
            print('Updating headers and finding non-linearity levels')
        for hdu, f, m in zip(self, self.datainfo['infile'], maxlist):

            """ Get the central wavelength of the filter being used """
            hdr = hdu.header
            defwave = self.inst.get_central_wavelength(hdr)

            """ Add the new wavelength header cards """
            hdr['effwave'] = defwave
            hdr['cenwave'] = defwave
            hdr['camname'] = 'narrow'

            """
            Add a header card for the count level where the detector goes
             non-linear
            The raw data level for non-linearity is 12,000 but we subtracted
             off a sky which changed this level. The sky level may have been
             scaled, so the level could be slightly different for every frame.
            """
            try:
                nonlinSky = hdr['skylev']
            except KeyError:
                nonlinSky = 0.
            coaddkey = self.inst.hdr_keys['coadds']
            try:
                coadds = hdr[coaddkey]
            except KeyError:
                coadds = 1
            sat_level_units = self.inst.get_saturation_level_units()
            if sat_level_units == 'DN':
                satLevel = coadds * self.inst.get_saturation_level(hdr) \
                           - nonlinSky
            elif sat_level_units == 'DN/coadd':
                satLevel = self.inst.get_saturation_level(hdr) - nonlinSky
            else:
                raise ValueError('Unexpected value for saturation level units')

            # satLevel = (coadds * self.inst.get_saturation_level()) - nonlinSky
            hdr['satlevel'] = satLevel

            """ Save the updated information """
            hdu.save(f)
            with open(m, 'w') as maxout:
                maxout.write(str(hdr['satlevel']))
            if verbose:
                print('Saved non-linearity level to %s.  Value = %.2f' %
                      (m, satLevel))

    #  ------------------------------------------------------------------------

    def make_mask(self, obsfilt, inpref='bp', crpref='crmask', maskpref='mask',
                  caldir='kaidefault', badColumns=None, verbose=True):
        """

        Make an individual mask for each image.  This mask is the combination
        of two separate masks:
            1. The "static mask", which is in turn a combination of two inputs:
                a. The "supermask" created from the dark and domeflat files
                b. A list of bad columns
            2. A cosmic ray mask, which is created by this make_mask method.

        The cosmic ray mask is created  by (eventually) calling the iraf
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

        """
        Make a static pixel mask, which is the supermask plus bad columns
        """
        if self.maskdir is None:
            self.set_caldirs(caldir=caldir)
        if caldir == 'kaidefault':
            relmaskdir = '../../%s' % self.maskdir
        else:
            relmaskdir = self.maskdir
        _supermask = os.path.join(relmaskdir, 'supermask.fits')
        _statmask = 'static_mask.fits'
        util.rmall([_statmask])
        data.clean_get_supermask(_statmask, _supermask, badColumns)
        staticMask = pf.getdata(_statmask)

        """ Make the input and output lists """
        inlist = self.datainfo['basename']
        crlist = self.make_outlist(inpref, crpref)
        masklist = self.make_outlist(inpref, maskpref)

        """ Loop through the files, creating cosmic ray masks """
        if verbose:
            print('')
            print('Creating masks')
            print('-------------------------')
        for i, cr, msk in zip(inlist, crlist, masklist):
            if verbose:
                print('Making cosmic-ray mask:   %s ---> %s' % (i, cr))
            if os.path.isfile(cr):
                os.remove(cr)
            data.clean_cosmicrays(i, cr, obsfilt.lower(), _supermask)

            """ Combine the new cosmic ray mask with the static mask """
            cosmicMask = pf.getdata(cr)
            mask = staticMask + cosmicMask

            """ Add the mid-IR mask if appropriate """
            if self.instrument == 'nirc2' and mask.shape[0] > 512 and \
                    ('lp' in obsfilt or 'ms' in obsfilt):
                module_dir = os.path.dirname(data.__file__)
                _lpmask = module_dir + '/masks/nirc2_lp_edgemask.fits'
                lpmask = pf.getdata(_lpmask)
                mask += lpmask

            """
            Make the final mask. 
            Note that drizzle expects bad pixels to have value 0 and good 
             pixels to have value 1, which is the opposite of the masks that
             have been made so far
            """
            outMask = np.zeros(mask.shape)
            outMask[mask == 0] = 1

            """
            Trim 12 rows from top and bottom for NIRC2 b/c the distortion 
            solution introduces a torque to the image.
            """
            if self.instrument == 'nirc2':
                outMask[1012:1024, 0:1024] = 0
                outMask[0:12, 0:1024] = 0

            """ Save the output mask """
            pf.PrimaryHDU(outMask).writeto(msk)
            if verbose:
                print('Created mask for drizzle: %s' % msk)

    #  ------------------------------------------------------------------------

    def dewarp(self, inpref='bp', outpref='ce', whtpref='wgt', logpref='driz',
               xdistmap=None, ydistmap=None, fixDAR=True,
               use_koa_weather=False, verbose=True):
        """

        Corrects for the optical distortion of the files, i.e., "dewarps" them,
        by drizzling them using the provided distortion map.

        Inputs:
        xdistmap - x distortions as a function of position.  The default value
                   (None), means use the KAI instrument class to get the
                   distortion map
        ydistmap - y distortions as a function of position.  The default value
                   (None), means use the KAI instrument class to get the
                   distortion map
        inpref   - prefix in the input filenames to be replaced in the output
                   filenames
        outpref  - prefix to be used for the output files
        whtpref  - prefix to be used for the output weight files
        logpref  - prefix to be used for the logfile that stores the text
                   generated by the drizzle process
        fixDAR   - Fix differential atmospheric distortion? Default is True
        use_koa_weather - ?

        """

        """
        Prep drizzle stuff, using the image size from the first image.
        The image size is needed just in case the image isn't 1024x1024
         (e.g., NIRC2 sub-arrays).
        Also, if it's rectangular, choose the larger dimension and make it
         square
        """
        hdr1 = self[0].header
        imgsizeX = int(hdr1['NAXIS1'])
        imgsizeY = int(hdr1['NAXIS2'])
        if imgsizeX >= imgsizeY:
            imgsize = imgsizeX
        else:
            imgsize = imgsizeY
        data.setup_drizzle(imgsize)

        """ Set up the distortion maps"""
        if xdistmap is not None and ydistmap is not None:
            distXgeoim = xdistmap
            distYgeoim = ydistmap
        else:
            distXgeoim, distYgeoim = self.inst.get_distortion_maps(hdr1)

        """ Set up the output file lists """
        inlist = self.datainfo['infile'].copy()
        outlist = self.make_outlist(inpref, outpref)
        whtlist = self.make_outlist(inpref, whtpref)
        loglist = self.make_outlist([inpref, '.fits'], [logpref, '.log'])

        """ Loop through the images, dewarping each one """
        if verbose:
            print('Dewarping images (using drizzle)')
            print('--------------------------------')
        for _bp, _ce, _wgt, _dlog in zip(inlist, outlist, whtlist, loglist):
            if verbose:
                print('Drizzing individual file: %s --> %s' % (_bp, _ce))
            data.clean_drizzle(distXgeoim, distYgeoim, _bp, _ce, _wgt, _dlog,
                               fixDAR=fixDAR, instrument=self.inst,
                               use_koa_weather=use_koa_weather)
            print('')
