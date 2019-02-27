"""

niresfuncs.py

"""

import os
from matplotlib import pyplot as plt
from astropy.modeling.blackbody import blackbody_lambda
from specim import specfuncs as ss
from keckcode.nires import nsxspec as nsx
from keckcode.nires import nsxset, nires1d

# -----------------------------------------------------------------------

def set_inputs(inroot, inframes):
    """

    Sets up the list of input files and frames.

    NOTE: This assumes that all of the data were obtained using an ABBA
    dither pattern.  If not then the files2 and frames2 variables that
    are returned by this function will have to be reset manually

    """

    """
    Make the list of second-object frames, assuming an ABBA dither pattern
    """
    frames2 = []
    for i in range(0, len(inframes), 2):
        j = i+1
        f1 = inframes[i]
        f2 = inframes[j]
        frames2.append(f2)
        frames2.append(f1)

    """ Make the list of filenames """
    infiles = []
    infiles2 = []
    for frame, frame2 in zip(inframes, frames2):
        file1 = '%s_%04d.fits' % (inroot, frame)
        file2 = '%s_%04d.fits' % (inroot, frame2)
        infiles.append(file1)
        infiles2.append(file2)

    """ Return the relevant information """
    return infiles, infiles2, frames2

# -----------------------------------------------------------------------

def copy_raw(infiles, rawdir='../../Raw'):
    """

    Copies the files over from the raw directory.  They are copied
    just in case something corrupts the files.

    """

    """ Copy over the raw data files """
    for f in infiles:
        rawfile = os.path.join(rawdir, f)
        os.system('cp -v %s .' % rawfile)

# -----------------------------------------------------------------------

def extract_nsx(infiles, nsxmode, obj=None, bkgd=None):
    """

    Calls Tom Barlow's nsx program to:
      1. Calibrate the raw data
      2. Rectifiy the orders
      3. Produce spatial profile information
      4. Extract the spectrum from each order, if requested

    Required inputs:
      infiles - one or two file names.  If two filenames are given, they
                are expected to represent an AB pair, for which the
                telescope was dithered between the exposures.  Both files
                will be processed, but only one spectrum will be extracted
      nsxmode - One of three choices for how to run nsx:
                 1. 'profonly' - run nsx so that it calibrates and rectifies
                      the spectra, and produces profile information, but
                      does not extract any 1d spectra
                 2. 'auto' - do all four steps above, where the object and
                      background apertures for the extraction are chosen
                      automatically
                 3. 'manual' - do all four steps above, where the object and
                      background apertures for the extraction are chosen
                      manually.  This aperture information is provided via
                      the obj and bkgd parameters
    """

    """ Make a string containing the file information """
    if isinstance(infiles, list):
        files = '%s %s' % (infiles[0], infiles[1])
    elif isintance(infiles, str):
        files = '%s'
    else:
        errstr = 'ERROR: infiles must be a string or a list of strings'
        raise ValueError(errstr)

    """ Run nsx in the requested mode """
    if nsxmode == 'profonly':
        os.system('nsx %s' % files)
    elif nsxmode == 'auto':
        os.system('nsx %s -autox' % files)
    elif nsxmode == 'manual':
        if obj is not None:
            objstr = 'sp=%s,%s' % (obj[0], obj[1])
        else:
            errstr = '\nERROR: manual nsx mode requested but '
            errstr += 'obj parameter is None\n\n'
            raise ValueError(errstr)
        if bkgd is not None:
            if isinstance(bkgd, list):
                bkgdstr = ''
                if isinstance(bkgd[0], list) or isinstance(bkgd[0], tuple):
                    for bk in bkgd:
                        bkgdstr += 'bk=%s,%s ' % (bk[0], bk[1])
                else:
                    bkgdstr += 'bk=%s,%s ' % (bkgd[0], bkgd[1])
                # for bk in bkgd:
                #     bkgdstr += 'bk=%s,%s ' % (bk[0], bk[1])
            elif isinstance(bkgd, tuple):
                bkgdstr = 'bk=%s,%s ' % (bkgd[0], bkgd[1])
            else:
                errstr = '\nERROR: manual nsx mode requested but '
                errstr += 'bkgd parameter is None\n\n'
                raise ValueError(errstr)
        execstr = 'nsx %s %s %s' % (files, objstr, bkgdstr)
        # print(execstr)
        os.system(execstr)

# -----------------------------------------------------------------------

def coadd_nsx(inroot, inframes, outroot, inframes2=None, outsuff='coaddv1',
              verbose=True, **kwargs):
    """

    Coadds the extracted spectra

    """

    if verbose:
        print('Coadding input frames')
        print('---------------------')

    """ Load the files into a NsxSet structure and coadd them """
    specset = nsxset.NsxSet(inroot, inframes, inframes2)
    coaddspec = specset.coadd(**kwargs)

    """ Save the output file """
    outbase = '%s_%s' % (outroot, outsuff)
    coaddspec.save_multi(outbase)
    return coaddspec

# -----------------------------------------------------------------------

def _mask_line_linfo(spec, lines, linelist, mode):
    """

    Actually calls the line masking code.
    This method is called from mask_lines_std and is not meant to be
     called interactively.

    """

    for line in linelist:
        linfo = lines[line]
        spec.mask_line([linfo[0], linfo[1]], linfo[2], mode=mode)
        if mode == 'atmcorr':
            mask = (spec['wav'] > linfo[0]) & (spec['wav'] < linfo[1])
            spec['flux'][mask] = spec.atmcorr['flux'][mask]

# -----------------------------------------------------------------------

def mask_lines_std(stdspec, mode='atmcorr'):
    """

    For each of the echelle orders in a telluric standard spectrum,
    replace the spectral features due to stellar absorption with
    an estimate of the continuum in that region.

    """

    lines = {
        'Pa-delta': [10006., 10100., 100.],
        'Pa-gamma': [10877., 11003., 100.],
        'Pa-beta': [12763., 12872., 100.],
        'Br12': [15520., 15584., [100., 70.]],
        'Br11': [15657., 15851., 70.],
        'Br10': [15846., 16943., 70.],
        'Br9': [16040., 16224., [70., 70.]],
        'Br8': [16318., 16495., [70.,100.]],
        'Br7': [16723., 16907., 100.],
        'Br6': [17303., 17450., 100.],
        'Br-delta': [19427., 19515., 100.],
        'Br-gamma': [21488., 21773., 300.]
        }

    for i, spec in enumerate(stdspec):
        if i == 0:
            linelist = ['Pa-delta',]
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 1:
            linelist = ['Pa-gamma', 'Pa-delta']
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 2:
            linelist = ['Pa-beta',]
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 3:
            linelist = ['Br6', 'Br7', 'Br8', 'Br9', 'Br10', 'Br11', 'Br12']
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 4:
            linelist = ['Br-delta', 'Br-gamma']
            _mask_line_linfo(spec, lines, linelist, mode)

# -----------------------------------------------------------------------

def telluric_corr(inspec, atmcorr, tellfile=None, airmass=1.0,
                  airmass_std=1.0, doplot=True, verbose=True, **kwargs):

    if verbose:
        print('Doing telluric correction')
        print('-------------------------')
    if atmcorr == 'model':
        for spec in inspec:
            spec.atm_corr()
    elif atmcorr == 'telluric' and tellfile is not None:
        telluric = nires1d.Nires1d(tellfile)
        for spec, tellspec in zip(inspec, telluric):
            spec.atm_corr(atm='telluric', model=tellspec)
    else:
        print('')
        print('ERROR: Either atmcorr was set to an invalid value or')
        print('  atmcorr was set to "telluric" but no valid file was given')
        print('  via the tellfile parameter')
        return
    if doplot:
        inspec.plot_all(mode='atmcorr', **kwargs)

# -----------------------------------------------------------------------

def do_respcorr(inspec, respcorr, respfile=None, Tstar=1.e4, doplot=True,
                **kwargs):
    """
    Does a response correction.

    The default behavior, which occurs when respcorr is 'model', involves
     generating a set of blackbody curves over the wavelength ranges of the
     input spectra and then multiplying the spectra by those curves.
     This should work because the assumption is that the input spectra have
     already been divided by the standard to correct for the atmosphere.
     
    """
    if respcorr == 'model':
        for spec in inspec:
            bbspec = blackbody_lambda(spec['wav'], Tstar).value
            spec.resp_corr(bbspec, mode='atmcorr')
    else:
        print('')
        print("NOTE: only respcorr='model' allowed at this time")
        print('')
        return

# -----------------------------------------------------------------------

def redux(inroot, inframes, copyraw, donsx, nsxmode, docoadd, outroot,
          rawdir='../../Raw', aplist=None, bkgd=None, echfile=None,
          atmcorr='model', tellfile=None, airmass=1.0,
          respcorr='model', respfile=None,
          smo=None, z=None, **kwargs):
    """

    Code to reduce NIRES data files associated with a given target.
    This code:
      1. Copies over the data files from the raw data directory
      2. Runs the nsx code to reduce the spectra
      3. Coadds the data files
      4. Applies a response correction and an atmospheric correction to the
         coadded data

    The steps that are actually run are selected by the user via the passed
     parameters
    """

    """ Set some defaults """
    mode = 'input'
    fignum = 0

    """ Set up the list of input files / frames """
    infiles, infiles2, sframes2 = set_inputs(inroot, inframes)

    """ Copy over the raw data files """
    if copyraw:
        copy_raw(infiles, rawdir)

    """ Run the nsx code on the raw data files """
    if donsx:
        if nsxmode == 'auto' or nsxmode == 'profonly':
            for file1, file2 in zip(infiles, infiles2):
                extract_nsx([file1, file2], nsxmode)
                extract_nsx([file2, file1], nsxmode)
        elif nsxmode == 'manual':
            if aplist is not None:
                for file1, file2, obj in zip(infiles, infiles2, aplist):
                    extract_nsx([file1, file2], nsxmode, obj, bkgd)
                    extract_nsx([file2, file1], nsxmode, obj, bkgd)
        else:
            raise NameError('Invalid choice for nsxmode')

    if nsxmode == 'profonly':
        return

    """ Coadd the data """
    if docoadd:
        fignum += 1
        coadd = coadd_nsx(inroot, inframes, outroot, sframes2, doplot=False,
                          **kwargs)
        mode = 'input'
    elif echfile is not None:
        coadd = nires1d.Nires1d(echfile)
        mode = 'input'
    else:
        print('')
        print('No coadded file provided.  Stopping before atmospheric'
              'correction')
        print('')
        return
    coadd.plot_all(mode=mode, z=z, smo=smo, **kwargs)

    """ Do an atmospheric absorption correction """
    if atmcorr is not None:
        print('')
        print('atmcorr = %s' % atmcorr)
        print('tellfile = %s' % tellfile)
        print('')
        fignum += 1
        plt.figure(fignum)
        telluric_corr(coadd, atmcorr, tellfile=tellfile, smo=smo,
                      title='Corrected for atmosphere')
        mode = 'atmcorr'
        coadd.save_multi('%s_atmcorr' % outroot, mode='atmcorr')

    """
    Do the response corrections
    This must be done after the atmospheric correction has been done.
    If the respcorr parameter is 'model' then just generate a blackbody
     curve for the wavelength range of the spectrum to be corrected
    """
    if respcorr is not None:
        # fignum += 1
        # plt.figure(fignum)
        do_respcorr(coadd, respcorr, respfile=respfile)
        mode = 'respcorr'
        mode = 'atmcorr' # Until a proper respcorr method can be written
        coadd.save_multi('%s_respcorr' % outroot, mode='respcorr')

    """ """

    """ Plot the spectra """
    fignum += 1
    plt.figure(fignum)
    coadd.plot_all(mode=mode, z=z, smo=smo, **kwargs)

    return coadd
