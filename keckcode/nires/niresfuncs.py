"""

niresfuncs.py

"""

import os
from matplotlib import pyplot as plt
from specim import specfuncs as ss
from specim.specfuncs import echelle1d as ech1d
from keckcode.nires import nsxspec as nsx
from keckcode.nires import nsxset

# -----------------------------------------------------------------------

def copy_raw(inroot, inframes, rawdir='../../Raw'):
    """

    Copies the files over from the raw directory.  They are copied
    just in case something corrupts the files.

    """

    """ Make the list of filenames """
    infiles = []
    for frame in inframes:
        datafile = '%s_%04d.fits' % (inroot, frame)
        infiles.append(datafile)

    """ Copy over the raw data files """
    if copyraw:
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
            linelist = ['Br-delta', 'Br-gamma']
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 1:
            linelist = ['Br6', 'Br7', 'Br8', 'Br9', 'Br10', 'Br11', 'Br12']
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 2:
            linelist = ['Pa-beta',]
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 3:
            linelist = ['Pa-gamma', 'Pa-delta']
            _mask_line_linfo(spec, lines, linelist, mode)
        if i == 4:
            linelist = ['Pa-delta',]
            _mask_line_linfo(spec, lines, linelist, mode)

# -----------------------------------------------------------------------

def telluric_corr(inspec, atmcorr, tellfile=None, doplot=True, **kwargs):
    if atmcorr == 'model':
        for spec in inspec:
            spec.atm_corr()
    elif atmcorr == 'telluric' and tellfile is not None:
        telluric = ech1d.Ech1d(echfile=tellfile)
        for spec, tellspec in zip(inspec, telluric):
            spec.atm_corr(atm='telluric', model=tellspec)
    else:
        print('')
        print('ERROR: Either atmcorr was set to an invalid value or')
        print('  atmcorr was set to "telluric" but no valid file was given')
        print('  via that tellfile parameter')
        return
    if doplot:
        inspec.plot_all(mode='atmcorr', **kwargs)

# -----------------------------------------------------------------------

def redux(inroot, inframes, copyraw, donsx, nsxmode, docoadd, outroot,
          rawdir='../../Raw', obj=None, bkgd=None, echfile=None,
          resp345file='../rspec_345.txt', respformat='345text', 
          atmcorr='model', tellfile=None, airmass=1.0, smo=None,
          z=None, **kwargs):
    """

    Code to reduce NIRES data files associated with a given targed.
    This code:
      1. Copies over the data files from the raw data directory
      2. Runs the nsx code to reduce the spectra
      3. Coadds the data files
      4. Applies a response correction and an atmospheric correction to the
         coadded data

    The steps that are actually run are selected by the user via the passed
     parameters
    """

    """ Make the list of filenames """
    infiles = []
    for frame in inframes:
        datafile = '%s_%04d.fits' % (inroot, frame)
        infiles.append(datafile)

    """ Copy over the raw data files """
    if copyraw:
        copy_raw(inroot, inframes, rawdir)

    """ Run the nsx code on the raw data files """
    if donsx:
        for i in range(0, len(inframes), 2):
            j = i+1
            file1 = infiles[i]
            file2 = infiles[j]
            extract_nsx([file1, file2], nsxmode, obj, bkgd)
            extract_nsx([file2, file1], nsxmode, obj, bkgd)
    if nsxmode == 'profonly':
        return

    """ Coadd the data """
    if docoadd:
        # sframes1 = []
        sframes2 = []
        for i in range(0, len(inframes), 2):
            j = i+1
            f1 = inframes[i]
            f2 = inframes[j]
            sframes2.append(f2)
            sframes2.append(f1)
        coadd = coadd_nsx(inroot, inframes, outroot, sframes2, **kwargs)
    elif echfile is not None:
        coadd = ech1d.Ech1d(echfile=echfile)
    else:
        print('')
        print('No coadded file provided.  Stopping before atmospheric'
              'correction')
        print('')
        return

    """ Do the response corrections (so far only 3 reddest orders) """
    # resp345 = ss.Spec1d(resp345file, informat='text')
    # for i in range(3):
    #     print(len(coadd[i]), len(resp345))
    #     coadd[i]['flux'] /= resp345['flux']
    #     coadd[i]['var'] /= resp345['flux']**2
    # coadd.save_multi('respcorr')

    """ Do an atmospheric absorption correction """
    if atmcorr is not None:
        plt.figure(2)
        telluric_corr(coadd, atmcorr, tellfile=tellfile)

    """ Smooth the spectra """
    plt.figure(3)
    coadd.plot_all(mode='atmcorr', z=z, smo=smo, **kwargs)
