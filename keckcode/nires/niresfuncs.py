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

def redux(inroot, inframes, copyraw, donsx, nsxmode, docoadd,
          rawdir='../../Raw', resp345file='../rspec_345.txt', 
          respformat='345text', atm='model', airmass=1.0, smo=None,
          z=None):
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
        for f in infiles:
            rawfile = os.path.join(rawdir, f)
            os.system('cp -v %s .' % rawfile)

    """ Run the nsx code on the raw data files """
    if donsx:
        for i in range(0, len(inframes), 2):
            j = i+1
            file1 = infiles[i]
            file2 = infiles[j]
            if nsxmode == 'profonly':
                os.system('nsx %s %s -autox' % (file1,file2))
                os.system('nsx %s %s -autox' % (file2,file1))
            elif nsxmode == 'auto':
                os.system('nsx %s %s -autox' % (file1,file2))
                os.system('nsx %s %s -autox' % (file2,file1))

    """ Coadd the data """
    if docoadd:
        print('Coadding input frames')
        sframes1 = []
        sframes2 = []
        for i in range(0, len(inframes), 2):
            j = i+1
            f1 = inframes[i]
            f2 = inframes[j]
            sframes1.append(f1)
            sframes1.append(f2)
            sframes2.append(f2)
            sframes2.append(f1)
        specset = nsxset.NsxSet(inroot, sframes1, sframes2)
        coadd = specset.coadd()
        coadd.save_multi('tmpcoadd')
    else:
        coadd = ech1d.Ech1d(echfile='tmpcoadd.fits')

    """ Do the response corrections (so far only 3 reddest orders) """
    resp345 = ss.Spec1d(resp345file, informat='text')
    for i in range(3):
        print(len(coadd[i]), len(resp345))
        coadd[i]['flux'] /= resp345['flux']
        coadd[i]['var'] /= resp345['flux']**2
    coadd.save_multi('respcorr')

    """ Do an atmospheric absorption correction """
    for spec in coadd:
        spec.atm_corr()
    plt.figure(2)
    coadd.plot_all(mode='atmcorr')

    """ Smooth the spectra """
    plt.figure(3)
    coadd.plot_all(mode='atmcorr', z=z, smo=smo, linetype='em')
