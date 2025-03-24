
def resample(inprefix,filename,outprefix,mswave=None,nobgsub=True,doCR=False,clobber=False):
    import pyfits,cPickle
    import scipy
    from mostools.MOSskysub import doskysub as skysub
    from math import ceil

    f = open(inprefix+"_wave.calib")
    obs = cPickle.load(f)
    f.close()
    del f

    if obs.side is not None:
        obs.set_side()

    """
    Resampling part of the code.
    """
    obs.inst = None
    oldFile = obs.sciencefiles[0]
    obs.add_files([filename])

    npix = obs.shape[1]
    nsci = len(obs.sciencefiles)

    obs.prep_science()
    obs.data[oldFile] = obs.data[filename]
    obs.sciencefiles.remove(filename)
    del obs.data[filename]

    npix = obs.shape[1]
    nsci = len(obs.sciencefiles)

    center = scipy.zeros(nsci)

    osize = 0
    csize = 0
    for i in range(len(obs.slits)):
        slit = obs.slits[i]
        a,b = slit.bottom,slit.top
        osize += b-a+5
        csize += b-a+5+ceil(center.max())
    osize -= 5
    csize -= 5

    outlength = int(npix*1.6)
    out = scipy.zeros((nsci,osize,outlength),scipy.float32)*scipy.nan
    if nobgsub is False:
        out2 = scipy.zeros((2,csize,outlength),scipy.float32)*scipy.nan
    else:
        varout = out.copy()
        bgout = out.copy()

    pos = 0
    pos2 = 0

    del obs.flatnorm

    from mostools import MOSskysub
    MOSskysub.RESAMPLE = nobgsub
    dichroic = obs.inst['dichroic']
    if obs.inst['instrument']=='LRISBLUE':
        if dichroic=='' or dichroic=='mirror':
            cutoff = 3500.
        elif dichroic=='680':
            cutoff = [3500.,6900.]
        elif dichroic=='560':
            cutoff = [3500.,5700.]
        elif dichroic=='500':
            cutoff = [3500.,5100.]
        elif dichroic=='460':
            cutoff = [3500.,4800.]
    else:
        if dichroic=='' or dichroic=='mirror':
            cutoff = 3900.
        elif dichroic=='680':
            cutoff = 6700.
        elif dichroic=='560':
            cutoff = 5500.
        elif dichroic=='500':
            cutoff = 4900.
        elif dichroic=='460':
            cutoff = 3900. 

    if mswave is None:
        mswave = obs.data[obs.sciencefiles[0]]['mswave']

    print "Determining X-Y-Lambda Solutions"
    for i in range(obs.nslits):
        slit = obs.slits[i]
        slit.set_xywsoln()

    print "Resampling"
    for i in range(obs.nslits):
        slit = obs.slits[i]
        a,b = slit.widebottom,slit.widetop

        sky2x = []
        sky2y = []
        ccd2wave = []

        sci = scipy.empty((nsci,b-a,npix))
        for j in range(nsci):
            f = obs.sciencefiles[j]

            sci[j] = obs.data[f]['data'][a:b].copy()

            sky2x.append(slit.xywsoln[f]['sky2x'])
            sky2y.append(slit.xywsoln[f]['sky2y'])
            ccd2wave.append(slit.xywsoln[f]['ccd2wave'])

            a,b = slit.widebottom,slit.widetop

        scale = obs.scale

        strt,bgsub,varimg,wave = skysub(slit.bottom,slit.top-slit.bottom,
                                outlength,sci,obs.ysoln['back'][a:b],sky2x,
                                sky2y,ccd2wave,scale,mswave,
                                center,cutoff,obs.readvar,doCR)
#                                center,cutoff,[obs.data[f]['airmass']])

#        varimg += 36. # include readnoise!
        varimg += obs.readvar
        out[:,pos:pos+slit.top-slit.bottom,:] = strt.copy()

        if nobgsub is False:
            h = bgsub.shape[0]
            out2[0,pos2:pos2+h,:] = bgsub.copy()
            out2[1,pos2:pos2+h,:] = varimg.copy()
            pos2 += h+5
        else:
            varout[:,pos:pos+slit.top-slit.bottom,:] = varimg.copy()
            bgout[:,pos:pos+slit.top-slit.bottom,:] = bgsub.copy()
        pos += (slit.top-slit.bottom)+5


    tmp = out[0].copy()
    tmp = scipy.where(scipy.isnan(tmp),0,1)
    mod = scipy.where(tmp.sum(axis=0)!=0)
    start = mod[0][0]
    end = mod[0][-1]+1
    wave = wave[start:end]

    outname = outprefix+"_straight.fits"
    outfile = pyfits.PrimaryHDU(out[:,:,start:end])
    outfile.header.update('CTYPE1','LINEAR')
    outfile.header.update('CRPIX1',1)
    outfile.header.update('CRVAL1',wave[0])
    outfile.header.update('CD1_1',wave[1]-wave[0])
    outfile.header.update('CRPIX2',1)
    outfile.header.update('CRVAL2',1)
    outfile.header.update('CD2_2',1)
    if nsci>1:
        outfile.header.update('CRPIX3',1)
        outfile.header.update('CRVAL3',1)
        outfile.header.update('CD3_3',1)
    outfile.writeto(outname,clobber=clobber)

    if nobgsub:
        hdr = outfile.header
        outname = outprefix+"_var.fits"
        outfile = pyfits.PrimaryHDU(varout[:,:,start:end])
        outfile.header = hdr
        outfile.writeto(outname,clobber=clobber)
        outname = outprefix+"_bgsub.fits"
        outfile = pyfits.PrimaryHDU(bgout[:,:,start:end])
        outfile.header = hdr
        outfile.writeto(outname,clobber=clobber)
        return

    tmp = out2[0].copy()
    tmp = scipy.where(scipy.isnan(tmp),0,1)
    mod = scipy.where(tmp.sum(axis=0)!=0)
    start = mod[0][0]
    end = mod[0][-1]+1
    del tmp

    # Background Subtracted
    outname = outprefix+"_bgsub.fits"
    outfile = pyfits.PrimaryHDU(out2[0,:,start:end])
    outfile.header.update('CTYPE1','LINEAR')
    outfile.header.update('CRPIX1',1)
    outfile.header.update('CRVAL1',wave[0])
    outfile.header.update('CD1_1',wave[1]-wave[0])
    outfile.header.update('CRPIX2',1)
    outfile.header.update('CRVAL2',1)
    outfile.header.update('CD2_2',1)
    outfile.writeto(outname,clobber=clobber)
    hdr = outfile.header.copy()

    # An estimate of the variance in each pixel
    outname = outprefix+"_var.fits"
    outfile = pyfits.PrimaryHDU(out2[1,:,start:end])
    outfile.header = hdr
    outfile.writeto(outname,clobber=clobber)
    del out2,hdr,outfile

