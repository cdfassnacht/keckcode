
def XSolve(outprefix,newPrefix=None,indir=None,prefix=None,files=None):
    import numpy,cPickle
    from mostools import spectools
    from LRIS import cross
    obs = numpy.load('%s_flat.calib'%outprefix)

    if files is not None:
        obs.sciencefiles = []
        filenames = []
        for f in files:
            filenames.append('%s/%s%04d.fits'%(indir,prefix,f))
        obs.add_files(filenames)

    print "Preparing Science Files"
    obs.set_ysoln_arrays()
    obs.prep_science()
    obs.set_centroids()
    print "Finding approximate line widths"
    obs.set_linewidth()

    skyspex = []
    arcspex = []
    offset = []
    for i in range(len(obs.sciencefiles)):
        offset.append([])

    print "Determining Science X Solutions"
    for i in range(obs.nslits):
        slit = obs.slits[i]
        wb,wt,b,t = slit.widebottom,slit.widetop,slit.bottom,slit.top

        # First get arc solution
        slit.set_xsoln()

        yforw = obs.ysoln['forw'][b:t].copy()
        arc2d = obs.tmp['arc'][b:t].copy()
        arc = spectools.resample1d(arc2d,slit.arcsoln['forw'],'x')
        arc = numpy.sort(arc,0)

        # I guess this is to avoid reflections....
        if obs.inst['instrument']=='LRISBLUE':
            arc = arc[2]
        else:
            arc = arc[arc.shape[0]/2]

        # For top and bottom slits we must account for curvature off the mask
        if i==0 or i==obs.nslits-1:
            arc = obs.tmp['arc'][b:t].copy()
            #arc[abs(arc)<1e-3] = -1e99
            arc[arc==0] = -1e99
            indx = numpy.where(arc.sum(1)>0)[0]
            arc = arc[indx[indx.size/2]]

        arcspex.append(arc)

        trace = []
        for f in slit.obs.sciencefiles:
            xsoln = slit.xsoln[f]['forw']
            sci = obs.data[f]['data'][wb:wt].copy()
            sci = spectools.resampley(sci,yforw,wb)
            sci = spectools.resample1d(sci,xsoln,'x')
            trace.append(sci.copy())
            sci = numpy.sort(sci,0)
            sky = sci[sci.shape[0]/4]
            if i==0 or i==obs.nslits-1:
                tmp = sci.copy()
                tmp[abs(tmp)>1e-3] = 1
                tmp[tmp!=1] = 0
                cnt = tmp.sum(0).astype(numpy.int16)
                sci[abs(sci)<1e-3] = 1e9
                sci = numpy.sort(sci,0)
                tmp = sci[cnt/4]
                sky = (numpy.eye(tmp.shape[0])*tmp).sum(1)
            trace[-1] -= sky
            trace[-1] = trace[-1].sum(1)
            skyspex.append(sky.copy())
            """ START RE-ARC """
            arc = spectools.resample1d(arc2d,xsoln,'x')
            arc = numpy.sort(arc,0)

            # I guess this is to avoid reflections....
            if obs.inst['instrument']=='LRISBLUE':
                arc = arc[2]
            else:
                arc = arc[arc.shape[0]/2]

            # For top and bottom slits we must account for curvature off the mask
            if i==0 or i==obs.nslits-1:
                arc = obs.tmp['arc'][b:t].copy()
                #arc[abs(arc)<1e-3] = -1e99
                arc[arc==0] = -1e99
                indx = numpy.where(arc.sum(1)>0)[0]
                arc = arc[indx[indx.size/2]]
            arcspex[-1] = arc
            """ END RE-ARC""" 

        offset[0].append(0.)
        for i in range(1,len(offset)):
            corr = signal.correlate(trace[0],trace[i],mode='full')
            corr = ndimage.zoom(corr,4)
            offset[i].append(float(corr.argmax()-trace[0].size*4)/4.)
    for i in range(len(offset)):
        offset[i] = numpy.array(offset[i]).mean()
    offset = numpy.array(offset)
    offset -= offset.min()
    c = {}
    i = 0
    for f in obs.sciencefiles:
        c[f] = offset[i]
        i += 1
    obs.center = c

    obs.make_skymodel()
    obs.set_arcmodel()
    obs.arcspex = arcspex
    obs.skyspex = skyspex

    if newPrefix is not None:
        outprefix = newPrefix

    for key in obs.tmp.keys():
        del obs.tmp[key]
    for key in obs.data.keys():
        del obs.data[key]['data']
    del obs.ysoln['forw']
    del obs.ysoln['back']
    obs.flatnorm = None

    f = open('%s_xsolve.calib'%outprefix,'wb')
    cPickle.dump(obs,f,2)
    f.close()


def SlitCross(outprefix,forceBlue=False,showFit=True):
    import numpy,cPickle
    from LRIS.cross import make_standard
    obs = numpy.load('%s_xsolve.calib'%outprefix)
    arcspex = obs.arcspex
    skyspex = obs.skyspex
    c = obs.inst['CENTRAL']
    print "Slit cross-identification"
    if obs.inst['instrument']=='LRISBLUE' or forceBlue==True:
        solutions = make_standard(arcspex,skyspex,showFit,usearc=True)
    else:
        solutions = make_standard(skyspex,skyspex,showFit)
    f = open('%s_slitcross.calib'%outprefix,'wb')
    cPickle.dump(solutions,f,2)
    f.close()


def WaveSolve(outprefix,pre=None,forceBlue=False,offset=None,showFit=False,order=3):
    import numpy,cPickle
    from LRIS.WaveSolve import getWaveSolution as gWS
    from mostools import specfit as spf
    obs = numpy.load('%s_xsolve.calib'%outprefix)
    solutions = numpy.load('%s_slitcross.calib'%outprefix)
    obs.set_ysoln_arrays()
    obs.prep_science()

    arcspex = obs.arcspex
    skyspex = obs.skyspex
    c = obs.inst['CENTRAL']
    print "Finding Wavelength Solution"
    if obs.inst['instrument']=='LRISBLUE' or forceBlue==True:
        wsoln = gWS(solutions,obs.arcmodel,c,arcspex,True,pre,showFit=showFit)
        obs.arcmodel = wsoln[2]
    else:
        wsoln = gWS(solutions,obs.skymodel,c,skyspex,pre=pre,showFit=showFit)
    obs.wsoln = wsoln[1:]
    wavesolutions = wsoln[0]

    for i in range(obs.nslits):
        slit = obs.slits[i]
        slit.wsoln = {}
        for f in obs.sciencefiles:
            slit.wsoln[f] = wavesolutions[0]
            del wavesolutions[0]

    print "Refining wavelength solution"
    for i in range(obs.nslits):
        slit = obs.slits[i]
        if obs.inst['instrument']=='LRISBLUE' or forceBlue==True:
            #slit.set_wsoln_lines(order=order,skylines=False,arclines=True,offset=offset)
            slit.wsoln[f] = spf.wave_arclines(arcspex[i],obs.arcmodel,skyspex[i],slit.wsoln[f],neworder=order)
        else:
            #slit.set_wsoln_lines(usesky=True)
            slit.wsoln[f] = spf.wave_skylines(skyspex[i],slit.wsoln[f])

    for key in obs.tmp.keys():
        del obs.tmp[key]
    for key in obs.data.keys():
        del obs.data[key]['data']
    del obs.ysoln['forw']
    del obs.ysoln['back']
    obs.flatnorm = None

    f = open('%s_wave.calib'%outprefix,'wb')
    cPickle.dump(obs,f,2)
    f.close()
