import scipy,numpy
import special_functions as sf


def new_xsoln(data,xord=3,yord=3,tol=2.):
    data = data.copy()

    height,width = data.shape

    indx = height/4
    mid = data[indx]

    xvals = scipy.arange(mid.size)
    middle_lines = get_lines(xvals,mid,False)

    ref = []
    for row in range(height):
        ref.append(get_lines(xvals,data[row],False))
    lines = []
    for l in middle_lines:
        """
        off is used to track the systematic shift of lines as we move
            up and down the slit. This is necessary for strongly
            curved lines that might be substantially shifted wrt
            the central row.
        """
        off = 0.
        for row in range(indx,height):    # The top of the slit
            if len(ref[row])==0:
                continue
            p = ref[row]+off
            diff = abs(l-p)
            if diff.min()<tol:
                pos = diff.argmin()
                lines.append([ref[row][pos],row,l])
                off = l-ref[row][pos]
        off = 0.
        for row in range(indx-1,-1,-1):    # The bottom of the slit
            if len(ref[row])==0:
                continue
            p = ref[row]+off
            diff = abs(l-p)
            if diff.min()<tol:
                pos = diff.argmin()
                lines.append([ref[row][pos],row,l])
                off = l-ref[row][pos]
    lines = scipy.asarray(lines)
    soln = sf.lsqfit(lines,'chebyshev',xord,yord)

    tmp = lines[:,0].copy()
    lines[:,0] = lines[:,2].copy()
    lines[:,2] = tmp.copy()
    soln2 = sf.lsqfit(lines,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def update_xsoln(data,soln,xord=3,yord=3,skip=False):
    data = data.copy()

    height,width = data.shape

    if skip:
        slice = data.mean(1)
        indx = slice.argsort()[slice.size/4]
        mid = data[indx]
        badrow = scipy.where(slice>slice[indx]+5.*(slice[indx])**0.5)[0]
        if badrow.size==0:
            badrow = scipy.array([-1])
    else:
        indx = height/2
        mid = data[indx]
        badrow = scipy.array([-1])

    xvals = scipy.arange(mid.size)
    middle_lines = get_lines(xvals,mid,False)
    straight = sf.genfunc(middle_lines,indx,soln['back'])
    ref = []
    lines = []

    for row in range(height):
        if abs(row-badrow).min()==0:
            continue
        current =  get_lines(xvals,data[row],False)
        if current.size==0:
            continue
        guess = sf.genfunc(current,row,soln['back'])
        for l in straight:
            diff = abs(l-guess)
            if diff.min()<2:
                pos = diff.argmin()
                lines.append([current[pos],row,l])
    lines = scipy.asarray(lines)
    newsoln = sf.lsqfit(lines,'chebyshev',xord,yord)

    tmp = lines[:,0].copy()
    lines[:,0] = lines[:,2].copy()
    lines[:,2] = tmp.copy()
    del tmp
    newsoln2 = sf.lsqfit(lines,'chebyshev',xord,yord)

    return {'back':newsoln,'forw':newsoln2}


def newModMatch(data,model,scale,order,center,showFit=False):
    import numpy,pylab
    import special_functions as sf
    from scipy import ndimage,interpolate
    x0 = numpy.arange(data.size).astype(numpy.float32)
    datalines = get_lines(x0[20:-20],data[20:-20]).astype(numpy.int32)
    datalines = datalines[data[datalines]>2*numpy.median(data)]
    Xstart = datalines.min()-15
    Xend = datalines.max()+15

    if Xend>data.size/2:
        Xend = data.size/2
    tdata = data[Xstart:Xend].copy()

    wave = numpy.linspace(center-scale*data.size,center,data.size)
    mod = interpolate.splev(wave,model['matched'])
    mod[wave>10300.] = 0.
    mod[wave<10300.] = specContinuumSub(mod[wave<10300.])
    tmp = specContinuumSub(tdata)
    ind = mod*0.
    ind[:tdata.size] = tmp#mod[100:100+tdata.size]
    fftd = numpy.fft.rfft(ndimage.gaussian_filter(ind,9))
    fftm = numpy.fft.rfft(ndimage.gaussian_filter(mod,9)[::-1])
    T = numpy.fft.fftshift(numpy.fft.irfft(fftd*fftm)).real
    Wstart = wave[mod.size/2-T.argmax()]

    if 1==2: # DEBUG
        pylab.plot(wave,mod)
        w0 = wave[wave>Wstart]
        if w0.size>tdata.size:
            pylab.plot(w0[:tdata.size],tdata)
        else:
            pylab.plot(w0,tdata[:w0.size])
        pylab.show()

    Xend = datalines.max()+15
    tdata = data[Xstart:Xend].copy()
    tdata = specContinuumSub(tdata)

    pars = numpy.zeros(order+2)
    pars[0] = Wstart-Xstart*scale
    pars[1] = scale
    pars[-1] = 10.
    cov = numpy.zeros(order+2)
    cov[0] = scale*5.
    cov[1] = scale/50.
    for i in range(2,order+1):
        cov[i] = (2./data.size)**i/10.
    cov[-1] = 0.1

    XOFF = 0.5*(Xstart+Xend)
    w0 = Wstart-Xstart*scale+XOFF*scale
    import pymc
    priors = []
    priors.append(pymc.Uniform('p0',w0-100*scale,w0+100*scale,value=w0))
    priors.append(pymc.TruncatedNormal('p1',scale,1./(scale/100.)**2,scale-scale/10.,scale+scale/10.,value=scale))
    for i in range(2,order+1):
        priors.append(pymc.Uniform('p%d'%(i+2),-10000*(2./data.size)**i,10000*(2./data.size)**i,value=0.))
    priors.append(pymc.Uniform('Norm',0.,1e5,value=100.))

    indata = []
    @pymc.observed
    def optfunc(value=0.,pars=priors):
        if len(indata)==0:
            return -1e200
        x,d,mod,sig = indata
        coeff = numpy.atleast_2d(numpy.array(pars[:-1])).T
        fit = {'coeff':coeff,'type':'polynomial'}
        x0 = sf.genfunc(x,0.,fit)
        try:
            model = interpolate.splev(x0,mod)
        except:
            return -1e200
        model[x0>10390.] = 0.
        model = model*pars[-1]
        resid = (model-d)/sig
        lp = -0.5*(resid**2).sum()
        return lp


    def clip(arr,nsig=3.5):
        a = arr.copy()
        m,s,l = a.mean(),a.std(),a.size
        while 1:
            a = a[abs(a-m)<nsig*s]
            if a.size==l:
                return m,s
            m,s,l = a.mean(),a.std(),a.size

    from SampleOpt import AMAOpt
    Wend = Wstart+(Xend-Xstart)*scale
    if Wend>10000:
        end = int((10000-Wend)/scale)
    else:
        end = -1
    print end
    niter = 2000*order*2
    while 1:
        xvals = x0[Xstart:Xend+end].copy()
        spec = ndimage.gaussian_filter(tdata[:end],9)
        sig = (abs(spec)+clip(tdata[:end])[1]**2)**0.5
        indata = [xvals-XOFF,spec,model['wide'],sig]
        fit = {'coeff':scipy.atleast_2d([p.value for p in priors[:-1]]).T,'type':'polynomial'}
        w = sf.genfunc(indata[0],0.,fit)
        MM = interpolate.splev(w,indata[2])
        norm = indata[1].mean()/MM.mean()
        MM *= norm
        priors[-1].value = norm
        print -0.5*(((indata[1]-MM)/indata[3])**2).sum()
        print w,Wstart
        Sampler = AMAOpt(priors,[optfunc],[],numpy.array(cov))
        Sampler.set_minprop(len(priors)*10)
        Sampler.sample(niter)
        logp,trace,det = Sampler.result()
        coeff = trace[-1].copy()
        print trace[-3:]
        fit = {'coeff':scipy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}
        print coeff
        w = sf.genfunc(indata[0],0.,fit)
        pylab.plot(w,spec)
        pylab.plot(w,interpolate.splev(w,model['wide'])*coeff[-1])
        MM = interpolate.splev(w,indata[2])*coeff[-1]
        print -0.5*(((indata[1]-MM)/indata[3])**2).sum()
        pylab.show()
        niter = 1000
        if w[-1]>10000.:
            end -= 50
        else:
            break

    Xend += end
    tdata = data[Xstart:Xend].copy()
    xvals = x0[Xstart:Xend].copy()
    tdata = specContinuumSub(tdata)

    spec = ndimage.gaussian_filter(tdata,9)
    sig = (abs(spec)+clip(tdata)[1]**2)**0.5
    indata = [xvals-XOFF,spec,model['wide'],sig]
    Sampler = AMAOpt(priors,[optfunc],[],numpy.array(cov))
    Sampler.set_minprop(len(priors)*5)
    Sampler.sample(4000)
    logp,trace,det = Sampler.result()
    coeff = trace[-1].copy()
    print coeff

    fit = {'coeff':scipy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}
    fit = sf.recenter(fit,-XOFF)
    w = sf.genfunc(x0[Xstart:Xend],0.,fit)
    m = interpolate.splev(w,model['wide'])*coeff[-1]
    if showFit==True:
        pylab.plot(w,m)
        pylab.plot(w,spec)
        pylab.show()


    spec = ndimage.gaussian_filter(tdata,5)
    sig = (abs(spec)+clip(tdata)[1]**2)**0.5
    indata = [xvals-XOFF,spec,model['wide'],sig]
    Sampler.sample(1000)

    spec = ndimage.gaussian_filter(tdata,3)
    sig = (abs(spec)+clip(tdata)[1]**2)**0.5
    indata = [xvals-XOFF,spec,model['wide'],sig]
    Sampler.sample(1000)

    spec = tdata.copy()
    sig = (abs(spec)+clip(tdata)[1]**2)**0.5
    indata = [xvals-XOFF,spec,model['matched'],sig]
    Sampler.sample(1000*order**2)
    logp,trace,det = Sampler.result()
    coeff = trace[-1].copy()

    fit = {'coeff':scipy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}
    fit = sf.recenter(fit,-XOFF)
    w = sf.genfunc(x0[Xstart:Xend],0.,fit)
    m = interpolate.splev(w,model['matched'])*coeff[-1]
    if showFit==True:
        pylab.plot(w,m)
        pylab.plot(w,spec)
        pylab.show()
    return fit

def arcwave(sky,arc,arcmodel,skymodel,scale,order):
    from scipy import ndimage,stats,interpolate,optimize
    sky = sky.copy()
    arc = arc.copy()

    sky = scipy.median(sky,0)

    x = scipy.arange(sky.size)
    x_orig = x.copy()

    wave = scipy.arange(3000.,10000.,scale)
    arc_wide = ndimage.gaussian_filter(arc,5)
    m = interpolate.splev(wave,arcmodel['norm'])

    a = arc.copy()
    aw = arc_wide.copy()
    arc = a[:a.size/2.]

    x = x_orig[:a.size/2.]
    arclines = get_lines(x,arc)
    fit = scipy.zeros(3*arclines.size+1)
    index = 1
    for i in range(arclines.size):
        fit[index] = 1.
        fit[index+1] = arclines[i]
        fit[index+2] = 15.*scale
        index += 3
    arc_wide = sf.ngauss(x,fit)
    """
    Do an approximate chi-square between the sky model and the data over a
        range of offsets using a broadened data and sky model.
    """
    max = 0.
    mid = 0

    delta = scale/10.
    s = scipy.arange(scale-delta,scale+delta,delta/10.)
    for stmp in s:
        wtmp = scipy.arange(2000.,10000.,stmp)
        m = interpolate.splev(wtmp,arcmodel['norm'])
        conv = scipy.empty(m.size-arc_wide.size+1)
        for i in range(conv.size):
            tmp = m[i:i+arc_wide.size].copy()
            if tmp.max()<0.1:
                conv[i] = 0.
                continue
            conv[i] = (tmp*arc_wide).sum()
            conv[i] = 1./((tmp-arc_wide)**2).sum()
        curr = conv.max()
        if curr>max:
            mid = conv.argmax()
            scale = stmp
            max = conv.max()
            wave = wtmp.copy()

    """
    Refine the starting wavelength position using the 'true' (ie narrow) model
        of the sky. Searches for a minimum around the minimum found in the
        previous optimization.
    """
    m = interpolate.splev(wave,arcmodel['matched'])
    conv = scipy.empty(m.size-arc.size+1)
    for i in range(conv.size):
        tmp = m[i:i+arc.size].copy()
        ratio = arc.max()/tmp.max()
        if tmp.max()<1.:
            conv[i] = 0.
            continue
        tmp *= ratio
        conv[i] = (tmp*arc).sum()
    pt = conv[mid-50:mid+51].argmax()+mid-50


    initial_pars = [wave[pt],scale]
    for i in range(order+1-len(initial_pars)):
        initial_pars.append(0.)
    modellines = get_lines(wave,m,std=10.)
    modellines = modellines[modellines>wave[pt]]
    modellines = arcmodel['lines']
    modellines = modellines[modellines>wave[pt]]

    fit = {'coeff':scipy.atleast_2d(initial_pars).T,'type':'polynomial'}

    for o in [1,2]:
        w = sf.genfunc(arclines,0.,fit)
        matches = []
        for j in range(w.size):
            diff = abs(w[j]-modellines)
            if diff.min()<5.*scale:
                matches.append([arclines[j],modellines[diff.argmin()]])
        fit = sf.lsqfit(scipy.asarray(matches),'polynomial',o)

    left_matches = [i for i in matches]    
    wmin = sf.genfunc(a.size*0.45,0.,fit)

    arc = a[a.size/2.:].copy()
    x = scipy.arange(arc.size).astype(scipy.float32)+a.size/2.
    arclines = get_lines(x,arc)
    fit = scipy.zeros(3*arclines.size+1)
    index = 1
    for i in range(arclines.size):
        fit[index] = 1.
        fit[index+1] = arclines[i]
        fit[index+2] = 10.*scale
        index += 3
    arc_wide = sf.ngauss(x,fit)
    """
    Do an approximate chi-square between the sky model and the data over a
        range of offsets using a broadened data and sky model.
    """
    max = 0.
    mid = 0
    delta = scale/10.
    s = scipy.arange(scale-delta,scale+delta,delta/10.)
    for stmp in s:
        wtmp = scipy.arange(wmin,10000.,stmp)
        m = interpolate.splev(wtmp,arcmodel['norm'])
        conv = scipy.empty(m.size-arc_wide.size+1)
        for i in range(conv.size):
            tmp = m[i:i+arc_wide.size].copy()
            if tmp.max()<0.1:
                conv[i] = 0.
                continue
            conv[i] = (tmp*arc_wide).sum()
        curr = conv.max()
        if curr>max:
            mid = conv.argmax()
            scale = stmp
            max = conv.max()
            wave = wtmp.copy()
    """
    Refine the starting wavelength position using the 'true' (ie narrow) model
        of the sky. Searches for a minimum around the minimum found in the
        previous optimization.
    """
    m = interpolate.splev(wave,arcmodel['matched'])
    conv = scipy.empty(m.size-arc.size+1)
    for i in range(conv.size):
        tmp = m[i:i+arc.size].copy()
        ratio = arc.max()/tmp.max()
        if tmp.max()<1.:
            conv[i] = 0.
            continue
        tmp *= ratio
        conv[i] = (tmp*arc).sum()
    pt = conv[mid-50:mid+51].argmax()+mid-50
    wavept = wave[pt]

    initial_pars = [wavept,scale]
    for i in range(order+1-len(initial_pars)):
        initial_pars.append(0.)
    modellines = get_lines(wave,m,std=10.)
    modellines = modellines[modellines>wavept]
    modellines = arcmodel['lines']
    modellines = modellines[modellines>wavept]

    fit = {'coeff':scipy.atleast_2d(initial_pars).T,'type':'polynomial'}
    for o in [1,2]:
        # The (o-2) bit is to correct the offset after the first loop
        w = sf.genfunc(arclines+(o-2)*a.size/2.,0.,fit)
        matches = []
        for j in range(w.size):
            diff = abs(w[j]-modellines)
            if diff.min()<5.*scale:
                matches.append([arclines[j],modellines[diff.argmin()]])
        fit = sf.lsqfit(scipy.asarray(matches),'polynomial',o)

    arc = a.copy()
    arc_wide = aw.copy()

    w = sf.genfunc(arclines,0.,fit)
    for i in range(w.size):
        diff = abs(w[i]-modellines)
        if diff.min()<5.*scale:
            left_matches.append([arclines[i],modellines[diff.argmin()]])

    fit = sf.lsqfit(scipy.asarray(left_matches),'polynomial',order)


    """ Optimization function for refining the wavelength solution. """
    def dofit(p,x,data,model):
        fit = {'coeff':scipy.atleast_2d(p).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        return (m-data)

    x = scipy.arange(arc.size).astype(scipy.float32)

    initial_pars = fit['coeff'][:,0].tolist()
    coeff,ier = optimize.leastsq(dofit,initial_pars,
                        (x,arc_wide,arcmodel['wide']),maxfev=100000)
    coeff,ier = optimize.leastsq(dofit,coeff,
                        (x,arc,arcmodel['matched']),maxfev=100000)
    outcoeff = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}

    def skycorrect(p,arc,sky,arcmodel,skymodel):
        fit = {'coeff':scipy.atleast_2d(p[:-1]).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        arcm = interpolate.splev(w+p[-1],arcmodel)
        chi_arc = (arcm-arc)
        s = sky[w>5100.]
        skym = interpolate.splev(w[w>5100.],skymodel)
        skym *= scipy.median(s/skym)
        chi_sky = 5.*(skym-s)#/abs(m)**0.5
        chi = scipy.concatenate((chi_arc,chi_sky))
        return chi

    newcoeff = coeff.tolist()
    newcoeff.append(0.)
    coeff,ier = optimize.leastsq(skycorrect,newcoeff,
                        (arc,sky,arcmodel['matched'],skymodel['matched']),
                        maxfev=100000)
    outcoeff = {'coeff':scipy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}

    """
    wave = sf.genfunc(x,0.,outcoeff)
    sky = sky[wave>5000.]
    wave = wave[wave>5000.]

    m = interpolate.splev(wave,wavemodel['matched'])
    ratio = scipy.median(sky/m)
    import pylab
    pylab.plot(wave,sky)
    pylab.plot(wave,m*ratio)
    pylab.show()

    offset,ier = optimize.leastsq(skycorrect,[0.],
                        (wave,sky,wavemodel['matched']),maxfev=100000)
    print offset
    outcoeff['coeff'][0] += offset
    """

    return outcoeff


def wave_skylines(sky,solution):
    STD_LINES = [5197.928,5200.286,5202.977,
                 5460.735,5577.345,5867.5522,5915.308,5932.864,6257.970,
                 6300.320,6363.810,
                 6533.040,6553.610,6863.971,6912.620,6923.210,6939.520,
                 7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,
                 7821.510,7841.266,7993.332,8310.719,8344.613,8399.160,8415.231,
                 8430.170,8791.186,8885.830,8943.395,8988.384,9038.059,9337.854,
                 9375.977,9419.746,9439.670,9458.524]
    x = scipy.arange(sky.size).astype(scipy.float32)
    lines = get_lines(x,sky)

    w = sf.genfunc(x,0.,solution)
    scale = sf.lsqfit(w,'polynomial',1)['coeff'][1]
    order = solution['coeff'].size - 1

    if scale>1.5:
        STD_LINES.insert(6,5891.)
    elif scale<0.95:
        STD_LINES.insert(6,5895.93)
        STD_LINES.insert(6,5889.96)
    w = sf.genfunc(lines,0.,solution)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    lines = get_lines(x,sky,nstd=7.)
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    return fit


def wave_arclines(arc,arcmodel,sky,solution,offset=None,neworder=None):
    from scipy import interpolate
    STD_LINES = scipy.sort(arcmodel['lines'])
    SKYCORR = [5577.338,5867.55,6235.6,6257.9,6300.304,6363.78,6553.617]
    SKYLINES = [5577.338,6300.304,6363.78,6498.736]#,6553.617]#,6912.62]
    if offset is None:
        offset = 0.
    x = scipy.arange(arc.size).astype(scipy.float32)
    lines = get_lines(x,arc)
    skycoords = get_lines(x,sky)

    arccoords = lines.copy()
    scale = arcmodel['scale']

    if scale>1.5:
        SKYLINES.insert(1,5891.)
    elif scale<0.6:
        SKYLINES.insert(1,5895.9)

    order = solution['coeff'].size - 1
    fit = solution
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for line in STD_LINES:
        diff = abs(w-line)
        if diff.min()<5.*scale:
            pos = diff.argmin()
            matches.append([lines[pos],line])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    if skycoords.size==0:
        print "Warning, no sky lines found; returning fit to arcs only"
        return fit

    if neworder is not None and neworder>order:
        order = neworder

    for i in range(7):
        fit['coeff'][0] += offset
        matched = [m for m in matches]
        if skycoords.size==0:
            offset = 0.
            break
        skyline = sf.genfunc(skycoords,0.,fit)
#        print skyline
        offsets = []
        for line in SKYCORR[:i*2+2]:
            diff = abs(line-skyline)
            if diff.min()<10.*scale:
                offsets.append(line-skyline[diff.argmin()])
        if len(offsets)==0:
            offset = 0.
            break
        offset = scipy.median(scipy.asarray(offsets))+offset
#        print offset
        for line in SKYCORR[:i*2+2]:
            diff = abs(line-skyline)
            if diff.min()<5.*scale:
                pos = diff.argmin()
                matched.append([skycoords[pos],line-offset])
        fit = sf.lsqfit(scipy.asarray(matched),'polynomial',order)

    def opt(p,x,w,weight,cond,getFit=False):
        if abs(p[2]-scale)/scale>0.05 or abs(p[0])>15:
            return w/weight
        fit = {'coeff':numpy.atleast_2d(p[1:]).T,'type':'polynomial'}
        m = sf.genfunc(x,0.,fit)
        m[cond] = sf.genfunc(x[cond]+p[0],0.,fit)
        res = (m-w)/weight
        if getFit==True:
            return fit
        return res

    fit['coeff'][0] += offset
    sm = []
    from scipy import optimize
    from numpy import linalg
    off = -1*offset/scale#fit['coeff'][1]
    off = numpy.array([off]+fit['coeff'].ravel().tolist())
    for S in SKYLINES:
        skyline = sf.genfunc(skycoords,0.,fit)
        diff = abs(S-skyline)
        if diff.min()<10.*scale:
            pos = diff.argmin()
            sm.append([skycoords[pos],S])
        else:
            break
        x0 = numpy.array([m[0] for m in matches]+[m[0] for m in sm])
        w = numpy.array([m[1] for m in matches]+[m[1] for m in sm])
        weight = numpy.array([scale/3. for m in matches]+[scale/10. for m in sm])
        cond = numpy.arange(x0.size)<len(matches)
        x1 = x0.mean()
        x0 -= x1
        A = numpy.zeros((x0.size,order+2))
        A[:len(matches),0] = 1.
        A[len(matches):,1] = 1.
        for i in range(1,order+1):
            A[:,i+1] = x0**i
        coeff = linalg.lstsq(A,w)[0]
        coeff[0] = off[0]
        off,ier = optimize.leastsq(opt,coeff,(x0,w,weight,cond),epsfcn=1e-5)
        fit = opt(off,x0,w,weight,cond,True)
        fit = sf.recenter(fit,-x1)

    print "Shifting arc by %4.2f pixels"%(off[0])
    return fit


def wave_arcsky(arc,arcmodel,sky,solution):
    """
    First find the best solution with the skylines, then apply this solution
        to all arclines within the bounds of the lowest/highest wavelength
        skylines, solving for the delta_pixel offset between the sky and the
        arc. Then find the solution for all lines (sky and delta_pixel-offset
        arcs).
    """
    def clip(arr):
        a = arr.copy()
        m,s,l = a.mean(),a.std(),a.size
        while 1:
            a = a[abs(a-m)<3.*s]
            if a.size==l:
                return m,s
            m,s,l = a.mean(),a.std(),a.size

    STD_LINES = [5197.928,5200.286,5202.977,
                 5460.735,5577.345,5867.5522,5915.308,5932.864,6257.970,
                 6300.320,6363.810,
                 6533.040,6553.610,6863.971,6912.620,6923.210,6939.520,
                 7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,
                 7821.510,7841.266,7993.332,8310.719,8344.613,8399.160,8415.231,
                 8430.170,8791.186,8885.830,8943.395,8988.384,9038.059,9337.854,
                 9375.977,9419.746,9439.670,9458.524]
    x = scipy.arange(sky.size).astype(scipy.float32)
    lines = get_lines(x,sky)

    global fit
    scale = solution['coeff'][1]
    order = solution['coeff'].size - 1

    if scale>1.5:
        STD_LINES.insert(6,5891.)

    w = sf.genfunc(lines,0.,solution)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    lines = get_lines(x,sky,nstd=7.)
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    matches = scipy.asarray(matches)
    fit = sf.lsqfit(matches,'polynomial',order)
    revfit = sf.lsqfit(matches[:,::-1],'polynomial',order)

    ARCS = scipy.sort(arcmodel['lines'])
    alines = get_lines(x,arc)
    xmin,xmax = matches[0,0],matches[-1,0]
    arc_x = sf.genfunc(ARCS,0.,revfit)
    offset = []
    for i in range(arc_x.size):
        if arc_x[i]<xmin-2 or arc_x[i]>xmax+2:
            continue
        diff = arc_x[i]-alines
        if abs(diff).min()<9.:
            offset.append(diff[abs(diff).argmin()])
    offset = scipy.asarray(offset)
    off,width = clip(offset)

    aw = sf.genfunc(alines+off,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    k = len(matches)
    for i in range(aw.size):
        diff = abs(aw[i]-ARCS)
        if diff.min()<3.*scale:
            matches.append([alines[i]+off,ARCS[diff.argmin()]])
    matches = scipy.asarray(matches)
    fit = sf.lsqfit(matches,'polynomial',order)
    return fit


def get_lines(xvals,data,dosat=True,std=None,nstd=15.):
    from scipy import ndimage
    def clip(arr,nsig):
        a = arr.copy()
        m,std,size = a.mean(),a.std(),a.size

        while 1:
            flag = a*0.
            cond = abs(a-m)>nsig*std
            flag[cond] = 1.
            flag = ndimage.maximum_filter(flag,5)
            a = a[flag==0]
            if a.size==size or a.size<2:
                return m,std
            m,std,size = a.mean(),a.std(),a.size

    bgsub = specContinuumSub(data)
    if std is None:
        m,std = clip(bgsub,3.)
    maxfilt = ndimage.maximum_filter(bgsub,9)
    peaks = numpy.where((maxfilt==bgsub)&(bgsub>nstd*std))[0]

    out = []
    fit = scipy.empty(4)
    for p in peaks:
        if dosat:
            if p-14<0 or p+15>xvals.size:
                continue
        else:
            if p-4<0 or p+5>xvals.size:
                continue

        if data[p]<54000.:
            fitdata = scipy.empty((9,2))
            fitdata[:,0] = xvals[p-4:p+5].copy()
            fitdata[:,1] = bgsub[p-4:p+5].copy()
        elif dosat:
            fitdata = scipy.empty((25,2))
            fitdata[:,0] = xvals[p-12:p+13].copy()
            fitdata[:,1] = bgsub[p-12:p+13].copy()
            fitdata = fitdata[fitdata[:,1]<54000.]
        else:
            continue

        answer = (fitdata[:,0]*fitdata[:,1]).sum()/fitdata[:,1].sum()
        out.append(answer)

    return scipy.asarray(out)


def specContinuumSub(spec,filtp=10.,filts=51):
    from scipy import ndimage
    cont = ndimage.percentile_filter(spec,filtp,filts)
    fft = numpy.fft.rfft(cont)
    k = numpy.ones(fft.size)
    k[:10] = 0.
    cont -= numpy.fft.irfft(fft*ndimage.gaussian_filter(k,5),n=cont.size)
    return spec-cont


def spexCorr(s1,s2,contsub=False):
    if contsub:
        fft1 = numpy.fft.rfft(specContinuumSub(s1))
        fft2 = numpy.fft.rfft(specContinuumSub(s2)[::-1])
    else:
        fft1 = numpy.fft.rfft(s1)
        fft2 = numpy.fft.rfft(s2[::-1])
    return numpy.fft.fftshift(numpy.fft.irfft(fft1*fft2)).real


def combine_xw(coords,xsoln,wsoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newx = sf.genfunc(x,y,xsoln['back'])
    wave = sf.genfunc(newx,0.,wsoln)

    data = scipy.empty((x.size,3))
    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = wave.copy()
    soln = sf.lsqfit(data,'chebyshev',xord,yord)

    tmp = data[:,0].copy()
    data[:,0] = data[:,2].copy()
    data[:,2] = tmp.copy()
    soln2 = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def combine_xy(offset,coords,xsoln,ysoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newy = ysoln['back'].flatten()
    newx = sf.genfunc(x,y-offset,xsoln['back'])

    data = scipy.empty((x.size,3))
    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = newx.copy()
    soln = sf.lsqfit(data,'chebyshev',xord,yord)

    tmp = data[:,0].copy()
    data[:,0] = data[:,2].copy()
    data[:,2] = tmp.copy()
    soln2 = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def combine_xyw(coords,xsoln,ysoln,wsoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newy = ysoln.flatten()
    from scipy import random
    k = random.random(x.size)
    args = k.argsort()
    x = x[args[:y.size/10]]
    newy = newy[args[:y.size/10]]
    y = y[args[:y.size/10]]

    newx = sf.genfunc(x,y,xsoln['back'])
    wave = sf.genfunc(newx,0.,wsoln)

    data = scipy.empty((x.size,3))
    data[:,0] = wave.copy()
    data[:,1] = y.copy()
    data[:,2] = x.copy()
    output_to_ccdx = sf.lsqfit(data,'chebyshev',xord,yord)

    data[:,2] = newy.copy()
    output_to_ccdy = sf.lsqfit(data,'chebyshev',xord,yord)

    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = wave.copy()
    ccdx_ycor_to_wave = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'sky2x':output_to_ccdx,'sky2y':output_to_ccdy,'ccd2wave':ccdx_ycor_to_wave}

