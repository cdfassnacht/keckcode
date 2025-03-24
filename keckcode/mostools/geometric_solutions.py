import scipy
import special_functions as sf


def new_xsoln(data,xord=3,yord=3,tol=2.):
    data = data.copy()

    height,width = data.shape

    indx = height/2
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


def skywave(data,skymodel,scale,order=3,cutoff=None,hicutoff=9100.,minwave=None):
    from scipy import ndimage,stats,interpolate,signal,optimize
    import pylab

    data = data.copy()

    if data.ndim==2:
        sky = stats.stats.nanmedian(data,0)
    else:
        sky = data.copy()
    sky[scipy.isnan(sky)] = 0.
    osky = sky.copy()
    sky = ndimage.gaussian_filter(sky,5.)

    if cutoff is None:
        outw = scipy.arange(3000.,9000.,1.)
    else:
        outw = scipy.arange(cutoff-0.5*osky.size*scale,9000.,1.)
    if minwave is not None:
        outw = outw[outw>minwave]

    out = outw*0.
    outc = out*0.

    x = scipy.arange(osky.size)
    l = get_lines(x,sky,nstd=15.)

    l = l[sky[l.astype(scipy.int16)]>scipy.median(sky)]
    offset = int(l.min()-50.)
    if offset<0:
        offset = 0
    sky = sky[offset:]

    print offset
    width = 250.
    if cutoff is None:
        w0 = 5150.-width
    else:
        w0 = cutoff-width

    while w0+width<hicutoff:
        w0 += width
        if cutoff is not None and w0<cutoff:
            continue
        w1 = scipy.arange(w0,w0+width,scale)
        t1 = interpolate.splev(w1,skymodel['wide'])
        corr = signal.correlate(t1,sky,'valid')

        chi = corr*0.
        for j in range(corr.size):
            med = scipy.median(sky[j:j+t1.size]/t1)
            chi[j] = ((t1*med-sky[j:j+t1.size])**2/sky[j:j+t1.size]).sum()

        ratio = corr.copy()/chi.copy()
        ratio /= ratio.mean()
        ratio = ndimage.gaussian_filter(ratio,5)
        mywave = w0-scipy.arange(ratio.size)*scale

        mywave = mywave[::-1]
        ratio = ratio[::-1]
        cond = (outw>=mywave[0])&(outw<=mywave[-1])

        mod = interpolate.splrep(mywave,ratio,s=0)
        out[cond] += interpolate.splev(outw[cond],mod)
        outc[cond] += 1.

    out[outc>0] /= outc[outc>0]
    start = outw[out.argmax()]-offset*scale
    print start
    sky = osky[offset:].copy()
    sky_wide = ndimage.gaussian_filter(sky,5)


    cutoff = 5500.
    """ Optimization function for refining the wavelength solution. """
    def dofit(p,x,data,model):
        if scipy.isnan(p).any():
            return x*0.+1e7
        fit = {'coeff':scipy.atleast_2d(p[1:]).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        m *= p[0]
        chi = (m-data)/abs(data)**0.5
        cond = ~scipy.isnan(chi)
        cond = cond&scipy.isfinite(chi)
        cond = cond&(w>cutoff)&(w<10400.)
        return chi[cond]/chi[cond].size

    m = interpolate.splev(scipy.arange(start+offset*scale,10400.,scale),skymodel['matched'])

    x = scipy.arange(sky.size).astype(scipy.float32)+offset
    initial_pars = [scipy.median(sky/m[:sky.size]),start,scale]
    diag = [1.,1e3,1.]
    for i in range(order+2-len(initial_pars)):
#        initial_pars.append(1e-4**(i+2))
        initial_pars.append(0.)
        diag.append(1e3**(i+2))

    cond = sky>1e-5
    sky = sky[cond]
    sky_wide = sky_wide[cond]
    x = x[cond]

    print initial_pars
    coeff,ier = optimize.leastsq(dofit,initial_pars,
                        (x,sky_wide,skymodel['wide']),maxfev=1e5,epsfcn=1e-15,diag=diag,xtol=1e-15,ftol=1e-15)
    coeff,ier = optimize.leastsq(dofit,coeff,
                        (x,sky,skymodel['matched']),maxfev=1e5,epsfcn=1e-15,diag=diag,xtol=1e-15,ftol=1e-15)


    return {'coeff':scipy.atleast_2d(coeff[1:]).T,'type':'polynomial'}


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

    import pylab
    pylab.plot(conv)
    pylab.show()

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
    scale = solution['coeff'][1]
    order = solution['coeff'].size - 1

    STD_LINES = [5460.735,5577.345,5915.308,5932.864,6257.970,6300.320,6363.810,
                 6533.040,6553.610,6863.971,6912.620,6923.210,6939.520,
                 7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,
                 7821.510,7841.266,7993.332,8310.719,8344.613,8399.160,8415.231,
                 8430.170,8791.186,8885.830,8943.395,8988.384,9038.059,9337.854,
                 9375.977,9419.746,9439.670,9458.524]
    if scale<1.5:
        STD_LINES.append(5889.959)
        STD_LINES.append(5895.932)
        STD_LINES.sort()

    x = scipy.arange(sky.size).astype(scipy.float32)
    lines = get_lines(x,sky,10.)

    w = sf.genfunc(lines,0.,solution)
    print w
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

    lines = get_lines(x,sky,nstd=5.)
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    return fit


def wave_arclines(arc,arcmodel,sky,solution):
    from scipy import interpolate
    STD_LINES = scipy.sort(arcmodel['lines'])
    SKYLINES = [5577.338,6300.304,6363.78,6553.617,6912.62]
    x = scipy.arange(arc.size).astype(scipy.float32)
    lines = get_lines(x,arc)

    scale = solution['coeff'][1]
    order = solution['coeff'].size - 1
    w = sf.genfunc(lines,0.,solution)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    sky = get_lines(x,sky)
    sky = sf.genfunc(sky,0.,fit)
    offsets = []
    for line in SKYLINES:
        diff = abs(line-sky)
        if diff.min()<5.*scale:
            offsets.append(line-sky[diff.argmin()])
    if len(offsets)==0:
        offset = 0.
    else:
        offset = scipy.median(scipy.asarray(offsets))

    fit['coeff'][0] += offset

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
            if a.size==size or a.size==0:
                return m,std
            m,std,size = a.mean(),a.std(),a.size

    bg = ndimage.percentile_filter(data,10,51)

    if std is None:
        m,std = clip(data-bg,3.)
    max = ndimage.maximum_filter(data,9)
    peaks = scipy.where((max==data)&((data-bg)>nstd*std))[0]

    out = []
    fit = scipy.empty(4)
    for p in peaks:
        if p-14<0 or p+15>xvals.size:
            continue
        fit[0] = bg[p]
        fit[1] = data[p]-fit[0]
        fit[2] = xvals[p]
        fit[3] = 1.

        if data[p]<54000.:
            fitdata = scipy.empty((9,2))
            fitdata[:,0] = xvals[p-4:p+5].copy()
            fitdata[:,1] = data[p-4:p+5].copy()-bg[p-4:p+5]
        elif dosat:
            fitdata = scipy.empty((25,2))
            fitdata[:,0] = xvals[p-12:p+13].copy()
            fitdata[:,1] = data[p-12:p+13].copy()-bg[p-12:p+13]
            fitdata = fitdata[fitdata[:,1]<54000.]
        else:
            continue

        answer = (fitdata[:,0]*fitdata[:,1]).sum()/fitdata[:,1].sum()
        out.append(answer)

    return scipy.asarray(out)


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


def combine_xyw(coords,xsoln,ysoln,wsoln,xord,yord,nskip=10):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newy = ysoln.flatten()
    from scipy import random
    k = random.random(x.size)
    args = k.argsort()
    x = x[args[:y.size/nskip]]
    newy = newy[args[:y.size/nskip]]
    y = y[args[:y.size/nskip]]

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


def arcwave2(arc,arcmodel,scale,order=3,bcutoff=2e3,rcutoff=1e4):
    from scipy import ndimage,stats,interpolate,signal,optimize
    import pylab

    arc = arc.copy()
    arc[scipy.isnan(arc)] = 1.
    arc[arc<=1.] = arc.mean()

    wave,model = arcmodel['orig']
    model[scipy.isnan(model)] = 0.
    wave = wave.copy()
    model = model.copy()
    cond = (wave>bcutoff)&(wave<rcutoff)
    corrmodel = model.copy()
    corrmodel[(wave<bcutoff)|(wave>rcutoff)] = 0.

    corr = signal.correlate(arc,corrmodel,mode='valid')
    offset = corr.argmax()

    lines = arcmodel['lines'].copy()
    bc = lines.min()-scale*20.
    rc = lines.max()+scale*20.

    x = scipy.arange(arc[offset:].size)
    w = wave[offset:offset+x.size].copy()
    cond = (w>bc)&(w<rc)
    fit = scipy.empty((x[cond].size,2))
    fit[:,0] = x[cond].copy()
    fit[:,1] = w[cond].copy()
    pars = sf.lsqfit(fit,'polynomial',1)

    pars['coeff'][0] = wave[offset]
    pars['coeff'][1] = scale
#    pars = [wave[offset],scale]
#    for i in range(2,order+1):
#        pars.append(1e-5**i)
    pylab.plot(wave,model)
    w = sf.genfunc(scipy.arange(x.size),0.,pars)
    pylab.plot(w,interpolate.splev(w,arcmodel['matched']))
    pylab.show()
    print sf.genfunc(x[cond],0.,sf.lsqfit(fit,'polynomial',1))
    print w[cond]
    print x[cond]
    print pars

    def arcfit(p,x,arc,mod):
        fit = {'coeff':scipy.atleast_2d(p).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        cond = (w>bcutoff)&(w<rcutoff)
        m = interpolate.splev(w[cond],mod)
        chi = (m-arc[cond])/abs(arc[cond])**0.5
        return chi

    widearc = ndimage.gaussian_filter(arc,7.)
    x = scipy.arange(arc[offset:].size)
    coeff,ier = optimize.leastsq(arcfit,pars,(x,widearc[offset:].copy(),
                        arcmodel['wide']),maxfev=100000)

    print coeff
    fit = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}


    x = scipy.arange(arc.size)
    l = get_lines(x,arc,nstd=15.)
    lw = sf.genfunc(l-offset,0.,fit)
    lines = []
    for i in range(l.size):
        diff = abs(lw[i]-arcmodel['lines'])
        if diff.min()>5.*scale:
            continue
        lines.append([l[i],arcmodel['lines'][diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(lines),'polynomial',order)

    pars = fit['coeff'].flatten()
    coeff,ier = optimize.leastsq(arcfit,pars,(x[offset:],arc[offset:],
                        arcmodel['matched']),maxfev=100000)

    fit = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}

    return fit
