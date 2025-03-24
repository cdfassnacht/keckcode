import scipy,pickle,numpy
from scipy import signal,ndimage,interpolate,optimize,stats
import special_functions as sf
import pylab
from mostools import specfit as spf


def make_standard(arc,sky,showFit=True,usearc=False):
    nspec = len(arc)

    # Create smoothed transforms
    oarc = []
    midarc = []
    widearc = []
    osky = []
    lines = []
    x = numpy.arange(arc[0].size)*1.
    for i in range(nspec):
        o = arc[i].copy()*0
        o[5:-5] = spf.specContinuumSub(arc[i][5:-5].copy())
        oarc.append(o)
        midarc.append(ndimage.gaussian_filter(o,5))
        widearc.append(ndimage.gaussian_filter(o,11))
        l = spf.get_lines(x,arc[i])
        coeff = numpy.ones(l.size*3+1)*3.
        coeff[0] = 0.
        coeff[2::3] = l
        lines.append(sf.ngauss(x,coeff))

    # Get order
    sols = numpy.empty((nspec,nspec))
    for i in range(nspec):
        order = []
        for j in range(nspec):
            T = spf.spexCorr(lines[i],lines[j])
            order.append(T.size/2-T.argmax())
        sols[i] = numpy.array(order).argsort()

    order = numpy.median(sols,0).astype(numpy.int32)

    # Get offsets from cross correlation
    offset = [0]
    for i in range(1,nspec):
        corr = spf.spexCorr(lines[order[i-1]],lines[order[i]])
        offset.append(corr.size/2-corr[:corr.size/2+10].argmax())

    # Reorder offsets
    offset = numpy.array(offset)
    offset[0:-1] = offset[1:]
    offset = offset[::-1]
    offset[0] = 0
    order = numpy.array(order)[::-1]

    # Prepare output arrays
    specsize = oarc[0].size
    outx = numpy.arange(specsize*2)
    out = numpy.empty((nspec,outx.size))*numpy.nan
    if usearc==True:
        out[0][:specsize] = oarc[order[0]].copy()
        m,s = clip(oarc[order[0]])
        out[0] /= numpy.median(oarc[order[0]][oarc[order[0]]>m+s])
    else:
        out[0][:specsize] = sky[order[0]].copy()
        out[0] /= stats.stats.nanmedian(out[0])

    x = numpy.arange(specsize).astype(numpy.float32)
    xsoln = [x.copy()]
    print "Matching slits"
    def dofit(p,x,data,model,off,var=0.):
        c = (p[1:]-1.)/rescale[1:]
        if numpy.isnan(p).any():
            return data
        if abs(c[1]-1.)>0.05:
            return data
        if abs(c[0]-off)>30:
            return data
        fit = {'coeff':numpy.atleast_2d(c).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        try:
            m = interpolate.splev(w,model)
        except:
            return data
        m *= p[0]
        chi = (m-data)/abs(data+m+var)**0.5
        cond = ~numpy.isnan(chi)
        cond = cond&numpy.isfinite(chi)
        badval = abs(chi[cond]).max()
        chi[~cond] = 2*badval
        return chi

    def show(p,x,d,mod):
        if type(p)==type({}):
            fit = p
        else:
            c = (p[1:]-1.)/rescale[1:]
            print c
            fit = {'coeff':scipy.atleast_2d(c).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,mod)
        if type(p)!=type({}):
            m *= p[0]
        pylab.plot(x,d)
        pylab.plot(x,m)
        pylab.show()

    porder = 3
    for i in range(1,order.size):
        off = offset[i]
        if off==0:
            off = 1.
        rescale = [1.,1./off,1.]
        for j in range(porder+2-len(rescale)):
            rescale.append((oarc[i].size/2.)**(j+2))
        rescale = numpy.asarray(rescale)
        diag = numpy.ones(rescale.size)
        initial_pars = numpy.ones(diag.size)
        initial_pars[0] = 1.
        initial_pars[1:3] += 1.
        coeff = numpy.asarray(initial_pars)

        wmod = interpolate.splrep(x[5:-5],widearc[order[i-1]][5:-5])
        mmod = interpolate.splrep(x[5:-5],midarc[order[i-1]][5:-5])
        mod = interpolate.splrep(x[5:-5],oarc[order[i-1]][5:-5])
        wide = widearc[order[i]][5:-off-15]
        mid = midarc[order[i]][5:-off-15]
        spec = oarc[order[i]][5:-off-15]
        var = clip(spec)[1]**2
        x0 = x[5:-off-15]
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,wide,wmod,off,var),
                        maxfev=100000,epsfcn=1e-3,diag=diag,ftol=1e-15,factor=99.)
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,wide,wmod,off,var),
                        maxfev=100000,epsfcn=1e-7,diag=diag,ftol=1e-15,factor=99.)
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,mid,mmod,off,var),
                        maxfev=100000,epsfcn=1e-5,diag=diag,ftol=1e-15,factor=99.)
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,mid,mmod,off,var),
                        maxfev=100000,epsfcn=1e-7,diag=diag,ftol=1e-15,factor=99.)
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,spec,mod,off,var),
                        maxfev=100000,epsfcn=1e-5,diag=diag,ftol=1e-15,factor=99.)
        coeff,ier = optimize.leastsq(dofit,coeff,(x0,spec,mod,off,var),
                        maxfev=100000,epsfcn=1e-9,diag=diag,ftol=1e-15,factor=99.)
        #show(coeff,x0,spec,mod)

        coeff = numpy.atleast_2d(numpy.array((coeff[1:]-1.)/rescale[1:])).T
        fit = {'coeff':coeff,'type':'polynomial'}
        xtmp = sf.genfunc(x,0.,fit)
        lines = spf.get_lines(x[5:-5],arc[order[i]][5:-5])
        LINES = spf.get_lines(x[5:-5],arc[order[i-1]][5:-5])
        xlines = sf.genfunc(lines,0.,fit)
        matches = []
        for j in range(xlines.size):
            diff = abs(LINES-xlines[j])
            if diff.min()<5.:
                matches.append([lines[j],LINES[diff.argmin()]])
        fit = sf.lsqfit(numpy.asarray(matches),'polynomial',porder)
        #show(fit,x0,spec,mod)

        x0 = sf.genfunc(x[5:-5],0.,fit)
        cond = (outx>x0[0])&(outx<x0[-1])
        if usearc==True:
            mod = interpolate.splrep(x0,oarc[order[i]][5:-5],s=0)
        else:
            mod = interpolate.splrep(x0,sky[order[i]][5:-5],s=0)
        out[i][cond] = interpolate.splev(outx[cond],mod)
        if usearc==True:
            m,s = clip(out[i][cond])
            out[i] /= numpy.median(out[i][cond][out[i][cond]>m+s])
        else:
            out[i] /= stats.stats.nanmedian(out[i])
#        mod = interpolate.splrep(x0,arc[order[i]][5:-5],s=0)
#        out2[i][cond] = interpolate.splev(outx[cond],mod)
#        out2[i] /= stats.stats.nanmedian(out2[i])

        #pylab.plot(out[i-1][cond]+1)
        #pylab.plot(out[i][cond])
        #pylab.show()
        x0 = sf.genfunc(x,0.,fit)
        x = x0.copy()
        xsoln.append(x0.copy())

    ""
    if out.shape[0]>1 and showFit==True:
        pylab.figure()
        for i in range(out.shape[0]):
            pylab.plot(outx,out[i])
        pylab.show()
    ""

    out = stats.stats.nanmedian(out,0)
    good = scipy.where(numpy.isfinite(out))[0]
    offset = good.min()

    spec = out[good.min():good.max()+1]
    outx = outx[good.min():good.max()+1]
    spec[~numpy.isfinite(spec)] = 0.

    if len(order)==1:
        spec = arc[0]
        outx = x

    return outx,spec,order,xsoln


def clip(arr,nsig=3.5):
    a = arr.copy()
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<nsig*s]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size

def lineMatch(line1,oline2,line2,tol,order):
    matches = []
    for line in line1:
        diff = abs(line-line2)
        if diff.min()<tol:
            matches.append([oline2[diff.argmin()],line])
    if len(matches)<order+1:
        return 0
    matches = numpy.asarray(matches)
    return sf.lsqfit(matches,'polynomial',order)
