import scipy,pickle,numpy
from scipy import signal,ndimage,interpolate,optimize,stats
import special_functions as sf
import pylab
from mostools import specfit as spf


def getWaveSolution(pars,model,centralwave,lines,isBlue=False,pre=None,showFit=False):
    outx,spec,order,xsoln = pars

    if pre is None:
        if isBlue==True:
            from LRIS import id_spec
            wtmp = numpy.arange(model['blue'],model['red'],model['scale']/2.)
            mtmp = interpolate.splev(wtmp,model['matched'])
            model['lines'] = spf.get_lines(wtmp,mtmp)
            model,soln = id_spec.id_spec(spec,model)
        else:
#            scale = model['scale']
#            soln = spf.newModMatch(spec,model,scale,3,centralwave)
            from LRIS import red_id_spec
            soln = red_id_spec.id_spec(spec,model)
    else:
        obs = numpy.load('%s_wave.calib'%pre)
        soln,model = obs.wsoln

    x = scipy.arange(spec.size)
    wave = sf.genfunc(x,0.,soln)
    wavemodel = interpolate.splrep(outx,wave,k=3,s=0)
    x0 = xsoln[0].copy()
    fitdata = scipy.empty((x0.size,2))
    fitdata[:,0] = x0.copy()
    solution = [None for i in range(order.size)]

    def refineFit(p,x,d,model):
        c = (p[1:]-1.)*rescale
        if scipy.isnan(p).any():
            return data
        fit = {'coeff':scipy.atleast_2d(c).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        m *= p[0]
        chi = (m-data)/abs(data)**0.5
        cond = ~scipy.isnan(chi)
        cond = cond&scipy.isfinite(chi)
        cond = cond&(w<10350.)
        badval = abs(chi[cond]).max()
        chi[~cond] = 2*badval
        return chi
    x0 = numpy.arange(lines[0].size).astype(numpy.float32)
    for i in range(order.size):
        x = xsoln[i]
        data = lines[order[i]]
        w = interpolate.splev(x,wavemodel)
        fitdata[:,1] = w.copy()
        fit = sf.lsqfit(fitdata,'polynomial',3)
        if order.size==1:
            solution = [fit]
            break
        w = sf.genfunc(x0,0.,fit)
        m = interpolate.splev(w,model['matched'])
        data /= numpy.median(data/m)

        rescale = fit['coeff'].flatten()
        coeff = numpy.ones(rescale.size+1)
        coeff[1:] += 1.
        coeff,ier = optimize.leastsq(refineFit,coeff,(x0,data,model['matched']),
                epsfcn=1e-7,factor=99.)
        coeff = (coeff[1:]-1.)*rescale
        fit = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}
        solution[order[i]] = fit

    return solution,soln,model

