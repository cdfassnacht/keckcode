import esi
import special_functions

import pyfits,scipy,pickle,numpy
from scipy import stats,io,ndimage

def clip(data,clip=3.):
    d = data.copy()
    while 1:
        m,s,l = d.mean(),d.std(),d.size
        d = d[abs(d-m)<clip*s]
        if d.size==l:
            return m,s

def solve(d,orders):
    path = esi.__path__[0]

    lines = {}
    #lines['cuar'] = io.read_array(open(path+"/data/cuar.lines"))
    #lines['hgne'] = io.read_array(open(path+"/data/hgne.lines"))
    #lines['xe'] = io.read_array(open(path+"/data/xe.lines"))
    lines['cuar'] = numpy.loadtxt(path+"/data/cuar.lines")
    lines['hgne'] = numpy.loadtxt(path+"/data/hgne.lines")
    lines['xe'] = numpy.loadtxt(path+"/data/xe.lines")

    #startsoln = pickle.load(open(path+"/data/esi_wavesolution.dat"))
    startsoln = numpy.load(path+"/data/esi_wavesolution.dat")

    arclist = d.keys()
    soln = []
    xvals = scipy.arange(4096).astype(scipy.float32)
    for i in range(10):
      solution = startsoln[i]
      start,end = orders[i]

      for converge in range(10):
        wave = 10**special_functions.genfunc(xvals,0.,solution)

        peaks = {}
        trace = {}
        WIDTH = 4
        import pylab
        for arc in arclist:
            data = stats.stats.nanmedian(d[arc][start:end],axis=0)
            data[scipy.isnan(data)] = 0.
            if i==0 and arc=='cuar':
                data[3860:3940] = scipy.median(data)
            trace[arc] = data.copy()
            bak = ndimage.percentile_filter(data,50.,75.)
            data -= bak
            p = ndimage.maximum_filter(data,9)
            std = clip(scipy.trim_zeros(data),3.)[1]
            peak = scipy.where((p>30.*std)&(p==data))[0]
            peaks[arc] = []
            for p in peak:
                if p-WIDTH<0 or p+WIDTH+1>xvals.size:
                    continue
                x = xvals[p-WIDTH:p+WIDTH+1].copy()#-xvals[p]
                f = data[p-WIDTH:p+WIDTH+1].copy()
                fitdata = scipy.array([x,f]).T
                fit = scipy.array([0.,f.max(),xvals[p],1.])
                fit,chi = special_functions.ngaussfit(fitdata,fit,weight=1)
                peaks[arc].append(fit[2])#+xvals[p])
            #for i in lines[arc]:
            #    if i>wave[0] and i<wave[-1]:
            #        pylab.axvline(i)
            #pylab.plot(wave,data)
            #print arc
            #pylab.show()
        refit = []
        corr = []
        err = wave[wave.size/2]-wave[wave.size/2-1]
        for arc in arclist:
            p = 10.**special_functions.genfunc(peaks[arc],0.,solution)
            for k in range(p.size):
                cent = p[k]
                diff = cent-lines[arc]
                corr.append(diff[abs(diff).argmin()])
        corr = numpy.array(corr)
        m,s = clip(corr)
        corr = numpy.median(corr[abs(corr-m)<5.*s])

        for arc in arclist:
            p = 10.**special_functions.genfunc(peaks[arc],0.,solution)
            for k in range(p.size):
                pos = peaks[arc][k]
                cent = p[k]
                diff = abs(cent-lines[arc]-corr)
                if diff.min()<2.*err:
                    refit.append([pos,lines[arc][diff.argmin()]])
        refit = scipy.asarray(refit)
        solution = special_functions.lsqfit(refit,'polynomial',3)
        refit = []
        err = solution['coeff'][1]
        for arc in arclist:
            data = trace[arc]
            for pos in peaks[arc]:
                cent = special_functions.genfunc(pos,0.,solution)
                delta = 1e9
                match = None
                for j in lines[arc]:
                    diff = abs(cent-j)
                    if diff<delta and diff<1.*err:
                        delta = diff
                        match = j
                if match is not None:
                    refit.append([pos,match])
        refit = scipy.asarray(refit)
        refit[:,1] = scipy.log10(refit[:,1])
        solution = special_functions.lsqfit(refit,'chebyshev',3)

        #refit[:,0],refit[:,1] = refit[:,1].copy(),refit[:,0].copy()
        refit = numpy.array([refit[:,1],refit[:,0]]).T
        solution2 = special_functions.lsqfit(refit,'chebyshev',3)
        #soln.append([solution,solution2])

        w  = 10**special_functions.genfunc(xvals,0.,solution)
        if (w==wave).all():
            print "Order %d converged in %d iterations"%(i,converge)
            soln.append([solution,solution2])
            break
            for arc in arclist:
                pylab.plot(w,trace[arc])
            for j in 10**refit[:,0]:
                if j>w[0] and j<w[-1]:
                    pylab.axvline(j)
            pylab.show()
            break

    return soln
