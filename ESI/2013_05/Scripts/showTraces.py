import numpy,pyfits,pylab
import special_functions as sf


blue = [0,1500,1400,1300,1200,1100,900,600,200,0,0,0]
red = [0,3000,3400,3700,-1,-1,-1,-1,-1,-1,-1]


d = pyfits.open('STR_2938832465_0028_bgsub.fits')
d = pyfits.open('Feige66_0054_bgsub.fits')

out = []
for i in range(1,len(d)):
    slit = d[i].data.copy()
    trace = numpy.median(slit[:,blue[i]:red[i]],1)
    x = numpy.linspace(0,100,trace.size)
    pylab.plot(x,trace)
    xx = numpy.arange(trace.size)
    fit = sf.ngaussfit(numpy.array([xx,trace]).T,numpy.array([0.,trace[trace.argmax()],trace.argmax(),1.]))[0]
    out.append(fit[2])

pylab.figure()
d = pyfits.open('Feige66_0055_bgsub.fits')

out2 = []
for i in range(1,len(d)):
    slit = d[i].data.copy()
    trace = numpy.median(slit[:,blue[i]:red[i]],1)
    x = numpy.linspace(0,100,trace.size)
    pylab.plot(x,trace)
    xx = numpy.arange(trace.size)
    fit = sf.ngaussfit(numpy.array([xx,trace]).T,numpy.array([0.,trace[trace.argmax()],trace.argmax(),1.]))[0]
    out2.append(fit[2]-41.)

pylab.figure()
d = pyfits.open('Feige66_0056_bgsub.fits')

out3 = []
for i in range(1,len(d)):
    slit = d[i].data.copy()
    trace = numpy.median(slit[:,blue[i]:red[i]],1)
    x = numpy.linspace(0,100,trace.size)
    pylab.plot(x,trace)
    xx = numpy.arange(trace.size)
    fit = sf.ngaussfit(numpy.array([xx,trace]).T,numpy.array([0.,trace[trace.argmax()],trace.argmax(),1.]))[0]
    out3.append(fit[2]-41.)


pylab.figure()
out = numpy.array(out)
out2 = numpy.array(out2)
out3 = numpy.array(out3)

pylab.plot(out-out.mean(),'ko')
pylab.plot(out2-out2.mean(),'ro')
pylab.plot(out3-out3.mean(),'bo')


out = (out+out2+out3-out.mean()-out2.mean()-out3.mean())/3.
print sf.lsqfit(out,'polynomial',1)


pylab.plot([0.,9.],[9.47933,9.47933-2.1065*9],'k')



pylab.show()
