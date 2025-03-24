import matplotlib as mpl
import pylab,numpy,pyfits,cPickle,sys
from scipy import interpolate,ndimage
from mostools import spectools as st


class myDraw:

    def __init__(self,w,d):
        self.w = w
        self.d = d
        self.var = d.var()
        print self.var
        self.s = pylab.plot(w,d)[0]
        self.m = pylab.plot(w,d)[0]
        self.f = pylab.plot(w,d*numpy.nan)[0]
        self.canvas = self.s.get_figure().canvas
        self.ax = self.s.get_axes()
        self.smooth = 0

#        self.m.set_xdata([])
#        self.m.set_ydata([])


        self.domotion = False
        self.gap = []
        self.keyid = self.canvas.mpl_connect('key_press_event',self.key_press)

        self.connect()
        self.xpoints = []
        self.ypoints = []
        pylab.show()

    def connect(self):
        """ Connect the mouse to the plot """
        self.pressid = self.canvas.mpl_connect('button_press_event',
                                                self.on_press)
        self.moveid = self.canvas.mpl_connect('motion_notify_event',
                                                self.on_motion)
        self.offid = self.canvas.mpl_connect('button_release_event',
                                                self.on_release)

    def on_press(self,event):
        """
        Deal with mouse button presses, including stretching, shifting,
            and line identification.
        """
        """ Turn off plot tools """
        if self.canvas.toolbar.mode!='':
            if event.button==2:
                self.canvas.toolbar.zoom()
                self.canvas.toolbar.pan()
                self.canvas.toolbar.pan()
            return

        if event.xdata==None or event.ydata==None: # Not in axes
            return

        self.domotion = True
        self.gap.append(len(self.xpoints))
        self.xpoints.append(event.xdata)
        self.ypoints.append(event.ydata)


    def on_motion(self,event):
        if self.domotion is False:
            return

        if event.button==1 and (event.xdata is not None)  \
                and (event.ydata is not None):
            #self.xpoints.append(event.xdata)
            #self.ypoints.append(event.ydata)
            xpoints = self.m.get_xdata()
            ypoints = self.m.get_ydata()
            arg = abs(xpoints-event.xdata).argmin()
            ypoints[arg] = event.ydata
            self.xpoints = xpoints
            self.ypoints = ypoints
            self.m.set_xdata(xpoints)
            self.m.set_ydata(ypoints)

            pylab.draw()

    def on_release(self, event):
        """ If the mouse button is released, reset! """
        if self.domotion:
            self.domotion = False
#            args = numpy.array(self.xpoints).argsort()
#            self.xpoints = self.xpoints[args]
#            self.ypoints = self.ypoints[args]
            pylab.draw()

    def key_press(self,event):
        if event.key.lower()=='b' and len(self.gap)>0:
            self.xpoints = self.xpoints[:self.gap[-1]]
            self.ypoints = self.ypoints[:self.gap[-1]]
            self.m.set_xdata(self.xpoints)
            self.m.set_ydata(self.ypoints)
            del self.gap[-1]
            pylab.draw()

        if event.key.lower()=='m':
            x = numpy.array(self.xpoints)
            x += numpy.random.random(x.size)*1e-7
            y = numpy.array(self.ypoints)
            args = x.argsort()
            x = x[args]
            self.model = interpolate.splrep(x,y[args],s=self.smooth*self.var,k=3)
            self.f.set_ydata(interpolate.splev(self.w,self.model))
            pylab.draw()
            while 1:
                resp = raw_input('Enter smoothing (current=%4.1f, q to quit): '%self.smooth)
                if resp=='q':
                    break
                self.smooth = float(resp)
                self.model = interpolate.splrep(x,y[args],s=self.smooth*self.var,k=3)
                self.f.set_ydata(interpolate.splev(self.w,self.model))
                pylab.draw()

        if event.key.lower()=='w':
            oname = raw_input('Enter output filename: ')
            f = open(oname,'wb')
            cPickle.dump(self.model,f,2)
            f.close()
            print "Successfully wrote spline model to %s"%oname

        if event.key.lower()=='read':
            oname = raw_input('Enter input filename: ')
            self.model = numpy.load(oname)
            self.f.set_ydata(interpolate.splev(self.w,self.model))
            pylab.draw()


def doResponse(modelName,specName):
    d = numpy.loadtxt(modelName)
    M = interpolate.splrep(d[:,0],d[:,1],s=0,k=3)

    w = st.wavelength(specName,1)
    d = pyfits.open(specName)[1].data.copy()

    m = interpolate.splev(w,M)
    d = myDraw(w,m/d)
