import scipy,pickle,numpy
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    import pylab as plt
import special_functions as sf
from scipy import ndimage
from mostools import specfit as spf
STANDARD = None

class IDSpectrum:
    """
    IDSpectrum class for identification of spectral lines starting with an
        initial model of the spectrum.
    """
    def __init__(self, arc,sky,wave,standardlines):
        """ Plot data """
        self.arc = plt.plot(wave,arc,c='b')[0]
        self.sky = plt.plot(wave,sky,c='g')[0]
        self.canvas = self.arc.get_figure().canvas
        self.ax = self.arc.get_axes()
        self.xdata = self.arc.get_xdata().copy()
        self.start = [arc.copy(),sky.copy(),wave.copy(),standardlines.copy()]
        
        self.standard = standardlines

        """ Get metadata (ie line locations) for arcs """
        arc = self.arc.get_ydata()
        self.arcpeaks = ndimage.maximum_filter(arc,9)
        tmp = scipy.sort(arc)
        thresh = tmp[tmp.size*0.95]
        cond = (arc==self.arcpeaks)&(arc>thresh)
        self.arcpeaks = scipy.where(cond)[0]
        #tmpX = numpy.arange(arc.size)*1.
        #self.arcpeaks = spf.get_lines(tmpX,arc)
        self.arcsel = self.arcpeaks*0
        self.arclines = []
        for peak in self.arcpeaks:
            l = plt.axvline(self.xdata[peak],c='k')
            self.arclines.append(l)

        sky = self.sky.get_ydata().copy()
        bg = ndimage.percentile_filter(sky,25,51)
        self.skypeaks = ndimage.maximum_filter(sky,9)
        tmp = scipy.sort(sky)
        thresh = tmp[tmp.size*0.99]
        cond = (sky==self.skypeaks)&(sky>thresh)
        self.skypeaks = scipy.where(cond)[0]
        self.skysel = self.skypeaks*0
        self.skylines = []
        for peak in self.skypeaks:
            l = plt.axvline(self.xdata[peak],c='k',ls=':')
            self.skylines.append(l)

        self.spec = self.arc
        self.peaks = self.arcpeaks
        self.selected = self.arcsel
        self.lines = self.arclines

        """ Set useful flags """
        self.domotion = False
        self.origx = None
        self.soln = None
        self.pick = False
        self.fitlines = None
        self.keyid = self.canvas.mpl_connect('key_press_event',self.key_press)

        self.connect()

        print """
Mouse Controls:
    - left button drags single lines (rescales spectrum!)
    - middle button drags all lines (or exits from pan/zoom modes)
    - right button selects/deselects lines
Keyboard Commands:
    g - add new line (use mouse to select the line)
    m - fit a polynomial solution to the ID'd lines
    w - write the current state to disk
    r - read a saved state
    n - reset to the initial state
"""
        plt.show()

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

        self.xdata = self.spec.get_xdata().copy()
        if event.xdata==None or event.ydata==None: # Not in axes
            return

        ind = abs(self.xdata-event.xdata).argmin()
        indx = abs(self.peaks-ind).argmin()

        if abs(self.peaks-ind).min()>4.: # Not near any peaks
            return

        """ Select/unselect lines """
        if event.button==3:
            if self.selected[indx]==1:
                self.selected[indx] = 0
                self.lines[indx].set_color('k')
            else:
                self.selected[indx] = 1
                self.lines[indx].set_color('r')
            plt.draw()
            return


        self.origx = self.xdata[self.peaks[indx]]
        self.lbound = self.xdata[0]
        self.lindx = 0
        for i in range(indx):
            if self.xdata[self.peaks[i]]>self.origx:
                break
            if self.selected[i]==1:
                self.lbound = self.xdata[self.peaks[i]]
                self.lindx = self.peaks[i]

        self.rbound = self.xdata[-1]
        self.rindx = -1
        for i in range(indx+1,self.peaks.size):
            if self.xdata[self.peaks[i]]<self.origx:
                continue
            if self.selected[i]==1:
                self.rbound = self.xdata[self.peaks[i]]
                self.rindx = self.peaks[i]
                break
        self.selected[indx] = 1
        self.indx = self.peaks[indx]
        self.domotion = True
        self.lines[indx].set_color('r')
        self.pick = False

    def on_motion(self, event):
        """ Controls the sliding/stretching of the spectra """

        """
        Ignore this if we aren't in slide/stretch mode (ie pressing the
            mouse button
        """
        if self.domotion is False:
            return

        xdata = self.xdata.copy()
        """ Left mouse button is for stretching """
        if event.button==1 and (event.xdata is not None)  \
                and (event.ydata is not None) and event.xdata>self.lbound \
                and event.xdata<self.rbound:
            leftpts = self.xdata[self.lindx+1:self.indx+1].copy()
            left = scipy.linspace(leftpts[0],event.xdata,leftpts.size)
            rightpts = self.xdata[self.indx:self.rindx].copy()
            right = scipy.linspace(event.xdata,rightpts[-1],rightpts.size)[1:]
            xd = scipy.concatenate((left,right))
            xdata[self.lindx+1:self.rindx] = xd.copy()
            self.arc.set_xdata(xdata)
            self.sky.set_xdata(xdata)
        """ Middle mouse button is for sliding """
        if event.button==2:
            offset = event.xdata-self.origx
            xdata = xdata + offset
            self.arc.set_xdata(xdata)
            self.sky.set_xdata(xdata)
        for i in range(self.arcpeaks.size):
            x = xdata[self.arcpeaks[i]]
            l = self.arclines[i]
            l.set_xdata([x,x])
        for i in range(self.skypeaks.size):
            x = xdata[self.skypeaks[i]]
            l = self.skylines[i]
            l.set_xdata([x,x])
        plt.draw()

    def on_release(self, event):
        """ If the mouse button is released, reset! """
        if self.domotion:
            self.domotion = False
            self.xdata = self.spec.get_xdata().copy()
            plt.draw()

    def key_press(self,event):
        """
        m is for fitting, p is for selecting new lines (or leaving select
            mode) 
        """
        if event.key.lower()=='m':
            self.do_fit()
        elif event.key.lower()=='a':
            if self.pick:
                self.pick = False
                self.canvas.mpl_disconnect(self.addid)
                self.connect()
                self.pick = False
            else:
                self.pick = True
                self.disconnect()
                print "Choose the line to add (a to exit)"
                self.addid = self.canvas.mpl_connect('button_press_event',
                                                        self.add_line)
        elif event.key.lower()=='w':
            self.write()
        elif event.key.lower()=='r':
            self.read()
        elif event.key.lower()=='n':
            self.arc.remove()
            self.sky.remove()
            for i in self.lines:
                i.remove()
            a,b,c,d = self.start
            self.disconnect()
            self.canvas.mpl_disconnect(self.keyid)
            self.__init__(a,b,c,d)
            plt.draw()

    def do_fit(self):
        xdata = self.spec.get_xdata()
        ydata = self.arc.get_ydata()
        STANDARD = self.standard

        """ Get the centroided peaks of the selected lines """
        print "Finding line centroids..."
        nlines = self.arcsel[self.arcsel==1].size
        fit = scipy.zeros(nlines*3+1)
        for i in range(nlines):
            fit[i*3+1] = ydata[self.arcpeaks[self.arcsel==1][i]]
            fit[i*3+2] = float(self.arcpeaks[self.arcsel==1][i])
            fit[i*3+3] = 1.
        ofit,chi = sf.ngaussfit(ydata,fit)
        xpeaks = []
        wpeaks = []
        x = xdata[self.arcpeaks[self.arcsel==1]]
        for i in range(nlines):
            if abs(fit[i*3+2]-ofit[i*3+2])>2.5:
                continue
            diff = abs(STANDARD-x[i])
            if diff.min()<5.:
                xpeaks.append(ofit[i*3+2])
                wpeaks.append(STANDARD[diff.argmin()])

        fitdata = scipy.empty((len(xpeaks),2))
        fitdata[:,0] = scipy.asarray(xpeaks)
        fitdata[:,1] = scipy.asarray(wpeaks)
        self.fitlines = fitdata[:,1].copy()
        ord = int(raw_input('Enter order of fit: '))
        fit = sf.lsqfit(fitdata,'polynomial',ord)
        self.soln = fit
        xdata = sf.genfunc(scipy.arange(xdata.size),0.,fit)
        self.sky.set_xdata(xdata)
        self.arc.set_xdata(xdata)
        for i in range(self.arcpeaks.size):
            x = xdata[self.arcpeaks[i]]
            l = self.arclines[i]
            l.set_xdata([x,x])
        for i in range(self.skypeaks.size):
            x = xdata[self.skypeaks[i]]
            l = self.skylines[i]
            l.set_xdata([x,x])
        plt.draw()
    def add_line(self,event):
        if self.canvas.toolbar.mode!='':
            if event.button==2:
                self.canvas.toolbar.zoom()
                self.canvas.toolbar.pan()
                self.canvas.toolbar.pan()
            return
        if event.xdata==None or event.ydata==None:
            print 'Invalid data'
            return
        xpos = event.xdata
        xdata = self.spec.get_xdata()
        ydata = self.spec.get_ydata()
        p = ndimage.maximum_filter(ydata,9)
        ind = abs(xdata-xpos).argmin()
        p = scipy.where((p==ydata))[0]
        if abs(p-ind).min()>5:
            print 'Not a line'
            return
        indx = p[abs(p-ind).argmin()]
        for i in self.peaks:
            if abs(indx-i)<9:
                print 'Too close to another line.'
                return
        if indx<5. or indx>xdata.size-6:
                print 'Too close to edge.'
                return
        peaks = scipy.arange(self.peaks.size+1)
        n = self.peaks[self.peaks<indx].size
        peaks[:n] = self.peaks[:n].copy()
        peaks[n] = indx
        peaks[n+1:] = self.peaks[n:].copy()
        sel = scipy.arange(peaks.size)
        sel[:n] = self.selected[:n].copy()
        sel[n] = 0
        sel[n+1:] = self.selected[n:].copy()
        self.peaks = peaks.copy()
        self.selected = sel.copy()
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.spec==self.arc:
            l = plt.axvline(xdata[indx],c='k')
            self.lines.insert(n,l)
            self.arcpeaks = self.peaks
            self.arcsel = self.selected
            self.arclines = self.lines
        else:
            l = plt.axvline(xdata[indx],c='k',ls=':')
            self.lines.insert(n,l)
            self.skypeaks = self.peaks
            self.skysel = self.selected
            self.skylines = self.lines
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        plt.draw()
        self.canvas.mpl_disconnect(self.addid)
        self.connect()
        self.pick = False
    def disconnect(self):
        self.canvas.mpl_disconnect(self.pressid)
        self.canvas.mpl_disconnect(self.moveid)
        self.canvas.mpl_disconnect(self.offid)
    def write(self):
        oname = raw_input('Name of output file: ')
        f = open(oname,'w')
        xdata = self.spec.get_xdata()
        selected = self.arcsel
        peaks = self.arcpeaks
        soln = self.soln
        pickle.dump([xdata,selected,peaks,soln],f)
        f.close()
        print "Writing Complete"
    def read(self):
        oname = raw_input('Name of input file: ')
        try:
            f = open(oname,'r')
            xdata,selected,peaks,soln = pickle.load(f)
            f.close()
        except:
            print "Could not read file: %s"%oname
            return
        self.arc.set_xdata(xdata)
        self.arcpeaks = peaks
        self.arcsel = selected
        self.soln = soln
        for i in self.lines:
            i.remove()
        for i in range(len(self.lines)):
            del self.lines[0]
        for i in range(len(peaks)):
            c = 'k'
            if selected[i]==1:
                c = 'r'
            l = plt.axvline(xdata[peaks[i]],c=c)
            self.arclines.append(l)
        self.lines = self.arclines
        self.peaks = self.arcpeaks
        self.selected = self.arcsel
        plt.draw()



def id_spec(spec,model):
    import numpy
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fmt_xdata = plt.FormatStrFormatter('%4.2f')
    ax.fmt_ydata = plt.FormatStrFormatter('%4.2f')
    from scipy import ndimage,interpolate

    STANDARD = model['lines']
    arcmodel = model['matched']
    arc = spec.copy()
    blue,red,scale = model['blue'],model['red'],model['scale']
    while (red-blue)/scale<arc.size:
        red += scale
        blue -= scale
    wave = numpy.arange(blue,red,scale)
    blue,red = model['blue'],model['red']
    model = interpolate.splev(wave,arcmodel)
    c = (wave<blue)|(wave>red)
    model[c] = model[~c].min()
    sci = arc*0. 

    model /= model.max()/arc[numpy.isfinite(arc)].max()
    plt.plot(wave,model,c='gray')

    lines = []

    dr = IDSpectrum(arc,sci,wave[:arc.size],STANDARD)

    if dr.soln is None:
        dr.do_fit()
    xdata = dr.arc.get_xdata()
    ydata = dr.arc.get_ydata()
    wave = sf.genfunc(numpy.arange(xdata.size),0.,dr.soln)
    wide = ndimage.gaussian_filter(ydata,3.)
    match = interpolate.splrep(wave,ydata)
    wide = interpolate.splrep(wave,wide)
    model = {'wide':wide,'matched':match,'orig':[wave,ydata],
            'lines':dr.fitlines,'red':red,'blue':blue,'scale':scale}
    return model,dr.soln
