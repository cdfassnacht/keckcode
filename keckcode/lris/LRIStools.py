import pyfits,sys,cPickle,numpy
from LRIS import LRIS_MOSobs,LRIS_MOSmask,LRIS_MOSslit


def waveSolve(outprefix,indir=None,prefix=None,files=None):
    obs = numpy.load('%s_flat.calib'%outprefix)

    if files is not None:
        obs.sciencefiles = []
        filenames = []
        for f in files:
            filenames.append('%s/%s%04d.fits'%(indir,prefix,f))
        obs.add_files(filenames)




def slitID(indir,prefix,files,outprefix,slits=None,side='top',minwidth=15):
    import pyfits,sys
    from LRIS import LRIS_MOSobs,LRIS_MOSmask,LRIS_MOSslit

    filenames = []
    for f in files:
        filenames.append('%s/%s%04d.fits'%(indir,prefix,f))
    maskname = pyfits.open(filenames[0])[0].header['SLITNAME']
    if maskname[:5]=='long_':
        isLongslit = True
    else:
        isLongslit = False

    mask = LRIS_MOSmask.LRIS_MOSmask(maskname)
    obs = LRIS_MOSobs.LRIS_MOSobs(mask)
    obs.add_files(filenames)
    obs.set_side(side)

    if slits is None or type(slits)==type('a'):
        print 'Determining Y solution'
        obs.set_ysoln()

        print 'Automatically Identifying Slits'
        obs.set_slits(useFlat=True)
        img = obs.tmp['flat']
        cut = img[:,img.shape[1]*3/8:img.shape[1]*5/8].mean(1)

        if obs.nslits == 0:
            print 'No slits found!'
            sys.exit()
        regions = []
        for i in range(obs.nslits,0,-1):
            slit = obs.slits[i-1]
            if obs.inst['instrument']=='LRIS_R_new' and slit.top>1005:
                slit.top = 1005
            width = slit.top-slit.bottom
            if width<minwidth:
                obs.nslits -= 1
                del obs.slits[i-1]
            else:
                regions.append([slit.bottom,slit.top])
        regions = regions[::-1]

        obj = IDLines(cut,regions)
        obs.slits = []
        obj.regions.sort()
        for i,j in obj.regions:
            obs.slits.append(obs.INST_MOSslit(i,j,obs))
        obs.nslits = len(obs.slits)

        for i in range(obs.nslits):
            b,t = obs.slits[i].bottom,obs.slits[i].top
            print "   %4d %4d"%(b,t)

        if type(slits)==type('a'):
            print "Only using slits:"
            use = slits.split()
            slits = []
            for i in range(len(use)):
                slits.append(obs.slits[int(use[i])])
                print "   %4d %4d"%(slits[-1].bottom,slits[-1].top)
            obs.nslits = len(slits)
            obs.slits = [s for s in slits]
    else:
        if isLongslit==False:
            obs.set_ysoln()
        else:
            obs.ysoln = {}
            obs.ysoln['ytrue'] = {'coeff':numpy.array([[0.,1.]]),'type':'polynomial'}
            obs.ysoln['ymap'] = {'coeff':numpy.array([[0.,1.]]),'type':'polynomial'}
            data = obs.get_flatavg()
            obs.tmp['flatavg'] = data.copy()
            obs.shape = data.shape
            coords = numpy.indices(data.shape)
            x = coords[1].flatten()
            y = coords[0].flatten()
            del coords
            from special_functions import genfunc
            obs.ysoln['forw'] = genfunc(x,y,obs.ysoln['ytrue']).reshape(data.shape)
            obs.ysoln['back'] = genfunc(x,y,obs.ysoln['ymap']).reshape(data.shape)
        obs.slits = []
        obs.nslits = 0
        for i,j in slits:
            obs.slits.append(obs.INST_MOSslit(i,j,obs))
            obs.nslits += 1

    arcXsoln(outprefix,obs)


def arcXsoln(outprefix,obs=None):
    if obs is None:
        obs = numpy.load('%s_flat.calib'%outprefix)
    print "Determining Arc X Solutions"
    for i in range(obs.nslits):
        slit = obs.slits[i]
        wb,wt,b,t = slit.widebottom,slit.widetop,slit.bottom,slit.top

        slit.set_arcsoln()

    for key in obs.tmp.keys():
        del obs.tmp[key]
    del obs.ysoln['forw']
    del obs.ysoln['back']
    obs.flatnorm = None

    f = open('%s_flat.calib'%outprefix,'wb')
    cPickle.dump(obs,f,2)
    f.close()


import pylab,numpy
class IDLines:
    def __init__(self,data,regions=[]):
        self.data = data
        self.orig = [r for r in regions]
        self.x = numpy.arange(data.size)
        self.plot = pylab.plot(self.x,self.data)[0]
        self.canvas = self.plot.get_figure().canvas

        self.regions = regions
        self.shades = []
        for x1,x2 in regions:
            self.shades.append(pylab.axvspan(x1,x2,ec='k',alpha=0.3,fc='g'))

        self.start = None
        self.end = None
        self.domotion = False
        self.span = None

        print """
        Click and drag the left mouse button to highlight slits.
        Click the right mouse button to remove a slit.
        Press 'w' to write the current fit to disk.
        Press 'r' to read a previous fit.
        Press 'n' to reset to the original input.

        When finished simply close the window!
        """
        self.keyid = self.canvas.mpl_connect('key_press_event',self.key_press)
        self.connect()
        pylab.show()


    def connect(self):
        self.pressid = self.canvas.mpl_connect('button_press_event',
                                                self.on_press)
        self.moveid = self.canvas.mpl_connect('motion_notify_event',
                                                self.on_motion)
        self.offid = self.canvas.mpl_connect('button_release_event',
                                                self.on_release)

    def on_press(self,event):
        if self.canvas.toolbar.mode!='':
            if event.button==2:
                self.canvas.toolbar.zoom()
                self.canvas.toolbar.pan()
                self.canvas.toolbar.pan()
            return

        if event.xdata==None or event.ydata==None: # Not in axes
            return

        indx = abs(self.x-event.xdata).argmin()
        if event.button==3:
            for i in range(len(self.regions)):
                x1,x2 = self.regions[i]
                if indx>=x1 and indx<=x2:
                    self.shades[i].set_visible(False)
                    del self.shades[i]
                    del self.regions[i]
                    pylab.draw()
                    return

        if event.button==1:
            self.domotion = True
            self.start = indx
            self.span = pylab.axvspan(indx,indx,ec='k',alpha=0.3,fc='g')

    def on_motion(self,event):
        if self.domotion==False:
            return

        if event.xdata is not None:
            self.span.set_xy([[self.start,0.],[self.start,1.],[event.xdata,1.],[event.xdata,0.],[self.start,0.]])
            self.end = event.xdata
            pylab.draw()

    def on_release(self,event):
        if self.domotion is False:
            return
        if event.xdata is None:
            end = self.end
        else:
            end = event.xdata
        end = int(end)
        start = self.start
        if start>end:
            start,end = end,start
        self.span.set_xy([[start,0.],[start,1.],[end,1.],[end,0.],[start,0.]])
        if end-start>5:
            self.regions.append([start,end+1])
            self.shades.append(self.span)
        else:
            self.span.set_visible(False)
        self.span = None
        self.start = None
        self.end = None
        self.domotion = False
        pylab.draw()

    def key_press(self,event):
        if event.key.lower()=='n':
            for i in range(len(self.regions)):
                self.shades[0].set_visible(False)
                del self.shades[0]
                del self.regions[0]
            # Add the original
            self.regions = [r for r in self.orig]
            for x1,x2 in self.regions:
                self.shades.append(pylab.axvspan(x1,x2,ec='k',alpha=0.3,fc='g'))
        elif event.key.lower()=='r':
            oname = raw_input('Name of input file: ')
            try:
                regions = numpy.load(oname)
                print "Reading slits from %s"%(oname)
            except:
                print "Could not read file: %s"%oname
                return
            # Delete the old
            for i in range(len(self.regions)):
                self.shades[0].set_visible(False)
                del self.shades[0]
                del self.regions[0]
            # Add the new
            self.regions = regions
            for x1,x2 in regions:
                self.shades.append(pylab.axvspan(x1,x2,ec='k',alpha=0.3,fc='g'))
        elif event.key.lower()=='w':
            oname = raw_input('Name of output file: ')
            try:
                f = open(oname,'wb')
                cPickle.dump(self.regions,f,2)
                f.close()
                print "File written successfully"
            except:
                'Could not write file %s'%(oname)
                return
        pylab.draw()
