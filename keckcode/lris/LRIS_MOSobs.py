import pyfits,scipy,sys
from MOSredux.MOSobs import MOSobs
from LRIS_MOSslit import LRIS_MOSslit

DispScale = {
    'LRIS_R_new':
        {'150/7500':3.0,
        '300/5000':1.59,
        '400/8500':1.16,
        '600/5000':0.80,
        '600/7500':0.80,
        '600/10000':0.80,
        '831/8200':0.58,
        '900/5500':0.53,
        '1200/7500':0.40},
    'LRIS_R_old':
        {'150/7500':4.8,
        '300/5000':2.45,
        '400/8500':1.85,
        '600/5000':1.25,
        '600/7500':1.25,
        '600/10000':1.25,
        '831/8200':0.915,
        '900/5500':0.85,
        '1200/7500':0.64},
    'LRISBLUE':
        {'300/5000':1.45,
        '400/3400':1.05,
        '600/4000':0.63,
        '1200/3400':0.24}}

Gain = {
    'LRIS_R_old':2.,
    'LRIS_R_new':0.95,
    'LRISBLUE':1.6}

ReadVar = {
    'LRIS_R_old':36.,
    'LRIS_R_new':22.,
    'LRISBLUE':17.}

class LRIS_MOSobs(MOSobs):

    def __init__(self,input,files=None):
        if type(input)==type('str'):
            o = self.explode(input)
            self.arcfile = o.arcfile
            self.flatfiles = o.flatfiles
            self.sciencefiles = o.sciencefiles
            self.biasfiles = o.biasfiles
            self.inst = o.inst
            self.ysoln = o.ysoln
            self.slits = o.slits
            self.nslits = o.nslits
            self.maskname = o.maskname
            self.offsets = o.offsets
            self.INST_MOSslit = o.INST_MOSslit
            self.bias = o.bias
            self.skymodel = o.skymodel
            self.arcmodel = o.arcmodel
            self.flatnorm = o.flatnorm
            self.data = o.data
            self.linewidth = o.linewidth
            self.center = o.center
            self.scale = o.scale
            self.tmp = o.tmp

            if o.side is not None:
                self.set_side(o.side)
            self.side = o.side

            return None
        else:
            mask = input

        self.arcfile = None
        self.flatfiles = []
        self.sciencefiles = []
        self.biasfiles = []
        self.inst = None
        self.ysoln = None
        self.slits = None
        self.nslits = None
        self.maskname = mask.maskname
        self.offsets = None
        self.INST_MOSslit = LRIS_MOSslit
        self.shape = None
        self.blueside = None

        self.bias = None

        self.skymodel = None
        self.arcmodel = None
        self.flatnorm = None
        self.data = {}
        self.linewidth = None
        self.center = None

        self.tmp = {}
        self.scale = None

        if files is not None:
            self.add_files(files)


    def add_files(self,files):
        """
        add_files(files)

        Adds science or calibration files to the observation.
            files is a list of filenames
        """
        for f in files:
            hdu = pyfits.open(f)
            instrument = self.get_instrument(hdu)
            if instrument['maskname']!=self.maskname:
                print 'ERROR: Masks did not match parent mask.'
                return
            if self.inst is None:
                self.inst = instrument
            else:
                for k in self.inst.keys():
                    if self.inst[k]!=instrument[k]:
                        print 'ERROR: Instruments do not match: ',self.inst[k],instrument[k]
                        sys.exit()
                        return

          #  lamps = instrument['lamps']
            lamps = [int(j) for j in hdu[0].header['LAMPS'].split(',')]

            # The ordering of the LAMPS keyword changed in Feb 2001
            date = hdu[0].header['DATE'].split('-')
            year = int(date[0])
            month = int(date[1])
            if year<2001 or (year==2001 and month<2):
                tmp = lamps[-1]
                lamps[-1] = lamps[3]
                lamps[3] = tmp
                del tmp

            if max(lamps)==0:   # Lamps off -> Science file
                if f not in self.sciencefiles:
                    self.sciencefiles.append(f)
            else:
                if lamps[-1]==max(lamps[:-1]):  # The flat is on with an arc
                    print 'WARNING: Flat lamps and arc lamps on simultaneously for file %s.' % f
                if lamps[-1]==1:    # Flat lamp is on
                    if f not in self.flatfiles:
                        self.flatfiles.append(f)
                else:
                    self.arcfile = f
                    self.lamps = lamps


    def get_instrument(self,hdulist):
        if len(hdulist)<2:
            return self.get_oldinst(hdulist)
        else:
            return self.get_newinst(hdulist)

    def get_oldinst(self,hdulist):
        inst = {}
        hdu = hdulist[0]
        try:
            inst['maskname'] = hdu.header['SLITNAME']
            instrument = hdu.header['INSTRUME']
            if instrument=='LRIS':
                inst['instrument'] = 'LRIS_R_old'
            else:
                inst['instrument'] = 'LRISBLUE'
            if inst['instrument']=='LRISBLUE':
                inst['disperser'] = hdu.header['GRISNAME']
            else:
                inst['disperser'] = hdu.header['GRANAME']
            try:
                inst['dichroic'] = hdu.header['DICHNAME']
                if inst['dichroic'].lower()=='clear':
                    inst['dichroic'] = ''
            except:
                inst['dichroic'] = ''
        except:
            return None
        """
        lamps = [int(j) for j in hdulist[0].header['LAMPS'].split(',')]

        # The ordering of the LAMPS keyword changed in Feb 2001
        date = hdulist[0].header['DATE'].split('-')
        year = int(date[0])
        month = int(date[1])
        if year<2001 or (year==2001 and month<2):
            tmp = lamps[-1]
            lamps[-1] = lamps[3]
            lamps[3] = tmp
            del tmp
        """
        self.scale = DispScale[inst['instrument']][inst['disperser']]
        return inst


    def get_newinst(self,hdulist):
        hdu = hdulist[0]
        inst = {}
        try:
            inst['maskname'] = hdu.header['SLITNAME']
            instrument = hdu.header['INSTRUME'].strip()
            if instrument=='LRIS':
                inst['instrument'] = 'LRIS_R_new'
            else:
                inst['instrument'] = 'LRISBLUE'
            if inst['instrument']=='LRISBLUE':
                inst['disperser'] = hdu.header['GRISNAME']
                bindex = 0
            else:
                inst['disperser'] = hdu.header['GRANAME']
                bindex = 1
            inst['dichroic'] = hdu.header['DICHNAME']
            if inst['dichroic'].lower()=='clear':
                inst['dichroic'] = ''
            inst['binning'] = [int(i) for i in hdu.header['CCDSUM'].split()] 
        except:
            return None
        self.scale = DispScale[inst['instrument']][inst['disperser']]
        self.scale *= inst['binning'][bindex]
        return inst


    def set_side(self,side=None):
        import lrisbiastrim
        if side is None:
            side = self.side
        lrisbiastrim.SIDE = side
        self.biastrim = lrisbiastrim.biastrim
        self.side = side

    def prep_science(self):
        from mostools.offset import findoffset
        from scipy import stats
        self.set_side(self.side)
        try:
            tmp = self.ysoln['yforw']
            del tmp
        except:
            self.set_ysoln_arrays()

        if self.inst['maskname'][:5]=='long_':
            centkey = 'WAVELEN'
        else:
            centkey = 'MSWAVE'

        for f in self.sciencefiles:
            tmp = pyfits.open(f)
            if len(tmp)==1:
                data = self.biastrim(tmp[0].data.copy())
            else:
                data = self.biastrim(tmp)

            if self.flatnorm is None:
                self.normalized_flat()

#            if self.inst['instrument']!='LRISBLUE':
            data /= self.flatnorm

            airmass = tmp[0].header['AIRMASS']
            if self.inst['instrument']!='LRISBLUE':
                mswave = tmp[0].header[centkey]
            else:
                mswave = 4200.

            self.shape = data.shape
            self.data[f] = {'data':data,'airmass':airmass,'mswave':mswave}
            self.inst['CENTRAL'] = mswave
            self.readvar = ReadVar[self.inst['instrument']]

    def make_skymodel(self):
        import pickle,os,numpy
        from scipy import interpolate,ndimage
        from mostools import specfit as spf

        path = os.path.dirname(__file__)
        skyfile = path+"/data/uves_sky.model"
        f = open(skyfile)
        wavecalmodel = pickle.load(f)
        f.close()
        wave = scipy.arange(3400.,10400.,0.1)

        scale = self.scale
        if self.linewidth is None:
            self.set_linewidth()
        lw = self.linewidth*scale

        dichroic = self.inst['dichroic']
        if dichroic!='' and dichroic!='mirror':
            file = "%s/data/dichroic_%s_t.dat" % (path,dichroic)
            input = numpy.loadtxt(file)
            spline = interpolate.splrep(input[:,0],input[:,1],s=0)
            dich = interpolate.splev(wave,spline)
            dich[wave<4500.] = 0.
            dich[wave>8800.] = 1.
            del input,spline
            if self.inst['instrument']=='LRISBLUE':
                dich = (dich*-1.)+1.
        else:
            dich = scipy.ones(wave.size)

        if dichroic=='' or dichroic=='mirror':
            cutoff = None
        elif dichroic=='680':
            cutoff = 6700.
        elif dichroic=='560':
            cutoff = 5500.
        elif dichroic=='500':
            cutoff = 4900.
        elif dichroic=='460':
            cutoff = None
        if self.inst['instrument']=='LRISBLUE':
            rcutoff = cutoff
            bcutoff = None
        else:
            bcutoff = cutoff
            rcutoff = self.inst['CENTRAL']+3000*scale
            if rcutoff>10400.:
                rcutoff = 10400.

        wavemodel = interpolate.splev(wave,wavecalmodel)*dich
        wavemodel_match = ndimage.gaussian_filter1d(wavemodel,0.93*lw/0.1)
        owavemodel_match = interpolate.splrep(wave,wavemodel_match,s=0)
        wavemodel_wide = ndimage.gaussian_filter1d(wavemodel,5*lw/0.1)
        owavemodel_wide = interpolate.splrep(wave,wavemodel_wide,s=0)
        wavemodel = spf.specContinuumSub(wavemodel)
        wavemodel_match = ndimage.gaussian_filter1d(wavemodel,0.93*lw/0.1)
        wavemodel_match = interpolate.splrep(wave,wavemodel_match,s=0)
        wavemodel_wide = ndimage.gaussian_filter1d(wavemodel,5*lw/0.1)
        wavemodel_wide = interpolate.splrep(wave,wavemodel_wide,s=0)

        self.skymodel = {'wide':wavemodel_wide,'matched':wavemodel_match,
                            'red':rcutoff,'blue':bcutoff,'scale':scale,
                            'wideC':owavemodel_wide,'matchedC':owavemodel_match}


    def get_skymodel(self):
        if self.skymodel is None:
            self.make_skymodel()
        return self.skymodel


    def set_skymodel(self,model=None):
        if model is None:
            self.skymodel = self.get_skymodel()
            return
        from scipy import ndimage,interpolate
        wave = model['wave']
        sky = model['sky']
        blue,red,scale = model['blue'],model['red'],model['scale']
        wave = wave[~scipy.isnan(sky)]
        sky = sky[~scipy.isnan(sky)]
        wide = ndimage.gaussian_filter(sky,5.)
        match = interpolate.splrep(wave,sky,s=0)
        wide = interpolate.splrep(wave,wide,s=0)
        self.skymodel = {'wide':wide,'matched':match,'red':red,'blue':blue,
                            'scale':scale}


    def make_arcmodel(self):
        import os,pickle
        from scipy import interpolate,signal,ndimage
        from special_functions import ngauss
        from math import ceil,floor
        path = os.path.dirname(__file__)
        arcfile = path+"/data/arclamps.dat"
        lampkey = ['Hg','Ne','Ar','Cd','Zn']

        f = open(arcfile)
        lamplines = pickle.load(f)
        f.close()

        normfile = path+"/data/blue_normlines.dat"
        f = open(normfile)
        normlines = pickle.load(f)
        f.close()

        if self.linewidth is None:
            self.set_linewidth()

        lw = self.linewidth
        scale = self.scale

        wave = scipy.arange(2000.,10400.,scale)
        finemodel = scipy.zeros(wave.size)
        widemodel = finemodel*0.
        norm = finemodel*0.
        used_lines = []
        for i in range(len(self.lamps)):
            if self.lamps[i]==0:
                continue
            fit = [0.]
            lines = lamplines[lampkey[i]]
            for l,a in lines:
                fit.append(a)
                fit.append(l)
                fit.append(lw)
                used_lines.append(l)
            fit = scipy.asarray(fit)
            finemodel += ngauss(wave,fit)
            fit[fit==lw] = lw*10.
            widemodel += ngauss(wave,fit)
            fit = [0.]
            for l in normlines[lampkey[i]]:
                fit.append(1.)
                fit.append(l)
                fit.append(lw*5.)
            fit = scipy.asarray(fit)
            norm += ngauss(wave,fit)

        finemodel = interpolate.splrep(wave,finemodel)
        widemodel = interpolate.splrep(wave,widemodel)
        norm = interpolate.splrep(wave,norm)

        lines = scipy.sort(scipy.asarray(used_lines))
        red = ceil((lines.max()+100.)/500.)*500.
        blue = floor((lines.min()-100.)/500.)*500.
        self.arcmodel = {'wide':widemodel,'matched':finemodel,'norm':norm,
                            'lines':lines,'red':red,'blue':blue,'scale':scale}


    def get_arcmodel(self):
        if self.arcmodel is None:
            self.make_arcmodel()
        return self.arcmodel


    def set_arcmodel(self,model=None):
        if model is None:
            self.arcmodel = self.get_arcmodel()
            return
        from scipy import ndimage,interpolate
        wave = model['wave']
        arc = model['arc']
        blue,red,scale = model['blue'],model['red'],model['scale']
        wave = wave[~scipy.isnan(arc)]
        arc = arc[~scipy.isnan(arc)]
        wide = ndimage.gaussian_filter(arc,5.)
        match = interpolate.splrep(wave,arc,s=0)
        wide = interpolate.splrep(wave,wide,s=0)

        self.arcmodel = {'wide':wide,'matched':match,'orig':model['orig'],
                            'red':red,'blue':blue,'scale':scale}

