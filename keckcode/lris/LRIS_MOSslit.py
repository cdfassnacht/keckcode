import scipy
from MOSredux.MOSslit import MOSslit


class LRIS_MOSslit(MOSslit):

    def __init__(self,bottom,top,obs):
#        super(LRIS_MOSslit,self).__init__(bottom,top,obs)
        MOSslit.__init__(self,bottom,top,obs)
        self.wsoln = None
        self.arcwsoln = None
        self.xywsoln = None


    def set_wsoln(self,order=3,usesky=True):
        wsoln = {}
        for f in self.obs.sciencefiles:
            wsoln[f] = self.get_wsoln(f,order,usesky)
        self.wsoln = wsoln


    def set_wsoln_lines(self,order=3,usesky=True,skylines=True,arclines=False,offset=None):
        from mostools import spectools
        from mostools import geometric_solutions as gs
        from mostools import specfit as spf
        if skylines==False:
            arclines = True

        if self.wsoln is None:
            self.set_wsoln(order,usesky)

        b,t = self.bottom,self.top
        wb,wt = self.widebottom,self.widetop
        yforw = self.obs.ysoln['forw'][b:t].copy()

        if arclines==True:
            from scipy import stats
            arcmodel = self.obs.get_arcmodel()
            try:
                arc = self.obs.tmp['arc'][b:t].copy()
            except:
                self.obs.tmp['arc'] = self.obs.straighten_arc()
                arc = self.obs.tmp['arc'][b:t].copy()
            
        for f in self.obs.sciencefiles:
            xsoln = self.xsoln[f]['forw']
            sci = self.obs.data[f]['data'][wb:wt].copy()
            sci = spectools.resampley(sci,yforw,wb)
            sci = scipy.sort(spectools.resample1d(sci,xsoln,'x'),0)
            sci = sci[sci.shape[0]/4].copy()
            if skylines==True and arclines==False:
                self.wsoln[f] = spf.wave_skylines(sci,self.wsoln[f])
            elif skylines==True and arclines==True:
                a = spectools.resample1d(arc,xsoln,'x')
                a = scipy.median(a,0)
                self.wsoln[f] = spf.wave_arcsky(a,arcmodel,sci,self.wsoln[f])
            else:
#                a = spectools.resample1d(arc,self.xsoln[f]['forw'],'x')
                a = spectools.resample1d(arc,self.arcsoln['forw'],'x')
                a = scipy.median(a,0)
                self.wsoln[f] = spf.wave_arclines(a,arcmodel,sci,self.wsoln[f],offset=offset,neworder=order)


    def get_wsoln(self,f,order=3,usesky=True,bcutoff=None,rcutoff=None):
        from mostools import spectools
        from mostools import geometric_solutions as gs
        from mostools import specfit as spf

        scale = self.obs.scale
        b,t = self.bottom,self.top
        wb,wt = self.widebottom,self.widetop

        yforw = self.obs.ysoln['forw'][b:t].copy()

        if self.xsoln is None:
            self.set_xsoln()
        xsoln = self.xsoln[f]['forw']

        sci = self.obs.data[f]['data'][wb:wt].copy()
        sci = spectools.resampley(sci,yforw,wb)
        sci = spectools.resample1d(sci,xsoln,'x')

        if usesky:
            wavemodel = self.obs.get_skymodel()
            dichroic = self.obs.inst['dichroic']
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
            return spf.modelmatch(sci,wavemodel,scale,order,[cutoff,9000.])
#            return gs.skywave(sci,wavemodel,scale,order,cutoff)
        else:
            from scipy import stats
            try:
                arc = self.obs.tmp['arc'][b:t].copy()
            except:
                self.obs.tmp['arc'] = self.obs.straighten_arc()
                arc = self.obs.tmp['arc'][b:t].copy()
            arcmodel = self.obs.get_arcmodel()
            a = spectools.resample1d(arc,self.xsoln[f]['forw'],'x')
            a = stats.stats.median(a,0)
            if bcutoff is None:
                bcutoff = arcmodel['blue']
            if rcutoff is None:
                rcutoff = arcmodel['red']
            return gs.arcwave2(a,arcmodel,scale,order,bcutoff,rcutoff)

    def set_xywsoln(self,out_xord=5,out_yord=3,xord=2,yord=2,usesky=True,
                        order=3):
        from mostools import geometric_solutions as gs
        if self.xsoln is None:
            self.set_xsoln(xord,yord,usesky)

        if self.wsoln is None:
            self.set_wsoln(order,usesky)

        b,t = self.bottom,self.top
        wb,wt = self.widebottom,self.widetop
        height = t-b
        width = self.obs.shape[1] 
        xywsoln = {}
        coords = scipy.indices((height,width)).astype(scipy.float32)
        for f in self.obs.sciencefiles:
            xywsoln[f] = gs.combine_xyw(coords,self.xsoln[f],
                        self.obs.ysoln['forw'][b:t]-wb,self.wsoln[f],
                        out_xord,out_yord)
        self.xywsoln = xywsoln
