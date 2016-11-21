from LRIS.resample import resample
from LRIS.LRIStools import *
from LRIS.XSolve import *

indir = 'data'
pref = 'r111130_'
flat = 102
arc = 100
outname = 'mask1_red_bottom'
slitID(indir,pref,[flat,arc,65],outname,side='bottom')
oldName = None

# Update to determine whether you want to see the fit for every exposure
sf = True
oldName = '%s_65'%(outname)
for img in [67,68,69,70,71,72,73,74]:
    newName = '%s_%2d'%(outname,img)
    XSolve(outname,newName,indir,pref,[flat,arc,img])
    SlitCross(newName,showFit=sf)
    WaveSolve(newName,oldName,showFit=sf)
    resample(newName,nobgsub=False,clobber=True)
    oldName = newName

    # Don't show the fit anymore...
    sf = False
