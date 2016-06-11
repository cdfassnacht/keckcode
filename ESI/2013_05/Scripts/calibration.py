import esi
from esi.biastrim import make_bias,biastrim

from esi.flat import *
from esi.straighten import startrace,straighten,fullSolution,getOrders
from esi import wavesolve

from spectra import spectools,offset,measure_width
from spectra.extract import extract
import special_functions as sf

from pickle import dump,load
from math import floor,ceil,fabs
import os,sys

import numpy,scipy,pyfits
from scipy import stats,interpolate,ndimage
from scipy import io as sio


def message(string):
    sys.stdout.write(string)
    sys.stdout.flush()

def clip(arr,nsig=3.5):
    a = arr.flatten()
    a = a[numpy.isfinite(a)]
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<nsig*s]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def prepare(dir,prefix,bias,stars,hgne,cuar,xe,flat,out_prefix,onearc=False,redoWave=False,arc=None):
    biasnums = bias.split(",")
    starnums = stars.split(",")
    flatnums = flat.split(",")

    for i in range(len(biasnums)):
        biasnums[i] = dir+prefix+biasnums[i]+".fits"
    for i in range(len(starnums)):
        starnums[i] = dir+prefix+starnums[i]+".fits"
    for i in range(len(flatnums)):
        flatnums[i] = dir+prefix+flatnums[i]+".fits"

    hgne = dir+prefix+hgne+".fits"
    cuar = dir+prefix+cuar+".fits"
    xe = dir+prefix+xe+".fits"

    message("Bias and bad pixel mask ... ")
    try:
        bias = pyfits.open(out_prefix+"_bias.fits")[0].data.copy()
        bpm = pyfits.open(out_prefix+"_bpm.fits")[0].data.copy()
        message("opened.\n")
    except:
        bias,bpm = make_bias(biasnums)
        pyfits.PrimaryHDU(bias).writeto(out_prefix+"_bias.fits")
        pyfits.PrimaryHDU(bpm).writeto(out_prefix+"_bpm.fits")
        message("created.\n")

    message("Flat field ... ")
    try:
        flat = pyfits.open(out_prefix+"_flat.fits")[0].data.copy()
        message("opened.\n")
    except:
        flat = make_flat(flatnums,out_prefix)
        pyfits.PrimaryHDU(flat).writeto(out_prefix+"_flat.fits")
        message("created.\n")

    message("Star traces ... ")
    try:
        star = pyfits.open(out_prefix+"_trace.fits")[0].data.astype(scipy.float32)
        message("opened.\n")
    except:
        from esi.straighten import starstack
        star = starstack(starnums,out_prefix)
        pyfits.PrimaryHDU(star).writeto(out_prefix+"_trace.fits")
        message("created.\n")

    message("Order solutions ... ")
    try:
        orders,solutions,wideorders = numpy.load(out_prefix+"_ycor.dat")
        message("opened.\n")
    except:
        orders = find_orders(flat)
        solutions,wideorders = startrace(star,orders)
        file = open(out_prefix+"_ycor.dat","wb")
        dump([orders,solutions,wideorders],file,2)
        file.close()
        message("created.\n")

    message("Normalized flat ... ")
    try:
        respflat = pyfits.open(out_prefix+"_norm.fits")
        message("opened.\n")
    except:
        cflat = response(flat,orders,solutions,wideorders)
        pyfits.PrimaryHDU(cflat).writeto(out_prefix+"_norm.fits")
        message("created.\n")

    if onearc is False:
      try:
          hgne = pyfits.open(out_prefix+"_hgne.fits")[0].data.astype(scipy.float32)
      except:
          print "Straightening HgNe arc"
          hdu = pyfits.open(hgne)
          hgne = hdu[0].data.copy()
          hgne = biastrim(hgne,bias,bpm)
          hgne = straighten(hgne,solutions,orders,wideorders)
          hdu[0].data = hgne
          hdu[0].writeto(out_prefix+"_hgne.fits")

      try:
          cuar = pyfits.open(out_prefix+"_cuar.fits")[0].data.astype(scipy.float32)
      except:
          print "Straightening CuAr arc"
          hdu = pyfits.open(cuar)
          cuar = hdu[0].data.copy()
          cuar = biastrim(cuar,bias,bpm)
          cuar = straighten(cuar,solutions,orders,wideorders)
          hdu[0].data = cuar
          hdu[0].writeto(out_prefix+"_cuar.fits")
      try:
          xe = pyfits.open(out_prefix+"_xe.fits")[0].data.astype(scipy.float32)
      except:
          print "Straightening Xe arc"
          hdu = pyfits.open(xe)
          xe = hdu[0].data.copy()
          xe = biastrim(xe,bias,bpm)
          xe = straighten(xe,solutions,orders,wideorders)
          hdu[0].data = xe 
          hdu[0].writeto(out_prefix+"_xe.fits")

      arcs = {}
      arcs['cuar'] = cuar.copy()
      arcs['hgne'] = hgne.copy()
      arcs['xe'] = xe.copy()
    else:
        try:
            archdu = pyfits.open(out_prefix+"_arc.fits")[0] 
            arc = archdu.data.copy()
            arcs = {}
            if archdu.header['LAMPCU1']=='on':
                arcs['cuar'] = True
            if archdu.header['LAMPNE1']=='on':
                arcs['hgne'] = True
            if archdu.header['LAMPAR1']=='on':
                arcs['xe'] = True
            arcs['arc'] = arc
        except:
            print "Straightening arc"
            arc = dir+prefix+arc+".fits"
            archdu = pyfits.open(arc)[0]
            arc = archdu.data.copy()
            arc = biastrim(arc,bias,bpm)
            arc = straighten(arc,solutions,orders,wideorders)
            archdu.data = arc
            archdu.writeto(out_prefix+"_arc.fits")
            arcs = {}
            if archdu.header['LAMPCU1']=='on':
                arcs['cuar'] = True
            if archdu.header['LAMPNE1']=='on':
                arcs['hgne'] = True
            if archdu.header['LAMPAR1']=='on':
                arcs['xe'] = True
            arcs['arc'] = arc


    if redoWave==False:
        try:
            wave_solution = numpy.load(out_prefix+"_wave.dat")
        except:
            redoWave = True
    if redoWave==True:
        print "Determining wavelength solutions"
        if onearc is False:
            wave_solution = wavesolve.solve(arcs,orders)
        else:
            wave_solution = wavesolve.jointSolve(arcs,orders)

        file = open(out_prefix+"_wave.dat","wb")
        dump(wave_solution,file,2)
        file.close()
    #print wave_solution

    print "Creating full slit solutions"
    soln = fullSolution(bias.shape,solutions,orders,wideorders,wave_solution)
    file = open(out_prefix+"_full.dat",'wb')
    dump(soln,file,2)
    file.close()

    print ""
    print "**********************"
    print "Preparations complete!"
    print "**********************"
    print ""


