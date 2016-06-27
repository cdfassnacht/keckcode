
# coding: utf-8

# In[1]:

import numpy as np
import ccdredux as ccd
import imfuncs as imf
from astropy.io import fits as pf
from scipy.ndimage import filters
import shutil
from math import fabs
import glob


# In[2]:

""" 
Define a function to create a temporary sky frame from
three input frames
"""

def makesky_3(infiles, medians, indices):
    """ Makes a temporary sky file from 3 input files """
    
    tmpdat = pf.getdata(infiles[indices[0]])
    tmpshape = tmpdat.shape
    tmpsky = np.zeros((3,tmpshape[0],tmpshape[1]))
    for i in range(3):
        tmpsky[i,:,:] = pf.getdata(infiles[indices[i]]) / medians[indices[i]]
    outsky = np.median(tmpsky,axis=0)
    del tmpdat,tmpsky
    return outsky


# In[3]:

"""
Define a function that does a quick coadd of a list of input files given
pixel offsets between them.  The quick-and-dirty aspect to this processing
is that the function will just do integer pixel shifts.
"""

def quick_coadd(filelist, offsets, outfile):
    
    """ Start by shifting the offsets to be centered on the mean offset """
    dx = offsets[:,0] - (offsets[:,0].mean())
    dy = offsets[:,1] - (offsets[:,1].mean())
    dxrange = (dx.min(),dx.max())
    dyrange = (dy.min(),dy.min())
    
    """ Make a blank image of the appropriate size """
    dat0 = pf.getdata(filelist[0])
    xsize,ysize = dat0.shape[1],dat0.shape[0]
    del dat0
    outxsize = int(xsize + fabs(dxrange[0]) + fabs(dxrange[1]))
    outysize = int(ysize + fabs(dyrange[0]) + fabs(dyrange[1]))
    outim = np.zeros((outysize,outxsize))
    
    """ Insert the data with the appropriate offsets """
    x0 = fabs(dxrange[0])
    y0 = fabs(dyrange[0])
    x1 = (x0 - dx).astype(int)
    x2 = x1 + int(xsize)
    y1 = (y0 - dy).astype(int)
    y2 = y1 + int(ysize)
    print x1,x2
    print y1,y2
    print x2 - x1
    print y2 - y1
    for i in range(len(filelist)):
        tmp = pf.getdata(filelist[i])
        print i, tmp.shape
        try:
            outim[y1[i]:y2[i],x1[i]:x2[i]] += tmp
        except:
            print 'Failed on image %i (%s) with x1=%d,x2=%d,y1=%d,y2=%d' % (i,filelist[i],x1[i],x2[i],y1[i],y2[i])
        del tmp
    pf.PrimaryHDU(outim).writeto(outfile,clobber=True)


# In[4]:

""" 
A function that makes copies of the raw data files that we can freely modify.
Unfortunately, something is screwed up with the header in the raw files,
so we are just going to copy over the data, which is non-optimal
"""
def get_raw(sciframes):
    """ 
    Inputs:
      sciframes - list or array containing the frame numbers associated with
                  the observations of the desired object
      lensroot  - rootname for output files
    """
    rawdir = '../Raw_scam/'
    rawroot = 'jun18i00'   
    infiles = []
    for i in sciframes:
        rawname = '%s%s%02d.fits' % (rawdir,rawroot,i)
        workname = 'work_%02d.fits' % i
        infiles.append(workname)
        #shutil.copyfile(rawname,workname)
        data = pf.getdata(rawname,ignore_missing_end=True)
        pf.PrimaryHDU(data).writeto(workname)
        del data
    return infiles


# In[17]:

""" Get the raw data (modify for each lens-filter combination) """
sciframes = np.arange(46,59)
lensroot = 'J1618_Kp'
infiles = get_raw(sciframes)
print ''
print infiles


# In[18]:

""" Make the first sky frame """
skyname = '%s_sky1.fits' % lensroot
ccd.median_combine(infiles,skyname,normalize=True)


# In[19]:

""" 
Use the sky frame to make an initial bad pixel mask and then
update the sky frame itself to mask out the bad pixels.
"""
sky = imf.Image(skyname)

""" 
Do a 3-sigma clipping on the data and use the resulting clipped
mean and clipped rms to set the criterion for determining bad 
pixels
"""
sky.sigma_clip()
diff = np.fabs((sky.hdu[0].data - sky.mean_clip) / sky.rms_clip)

""" Create the bad pixel mask and write it out """
bpmask = diff>5.
tmp = bpmask.astype(int)
maskname = '%s_mask_sky1.fits' % lensroot
pf.PrimaryHDU(tmp).writeto(maskname,clobber=True)
del tmp

""" Replace the bad pixels with the median value in the image """
skydat = sky.hdu[0].data.copy()
skymed = np.median(skydat)
skydat[bpmask] = skymed

""" Save the result """
sky1v2name = '%s_sky1_v2.fits' % lensroot
pf.PrimaryHDU(skydat).writeto(sky1v2name,clobber=True)


# In[20]:

""" 
Replace the bad pixels in the input files using a 2-step process:
  1. Replace each of the bad pixels with the overall median value
  2. WAIT FOR NOW ON THIS

"""
for i in infiles:
    hdu = pf.open(i,mode='update')
    data = hdu[0].data
    datmed = np.median(data)
    data[bpmask] = datmed
    hdu.flush()
    del(hdu)


# In[21]:

""" Do a running sky subtraction """

""" Start by reading in all the files and calculating their median values """
alldat = np.zeros((len(infiles),bpmask.shape[0],bpmask.shape[1]))
allmed = np.zeros(len(infiles))
index_list = []
for i in range(len(infiles)):
    tmp = pf.getdata(infiles[i])
    allmed[i] = np.median(tmp)
    if i==0:
        index_list.append([0,1,2])
    elif i==(len(infiles)-1):
        index_list.append([i-2,i-1,i])
    else:
        index_list.append([i-1,i,i+1])
for i in range(len(infiles)):
    tmp = pf.getdata(infiles[i])
    tmpsub = tmp - allmed[i] * makesky_3(infiles,allmed,index_list[i])
    outname = 'tmpsub_%02d.fits' % sciframes[i]
    pf.PrimaryHDU(tmpsub).writeto(outname,clobber=True)
        


# In[22]:

"""
Create an updated sky frame by:
  1. Replacing the bad pixels in the sky frame with the median level of the image
"""

""" Replace the bad pixels with the median value in the image """
skydat = sky.hdu[0].data.copy()
skymed = np.median(skydat)
skydat[bpmask] = skymed

""" Smooth by quadrant MAY BE A MISTAKE"""
skysmo = np.ones(skydat.shape)
x = [0, 0, 128, 128]
y = [0, 128, 0, 128]
quad = np.ones((4,128,128))
for i in range(4):
    x1 = x[i]
    x2 = x1 + 128
    y1 = y[i]
    y2 = y1 + 128
    tmp = skydat[y1:y2,x1:x2]
    tmpsmo = filters.uniform_filter(tmp,3,cval=1.)
    skysmo[y1:y2,x1:x2] = tmpsmo.copy()

""" Save the result """
pf.PrimaryHDU(skydat).writeto('sky1_v2.fits')


# In[23]:

""" Do the initial sky subtraction using the updated sky frame """
sky1 = pf.getdata(sky1v2name)
skymed = np.median(sky1)
for i in infiles:
    scidat = pf.getdata(i).astype(float)
    scimed = np.median(scidat)
    scidat[bpmask] = scimed
    scimed2 = np.median(scidat)
    diff = scidat - (scimed2 / skymed) * sky1
    outfile = i.replace('work','sub1')
    pf.PrimaryHDU(diff).writeto(outfile)
    print 'Wrote sky-subtracted image (pass 1) to %s' % outfile
    del scidat,diff


# In[16]:

""" Do a quick-and-dirty coadd """
offsets = np.loadtxt('%s_offsets.txt' % lensroot)
subfiles = []
for i in infiles:
    subfiles.append(i.replace('work','sub1'))
quick_coadd(subfiles,offsets,'%s_coadd_quick.fits' % lensroot)


# In[ ]:



