"""
A set of functions related to preparing files related to using the DSIMULATOR
package for preparing DEIMOS slitmasks
"""

import numpy as np
import astropy
if astropy.__version__[:3] == '0.3':
   from astropy.coordinates import ICRS as SkyCoord
else:
   from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import catfuncs as cf
import sys

#---------------------------------------------------------------------------

class dsimCat(cf.Secat):
   """
   This class is essentially just the Secat class from catfuncs.py, but
   with a few extra components that will make it easier to interface with
   the files needed for and produced by DSIMULATOR
   """

   def __init__(self, catfile, catformat, verbose=True, namecol=None,
                racol=None, deccol=None, usecols=False):
      """
      Creates an instance of the dsimCat class and loads it with data from
      the input file (catfile)
      """

      """ Initialize some variables """
      self.id        = None
      self.pri       = None
      self.theta     = None
      self.theta_out = None
      self.selband   = None
      self.bandname  = None
      self.starband  = None
      self.selmask   = None
      self.dstab     = None

      """ 
      Read in the catalog and call the superclass initialization (to cf.Secat)
      for useful attributes
      """
      cf.Secat.__init__(self,catfile,catformat=catformat)

      """ 
      Make sure that there are RA and Dec coordinates, otherwise the catalog
      is useless for DSIMULATOR work
      """

      self.get_radec()
      if self.radec is None:
         print ''
         print 'ERROR: Input file %s does not have recognized RA and Dec data' \
             % self.infile
         print ''
         sys.exit()

      """ 
      Set up a second table that, for now, has dummy values but which will
      contain DSIMULATOR-specific information
      """
      #name  = np.zeros(self.nrows,dtype='S16')
      #pri   = np.zeros(self.nrows)
      #theta = np.zeros(self.nrows)
      #self.dstab = Table([name,pri,theta], names=('name','pri','theta'))

   #----------------------------------------------------------------------

   def write_gal(self, outfile, verbose=False):
      """
      Writes out the selected galaxies in the format needed for DSIMULATOR
      input. 

      NOTE: Many of the parameters that used to be passed to the stand-alone
      function, e.g., the priority vector, are now contained within the
      dsimCat structure itself and so no longer have to be explicitly
      passed in the call.
      """

      """ 
      Make sure that the priorities, etc. have been put into the dstab
      array within the dsimCat structure
      """
      if self.dstab is None:
         print ''
         print 'ERROR: Cannot write out galaxies without dstab (priorities,etc.)'
         print ' being set'
         print ''
         return

      """ 
      Check to make sure that the selection mask has been made.  If not, make
      a mask that selects all of the catalog members.
      """
      if self.selmask is None:
         self.selmask = np.ones(self.nrows,dtype=bool)

      """ Use the selection mask to pull out the selected objects """
      seldata  = self.data[self.selmask]
      selradec = self.radec[self.selmask]
      #selinfo  = self.dstab[self.selmask]

      """ Convert the radec info into the appropriate format """
      tmpra = selradec.ra.to_string(unit=u.hourangle,decimal=False,sep=':',
                                    precision=3,pad=True)
      tmpdec = selradec.dec.to_string(decimal=False,sep=':',precision=3,
                                      alwayssign=True,pad=True)

      """ 
      Store the output information in a numpy array for improved (hopefully)
      speed in writing out the data.
      """
      nsel = self.selmask.sum()
      dfmt = ['S16','S12','S13',float,float,'S2',int,int,int,float]
      dnames = ['name','ra','dec','equinox','mag','band','pri','samp','sel',
                'pa']
      outarr = np.zeros(nsel,dtype={'names':dnames,'formats':dfmt})
      outarr['name'] = self.dstab['name']
      outarr['ra'] = tmpra
      outarr['dec'] = tmpdec
      outarr['equinox'] += 2000.
      outarr['mag'] = seldata[self.selband]
      outarr['band'] = self.selband
      outarr['pri'] = self.dstab['pri']
      outarr['samp'] += 1
      outarr['pa'] = self.dstab['pa']
      #print outarr

      """ Write to the output file """
      outfmt = '%-16s %s %s %.1f %5.2f %s %d %d %d %.1f'
      np.savetxt(outfile,outarr,fmt=outfmt)

#---------------------------------------------------------------------------

def write_dsim(outcat, outradec, outid, pri, theta, selband, magname, outfile,
               verbose=False):
    """
    Writes the selected objects in a catalog out in the proper format to
    be used as inputs to DSIMULATOR
    """

    ngal = len(outradec)
    f = open(outfile,'w')
    tmpra = outradec.ra.to_string(unit=u.hourangle,decimal=False,sep=':',
                                  precision=3,pad=True)
    tmpdec = outradec.dec.to_string(decimal=False,sep=':',precision=3,
                                    alwayssign=True,pad=True)
    mag = outcat[magname]
    for i in range(ngal):
        #g = outcat[i]
        tmpstr = '%-16s %s %s 2000.0 %5.2f %s %d 1 0 %.1f' % \
            (outid[i],tmpra[i],tmpdec[i],mag[i],selband,pri[i],theta[i])
        f.write('%s\n' % tmpstr)
        if verbose:
            print tmpstr

    f.close()
