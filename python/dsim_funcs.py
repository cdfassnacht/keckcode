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
                racol=None, deccol=None, usecols=False, rafield='ra',
                decfield='dec'):
      """
      Creates an instance of the dsimCat class and loads it with data from
      the input file (catfile)
      """

      """ Initialize some variables """
      self.id        = None
      self.pri       = None
      self.theta     = None
      self.pa        = None
      self.selband   = None
      self.magname   = None
      self.dstab     = None
      self.selmask   = None

      """ 
      Read in the catalog and call the superclass initialization (to cf.Secat)
      for useful attributes
      """
      cf.Secat.__init__(self,catfile,catformat=catformat,namecol=namecol,
                        racol=racol,deccol=deccol,rafield=rafield,
                        decfield=decfield,usecols=usecols)

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

   def make_ids(self, root, ndigits=4):
      """
      Creates an array of object ids based on the position within the array
      """

      self.id = np.zeros(self.nrows,dtype='S16')
      for i in range(self.nrows):
         self.id[i] = '%s%04d' % (root,i+1)

   #----------------------------------------------------------------------

   def make_dstab(self, pri=None):
      """
      Creates a table of additional information about the selected objects.
      There are several assumptions that go into this:
        1. The selection mask, represented by the selmask array within the
           dsimCat structure, has been set.  If not, then dstab will be
           set with all of the members of the catalog.
        2. The make_ids method has previously run.  If self->id is None,
           as it is initialized to be during the creation of the structure,
           then make_ids is called with a root of 'obj'
        3. The galaxy position angles, represented by the theta array within
           the dsimCat structure have been set.  If self->theta is None
           (as it is initialized to be when the structure is created) then
           the 'pa' values within dstab will all be set to zero.
      With these assumptions, the size of dstab is set and the 'id' and 'pa'
      columns can be filled with initial values.  The 'pri' column will be
      set to its passed value if that is not None, otherwise it will be set
      to be identically zero.
      """

      """ 
      Test for existence of arrays that are needed to create the dstab 
      """
      if self.selmask is None:
         self.selmask = np.ones(self.nrows,dtype=bool)
      nsel = self.selmask.sum()
      if self.id is None:
         self.make_ids('obj')
      if self.theta is None:
         self.theta = np.zeros(self.nrows)

      """ Set the colums of dstab and create the table """
      selid = self.id[self.selmask]
      pa = self.theta[self.selmask]
      if pri is None:
         pri = np.zeros(nsel)
      self.dstab = Table([selid,pri,pa], names=('id','pri','pa'))

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
      dnames = ['id','ra','dec','equinox','mag','band','pri','samp','sel',
                'pa']
      outarr = np.zeros(nsel,dtype={'names':dnames,'formats':dfmt})
      outarr['id'] = self.dstab['id']
      outarr['ra'] = tmpra
      outarr['dec'] = tmpdec
      outarr['equinox'] += 2000.
      outarr['mag'] = seldata[self.magname]
      outarr['band'] = self.selband
      outarr['pri'] = self.dstab['pri']
      outarr['samp'] += 1
      outarr['pa'] = self.dstab['pa']
      #print outarr

      """ Write to the output file """
      outfmt = '%-16s %s %s %.1f %5.2f %s %4d %d %d %.1f'
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
