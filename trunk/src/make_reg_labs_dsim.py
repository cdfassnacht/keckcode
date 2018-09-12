"""

Takes the fits file that is one of the main outputs of the DSIMULATOR program
and makes a ds9 region file that contains the object names at the correct
positions.  This region file is complementary to the slitmask region file
that DSIMULATOR produces, since the slitmask file does not have labels.

"""

import sys
from astropy.io import fits as pf

""" Check command line format """
if len(sys.argv) < 3:
    print('')
    print('Usage: python make_reg_labs_dsim.py [infits] [outreg]')
    print('')
    print('  where [infits] is the name of the fits file produced by'
          ' DSIMULATOR')
    print('  and [outreg] is the desired name of the output ds9 region file')
    print('')
    exit()

infits = sys.argv[1]
outreg = sys.argv[2]

hdu = pf.open(infits)
hdu.info()

tdat1 = hdu[1].data
f = open(outreg, 'w')
for i in tdat1:
    if i['ObjClass'] == 'Program_Target':
        f.write('fk5;text(%f,%f) # text={%s}\n' % \
                    (i['RA_OBJ'], i['DEC_OBJ'], i['OBJECT']))
f.close()
