"""

Code to take the object information in the fits file that is produced by
DSIMULATOR and create an ascii file that can be the basis for a redshift
summary file.

"""

import sys
from astropy.io import fits as pf

""" Check the command-line format """
if len(sys.argv) < 2:
    print('')
    print('Usage: python make_dsim_summary.py [mask_name]')
    print('')
    print(' Example: python amke_dsim_summary.py 1206m4')
    print('')
    exit()

""" Set up the file names """
maskname = sys.argv[1]
infile = '%s.fits' % maskname
sumfile = '%s_summary.txt' % maskname

""" Open the input file """
try:
    hdu = pf.open(infile)
except IOError:
    print('')
    print('ERROR: %s was not found' % infile)
    print('')
    exit()

""" Get the object information """
objcat = hdu[1].data

""" Transfer the desired information into the output file """
f = open(sumfile, 'w')
f.write('# Object_ID       RA2000     Dec2000    Mask'
        '   Inst   ObsDate\n')
f.write('#-------------- ---------- ---------- -------'
        ' ------- -------\n')
for i in objcat:
    if i['objclass'] == 'Program_Target':
        f.write('%-15s %10.6f %+10.6f %-7s DEIMOS  20XX_XX\n' \
                    % (i['object'], i['ra_obj'], i['dec_obj'], maskname))
f.close()
print('')
print('Wrote information to %s' % sumfile)
print('')
