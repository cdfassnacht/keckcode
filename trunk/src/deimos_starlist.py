#!/usr/bin/env python

"""Given a set of DEIMOS DSIMULATOR ASCII output files, Generates 
a Keck-formatted starlist.

Usage:
        Change to directory with your mask designs
        deimos_starlist [-D] [file1 .. fileN]

Switches:
        -D = debug mode; print each input line as it is processed

Args:
        fileN = DSIMULATOR output file list
       
Output: 
	Keck-formatted starlist file, "starlist"

External modules needed:
        None

Version:
	Tested in Python 2.4.6 on Solaris10

Examples:
        1) Generate a starlist for DSIMULATOR files mask1.out and
        mask2.out:
                python deimos_starlist.py mask1.out mask2.out

        2) Generate a starlist for DSIMULATOR files mask1.out and
        mask2.out with debugging output:
                python deimos_starlist.py -D mask1.out mask2.out

        3) Generate a starlist for all DSIMULATOR files *.out:
                python deimos_starlist.py *.out

Modification History:
        2003-Mar-28     DKM     Original version
        2012-Nov-11     GDW     Fixed print format to handle -10<Dec<0
        2013-Mar-10     GDW     - Reverted mods from Nov, which broke
                                behavior at -1 < Dec < 0;
                                - Removed code which turned all spaces in
                                target name to underscores in starlist;
                                - Round all guider coords to nearest int
                                - Move guide star labels closer to star
        2013-Aug-27     GDW     - Add sanity check for undefined values on
                                tv coords
                                - add offsets in tv x and y
        2013-Sep-08     GDW     - Fixed problem with double-minus in DEC
                                - Added offsets to printout
        2013-Dec-20     GDW     - Added debug mode switch
	2015-Feb-19	jlyke	-Adapted for just starlist output
                                """

import math
import re

class Mask:
    def __init__(self):
        self.selectedObjects = []
        self.alignmentStars = []
        self.guideStars = []
        self.debug = False

    def readMaskFile(self, maskfile):

        print "Processing file %s" % maskfile

        self.file = maskfile
        f = open(self.file,"r")
        l = f.readlines()
        f.close()

        # parse the second line of the file, which contains the field
        # name, field center, and equinox...
        if self.debug:
            print l[1]
        pattern = "^(.+)\s+(\d+:\d+:\d+\.\d+)\s+(\S*\d+:\d+:\d+\.\d+)\s+(\S+)\s+PA=\s*(\S*)\s+##"
        mobj = re.match( pattern, l[1])
        if mobj:
            self.name = mobj.group(1)
            self.ra = mobj.group(2)
            self.dec = mobj.group(3)
            self.equinox = float(mobj.group(4))
            self.pa = float(mobj.group(5))
        else:
            print '  ERROR: line 2 does not match expected format -- abort!'
            sys.exit(1)


if __name__ == '__main__':
    import os, sys, glob
    import math
    import re
    import getopt

    usage = "Usage: "+sys.argv[0]+" [-h] [-d] filename .. filenameN"

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'hD')
    except getopt.GetoptError,err:
        print(err)
        print usage
        sys.exit(2)

    debug = False
    for o,a in optlist:
        if o == "-h":
            print usage
            sys.exit(1)
        elif o in "-D":
            print "DEBUG mode enabled"
            debug = True
        else:
            assert False, "unhandled option"
        
    mlist = args
    slf = open('starlist','w')

    for ml in mlist:
        input = ml
        # make mask instanace
        m = Mask()
        m.debug = debug
        try:
            m.readMaskFile(input)
        except Exception, err:
            print 'ERROR: Unable to read mask file: '+ml
            print err
            sys.exit(1)

        # build starlist
        # Modified by DKM 2009-03-30: Added PA info
        sra = m.ra.split(':')
        sdec = m.dec.split(':')
        srah = int(sra[0])
        sram = sra[1]
        sras = sra[2]
        sdech = sdec[0]
        sdecm = sdec[1]
        sdecs = sdec[2]
        spa = float(m.pa)

        # determine whether the DEC is + or -.  NOTE: this is
        # necessary because simply printing the DEC in %+03i format
        # will yield "+00" for -1 < Dec < 0!
        mobj = re.match( "^-", sdech)
        if mobj:
            sign = '-'
        else:
            sign = '+'
        adech = abs(int(sdech))

        # Sky's got no love for a negative PA
        if spa < 0:
            spa += 360
        seqnx = str(m.equinox).split('.')[0]
        targname = m.name[:16]
        slist = (targname,srah,sram,sras,sign,adech,sdecm,sdecs,seqnx,spa)
        sformat = '%-16s  %02i %s %s %s%02i %s %s %s rotdest=%.2f rotmode=pa\n'
        print "  "+sformat % slist
        slf.write(sformat % slist)

