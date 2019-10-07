"""

cubeset - Defines a CubeSet class that contains code to handle operations
          on several IFU data cubes, e.g., coaddition

"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii
from .oscube import OsCube


class CubeSet(list):
    """

    The Cubeset class is effectively just a list of Spec1d instances that
    includes operations that act on the whole list (e.g., coadds)

    """

    # ------------------------------------------------------------------------

    def __init__(self, inlist, informat=None, verbose=True):
        """

        There are three ways to create a Cubeset instance

        Option 1: List of filenames
        Option 2: List of already existing OsCubeinstances
        Option 3: A filename, where the file contains within it a list
         of input files

        """

        """
        First check to see if inlist is a list, of either filenames (strings)
         or OsCube instances (options 1 and 2).
        If it's not a list, check to see if it is a single string, in which
         case it may be a file that contains a list of filenames (option 3)
        """

        if isinstance(inlist, list):
            """ Option 1 or 2: list of files or of OsCube instances """
            if isinstance(inlist[0], str):
                for f in inlist:
                    try:
                        cube = OsCube(f)
                    except IOError:
                        print('')
                        print('Could not open requested file: %s' % f)
                        return
                    self.append(cube)

            elif isinstance(inlist[0], OsCube):
                for cube in inlist:
                    self.append(cube)

            else:
                raise TypeError('input list must be a list of filenames or'
                                ' OsCube instances')

        elif isinstance(inlist, str):
            """ Option 3: a file containing a list of filenames """
            try:
                if informat is not None:
                    intab = ascii.read(inlist, guess=False, format=informat)
                else:
                    intab = ascii.read(inlist)
            except IOError:
                print('')
                print('Could not open requested file: %s' % inlist)
                print('')
                return

            for f in intab.columns[0]:
                try:
                    cube = OsCube(f)
                except IOError:
                    print('')
                    print('Could not open requested file: %s' % f)
                    return
                self.append(cube)

    # ------------------------------------------------------------------------

    def coadd(self, doplot=True, outfile=None, verbose=True, **kwargs):
        """

        Do a variance-weighted sum of the spectra

        """

        """ Initialize """
        nx = self[0].wsize
        print(nx)
