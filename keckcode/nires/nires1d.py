"""

nires1d.py 

Code to handle NIRES echelle data, once it has been converted from the
NSX format into the Ech1d format
"""

import numpy as np
from astropy.table import Table
from specim.specfuncs import echelle1d

# ===========================================================================


class Nires1d(echelle1d.Ech1d):
    """

    A class used to visualize and analyze NIRES data, once they have been
    converted from the NSX format to Ech1d format

    """

    def __init__(self, inspec, informat='text', summary=True, verbose=True):
        """

        Initializes a Nires1d instance, essentially by initializing an
        Ech1d instance with the NIRES order information

        """

        """
        Define the information pertaining to the NIRES echelle orders

        Note that the pixmin and pixmax are not used at this point since
        any trimming of the orders should have been done in previous
        steps.
        """
        dtype = [('order', int), ('pixmin', int), ('pixmax', int)]
        oinfo = np.array([
                (7, 0, -1),
                (6, 0, -1),
                (5, 0, -1),
                (4, 0, -1),
                (3, 0, -1),
                ], dtype=dtype)
        ordinfo = Table(oinfo)

        """ Initialize by calling the parent class """
        super(Nires1d, self).__init__(inspec, informat=informat, 
                                      ordinfo=ordinfo, summary=summary,
                                      verbose=verbose)

